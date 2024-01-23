#' revise_ranges
#'
#' Function to generate a consensus age for assemblages
#' of fossil data in x, given a table of taxonomic
#' ranges. The need for error-checking is informed by
#' the error codes for the individual fossil
#' occurrences within each collection - if there is no
#' error, then the consensus age is unchanged. If
#' errors are present, then a consensus age for a
#' threshold proportion of taxa is searched for using
#' the overlap of the ranges for those taxa, as given
#' in range table y. Taxa whose occurrences lie outside
#' this consensus age are flagged as potential taxonomic
#' errors. If the threshold consensus partially overlaps
#' with the assemblage age, this overlap is returned to
#' present overzealous alteration of the age - otherwise
#' the complete consensus age is returned. If a
#' consensus age cannot be found, the original assemblage
#' age is returned, and each occurrence in the collection
#' flagged as potential taxonomic errors.
#' @param x Fossil occurrence data grouped into
#' spatiotemporally distinct assemblages
#' @param y A stratigraphic range dataset from which
#' consensus assemblage ages will be derived
#' @param assemblage The column name of the assemblage
#' groups in x
#' @param srt The column name of stratigraphic bases for
#' each element in both x and y - i.e. x and y must have
#' this same name for that column
#' @param end The column name of stratigraphic tops for
#' each element in both x and y - i.e. x and y must have
#' this same name for that column
#' @param taxon The column name denoting the taxon names
#' in both x and y - i.e. x and y must have
#' this same name for that column
#' @param err The column name flagging age errors for
#' occurrences in x. This allows 100$ valid assemblages to
#' be skipped. Age errors can be derived using @seealso
#' flag_ranges. All error codes must be one of: "000" -
#' unchecked, "R1R" - valid, "0R0" - both FAD and LAD exceeded,
#' "00R" - totally older than range, "R00" - totally younger
#' than range, "01R" - FAD exceeded, "1R0" - LAD exceeded.
#' If not supplied, all assemblages will be checked, even if
#' they are already valid a priori.
#' @param do.flag Rather than supplying error codes, should
#' flag_ranges be called internally to generate error codes
#' for supply to the rest of revise_ranges? As with err, this
#' is useful to prefilter individual occurrences, allowing
#' assemblages contain all valid, all unchecked or a mixture
#' of such error codes to be skipped. This can massively speed
#' up processing time for large datasets.
#' @param prop A numeric, between 0 and 1, denoting the
#' threshold percentage of taxa in the assemblage for
#' which a consensus age must be found
#' @param allow.zero A logical determining if, in the
#' case of a collection LAD being equal to the consensus
#' age FAD (i.e. a pointwise overlap), that pointwise
#' age will be taken as the revised age. The resultant
#' collection age will have no uncertainty as a result,
#' which may be unrealistic. The default behaviour is
#' FALSE, in which case pointwise overlaps will be
#' ignored and the revised age taken instead
#' @param verbose A logical determining if the progress
#' of the redating procedure should be reported
#' @return A list of two dataframes, the first recording
#' the results of the consensus redating procedure for
#' each assemblage in x, the second recording any flags
#' (if any) for each occurrence in x
#' @import data.table
#' @importFrom stats complete.cases
#' @export
#' @examples
#' # load datasets
#' data("brachios")
#' data("sepkoski")
#' # subsample brachios to make for a short example runtime
#' set.seed(1)
#' brachios <- brachios[sample(1:nrow(brachios), 1000),]
#' # rename columns in Sepkoski to match brachios
#' colnames(sepkoski)[4:6] <- c("genus", "max_ma", "min_ma")
#' # flag and resolve against the Sepkoski Compendium, collection-wise
#' revrng <- revise_ranges(x = brachios, y = sepkoski, do.flag = TRUE, verbose = TRUE,
#'                         taxon = "genus", assemblage = "collection_no",
#'                         srt = "max_ma", end = "min_ma")
#' # append the revised occurrence ages and error codes to the dataset
#' brachios$newfad <- revrng$occurrence$FAD
#' brachios$newlad <- revrng$occurrence$LAD
#' brachios$errcode <- revrng$occurence$status

revise_ranges <- function(x, y, assemblage = "collection_no", srt = "max_ma", end = "min_ma", taxon = "genus",
                           err = NULL, do.flag = FALSE, prop = 0.75, allow.zero = TRUE, verbose = TRUE) {

  # x = brachios_c
  # y = sepkoski_c
  # assemblage = "collection_no"
  # srt = "max_ma"
  # end = "min_ma"
  # taxon = "genus"
  # err = NULL
  # do.flag = TRUE
  # prop = 0.75
  # allow.zero = TRUE
  # verbose = TRUE


  if(!exists("x") | !exists("y")) {
    stop("Both x and y must be supplied")
  }
  if(!is.data.frame(x) | !is.data.frame(y)) {
    stop("Both x and y must be dataframes")
  }
  if(!assemblage %in% colnames(x)) {
    stop("assemblage must be a column name in x")
  }
  if(!all(c(taxon, srt, end) %in% colnames(x))) {
    stop("Arguments taxon, srt and end must all be the same column names in x and y")
  }
  if(!all(c(taxon, srt, end) %in% colnames(y))) {
    stop("Arguments taxon, srt and end must all be the same column names in x and y")
  }
  if(!all(is.numeric(x[,srt]), is.numeric(x[,end]), is.numeric(y[,srt]), is.numeric(y[,end]))) {
    stop("srt and end columns in x and y must all be numeric")
  }
  if(any(x[,srt] < x[,end])) {
    stop("One or more maximum ages in x are smaller than their corresponding minimum ages")
  }
  if(any(y[,srt] < y[,end])) {
    stop("One or more maximum ages in y are smaller than their corresponding minimum ages")
  }
  if(length(prop) != 1 | !is.numeric(prop)) {
    stop("prop must be a numeric of length 1")
  }
  if(prop > 1 | prop < 0) {
    stop("prop must be greater than 0 and less than or equal to 1")
  }
  if(prop < 0.5 & verbose) {
    warning("Prop is quite a low value - a minimum of 0.6 may be desirable")
  }
  if(!is.null(err) & do.flag) {
    warning("If err is supplied, then do.flag will be ignored")
    do.flag <- FALSE
  }
  if(!is.null(err) & isFALSE(do.flag)) {
    if(!err %in% colnames(x)) {
      stop("err must be a column name in x")
    }
    if(!all(unique(x[,err]) %in% c("000", "R1R", "0R0", "00R", "R00", "0R1", "1R0"))) {
      stop("err contains a non-standard error code (see documentation")
    }
  } else {
    err <- "age_flag"
    x$age_flag <- rep("0R0", times = nrow(x))
  }
  if(do.flag) {
    if(verbose) {cat("Performing pre-revision taxon flagging", "\n")}
    range_comp <- flag_ranges(x = x, y = y, xcols <- c(taxon, srt, end), verbose = verbose)
    x$age_flag <- range_comp$occurrence$status
    err <- "age_flag"
  }
  # global variable workaround
  . <- lb <- ub <- b <- N <- .N <- .SD <- .EACHI <- NULL
  if(verbose) {cat("Beginning age revision", "\n")}

  # tabulate error codes in each collection (ignoring uncoded and correct ages R1R, 000)
  codes <- c("000", "R1R", "0R0", "00R", "R00", "0R1", "1R0")
  tabs <- data.frame(unique(x[,assemblage]), NA, NA, NA, NA, NA, NA, NA)
  for(i in 1:length(codes)) {
    vals <- tapply(x[,err], x[,assemblage], function(x) {length(which(x == codes[i]))})
    tabs[,i + 1] <- vals[match(tabs$unique.x...assemblage.., as.numeric(names(vals)))]
  }
  colnames(tabs) <- c("assemblage", codes)

  # objects to store output
  z <- data.frame(unique(x[,assemblage]), NA, NA, "Unchecked", "000", NA)
  colnames(z) <- c("assemblage", "FAD", "LAD", "revision", "status", "prop")
  z$FAD <- x[match(z$assemblage, x[,assemblage]), srt]
  z$LAD <- x[match(z$assemblage, x[,assemblage]), end]
  FAD1 <- FAD <- x[,(srt)]
  LAD1 <- LAD <- x[,(end)]
  tax_flag <- FAD_diff <- LAD_diff <- FAD
  tax_flag[] <- "000"
  FAD_diff[] <- LAD_diff[] <- NA

  # find assemblages with errors
  to_do <- apply(tabs[,4:ncol(tabs)], 1, sum)
  to_do <- tabs$assemblage[which(to_do != 0)]
  # convert to data.table for rapid indexing
  x$rnum <- 1:nrow(x)
  x1 <- data.table::data.table(x)
  data.table::setkeyv(x1, assemblage)

  # if there are range chart taxa, test for an overlapping range, assigning Unchecked > Unresolved > Retained > Revised
  for(i in 1:length(to_do)) {

    # get assemblage occs
    occs <- as.data.frame(x1[.(to_do[i]), ])
    # trim to unique, named taxa and get (if any) matches in the comparative database
    occs2 <- unique(occs[, c(taxon, srt, end)])
    occs2 <- occs2[stats::complete.cases(occs2), ]
    foo <- y[match(occs2[, taxon], y[, taxon]), c(taxon, srt, end)]
    foo <- foo[stats::complete.cases(foo), ]
    zpos <- match(to_do[i], z$assemblage)

    # set defaults (default age solution is current assemblage age)
    sol_out <- ores <- unlist(z[i, c("FAD", "LAD")])
    tprop <- NA
    revise <- TRUE
    revision <- "Revised"

    # if there are taxa to use for checking
    if (nrow(foo) != 0) {

      # find the overlaps
      dt <- data.table(lb = foo[[3]], ub = foo[[2]])
      mdt <- dt[, .(b = unique(unlist(.SD)))]
      sol <- dt[mdt, on = .(lb <= b, ub >= b), .N, by = .EACHI]
      sol$N <- sol$N / nrow(foo)

      # if there are any solutions above the threshold
      if(any(sol$N >= prop)) {

        # if zero-length spans are allowed
        if(allow.zero) {

          # get the most inclusive solution
          res <- rev(sol[N == max(N), range(lb)])
          # if the solution is no better than the original age
          if(res[1] >= ores[1] & res[2] <= ores[2]) {
            revision <- "Retained"
            tprop <- 1
          # otherwise check for overlaps, taking the overlap if present
          } else {
            sol2 <- rbind(ores, res)
            sol2 <- sol2[order(sol2[, 1], decreasing = TRUE),]
            if (sol2[2, 1] >= sol2[1, 2]) {
              sol_out <- c(min(sol2[, 1]), max(sol2[, 2]))
            } else {
              sol_out <- res
            }
            tprop <- max(sol$N)
          }

        # if zero-length spans are not allowed
        } else {

          # check non-zero combinations where each boundary exceeds the proportion
          nz_sols <- list()
          sol <- as.matrix(sol)
          sol <- sol[which(sol$N >= prop),]
          for(j in 1:nrow(sol)) {
            ob <- cbind(rep(sol[j,1], nrow(sol)), sol[,2], (rep(sol[j,3], nrow(sol)) + sol[,3]) / 2)
            ob <- ob[!ob[,1] == ob[,2],]
            if(nrow(ob) != 0) {nz_sols[[j]] <- ob[which.max(ob[,3]), ,drop = FALSE]}
          }

          # if there is still a solution
          nz_sols <- do.call(rbind, nz_sols)
          if(!is.null(nz_sols)) {
            res <- as.vector(nz_sols[which.max(nz_sols[,3]),1:2])

            # if the solution is no better than the original age
            if(res[1] >= ores[1] & res[2] <= ores[2]) {
              revision <- "Retained"
              tprop <- 1
              # otherwise check for overlaps, taking the overlap if present
            } else {
              sol2 <- rbind(ores, res)
              sol2 <- sol2[order(sol2[, 1], decreasing = TRUE),]
              if (sol2[2, 1] > sol2[1, 2]) {
                sol_out <- c(min(sol2[, 1]), max(sol2[, 2]))
              } else {
                sol_out <- res
              }
              tprop <- max(nz_sols$N)
            }

          # otherwise the age goes unresolved (no non-zero solutions)
          } else {
            revise <- FALSE
            revision <- "Unresolved"
          }
        }

        # otherwise the age goes unresolved (no solutions above the threshold)
      } else {
        revise <- FALSE
        revision <- "Unresolved"
      }

      # otherwise the age goes unchecked (no taxa in the database to check with)
    } else {
      revise <- FALSE
      revision <- "Unchecked"
    }

    # if a solution has been found, revise accordingly
    if(revise) {

      # update assemblage and occurrence ages and revision status
      z[zpos, c("FAD", "LAD")] <- sol_out
      z[zpos, "revision"] <- revision
      z[zpos, "prop"] <- tprop
      FAD[occs$rnum] <- sol_out[1]
      LAD[occs$rnum] <- sol_out[2]

      # flag occurrences if the consensus is not 100%
      if(tprop != 1) {

        # occurrences exceeding FAD or LAD
        foo2 <- foo[which((foo$max_ma > sol_out[1] & foo$min_ma >= sol_out[1])),taxon]
        foo3 <- foo[which((foo$max_ma <= sol_out[2] & foo$min_ma < sol_out[2])),taxon]
        rnum1 <- occs$rnum[which(occs[,taxon] %in% foo2)]
        rnum2 <- occs$rnum[which(occs[,taxon] %in% foo3)]

        # taxon-wise flagging
        tax_flag[rnum1] <- "00R"
        tax_flag[rnum2] <- "R00"
        FAD_diff[rnum1] <- abs(FAD[rnum1] - FAD1[rnum1])
        LAD_diff[rnum2] <- abs(LAD[rnum2] - LAD1[rnum2])

        # collection-wise flagging
        if(length(foo2) != 0 & length(foo3) != 0) {
          coll_flag <- "0R0"
          z[zpos, "status"] <- "0R0"
        } else {
          if(length(foo2) != 0) {
            coll_flag <- "00R"
            z[zpos, "status"] <- "00R"
          }
          if(length(foo3) != 0) {
            coll_flag <- "R00"
            z[zpos, "status"] <- "R00"
          }
        }

      # otherwise the revision is completely consistent
      } else {
        tax_flag[occs$rnum] <- "R1R"
        tax_flag[occs$rnum] <- "R1R"
        z[zpos, "status"] <- "R1R"
      }

    # otherwise assign unrevised/unchecked
    } else {
      z[zpos, c("FAD", "LAD")] <- sol_out
      z[zpos, "revision"] <- revision
      z[zpos, "prop"] <- tprop
    }

    # notify R
    if(verbose) {
      if(i != 1) {cat(paste0("\r"))}
      cat(paste0("Assemblage ", i, "/", length(to_do), " checked"))
      if(i == length(to_do)) {cat("\n")}
    }
  }

  # return
  if(verbose) {message("See $occurrence in output for the revised ages of individual occurrences\n")}
  per_occ <- cbind.data.frame(FAD, LAD, FAD_diff, LAD_diff, tax_flag)
  colnames(per_occ) <- c("FAD", "LAD", "FAD_diff", "LAD_diff", "tax_flag")
  out <- list()
  out[[1]] <- z
  out[[2]] <- per_occ
  names(out) <- c("collection", "occurrence")
  return(out)
}
