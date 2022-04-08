#' check_taxonomy
#'
#' Wrapper functions to implement a multi-step cleaning routine
#' for hierarchically structured taxonomic data. The first part
#' of the routine calls @seealso format_check to perform a few
#' presumptive checks on all columns, scanning for non-letter
#' characters and checking the number of words in each string.
#' By default, @seealso clean_name is called to ensure correct
#' formatting as this improves downstream checking. The second
#' part of the routine calls @seealso spell_check to flag spelling
#' discrepancies between names within a given taxonomic group. If
#' chosen, the function can automatically impose the more frequent
#' spelling. The third part of the routine calls @seealso discrete_ranks
#' to flag name re-use at different taxonomic levels. Some of these
#' cases may arise when a name has been unfortunately, (although
#' permissibly) used to refer to groups at different taxonomic levels,
#' or where a higher classification may have been inserted as a
#' placeholder for a missing lower classification. The fourth part of
#' the routine calls @seealso find_duplicates to flag variable higher
#' classifications for a given taxon, including cases where a higher
#' classification is missing for one instance of a taxon, but present
#' for the others. If chosen, @seealso resolve_duplicates is called
#' to ensure a consistent classification is imposed. For cases where
#' a name has been re-used at the same rank for genuinely different
#' taxa (not permissible, unlike name re-use at different ranks)
#' suffixes are added as capital letters, e.g. TaxonA, TaxonB. If
#' any of the automatic cleaning routines are employed (again the
#' default behaviour as clean_name is TRUE by default), the function
#' will return are a cleaned version of the dataset. If the use of
#' suffixes from @seealso resolve_duplicates is not desirable,
#' the function behaviour can be altered so that any suffixes are
#' dropped before returning.
#'
#' * Data supply arguments *
#' @param x A dataframe with hierarchically organised
#' taxonomic information. If x only comprises the taxonomic
#' information, @param ranks does not need to be specified, but the
#' columns must be in order of decreasing taxonomic rank
#' @param ranks The column names of the taxonomic data fields
#' in x. These must be provided in order of decreasing taxonomic
#' rank
#' @param species A logical indicating if x contains a species
#' column. As the data must be supplied in hierarchical order,
#' this column will naturally be the last column in x and
#' species-specific spell checks will be performed on this column.
#' @param species_sep A character vector of length one specifying
#' the genus name and specific epithet in the species column, if
#' present
#'
#' * Flagging routine arguments *
#' @param routine A character vector determining the flagging
#' and cleaning routines to employ. Valid values are format_check (check for
#' non letter characters and the number of words in names), spell_check
#' (flag potential spelling errors), discrete_rank (check that
#' taxonomic names are unique to their rank), duplicate_tax (flag
#' conflicting higher classifications of a given taxon)
#' @param report A logical of length one determining if the flagging outputs
#' of each cleaning routine should be returned to the user for inspection.
#' This is different to @param verbose, which controls whether flagging
#' should additionally be reported to the user on the console
#' @param verbose A logical determining if function progress and flagged
#' errors should be reported to the console
#'
#' * Cleaning routine arguments *
#' @param clean_name If TRUE, the function will return cleaned versions of
#' the columns in x using the routines in @seealso clean_name. These routines
#' can be altered using the 'term_set' and 'collapse_set' arguments.
#' @param clean_spell If TRUE, the function will return a cleaned version
#' of the supplied taxonomic dataframe, using the supplied threshold
#' for the similarity method given by method2, to automatically update
#' any names in pairs of flagged synonyms to the more frequent spelling.
#' This is not recommended, however, so the argument is FALSE by default
#' and the threshold left as NULL
#' @param thresh The threshold for the similarity method given by method2,
#' below which flagged pairs of names will be considered synonyms and
#' resolved automatically. See @seealso spell_check for details on method2
#' @param resolve_duplicates If TRUE, the function will return a cleaned version
#' of the supplied taxonomic dataframe, using @seealso resolve_duplicates
#' to resolve conflicts in the way documented by the function. Both
#' spell_clean and tax_clean can both be TRUE to return a dataset cleaned
#' by both methods
#' @param append If TRUE, any suffixes used during cleaning will
#' be retained in the cleaned version of the data. This is preferable
#' as it ensures that all taxonomic names are rank-discrete and
#' uniquely classified
#'
#' * Routine specific arguments *
#' @param term_set A character vector of terms (to be
#' used at all ranks) or a list of rank-specific terms
#' which will be supplied, element-wise as the @param collapse argument
#' called by @seealso clean_name. If a list, this
#' @param collapse_set A character vector of character strings (to be
#' used at all ranks) or a list of rank-specific strings
#' which will be supplied, element-wise as the @param collapse argument
#' called by @seealso clean_name. If a list, this
#' should be given in descending rank order
#' @param jw Called by @seealso spell_check
#' @param str Called by @seealso spell_check
#' @param str2 Called by @seealso spell_check
#' @param alternative Called by @seealso spell_check
#' @param q Called by @seealso spell_check
#' @param pref_set A character vector of prefixes (which
#' will be used at all ranks) or a list of rank-specific prefixes,
#' which will be supplied, element-wise as the @param pref argument
#' called by @seealso spell_check. If a list, this
#' should be given in descending rank order.
#' @param suff_set A character vector of suffixes (which
#' will be used at all ranks) or a list of rank-specific suffixes,
#' which will be supplied, element-wise as the @param suff argument
#' called by @seealso spell_check. If a list, this
#' should be given in descending rank order.
#' @param exclude_set A character vector of terms to exclude (which
#' will be used at all ranks) or a list of rank-specific exclusion terms,
#' which will be supplied, element-wise as the @param exclude argument
#' called by @seealso spell_check. If a list, this
#' should be given in descending rank order.
#' @param jump Called by @seealso resolve_duplicates
#' @param plot Called by @seealso resolve_duplicates
#' @return A list with elements corresponding to the outputs of the chosen
#' flagging routines (four by default: $formatting, $synonyms, $ranks, $duplicates),
#' plus a cleaned verison of the data ($data) if any of clean_name, clean_spell
#' or resolve_duplicates are TRUE. See @seealso format_check, @seealso spell_clean,
#' @seealso discrete_ranks and @seealso find_duplicates for details of the structure
#' of the flagging outputs
#' @importFrom stats na.omit
#' @export

check_taxonomy <- function(x, ranks = c("phylum", "class", "order", "family", "genus"), species = FALSE, species_sep = NULL,
                           routine = c("format_check", "spell_check", "discrete_ranks", "find_duplicates"), report = TRUE, verbose = TRUE,
                           clean_name = FALSE, clean_spell = FALSE, thresh = NULL, resolve_duplicates = FALSE, append = TRUE,
                           term_set = NULL, collapse_set = NULL,
                           jw = 0.1, str = 1, str2 = NULL, alternative = "jaccard", q = 1, pref_set = NULL, suff_set = NULL, exclude_set = NULL,
                           jump = 3, plot = FALSE) {

  #x = occs
  #ranks = c("clade", "genus")
  #species = FALSE
  #species_sep = NULL
  #routine = c("format_check", "spell_check", "discrete_ranks", "find_duplicates")
  #report = TRUE
  #verbose = TRUE
  #clean_name = FALSE
  #clean_spell = FALSE
  #thresh = NULL
  #resolve_duplicates = FALSE
  #append = TRUE
  #term_set = NULL
  #collapse_set = NULL
  #jw = 0.1
  #str = 1
  #str2 = NULL
  #alternative = "jaccard"
  #q = 1
  #pref_set = NULL
  #suff_set = NULL
  #exclude_set = NULL
  #jump = 3
  #plot = FALSE

  ######## ARG CHECKS ########

  # check that data has minimally been supplied
  if(!exists("x")) {
    stop("Please supply x as a dataframe of taxonomic assignments")
  }
  # coerce to dataframe with column names to be safe
  if(!is.data.frame(x)) {x <- as.data.frame(x)}
  if(is.null(colnames(x))) {colnames(x) <- as.character(1:ncol(x))}

  # check that ranks are column names of x
  if(is.null(ranks)) {ranks <- colnames(x)}
  if(!all(ranks %in% colnames(x))) {
    stop("Not all elements of argument ranks are column names in x")
  }
  # check that ranks are in hierarchical order
  if(length(unique(x[,ranks[length(ranks)]])) < length(unique(x[,ranks[(length(ranks) - 1)]]))) {
    warning("Higher taxonomy is more diverse than lower taxonomy. Are the columns in x
            or the column names specified in 'ranks' supplied in descending hierarchical order?")
  }
  # subset to columns
  x1 <- x
  x <- x[,ranks]

  # check that the data is character
  if(!all(apply(x, 2, class) == "character")) {
    stop("Not all columns in x are of class character")
  }
  # check species designator
  if(!is.logical(species) & length(species) != 1) {
    stop("Species should be a logical of length one, indicating whether species-level designations are present in x")
  }

  # check cleaning routines have been correctly supplied
  if(!is.character(routine)) {
    stop("Routine should be a character vector containing one or more of the following:
         clean_name, spell_check, discrete_ranks, find_duplicates")
  }
  if(!any(routine %in% c("format_check", "spell_check", "discrete_ranks", "find_duplicates"))) {
    stop("All elements of argument routine are invalid.
         Valid elements are format_check, spell_check, discrete_ranks, find_duplicates")
  }
  if(!all(routine %in% c("format_check", "spell_check", "discrete_ranks", "find_duplicates"))) {
    warning("Some elements of argument routine are invalid and will be ignored.
            Valid elements are format_check, spell_check, discrete_ranks, find_duplicates")
  }
  # ensure routine vector is clean and correctly ordered
  routine <- unique(routine)
  routine <- as.vector(na.omit(routine[match(c("format_check", "spell_check", "discrete_ranks", "find_duplicates"), routine)]))

  # check additional flags
  if(clean_name) {

    if(is.null(term_set)) {
      term_set <- as.list(1:ncol(x))
      terms_list <- lapply(term_set, function(x) {x <- NULL})
    } else {
      if(is.atomic(term_set)) {
        terms_list <- lapply(1:ncol(x), function(x) {x <- term_set})
      }
      if(is.list(term_set)) {terms_list <- term_set}
      if(length(terms_list) != length(ranks)) {
        stop("term_set should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if (!all(unlist(lapply(terms_list, class)) %in% c("NULL", "character"))) {
        stop("Not all elements of term_set are of class character")
      }
      terms_list <- lapply(terms_list, function(x) {as.vector(na.omit(x))})
    }

    if(is.null(collapse_set)) {
      collapse_set <- as.list(1:ncol(x))
      collapse_list <- lapply(collapse_set, function(x) {x <- NULL})
    } else {
      if(is.atomic(collapse_set)) {
        collapse_list <- lapply(1:ncol(x), function(x) {x <- collapse_set})
      }
      if(is.list(collapse_set)) {collapse_list <- collapse_set}
      if(length(collapse_list) != length(ranks)) {
        stop("Collapse should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if (!all(unlist(lapply(collapse_list, class)) %in% c("NULL", "character"))) {
        stop("Not all elements of collapse are of class character. Additionally, any regex special characters must be escaped using backslashes")
      }
      collapse_list <- lapply(collapse_list, function(x) {as.vector(na.omit(x))})
    }
  }

  # check additional flags
  if("spell_check" %in% routine) {

    # check any supplied prefixes
    if(is.null(pref_set)) {
      pref_set <- as.list(1:ncol(x))
      pref_list <- lapply(pref_set, function(x) {x <- NULL})
    } else {
      if(is.atomic(pref_set)) {
        pref_list <- lapply(1:ncol(x), function(x) {x <- pref_set})
      }
      if(is.list(pref_set)) {pref_list <- pref_set}
      if(length(pref_list) != length(ranks)) {
        stop("Prefixes should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if (!all(unlist(lapply(pref_list, class)) %in% c("NULL", "character"))) {
        stop("Not all elements of pref are of class character")
      }
      pref_list <- lapply(pref_list, function(x) {as.vector(na.omit(x))})
    }

    # check any supplied suffixes
    if(is.null(suff_set)) {
      suff_set <- as.list(1:ncol(x))
      suff_list <- lapply(suff_set, function(x) {x <- NULL})
    } else {
      if(is.atomic(suff_set)) {
        suff_list <- lapply(1:ncol(x), function(x) {x <- suff_set})
      }
      if(is.list(suff_set)) {suff_list <- suff_set}
      if(length(suff_list) != length(ranks)) {
        stop("Suffixes should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if (!all(unlist(lapply(suff_list, class)) %in% c("NULL", "character"))) {
        stop("Not all elements of suff are of class character")
      }
      suff_list <- lapply(suff_list, function(x) {as.vector(na.omit(x))})
    }

    # check any supplied exclusions
    if(is.null(exclude_set)) {
      exclude_set <- as.list(1:ncol(x))
      exclude_list <- lapply(exclude_set, function(x) {x <- NULL})
    } else {
      if(is.atomic(exclude_set)) {
        exclude_list <- lapply(1:ncol(x), function(x) {x <- exclude_set})
      }
      if(is.list(exclude_set)) {exclude_list <- exclude_set}
      if(length(exclude_list) != length(ranks)) {
        stop("Exclusions should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if (!all(unlist(lapply(exclude_list, class)) %in% c("NULL", "character"))) {
        stop("Not all elements of exc are of class character")
      }
      exclude_list <- lapply(exclude_list, function(x) {as.vector(na.omit(x))})
    }

    if(clean_spell) {
      if(verbose) {warning("As spell checking is approximate, automatic resolution is not recommended. Any flagged names should be properly checked")}
      if(is.null(thresh)) {
        stop("If spell_clean has been requested, thresh must be specified. This threshold is specific to method2, see spell_check documentation")
      }
      if(!is.numeric(thresh) | any(thresh < 0)) {
        stop("Thresh must be a positive numeric to be used at all ranks, or vector of values to be used individually at each rank. See spell_check for details on value choice")
      }
      if(length(thresh == 1)) {
        thresh <- rep(thresh, length(ranks))
      }
      if(length(thresh) != length(ranks)) {
        stop("Thresh must be a positive numeric to be used at all ranks, or vector of values to be used individually at each rank. See spell_check for details on value choice")
      }
    }
  }
  # set up output list and timer
  stage <- 1
  out <- list()

  ######## FORMATTING ########

  # check formatting
  if("format_check" %in% routine) {

    message(paste0("Checking formatting [", stage, "/", length(routine), "]"))
    stage <- stage + 1

    out[[2]] <- format_check(x = x, ranks = ranks, species = species, species_sep = species_sep, verbose = FALSE)
    names(out)[2] <- "formatting"
    if(verbose) {
      if(max(lengths(out[[2]][[1]])) > 0 | max(lengths(out[[2]][[2]])) > 0) {
        cat(" - formatting errors detected (see $formatting in output)\n")
      } else {
        cat(" - no formatting errors detected\n")
      }
    }
    if(clean_name) {
      for(i in 1:length(ranks)) {
        if(i != 1) {cat(paste("\r"))}
        cat(paste0(" + cleaning names at rank ", ranks[i], "        "))
        if(i == length(ranks)) {cat(paste("\n"))}
        x[,i] <- clean_name(x[,i], terms = terms_list[[i]], collapse = collapse_list[[i]], verbose = FALSE)
      }
    }
  }

  ######## SYNONYMS ########

  # check spelling
  if("spell_check" %in% routine) {

    message(paste0("Checking spelling   [", stage, "/", length(routine), "]"))
    stage <- stage + 1

    out[[3]] <- 1
    names(out)[3] <- "synonyms"
    spell_list <- list()
    for(i in 1:(length(ranks) - 1)) {

      foo <- spell_check(x = x, terms = ranks[i + 1], groups = ranks[i],
                         jw = jw, str = str, str2 = str2, alternative = alternative,
                         q = q, pref = pref_list[[i + 1]], suff = suff_list[[i + 1]], exclude = exclude_list[[i + 1]], verbose = FALSE)
      spell_list[[i]] <- foo
      if(!is.null(foo)) {
        spell_list[[i]] <- cbind.data.frame(level = rep(ranks[i + 1], nrow(spell_list[[i]])), spell_list[[i]])
      }
    }
    if(verbose) {
      if(length(spell_list) != 0) {
        cat(" - potential synonyms detected (see $synonyms in output)\n")
      } else {
        cat(" - no potential synonyms detected\n")
      }
    }
    out[[3]] <- do.call(rbind, spell_list)

    if(clean_spell) {
      cat(paste0(" + resolving synonyms by frequency\n"))

      ob <- out[[3]]
      if(nrow(ob) != 0) {
        for(j in 1:nrow(ob)) {
          if(ob$m2[j] < thresh[match(ob$level[j], ranks)]) {
            x[which(x[,ob$level[j]] == c(ob$t1[j], ob$t2[j])[which.min(c(ob$freq1[j], ob$freq2[j]))]), ob$level[j]] <-
              c(ob$t1[j], ob$t2[j])[which.max(c(ob$freq1[j], ob$freq2[j]))]
          }
        }
      }
    }
  }

  ######## RANKS ########

  # check rank discretion (only reporting if requested - done anyway as a prerequisite for resolve_duplicates)
  rank_check <- discrete_ranks(x, ranks = rev(ranks))
  if("discrete_ranks" %in% routine) {

    message(paste0("Checking ranks      [", stage, "/", length(routine), "]"))
    stage <- stage + 1

    out[[4]] <- rank_check
    names(out)[4] <- "ranks"
    if(verbose) {
      if(any(unlist(lapply(out[[4]][[1]], length))) > 0 | any(unlist(lapply(out[[4]][[2]], length))) > 0) {
        cat(" - cross-rank names detected (see $ranks in output)\n")
      } else {
        cat(" - no cross-rank names detected\n")
      }
    }


  }

  ######## CONFLICTS ########

  # check classification
  if("find_duplicates" %in% routine) {

    message(paste0("Checking taxonomy   [", stage, "/", length(routine), "]"))
    stage <- stage + 1

    out[[5]] <- find_duplicates(x = x, ranks = ranks)
    names(out)[5] <- "duplicates"

    if(verbose) {
      if(!is.null(out[[5]])) {
        cat(" - conflicting classifications detected (see $duplicates in output)\n")
      } else {
        cat(" - no conflicting classifications detected\n")
      }
    }

    if(resolve_duplicates) {

      # paste a numeric at the start of each column to ensure rank discretion
      for(i in 1:ncol(x)) {
        blank <- which(is.na(x[,i]))
        x[,i] <- paste0(i, x[,i])
        x[blank,i] <- NA
      }
      x <- resolve_duplicates(x = x, ranks = ranks, jump = jump, verbose = verbose)
      # remove numeric
      for(i in 1:length(ranks)) {
        x[,i] <- gsub("[0-9]", "", x[,i])
      }
    }
  }

  ######## RETURN ########

  # format output
  cleaned <- FALSE
  if(clean_name | clean_spell | resolve_duplicates) {
    cleaned = TRUE
    # remove the duplicate classification suffixes if requested
    if(!append) {x <- apply(x, 2, function(y) {gsub("[A-Z]$", "", x = y)})}
    x1[,ranks] <- x
    out[[1]] <- x1
    names(out)[1] <- "data"
    if(!report) {
      out <- x1
    }
    if(verbose & report & any(clean_name, clean_spell, resolve_duplicates)) {message("See $data in output for the cleaned dataset")}
  }
  # remove null elements of list if needed
  if(is.list(out)) {
    to_remove <- unlist(lapply(out, is.null))
    if(sum(to_remove) > 0) {out <- out[!to_remove]}
  }
  return(out)
}
