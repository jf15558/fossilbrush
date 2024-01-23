#' spell_check
#'
#' Function for checking for potential synonyms with alternate
#' spellings. Synonyms are checked for within group using using
#' a Jaro Winkler string distance matrix. Potential synonyms are
#' selected using the jw threshold. These can then be further
#' filtered by the number of shared letters at the beginning and
#' end of the a synonym pair, and by prefixes or suffixes which
#' may give erroneously high similarities.
#'
#' @param x a dataframe containing a column with terms, and a
#' further column denoting the groups within which terms will
#' be checked against one another. If supplying a dataframe with
#' just these columns, terms should be column 1
#' @param terms a character vector of length 1, specifying the
#' terms column in x. This is required if x contains more than
#' two columns. Alternatively, if x is not provided, terms can
#' be a character vector. If groups are not specified, all
#' elements of terms will be treated as part of the same group
#' @param groups a character vector of length 1, specifying the
#' groups column in x. This is required if x contains more than
#' two columns. Alternatively, if terms is supplied as a
#' character vector, groups can also be supplied in the same way
#' to denote their groups
#' @param jw a numeric greater than 0 and less than 1. This is
#' the distance threshold below which potential synonyms will be
#' considered
#' @param str A positive integer specifying the
#' number of matching characters at the beginning of synonym
#' pairs. By default 1, i.e. the first letters must match
#' @param str2 If not NULL, a positive integer specifying the
#' number of matching characters at the end of synonym pairs
#' @param alternative A character string of length one corresponding
#' to one of the methods used by @seealso afind. One of "osa",
#' "lv", "dl", "hamming", "lcs", "qgram", "cosine",
#' "running_cosine", "jaccard", or "soundex".
#' @param q q-gram size. Only used when alternative is "qgram",
#' "cosine" or "Jaccard".
#' @param pref If not NULL, a character vector of prefixes which
#' may result in erroneously low JW distances. Synonyms will only
#' be considered if both terms share the same prefix
#' @param suff If not NULL, a character vector of suffices which
#' may result in erroneously low JW distances. Synonyms will only
#' be considered if both terms share the same suffix
#' @param exclude If not NULL, a character vector of group names
#' which should be skipped - useful for groups which are known
#' to contain potentially similar terms
#' @param verbose A logical determining if function progress be reported using the
#' pbapply progress bar
#' @return a dataframe of synonyms (cols 1 and 2), the group in
#' which they occur, the frequencies of each synonym in the dataset
#' and finally the q-gram difference between the synonyms
#' @importFrom stats na.omit
#' @import stringdist
#' @import pbapply
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # define suffixes
#' b_suff <- c("ina", "ella", "etta")
#' # run function
#' spl <- spell_check(brachios, terms = "genus", groups = "family", suff = b_suff)

spell_check <- function(x, terms = NULL, groups = NULL, jw = 0.1, str = 1, str2 = NULL, alternative = "jaccard", q = 1,
                        pref = NULL, suff = NULL, exclude = NULL, verbose = TRUE) {

  # check arguments
  if(!exists("x")) {
    if(is.null(terms)) {
      stop("Please supply a dataframe with a column of names and a column denoting their groups,
         or terms (and optionally groups) as a vector")
    }
    if(!is.vector(terms)) {
      stop("If x is not specified, then terms must be a character vector")
    }
    if(is.character(terms)) {
      stop("Terms must be of class character")
    }
    if(is.null(groups)) {
      groups <- rep("group", length(terms))
    } else {
      if(length(groups) != length(terms)) {
        stop("Groups must be the same length as terms")
      }
    }
    x <- cbind.data.frame(terms, groups)
  }
  if(exists("x")) {
    if(ncol(x) > 2) {
      if(is.null(terms) || is.null(groups)) {
        stop("If x contains more than two columns, terms and groups must be specified")
      }
      if(!is.character(groups) && !is.character(terms)) {
        stop("terms and groups should both be character vectors of length 1")
      }
      if(length(terms) > 1) {
        warning("terms is of length > 1. Only the first element will be used")
        terms <- terms[1]
      }
      if(length(groups) > 1) {
        warning("groups is of length > 1. Only the first element will be used")
        groups <- groups[1]
      }
      if(!all(c(terms, groups) %in% colnames(x))) {
        stop("terms and groups must be column names of x")
      }
      x <- x[,c(terms, groups)]
    }
  }
  if(!all(is.character(x[,1]), is.character(x[,2]))) {
    stop("term and group columns must be of mode character")
  }
  if(!is.numeric(jw)) {
    stop("jw must be a numeric, greater than 0, less than 1")
  }
  if(!is.null(str)) {
    if(!is.numeric(str)) {
      stop("str must be a positive integer, or NULL")
    }
  }
  if(!is.null(str2)) {
    if(!is.numeric(str)) {
      stop("str must be a positive integer, or NULL")
    }
  }
  if(length(alternative) > 1) {
    warning("alternative contains more than one element. Only the first will be used")
    alternative <- alternative[1]
  }
  if(!alternative %in% c("osa", "lv", "dl", "hamming", "lcs", "qgram",
                     "cosine", "running_cosine", "jaccard", "soundex")) {
    stop("alternative must be one of 'osa', 'lv', 'dl', 'hamming', 'lc', 'qgram',
                     'cosine', 'running_cosine', 'jaccard', 'soundex'")
  }
  if(alternative %in% c("qgram", "cosine", "jaccard")) {
    if(length(q) > 1) {
      warning("q contains more than one element. Only the first will be used")
      q <- q[1]
    }
    if(!q >= 1) {
      stop("q must be an integer greater than or equal to 1")
    }
    q <- round(q)
  }
  if(!is.null(pref)) {
    if(!is.character(pref)) {
      stop("pref must be of class character")
    }
  }
  if(!is.null(suff)) {
    if(!is.character(suff)) {
      stop("suff must be of class character")
    }
  }
  if(!is.null(exclude)) {
    if(!is.character(exclude)) {
      stop("exclude must be of class character")
    }
  }
  if(!verbose) {
    baseopt <- getOption("pboptions")
    opb <- pboptions(type = "none")
  }

  # apply the comparison algorithm groupwise
  gps <- unique(na.omit(x[,2]))
  gps <- gps[!gps %in% exclude]
  sp <- pbsapply(gps, simplify = FALSE, function(y) {

    # all elements of x which belong to group y in gps
    ob <- x[x[,2] == y,]
    ob <- unique(na.omit(ob[,1]))
    # if there is are not multiple names in the group, skip
    if(length(ob) < 2) {
      flag <- NULL
    }

    else {
      # else get the JW distance matrix for the elements in the group
      test <- stringdistmatrix(a = ob, b = ob, method = "jw")
      colnames(test) <- rownames(test) <- ob
      # subset to those which exceed the threshold for matching
      flag <- which(test < jw, arr.ind = TRUE)
      # remove self matches (dist = 0)
      flag <- flag[which(flag[,1] != flag[,2]),]
      # if there are no remaining flagged names
      if(length(flag) == 0) {
        flag <- NULL
      }

      else {
        # retrieve names
        flag <- cbind(ob[flag[,1]], ob[flag[,2]], y)
        # drop equivalent rows (xy, yx pairs)
        eq <- duplicated(t(apply(flag, 1, function(z) {paste0(z[order(z)])})))
        flag <- flag[!eq,,drop = FALSE]
        # cull by first y letter non-matches
        if(!is.null(str)) {
          c1 <- substr(flag[,1], start = 1, stop = str)
          c2 <- substr(flag[,2], start = 1, stop = str)
          flag <- flag[which(c1 == c2), , drop = FALSE]
        }
        # cull by last y letter non-matches
        if(!is.null(str2)) {
          c1 <- substr(flag[,1], start = (nchar(flag[,1]) - (str2 - 1)), stop = nchar(flag[,1]))
          c2 <- substr(flag[,2], start = (nchar(flag[,2]) - (str2 - 1)), stop = nchar(flag[,2]))
          flag <- flag[which(c1 == c2), , drop = FALSE]
        }
        if(length(flag) == 0) {
          flag <- NULL
        }

        else {
          # cull by some common prefixes and suffixes which give high similarity, but are different
          common <- c(paste0("^", pref), paste0(suff, "$"))
          for(j in 1:length(common)) {
            # grep the common suffix in the first column
            g1 <- grepl(common[j], flag[,1])
            # grep that same suffix for the second column
            g2 <- grepl(common[j], flag[,2])
            # remove the non matching elements from the dataframe
            to_drop <- which(g1 != g2)
            if(length(to_drop) != 0) {
              flag <- flag[-to_drop, , drop = FALSE]
              if(length(flag) == 0) {
                flag <- NULL
              }
            }
          }
        }
      }
    }
    out <- flag
  })

  # reformat from list
  err <- sp[!unlist(lapply(sp, is.null))]
  err <- as.data.frame(do.call(rbind, err))
  err$f1 <- as.vector(table(x[,terms])[match(err$V1, names(table(x[,terms])))])
  err$f2 <- as.vector(table(x[,terms])[match(err$V2, names(table(x[,terms])))])

  # do qgram score
  val <- apply(err, 1, function(y) {afind(y[1], y[2], method = alternative, q = q)$distance})
  err2 <- cbind.data.frame(err, val)

  # return
  row.names(err2) <- NULL
  if(nrow(err2) == 0) {
    if(verbose) {message("No synonyms flagged, returning NULL")}
    err2 <- NULL
  } else {
    colnames(err2) <- c("t1", "t2", "group",
                        "freq1", "freq2", "m2")
  }
  if(!verbose) {opt <- pboptions(baseopt)}
  return(err2)
}
