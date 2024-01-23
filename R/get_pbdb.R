#' get_pbdb
#'
#' Function for downloading Paleobiology Database (PBDB) data
#' (saved to disk and/or imported into R) or generating PBDB
#' API-compatible URLs. If downloading data over timespans
#' greater than 100 Ma, the download is performed in 100 Ma
#' chunks to better track the download progress.
#'
#' @param taxon A character vector of taxon names. Prepending
#' a taxon name with ^ will exclude it from the PBDB search.
#' Alternatively @param ex_taxon can be used to do this
#' @param interval A numeric vector of length two with positive
#' ages in Ma, or a character vector containing one or two ICS
#' chronostratigraphic interval names
#' @param mode A character vector of length one specifying the
#' type of data to return: one of occurrence, collection, taxa,
#' specimen, measurement, strata, diversity, opinion or reference
#' @param res A character vector of length one specifying the
#' taxonomic resolution of the dataset: one of all, family, genus
#' species, lump_genus or lump_subgen. The latter two lump
#' multiple occurrences of genera or subgenera within collections
#' into a single representative occurrence
#' @param fields A character vector of PBDB vocabulary for
#' additional data fields to download:
#' see https://paleobiodb.org/data1.2/occs/list_doc.html
#' @param ex_taxon A character vector of taxon names to exclude
#' from the PBDB search
#' @param area If not NULL, then a numeric vector of length four
#' specifying, in order, the min lng, max lng, min lat and max
#' lat of the area from which occurrences will be returned, in
#' decimal degrees (equator = 0 lat, prime meridian = 0 lng).
#' Alternatively, a character vector of regions from which
#' occurrences will be returned: any valid country name or ISO2
#' code. Continent names and codes are also supported as follows:
#' ATA Antarctica, AFR Africa, ASI Asia, AUS Australia, EUR
#' Europe, IOC Indian Ocean, NOA North America, OCE Oceania,SOA
#' South America
#' @param ex_area If not NULL, then a character vector of
#' valid country names or ISO2 codes, as in @param area (), from
#' which occurrences will be excluded from a PBDB search
#' @param invert_area If TRUE, then regions specified in area
#' will be excluded from a PBDB search, except for the regions
#' specified in ex_area
#' @param litho If not NULL, a character vector of PBDB vocabulary
#' corresponding to which lithologies a PBDB search should return
#' @param invert_litho If TRUE, a character vector of PBDB
#' vocabulary corresponding to which lithologies a PBDB search
#' should exclude
#' @param env If not NULL, a character vector of PBDB vocabulary
#' corresponding to which environments a PBDB search should return
#' @param ex_env If not NULL, a character vector of PBDB vocabulary
#' corresponding to which environments a PBDB search should exclude
#' @param invert_env If TRUE, then environments specified in env
#' will be excluded from a PBDB search, except for the environments
#' specified in ex_env
#' @param pres A character vector of length one specifying the
#' preservation mode of the occurrences to return: one of regular,
#' form, ichno, or 'form,ichno'
#' @param idqual A character vector of length one specifying the
#' taxonomic certainty of the occurrences to return: one of certain,
#' genus_certain, uncertain, new"
#' @param return_url If TRUE, the function will return a correctly
#' formatted url suitable for use with curl or similar API functions,
#' comprising the search parameters set by the user
#' @param return_data If TRUE (default), the downloaded csv will
#' automatically be read into R (this must be assigned to an object)
#' @param save_as If not NULL, the file name to which the downloaded
#' data will be saved on the disk as a .csv
#' @param tscale A character vector of length one determining what
#' chronostratigraphic timescale will be applied to the data. "ICS2013"
#' will retain the PBDB ICS 2013 standard. "GTS2020" will update all
#' early and late interval ages to the GTS2020 standard, using a
#' lookup table supplied with the function. Alternatively, the
#' pathway to a custom .csv file  with columns Interval, FAD and
#' LAD where Interval are the names of the early and late intervals
#' in the PBDB, and FAD and LAD are the numeric lower and upper
#' boundaries of those intervals
#' @param wait The maximum wait time for the download in milliseconds,
#' as used by curl. This is set to no wait time by default
#' @return either a PBDB API compatible URL or a PBDB dataset
#' @importFrom stats na.omit
#' @import curl
#' @importFrom utils data
#' @importFrom data.table fread fwrite
#' @examples
#' # download Triassic dinosaurs (wait time set to meet CRAN example requirement)
#' tdinos <- fossilbrush:::get_pbdb(taxon = "Dinosauria", interval = "Triassic", wait = 499)

get_pbdb <- function(taxon = NULL, interval = NULL, mode = "occurrence", res = "all", fields = c("ident", "coords", "class"),
                     ex_taxon = NULL, area = NULL, ex_area = NULL, invert_area = FALSE, litho = NULL, invert_litho = FALSE,
                     env = NULL, ex_env = NULL, invert_env = NULL, pres = NULL, idqual = NULL,
                     return_url = FALSE, return_data = TRUE, save_as = NULL, tscale = "ICS2013", wait = Inf) {


  #taxon = "Dinosauria"
  #interval = "Triassic"
  #mode = "occurrence"
  #res = "all"
  #fields = c("ident", "coords", "class")
  #ex_taxon = NULL
  #area = NULL
  #ex_area = NULL
  #invert_area = FALSE
  #litho = NULL
  #invert_litho = FALSE
  #env = NULL
  #ex_env = NULL
  #invert_env = NULL
  #pres = NULL
  #idqual = NULL
  #return_url = FALSE
  #return_data = TRUE
  #save_as = NULL
  #tscale = "ICS2013"

  # check fundamental arguments (url vs download)
  if(is.null(taxon) & is.null(interval)) {
    stop("Taxon or interval must be specified at the very least")
  }
  if(isTRUE(return_url) & isTRUE(return_data)) {
    stop("return_url and auto_read cannot both be TRUE")
  }
  if(isFALSE(return_url) & isFALSE(return_data) & is.null(save_as)) {
    stop("If both return_url and return_data are FALSE, then the file location to which the data will be saved must be supplied")
  }

  if(!exists("taxon")) {
    taxon <- NULL
  }

  #######################
  # set up taxon argument
  if(!is.null(taxon)) {
    taxon <- unique(na.omit(taxon))
    if(!is.character(taxon) | length(taxon) == 0) {
      stop("Taxon must be a character vector of one or more names")
    }
    if(!is.null(ex_taxon)) {
      ex_taxon <- unique(na.omit(ex_taxon))
      if(!is.character(ex_taxon) | length(ex_taxon) == 0) {
        stop("Taxa to exclude must be provided as a character vector of one or more names")
      }
      if(any(ex_taxon %in% taxon)) {
        stop("One or more taxa are present in both the search and the exclusion")
      }
      ex_taxon <- paste0("^", ex_taxon, collapse = " ")
      taxon <- paste0(taxon, collapse = ",")
      taxon <- paste0(taxon, ex_taxon, collapse = " ")
    } else {
      taxon <- paste0(taxon, collapse = ",")
    }
  }
  if(!is.null(pres)) {
    pres <- na.omit(unique(pres))
    if(!is.character(pres)) {
      stop("Pres should be one of the following: regular, form, ichno, 'form,ichno'")
    }
    if(length(pres) != 0) {
      warning("More than one element was provided in pres - only the first will be used")
      pres <- pres[1]
    }
    if(!is.null(taxon)) {
      pres <- paste0("&pres=", pres)
    } else {
      pres <- paste0("pres=", pres)
    }
    taxon <- paste0(taxon, pres)
  }
  if(!is.null(idqual)) {
    idqual <- na.omit(unique(idqual))
    if(!is.character(idqual)) {
      stop("Idqual should be one of the following: certain, genus_certain, uncertain, new")
    }
    if(length(idqual) != 0) {
      warning("More than one element was provided in pres - only the first will be used")
      idqual <- idqual[1]
    }
    if(!is.null(taxon)) {
      idqual <- paste0("&idqual=", idqual)
    } else {
      idqual <- paste0("idqual=", idqual)
    }
    taxon <- paste0(taxon, idqual)
  }

  #######################
  # set up time argument if provided
  if(!is.null(interval)) {
    if((!class(interval) %in% c("character", "numeric")) | (is.character(interval) & length(grep("[0-9]", interval)) != 0)) {
      stop("Interval must be either a single interval name, a vector of two interval names, or numeric vector of length 2")
    }
    interval <- unique(na.omit(interval))
    if(is.character(interval)) {
      if(length(interval) > 2) {
        warning("Three or more interval names were provided - only the first and last will be used")
      }
      interval <- interval[c(1, length(interval))]
      GTS2020 <- get("GTS2020")
      interval <- c(GTS2020$FAD[match(interval[1], GTS2020$Interval)], GTS2020$LAD[match(interval[2], GTS2020$Interval)])
    }
    if(is.numeric(interval)) {
      if(length(interval) > 2) {
        warning("Three or more numerics were provided - only the largest and smallest values will be used")
      }
      interval <- c(max(na.omit(interval)), min(na.omit(interval)))
      if(abs(diff(interval)) > 100) {
        interval <- seq(interval[1], interval[2], by = -100)
        if(max(interval > 601)) {
          interval <- c(interval[1], interval[which(interval < 601)])
        }
      }
    }
  }

  #######################
  # set up mode
  mode <- na.omit(unique(mode))
  if(length(mode) != 1) {
    warning("Mode contains multiple elements. Only the first will be used")
    mode <- mode[1]
  }
  mode <- switch(mode,
                 "occurrence" = "https://paleobiodb.org/data1.2/occs/list.csv?",
                 "collection" = "https://paleobiodb.org/data1.2/colls/list.csv?",
                 "taxa" = "https://paleobiodb.org/data1.2/occs/taxa.csv?",
                 "specimen" = "https://paleobiodb.org/data1.2/specs/list.csv?",
                 "measurement" = "https://paleobiodb.org/data1.2/specs/measurements.csv?",
                 "strata" = "https://paleobiodb.org/data1.2/occs/strata.csv?",
                 "diversity" = "https://paleobiodb.org/data1.2/diversity.csv?",
                 "opinion" = "https://paleobiodb.org/data1.2/occs/opinions.csv?",
                 "reference" = "https://paleobiodb.org/data1.2/occs/refs.csv?")
  if(is.null(mode)) {
    stop("Mode must be one of the following: occurrence, collection, taxa, specimen,
         measurement, strata, diversity, opinion or reference")
  }

  #######################
  # check res argument
  res <- na.omit(unique(res))
  if(length(res) != 1) {
    warning("Only the first element of res will be used")
    res <- res[1]
  }
  if(!res %in% c("all", "family", "genus", "species", "lump_genus", "lump_gensub")) {
    stop("res must be one of all, family, genus, species, lump_genus or lump_subgen")
  }
  if(res == "all") {
    res <- NULL
  } else {
    res <- paste0("taxon_reso=", res, "&")
  }

  #######################
  # set up fields argument if provided
  if(!is.null(fields)) {
    if(!is.character(fields)) {
      stop("Fields must be a character vector with one or more elements corresponding to PBDB vocabulary")
    }
    fields <- na.omit(unique(fields))
    pbdb_fields <- get("pbdb_fields")
    inv <- fields[!fields %in% as.vector(unlist(pbdb_fields))]
    if(length(inv) != 0) {
      fields <- fields[!fields %in% inv]
      warning(paste0("The following elements in fields (", paste0(inv, collapse = ", "), ")
                     are not part of the PBDB vocabulary and have been dropped"))
    }
    bulk2 <- pbdb_fields$bulk[pbdb_fields$bulk %in% fields]
    if(length(bulk2) != 0) {
      for(i in 1:length(bulk2)) {
        rem <- as.vector(unlist(pbdb_fields[which(names(pbdb_fields) == bulk2[i])]))
        fields <- fields[!fields %in% rem]
      }
    }
    fields <- paste0(fields, collapse = ",")
    fields <- paste0("&show=", fields)
  }

  #######################
  # set up area argument if provided
  if(!is.null(area)) {
    area <- na.omit(area)
    if(!class(area) %in% c("character", "numeric") | length(grep("[0-9]", interval)) != 0) {
      stop("Area must be a numeric of length four with decimal degree values in order of min longitude, max longitude,
           min latitude, max latitude (prime meridian = 0 longitude, equator = 0 latitude), or a character vector of one
           or more valid country names, their ISO2 abbreviations, and/or one or more of the following continental regions:
           ATA Antarctica, AFR Africa, ASI Asia, AUS Australia, EUR Europe, IOC Indian Ocean, NOA North America, OCE Oceania,
           SOA South America")
    }
    if(is.numeric(area)) {
      if(length(area) != 4) {
        stop("If supplied as coordinates, Area must be a numeric of length four with decimal
        degree values in order of min longitude, max longitude, min latitude, max latitude (prime meridian = 0 longitude)")
      }
      if(any(area[1:2]) < -180 | any(area[1:2]) > 180 | any(area[3:4]) < -90 | any(area[3:4]) > 90) {
        stop("If supplied as coordinates, Area must be a numeric of length four with decimal degree values in order of
             min longitude, max longitude, min latitude, max latitude (prime meridian = 0 longitude, equator = 0 latitude)")
      }
      if(area[1] > area[2]) {
        stop("The minimum longitude cannot be greater than the maximum longitude")
      }
      if(area[3] > area[3]) {
        stop("The minimum latitude cannot be greater than the maximum longitude")
      }
      if(!is.null(ex_area)) {
        warning("As area has been supplied as a vector of decimal degrees, ex_area will be ignored")
      }
      area <- paste0("&lngmin=", area[1], "&lngmax=", area[2], "&latmin=", area[3], "&latmax=", area[4])
    }
    if(is.character(area)) {
      area <- unique(area)
      geog_lookup <- get("geog_lookup")
      test <- match(area, geog_lookup$Name)
      if(any(is.na(test))) {
        if(length(test) == length(area)) {
          stop("None of the area names are valid")
        }
        inv <- area[which(is.na(test))]
        warning(paste0("The following elements in area (", paste0(inv, collapse = ", "), ")
                     are not valid area names and have been dropped"))
        area <- area[-which(is.na(test))]
      }
      if(is.null(ex_area)) {
        area <- paste0(area, collapse = ",")
      }
    }
  }

  if(!is.null(ex_area)) {
    if(is.numeric(area)) {
      warning("As area has been supplied as a vector of decimal degrees, ex_area will be ignored")

    } else {
      ex_area <- na.omit(ex_area)
      if(!is.null(area)) {
        if(any(ex_area %in% area)) {
          stop("One or more areas are present in both the search and the exclusion")
        }
        geog_lookup <- get("geog_lookup")
        geog_lookup <- geog_lookup[10:nrow(geog_lookup),]
        test2 <- match(ex_area, geog_lookup$Name)
        if(any(is.na(test2))) {
          if(length(test2) == length(ex_area)) {
            stop("None of the area names are valid")
          }
          inv <- ex_area[which(is.na(test2))]
          warning(paste0("The following elements in area (", paste0(inv, collapse = ", "), ")
                     are not valid area names and have been dropped"))
          ex_area <- ex_area[-which(is.na(test2))]
        }
        ex_area <- paste0("^", ex_area, collapse = ",")
        if(invert_area) {
          area <- paste0("!", area, collapse = ",")
        }
        area <- paste0(area, ex_area, collapse = ",")
      } else {
        area <- paste0("!", paste0(ex_area, collapse = ","))
      }
    }
  }
  if(!is.null(area)) {
    area <- paste0("&cc=", area)
  }

  #######################
  # set up lithology argument if provided
  if(!is.null(litho)) {
    litho <- na.omit(unique(litho))
    if(!is.character(litho)) {
      stop("Lithology must be a character vector of one or more of the following elements:
           siliclastic, mixed, carbonate, evaporite, organic, chemical, volcanic, metasedimentary,
           metamorphic, other or unknown")
    }
    inv <- litho[!litho %in% c("siliclastic", "mixed", "carbonate", "evaporite", "organic", "chemical",
                               "volcanic", "metasedimentary", "metamorphic", "other", "unknown")]
    if(length(inv) != 0) {
      lith <- lith[!lith %in% inv]
      warning(paste0("The following elements in fields (", paste0(inv, collapse = ", "), ")
                     are not part of the PBDB lithology vocabulary and have been dropped"))
    }
    if(invert_litho) {
      litho <- paste0("!", paste0(litho, collapse = ","))
    }
    litho <- paste0("&lithology=", litho)
  }

  #######################
  # set up environment argument if provided
  if(!is.null(env)) {
    env <- na.omit(unique(env))
    if(!is.character(env)) {
      stop("Environment must be a character vector of one or more of the following elements (ignoring bracketed phrases):
           terrestrial, any marine, carbonate, siliciclastic, unknown, lacustrine, fluvial, karst,
           terrother (other terrestrial), marginal (marginal marine), reef, stshallow (shallow subtidal),
           stdeep (deep subtidal), offshore, slope (slope or basin), or marindet (indeterminate marine)")
    }
    inv <- env[!env %in% c("siliclastic", "mixed", "carbonate", "evaporite", "organic", "chemical",
                           "volcanic", "metasedimentary", "metamorphic", "other", "unknown",
                           "lacustrine", "fluvial", "karst", "terrother", "marginal", "reef", "stshallow",
                           "stdeep", "offshore", "slope", "marindet")]
    if(length(inv) != 0) {
      env <- env[!env %in% inv]
      warning(paste0("The following elements in environment (", paste0(inv, collapse = ", "), ")
                     are not part of the PBDB environment vocabulary and have been dropped"))
    }
    if(is.null(ex_env)) {
      env <- paste0(env, collapse = ",")
    }
  }

  if(!is.null(ex_env)) {
    ex_env <- na.omit(ex_env)
    if(!is.null(env)) {
      inv <- ex_env[!ex_env %in% c("lacustrine", "fluvial", "karst", "terrother", "marginal", "reef", "stshallow",
                                   "stdeep", "offshore", "slope", "marindet")]
      if(length(inv) != 0) {
        ex_env <- ex_env[!ex_env %in% inv]
        warning(paste0("The following elements in exclude environment (", paste0(inv, collapse = ", "), ")
                     are not part of the PBDB environment vocabulary and have been dropped"))
      }
      if(any(ex_env %in% env)) {
        stop("One or more environments are present in both the search and the exclusion")
      }
      ex_env <- paste0("^", ex_env, collapse = ",")
      if(invert_env) {
        env <- paste0("!", env, collapse = ",")
      }
      env <- paste0(env, ex_env, collapse = ",")
    } else {
      env <- paste0("!", paste0(ex_env, collapse = ","))
    }
  }
  if(!is.null(env)) {
    env <- paste0("&envtype=", env)
  }

  #######################
  # build url
  if(!is.null(interval)) {

    if(is.character(interval)) {
      interval <- paste0("interval=", interval)
    }
    if(is.numeric(interval)) {
      intnum <- paste0("max_ma=", interval[1], "&min_ma=", interval[length(interval)])
      if(length(interval) == 2) {
        interval <- intnum
      } else {
        interval1 <- interval
        interval <- vector()
        for(i in 1:(length(interval1) - 1)) {
          interval[i] <- paste0("max_ma=", interval1[i], "&min_ma=", interval1[i + 1])
        }
      }
    }
    if(is.null(taxon)) {
      if(is.character(interval)) {
        purl <- paste0(mode, res, interval, fields, area, litho)
        testurl <- purl
      }
      if(is.numeric(interval)) {
        purl <- paste0(mode, res, intnum, fields, area, litho)
        testurl <- vector()
        for(i in 1:length(interval)) {
          testurl[i] <- paste0(mode, res, interval[i], fields, area, litho)
        }
      }
    } else {
      if(is.character(interval)) {
        purl <- paste0(mode, paste0("base_name=", taxon, "&"), res, interval, fields, area, litho)
        testurl <- purl
      }
      if(is.numeric(interval)) {
        purl <- paste0(mode, paste0("base_name=", taxon, "&"), res, intnum, fields, area, litho)
        testurl <- vector()
        for(i in 1:length(interval)) {
          testurl[i] <- paste0(mode, paste0("base_name=", taxon, "&"), res, interval[i], fields, area, litho)
        }
      }
    }
  } else {
    purl <- paste0(mode, taxon, "&", res, fields, area, litho)
    testurl <- purl
  }

  #######################
  # do output
  if(return_url) {
    return(purl)

  } else {
    # clean tempdir
    del <- list.files(tempdir(), pattern = "_pbdb", full.names = TRUE)
    unlink(del)
    # download data (short loop used to allow messaging)
    for(i in 1:length(testurl)) {
      cat(paste0("Chunk ", i, "/", length(testurl), " downloading"), "\n")
      tfile <- tempfile(pattern = paste0(i, "_pbdb"), fileext = ".csv")
      cob <- new_handle()
      handle_setopt(cob, timeout = wait)
      foo <- try(curl_download(url = testurl[i], destfile = tfile, quiet = TRUE, handle = cob), silent = TRUE)
      if(inherits(foo, "try-error")) {break()}
    }
    if(inherits(foo, "try-error")) {
      message("Download wait time exceeded")
      return(NULL)
    } else {
      cat("Reading data", "\n")
      dload <- as.list(list.files(tempdir(), pattern = "_pbdb", full.names = TRUE))
      dload <- lapply(dload, data.table::fread, encoding = "UTF-8")
      dload <- do.call(rbind, dload)
      dload <- as.data.frame(dload[!duplicated(dload$occurrence_no), ])
      dload[dload == ""] <- NA

      # clean tempdir
      del <- list.files(tempdir(), pattern = "_pbdb", full.names = TRUE)
      unlink(del)
      # perform redating if specified
      do_redate = FALSE
      if(tscale == "ICS2013") {
        do_redate <- FALSE
      }
      if(tscale == "GTS2020") {
        tscale <- get("GTS2020")
        tscale <- tscale[,c("Interval", "FAD", "LAD")]
        do_redate <- TRUE
      }
      if(do_redate) {

        # this is a deconstruction of the redate_intervals function. The latter function exists
        # separately so it can be used with any dataset. Its functionality is deconstructed here so that
        # get_pbdb can stand alone as a single function
        xfad <- "early_interval"
        xlad <- "late_interval"
        xerl <- "max_ma"
        xlte <- "min_ma"
        cinterval <- "Interval"
        cfad <- "FAD"
        clad <- "LAD"
        # if the lad is missing, overwrite with the FAD
        dload[is.na(dload[,xlad]), xlad] <- dload[is.na(dload[,xlad]), xfad]
        # get new FADs and LADs
        new_fad <- tscale[match(dload[,xfad], tscale[,cinterval]), cfad]
        new_lad <- tscale[match(dload[,xlad], tscale[,cinterval]), clad]
        # in case of new intervals (which give NA), retain original date
        new_fad[is.na(new_fad)] <- dload[is.na(new_fad),xerl]
        new_lad[is.na(new_lad)] <- dload[is.na(new_lad),xlte]
        # overwrite
        dload[,xerl] <- new_fad
        dload[,xlte] <- new_lad
      }

      # build kingdom field if classification was requested
      if(grepl("class", fields)) {
        pbdb_kingdoms <- get("pbdb_kingdoms")
        animals <- pbdb_kingdoms$animals
        plants <- pbdb_kingdoms$plants
        protists <- pbdb_kingdoms$protists
        dload$kingdom <- NA
        dload$kingdom[dload$phylum %in% animals] <- "Animalia"
        dload$kingdom[dload$phylum %in% plants] <- "Plantae"
        dload$kingdom[dload$phylum %in% protists] <- "Protista"
        dload$kingdom[dload$phylum == "Pezizomycotina"] <- "Fungi"
        dload$kingdom[dload$phylum == "Cyanobacteria"] <- "Bacteria"
        if(which(colnames(dload) == "phylum") == 1) {
          dload <- cbind.data.frame(kingdom = dload$kingdom, dload[,-which(colnames(dload) == "kingdom"), drop = FALSE])
        } else {
          dload <- cbind.data.frame(dload[,1:(which(colnames(dload) == "phylum") - 1)], dload$kingdom, dload[,(which(colnames(dload) == "phylum")):(ncol(dload) - 1)])
        }
        colnames(dload)[which(colnames(dload) == "dload$kingdom")] <- "kingdom"
      }
      if(!is.null(save_as)) {
        save_as <- paste0(save_as, ".csv")
        cat("Writing data", "\n")
        data.table::fwrite(dload, save_as, bom = TRUE)
      }
      return(dload)
    }
  }
}
