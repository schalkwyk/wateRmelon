#### readEPIC

 # Author: Tyler Gorrie-Stone
 # LS add in support for epicv2 01-09-2023
 # Last Modified: 19-01-2016
 # Creation Date: 20-11-2015

 # Based heavily on 'methylumIDAT' and modified to read Epic arrays.
 # Functionality for 450k and 27k arrays remains (relatively) unchanged.
 # Although certain stops have been removed/moved to a different
 # place due to the implausibility of determining sample based on name
 # As differing between the three chips requires reading in .idats.

###
 # Usage:
 # readEPIC    : Instead of inputting the barcodes, it is possible to input filepath
 #               of a folder (containing idats or folders containing idats).
 #               This constructs a vector of the barcodes which is
 #               passed to methylumIDAT.
 #
 # methylumIDATepic: Usage remains the same as methylumIDAT. Although might be slightly changed
 #               due to the introduction of recursive=T argument in the list.files()
###
 # Leaving 27k and 450k array functionality in function for the purpose of convenience
 # however may be removed in the future when such arrays are not used anymore.
###

# {{{ # Utility function for dealing with single samples (still not 100% perfect...)
columnMatrix <- function(x, row.names = NULL) {
    if (is.null(dim(x)[2]))
        dim(x) = c(length(x), 1)
    if (!is.null(row.names))
        rownames(x) = row.names
    return(x)
}  
# }}}

# {{{ getMethylationBeadMappers2 
## require()s the appropriate package for annotating a chip & sets up mappings
## Potentially may temporarily require more indepth ordering file for 450k
## arrays until the 450k is completely replaced by epic chips.
getMethylationBeadMappers2 <- function(chipType = c("450k", "27k", "Epic", "Epicv2"),
    genome = c("hg19", "hg18", "hg38")) {
    genome <- match.arg(genome)  ## default to FDb.InfiniumMethylation.hg19
    pkg <- paste0("FDb.InfiniumMethylation.", genome)
    require(pkg, character.only = TRUE)  ## and




    chipType <- sub("^IlluminaHumanMethylation", "", chipType)
    if (class(chipType) %in% c("NChannelSet", "MethyLumiSet", "MethyLumiM")) {
        chipType <- sub("^IlluminaHumanMethylation", "", annotation(chipType))
        chipType <- sub(".db$", "", annotation(chipType))
    }
    chipType <- match.arg(chipType)  # backwards compatibility purposes

    ## addressA=U, addressB=M
    getProbes <- switch(chipType, `27k` = function(color = NULL, ...) {
        data(hm27.ordering)
        what <- c("Probe_ID", "M", "U")
        r <- split(hm27.ordering[, what], hm27.ordering$col)
        if (is.null(color)) return(r) else return(r[[substr(color, 1, 1)]])
    }, `450k` = function(design = NULL, color = NULL, ...) {
        data(hm450.ordering)
        what <- c("Probe_ID", "M", "U")
        r <- split(hm450.ordering[, c(what, "col")], hm450.ordering$DESIGN)
        r$I <- split(r$I[, what], r$I$col)
        r$II <- r$II[, what]
        if (is.null(design)) {
            return(r)
        } else if (design == "I" && !is.null(color)) {
            return(r[[design]][[substr(color, 1, 1)]])
        } else {
            return(r[[design]])
        }
    }, Epic = { function(design = NULL, color = NULL, ...) {
        # data(epic.ordering)
        epic.ordering <<- generateManifest("EPIC")
        what <- c("Probe_ID", "M", "U")
        r <- split(epic.ordering[, c(what, "col")], epic.ordering$DESIGN)
        r$I <- split(r$I[, what], r$I$col)
        r$II <- r$II[, what]
        if (is.null(design)) {
            r <- return(r)
        } else if (design == "I" && !is.null(color)) {
            r <- return(r[[design]][[substr(color, 1, 1)]])
        } else {
            r <- return(r[[design]])
        }}
    }, Epicv2 = function(design = NULL, color = NULL, ...) {
        # data(epic.ordering)
        epicV2.ordering <<- generateManifest("EPICv2")
        what <- c("Probe_ID", "M", "U")
        r <- split(epicV2.ordering[, c(what, "col")], epicV2.ordering$DESIGN)
        r$I <- split(r$I[, what], r$I$col)
        r$II <- r$II[, what]
        if (is.null(design)) {
            r <- return(r)
        } else if (design == "I" && !is.null(color)) {
            r <- return(r[[design]][[substr(color, 1, 1)]])
        } else {
            r <- return(r[[design]])
        }
    })
#}}}

    getControls <- switch(chipType, `27k` = function() { #{{{
        data(hm27.controls)
        return(hm27.controls)
    }, `450k` = function() {
        data(hm450.controls)
        return(hm450.controls)
    }, Epic = function() {
        data(epic.controls)
        return(epic.controls)
    }, Epicv2 = function() {
        data(epicV2.controls)
        return(epicV2.controls)
    })

    getOrdering <- switch(chipType, `27k` = function() {
        ord <- c("Probe_ID", "DESIGN", "COLOR_CHANNEL")
        return(hm27.ordering[, ord])
    }, `450k` = function() {
        ord <- c("Probe_ID", "DESIGN", "COLOR_CHANNEL")
        return(hm450.ordering[, ord])
    }, Epic = function() {
        ord <- c("Probe_ID", "DESIGN", "COLOR_CHANNEL")
        return(epic.ordering[, ord])
    }, Epicv2 = function() {
        ord <- c("Probe_ID", "DESIGN", "COLOR_CHANNEL")
        return(epicV2.ordering[, ord])
    })

    mapper <- list(probes = getProbes, controls = getControls, ordering = getOrdering)
    return(mapper)
}  # }}}

## IDATtoMatrix2: process a single IDAT (just the mean intensities) {{{
IDATtoMatrix2 <- function(x, fileExts = list(Cy3 = "Grn", Cy5 = "Red"), idatPath = ".") {
    chs = names(fileExts)
    names(chs) = fileExts
    processed = lapply(fileExts, function(ch) {
        ext = paste(ch, "idat", sep = ".")
        message(sprintf("Reading in: %s", paste(x, ext, sep = "_")))
        dat = readIDAT(file.path(idatPath, paste(x, ext, sep = "_")))
        Quants = data.matrix(dat$Quants)
        colnames(Quants) = paste(chs[ch], colnames(Quants), sep = ".")
        return(list(Quants = Quants, RunInfo = dat$RunInfo, ChipType = dat$ChipType))
    })
    probe.data = do.call(cbind, lapply(processed, function(x) x[["Quants"]]))
    probe.data <- probe.data[, c("Cy3.Mean", "Cy5.Mean", "Cy3.NBeads", "Cy5.NBeads")]
    # Nbeads for beadcount switch TGS, at this point add something that takes
    # the beadcount with the lower number of beads.
    attr(probe.data, "RunInfo") = processed[[1]][["RunInfo"]]
    attr(probe.data, "ChipType") = processed[[1]][["ChipType"]]
    return(probe.data)
}  # }}}

# {{{ IDATsToMatrices2
IDATsToMatrices2 <- function(barcodes, fileExts = list(Cy3 = "Grn", Cy5 = "Red"),
    parallel = F, idatPath = ".") {
    names(barcodes) = as.character(barcodes)
    if (parallel) {
        mats = .mclapply(barcodes, IDATtoMatrix2, fileExts = fileExts, idatPath = idatPath)
    } else {
        mats = lapply(barcodes, IDATtoMatrix2, fileExts = fileExts, idatPath = idatPath)
    }
    names(mats) = as.character(barcodes)
    return(mats)
}  # }}}

#{{{ extractAssayDataFromList2
## might consider doing the following in C++ via Rcpp to avoid copying data 
extractAssayDataFromList2 <- function(assay, mats, fnames) {
    d <- do.call("cbind", lapply(mats, function(x) x[fnames, assay]))  # Select Common Probes
    if (!is.matrix(d))
        d <- data.matrix(d)
    colnames(d) <- names(mats)
    rownames(d) <- fnames
    return(d)
}  # }}}

#{{{  DataToNChannelSet2
## a faster rewrite of DFsToNChannelSet() so that I can decommission it... 
## Stop-check to ensure different platforms are read simulatenously.  note EPIC
## and EPICv2 have the same chip type and there are several different EPIC
## manifests
DataToNChannelSet2 <- function(mats, chans = c(Cy3 = "GRN", Cy5 = "RED"), parallel = F,
    protocol.data = F, IDAT = TRUE, force = F) {
    epic = hm27 = hm450 = 0
    qw <- unlist(lapply(mats, function(x) attr(x, "ChipType")))
    epic = sum(grepl("BeadChip 8x5", qw))
    message(paste(epic, "HumanMethylationEpic / Epicv2 samples found"))
    hm450 = sum(grepl("BeadChip 12x8", qw))
    message(paste(hm450, "HumanMethylation450 samples found"))
    hm27 = sum(grepl("BeadChip 12x1", qw))
    message(paste(hm27, "HumanMethylation27 samples found"))
    if (hm27 > 0 && hm450 > 0 | hm27 > 0 && epic > 0 | hm450 > 0 && epic > 0) {
        stop("Cannot process multiple platforms simultaneously; please run separately.")
    }

    stopifnot(is(mats, "list"))
    assayNames = paste0(names(chans), c(".Mean"))
    assayNames = c(assayNames, paste0(names(chans), c(".NBeads")))
    names(assayNames) = assayNames
    assaylengths <- sapply(mats, nrow)
    names(assaylengths) <- names(mats)
    # Minfi Method for handling early version epic IDATS. from read.meth.R by
    # Kaspar Hansen.
    if (length(unique(assaylengths)) > 1) {
        commonAddresses <- as.character(Reduce("intersect", lapply(mats, function(x) rownames(x))))
        if (!force) {
            if (!all(assaylengths == length(commonAddresses))) {
                print(as.matrix(assaylengths))
                stop("Cannot combine IDATs of differing lengths, try force = T")
            }
        }
        fnames <- commonAddresses
    } else {
        fnames <- rownames(mats[[1]])
    }

    extract <- function(assay) extractAssayDataFromList2(assay, mats, fnames)
    assays = lapply(assayNames, extract)
    nb <- pmin(assays[["Cy3.NBeads"]], assays[["Cy5.NBeads"]])
    obj = new("NChannelSet", assayData = assayDataNew(R = assays[["Cy5.Mean"]], G = assays[["Cy3.Mean"]],
        N = nb))
    featureNames(obj) = fnames
    if (IDAT)
        {
            message("Attempting to extract protocolData() from list...")
            ChipType = attr(mats[[1]], "ChipType")
            RunInfo = lapply(mats, function(d) attr(d, "RunInfo"))
            if (protocol.data)
                {
                  scanDates = data.frame(
                     DecodeDate = rep(NA, length(mats)), 
                     ScanDate = rep(NA, length(mats))
                  )
                  rownames(scanDates) = names(mats)
                  for (i in seq_along(mats)) {
                    cat("decoding protocolData for", names(mats)[i], "...\n")
                    if (nrow(RunInfo[[i]]) >= 2) {
                        scanDates$DecodeDate[i] = RunInfo[[i]][1, 1]
                        scanDates$ScanDate[i] = RunInfo[[i]][2, 1]
                    }
                }
                protocoldata = new("AnnotatedDataFrame",
                    data = scanDates,
                    varMetadata = data.frame(
                        labelDescription = colnames(scanDates),
                        row.names = colnames(scanDates)
                    )
                )
                protocolData(obj) = protocoldata
                }

            message("Determining chip type from IDAT protocolData...")
            if (ChipType == "BeadChip 12x1") {
                annotation(obj) = "IlluminaHumanMethylation27k"
            } else if (ChipType == "BeadChip 12x8") {
                annotation(obj) = "IlluminaHumanMethylation450k"
            } else if (ChipType == "BeadChip 8x5") {
                annotation(obj) = "IlluminaHumanMethylationEpic"
            }

            if (ChipType == "BeadChip 8x5" && dim(obj)[1] > 1.1e+06) {
                annotation(obj) = "IlluminaHumanMethylationEpicv2"
            }
        }  
    return(obj)
}  # }}}

# {{{ getControlProbes2
getControlProbes2 <- function(NChannelSet) {
    fD <- getMethylationBeadMappers2(annotation(NChannelSet))$controls()
    ctls <- match(fD[["Address"]], featureNames(NChannelSet))

    ## FIXME: make this happen in the annotations, to avoid redundancy in
    ## names!
    rownames(fD) <- ctlnames <- make.names(fD[, "Name"], unique = T)
    fvD <- data.frame(labelDescription = c("Address of this control bead", "Purpose of this control bead",
        "Color channel for this bead", "Reporter group ID for this bead"))
    fDat <- new("AnnotatedDataFrame", data = fD, varMetadata = fvD)
    methylated <- assayDataElement(NChannelSet, "G")[ctls, , drop = FALSE]  # Cy3
    unmethylated <- assayDataElement(NChannelSet, "R")[ctls, , drop = FALSE]  # Cy5
    beadcount <- assayDataElement(NChannelSet, "N")[ctls, , drop = FALSE]  ## Testing TGS
    rownames(beadcount) <- rownames(methylated) <- rownames(unmethylated) <- ctlnames

    aDat <- assayDataNew(methylated = methylated, unmethylated = unmethylated, beadcount = beadcount)  #TGS
    new("MethyLumiQC", assayData = aDat, featureData = fDat, annotation = annotation(NChannelSet))
}  # }}}

## {{{ designItoMandU2 
## 27k design, both probes same channel; ~100,000 of the 450k probes as well
designItoMandU2 <- function(NChannelSet, parallel = F, n = F, n.sd = F, oob = T) {
    mapper <- getMethylationBeadMappers2(annotation(NChannelSet))
    probes <- mapper$probes(design = "I")  # as list(G=..., R=...)
    channels <- c("G", "R")
    names(channels) <- channels

    getIntCh <- function(NChannelSet, ch, al) {
        # {{{
        newprobes <- lapply(probes, function(y) {
            lapply(y, function(x) {
                x[probes[[ch]][[al]] %in% rownames(assayDataElement(NChannelSet,
                  ch))]
            })
        })
        a = assayDataElement(NChannelSet, ch)[as.character(newprobes[[ch]][[al]]),
            , drop = FALSE]
        rownames(a) = as.character(newprobes[[ch]][["Probe_ID"]])
        return(a)
    }  # }}}

    getOOBCh <- function(NChannelSet, ch, al) {
        # {{{
        ch.oob <- ifelse(ch == "R", "G", "R")
        newprobes <- lapply(probes, function(y) {
            lapply(y, function(x) {
                x[probes[[ch]][[al]] %in% rownames(assayDataElement(NChannelSet,
                  ch))]
            })
        })
        a = assayDataElement(NChannelSet, ch.oob)[as.character(newprobes[[ch]][[al]]),
            , drop = FALSE]
        rownames(a) = as.character(newprobes[[ch]][["Probe_ID"]])
        return(a)
    }  # }}}

    getNbeadCh <- function(NChannelSet, ch, al) {
        # {{{ #tgs
        newprobes <- lapply(probes, function(y) {
            lapply(y, function(x) {
                x[probes[[ch]][[al]] %in% rownames(assayDataElement(NChannelSet,
                  ch))]
            })
        })
        n = assayDataElement(NChannelSet, "N")[as.character(newprobes[[ch]][[al]]),
            , drop = FALSE]
        rownames(n) = as.character(newprobes[[ch]][["Probe_ID"]])
        return(n)
    }  # }}}

    getAllele <- function(NChannelSet, al, parallel = F, n = n, n.sd = T, oob = T) {
        # {{{
        fluor = lapply(channels, function(ch) getIntCh(NChannelSet, ch, al))
        fluor.oob = lapply(channels, function(ch) getOOBCh(NChannelSet, ch, al))
        bc = lapply(channels, function(ch) getNbeadCh(NChannelSet, ch, al))  #tgs
        res = list()
        res[["I"]] = fluor
        if (n)
            res[["BC"]] = bc
        if (oob)
            res[["OOB"]] = fluor.oob
        lapply(res, function(r) {
            names(r) = channels
            return(r)
        })
    }  # }}}

    signal <- lapply(c(M = "M", U = "U"), function(al) {
        getAllele(NChannelSet, al, parallel = F, n = n, n.sd = n.sd, oob = oob)
    })

    retval = list(methylated = rbind(signal$M$I$R, signal$M$I$G), unmethylated = rbind(signal$U$I$R,
        signal$U$I$G))

    if (n) {
        # tgs
        retval[["m.beadcount"]] = rbind(signal$M$BC$R, signal$M$BC$G)
        retval[["u.beadcount"]] = rbind(signal$U$BC$R, signal$U$BC$G)
        # NBeads takes the lowest bead count for each type I probe.
        retval[["NBeads"]] = pmin(retval[["m.beadcount"]], retval[["u.beadcount"]])
    }
    if (oob) {
        retval[["methylated.OOB"]] = rbind(signal$M$OOB$R, signal$M$OOB$G)
        retval[["unmethylated.OOB"]] = rbind(signal$U$OOB$R, signal$U$OOB$G)
    }
    return(retval)
}  # }}}

## {{{ designIItoMandU2 
## 450k/GoldenGate design (green=methylated, red=unmethylated, single address)
# loads the annotation DB so we can run SQL queries
designIItoMandU2 <- function(NChannelSet, parallel = F, n = F, n.sd = F, oob = T) {
    mapper <- getMethylationBeadMappers2(annotation(NChannelSet))
    probes2 <- mapper$probes(design = "II")
    probes2$M <- probes2$U  ## horrid kludge

    getIntCh <- function(NChannelSet, ch = NULL, al) {
        # {{{
        ch <- ifelse(al == "M", "G", "R")
        newprobes <- lapply(probes2, function(x) {
            x[probes2[[al]] %in% rownames(assayDataElement(NChannelSet, ch))]
        })
        a <- assayDataElement(NChannelSet, ch)[as.character(newprobes[[al]]), , drop = F]
        rownames(a) <- as.character(newprobes[["Probe_ID"]])
        return(a)
    }  # }}}

    getNbeadCh <- function(NChannelSet, ch = NULL, al) {
        # {{{ tgs
        ch <- ifelse(al == "M", "G", "R")
        newprobes <- lapply(probes2, function(x) {
            x[probes2[[al]] %in% rownames(assayDataElement(NChannelSet, ch))]
        })
        n <- assayDataElement(NChannelSet, "N")[as.character(newprobes[[al]]), ,
            drop = F]
        rownames(n) <- as.character(newprobes[["Probe_ID"]])
        return(n)
    }  # }}}

    getAllele <- function(NChannelSet, al, n = F, n.sd = F, oob = F) {
        # {{{
        ch <- ifelse(al == "M", "G", "R")
        res <- list()
        res[["I"]] <- getIntCh(NChannelSet, ch, al)
        if (n)
            res[["BC"]] <- getNbeadCh(NChannelSet, ch, al)
        if (oob) {
            res[["OOB"]] <- res[["I"]]
            is.na(res[["OOB"]]) <- TRUE
        }
        return(res)
    }  # }}}

    ## M == Grn/Cy3 and U == Red/Cy5, same address
    alleles = c(M = "M", U = "U")
    signal = lapply(alleles, function(a) getAllele(NChannelSet, a, n, n.sd, oob))

    retval = list(methylated = signal$M$I, unmethylated = signal$U$I)
    if (n) {
        # tgs
        retval[["m.beadcount"]] = signal$M$BC
        retval[["u.beadcount"]] = signal$U$BC
        retval[["NBeads"]] = signal$M$BC  #TGSbc
    }

    if (oob) {
        retval[["methylated.OOB"]] = signal$M$OOB
        retval[["unmethylated.OOB"]] = signal$U$OOB
    }
    return(retval)
}  # }}}

#{{{ mergeProbeDesigns2
## TGS: possible to do this better?
mergeProbeDesigns2 <- function(NChannelSet, parallel = F, n = F, n.sd = F, oob = T) {

    if (annotation(NChannelSet) == "IlluminaHumanMethylationEpicv2") {
        design1 = designItoMandU2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob)
        ## this is the source of the problem currently:
        design2 = designIItoMandU2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob)
        res <- list()
        for (i in names(design1)) {
            res[[i]] <- rbind(design1[[i]], design2[[i]])
            rownames(res[[i]]) <- c(rownames(design1[[i]]), rownames(design2[[i]]))
        }

    } else if (annotation(NChannelSet) == "IlluminaHumanMethylationEpic") {
        design1 = designItoMandU2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob)
        ## this is the source of the problem currently:
        design2 = designIItoMandU2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob)
        res <- list()
        for (i in names(design1)) {
            res[[i]] <- rbind(design1[[i]], design2[[i]])
            rownames(res[[i]]) <- c(rownames(design1[[i]]), rownames(design2[[i]]))
        }

    } else if (annotation(NChannelSet) == "IlluminaHumanMethylation450k") {
        design1 = designItoMandU2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob)
        ## this is the source of the problem currently:
        design2 = designIItoMandU2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob)
        res <- list()
        for (i in names(design1)) {
            res[[i]] <- rbind(design1[[i]], design2[[i]])
            rownames(res[[i]]) <- c(rownames(design1[[i]]), rownames(design2[[i]]))
        }
    } else if (annotation(NChannelSet) == "IlluminaHumanMethylation27k") {
        res <- designItoMandU2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob)
    } else {
        stop("don't know how to process chips of type", annotation(NChannelSet))
    }
    return(res)  # reorder on the way out...
}  # }}}

# {{{ NChannelSetToMethyLumiSet2
NChannelSetToMethyLumiSet2 <- function(NChannelSet, parallel = F, pval = 0.05, n = F,
    n.sd = F, oob = T) {
    history.submitted = as.character(Sys.time())
    results = mergeProbeDesigns2(NChannelSet, parallel = parallel, n = n, n.sd = n.sd,
        oob = oob)
    # The next part is somewhat messy, I will think of an better way to do this
    if (oob && n) {
        aDat <- with(results, assayDataNew(methylated = methylated, unmethylated = unmethylated,
            methylated.N = m.beadcount, unmethylated.N = u.beadcount, methylated.OOB = methylated.OOB,
            unmethylated.OOB = unmethylated.OOB, betas = methylated/(methylated +
                unmethylated), pvals = methylated/(methylated + unmethylated), NBeads = NBeads))
    } else if (oob) {
        aDat <- with(results, assayDataNew(methylated = methylated, unmethylated = unmethylated,
            methylated.OOB = methylated.OOB, unmethylated.OOB = unmethylated.OOB,
            betas = methylated/(methylated + unmethylated), pvals = methylated/(methylated +
                unmethylated)))
    } else if (n) {
        aDat <- with(results, assayDataNew(methylated = methylated, unmethylated = unmethylated,
            methylated.N = m.beadcount, unmethylated.N = u.beadcount, betas = methylated/(methylated +
                unmethylated), pvals = methylated/(methylated + unmethylated), NBeads = NBeads))

    } else {
        aDat <- with(results, assayDataNew(methylated = methylated, unmethylated = unmethylated,
            betas = methylated/(methylated + unmethylated), pvals = methylated/(methylated +
                unmethylated)))
    }
    rm(results)
    gc()

    ## now return the MethyLumiSet (which can be directly coerced to
    ## MethyLumiM)
    x.lumi = new("MethyLumiSet", assayData = aDat)
    x.lumi@QC <- getControlProbes2(NChannelSet)
    x.lumi@protocolData <- protocolData(NChannelSet)
    x.lumi@annotation <- annotation(NChannelSet)
    x.lumi@QC@annotation <- annotation(NChannelSet)
    pdat <- data.frame(barcode = sampleNames(NChannelSet))
    rownames(pdat) <- sampleNames(NChannelSet)
    pData(x.lumi) <- pdat
    varLabels(x.lumi) <- c("barcode")
    varMetadata(x.lumi)[, 1] <- c("Illumina BeadChip barcode")
    mapper <- getMethylationBeadMappers2(annotation(NChannelSet))
    fdat <- mapper$ordering()
    rownames(fdat) <- fdat$Probe_ID
    x.fnames <- rownames(betas(x.lumi))
    fdat <- fdat[x.fnames, ]

    ## Regression tests: fail noisily if there is an ordering issue
    stopifnot(identical(rownames(methylated(x.lumi)), rownames(unmethylated(x.lumi))))
    stopifnot(identical(rownames(betas(x.lumi)), rownames(unmethylated(x.lumi))))
    stopifnot(identical(rownames(fdat), rownames(betas(x.lumi))))
    # The only way to feasibly clean up this section would be to outsource the
    # possible metadatas to another script.
    fData(x.lumi) <- fdat
    possibleLabels <- c("Probe_ID", "DESIGN", "COLOR_CHANNEL", "PROBE_TYPE", "SNP10",
        "SYMBOL", "CHR36", "CPG36", "CPGS")
    fvarLabels(x.lumi) <- possibleLabels[1:ncol(fdat)]
    possibleMetadata <- c("Illumina probe ID from manifest", "Infinium design type (I or II)",
        "Color channel (for type I probes)", "Probe locus type (CpG, CpH, or SNP)",
        "SNP (dbSNP build 128) within 10bp of target?", "Gene symbol (if probe is annotated to a gene)",
        "Chromosome mapping for probe in hg18 assembly", "Coordinates of interrogated cytosine in hg18",
        "Number of CpG dinucleotides in probe sequence")
    fvarMetadata(x.lumi)[, 1] <- possibleMetadata[1:ncol(fdat)]
    pval.detect(x.lumi) <- pval  # default value
    history.finished <- as.character(Sys.time())
    history.command <- deparse(match.call())
    x.lumi@history <- rbind(x.lumi@history, data.frame(submitted = history.submitted,
        finished = history.finished, command = history.command))
    # if(normalize) x.lumi = normalizeMethyLumiSet(x.lumi)
    return(x.lumi)
}  # }}}

# {{{ methylumIDATepic
methylumIDATepic <- function(barcodes = NULL, pdat = NULL, parallel = F, n = F, n.sd = F,
    oob = T, idatPath = getwd(), force = F, ...) {
    if (is(barcodes, "data.frame"))
        pdat = barcodes
    if ((is.null(barcodes)) & (is.null(pdat) | (!("barcode" %in% names(pdat)))))
        {
            stop("\"barcodes\" or \"pdat\" (with pdat$barcode defined) must be supplied.")
        }
    if (!is.null(pdat) && "barcode" %in% tolower(names(pdat))) {
        names(pdat)[which(tolower(names(pdat)) == "barcode")] = "barcode"
        barcodes = pdat$barcode
        if (any(grepl("idat", ignore.case = TRUE, barcodes))) {
            message("Warning: filtering out raw filenames")
            barcodes = gsub("_(Red|Grn)", "", barcodes, ignore.case = TRUE)
            barcodes = gsub(".idat", "", barcodes, ignore.case = TRUE)
        }
        if (any(duplicated(barcodes)))
            {
                message("Warning: filtering out duplicates")
                pdat = pdat[-which(duplicated(barcodes)), ]
                barcodes = pdat$barcode
            } 
    } else {
        if (any(grepl("idat", ignore.case = TRUE, barcodes))) {
            message("Warning: filtering out raw filenames")
            barcodes = unique(gsub("_(Red|Grn)", "", barcodes, ignore.case = TRUE))
            barcodes = unique(gsub(".idat", "", barcodes, ignore.case = TRUE))
        }
        if (any(duplicated(barcodes))) {
            message("Warning: filtering out duplicate barcodes")
            barcodes = barcodes[which(!duplicated(barcodes))]
        }
    } 
    files.present = rep(TRUE, length(barcodes))  # {{{
    idats = sapply(barcodes, function(b) paste(b, c("_Red", "_Grn"), ".idat", sep = ""))
    for (i in colnames(idats)) for (j in idats[, i]) if (!j %in% list.files(idatPath,
        recursive = TRUE)) {
        message(paste("Error: file", j, "is missing for sample", i))
        files.present = FALSE
    }
    stopifnot(all(files.present))  # }}}

    mats <- IDATsToMatrices2(barcodes, parallel = parallel, idatPath = idatPath)
    dats <- DataToNChannelSet2(mats, IDAT = T, parallel = parallel, force = force)
    mlumi <- NChannelSetToMethyLumiSet2(dats, parallel = parallel, oob = oob, n = n)

    if (is.null(pdat)) {
        # {{{
        pdat = data.frame(barcode = as.character(barcodes))
        rownames(pdat) = pdat$barcode
        pData(mlumi) = pdat  # }}}
    } else {
        pData(mlumi) = pdat
    }
    if (!is.null(mlumi@QC))
        {
            sampleNames(mlumi@QC) = sampleNames(mlumi)
        } 

    # finally
    colnames(mlumi) <- as.character(colnames(mlumi))  # Fixing _that_ one bug.
    return(mlumi[sort(featureNames(mlumi)), ])
}  # }}}

#{{{ readEPIC
readEPIC <- function(idatPath, barcodes = NULL, pdat = NULL, parallel = F, n = T,
    oob = F, force = F, ...) {
    # path: file path to folder containing idat files, if working directory
    # is folder containing idats, use path <- './' the gsub will catch all
    # types of arrays as it is technically impossible to distinguish chiptype
    # through barcode.
    if (is.null(barcodes)) {
        barcodes <- idatPath
        methylumIDATepic(bfp(barcodes), pdat = pdat, parallel = parallel, n = n,
            n.sd = n.sd, oob = oob, idatPath = idatPath, force = force, ...)
    } else {
        methylumIDATepic(barcodes, pdat = pdat, parallel = parallel, n = n, n.sd = n.sd,
            oob = oob, idatPath = idatPath, force = force, ...)
    }

}  # }}}

bfp <- function(path) { #{{{
    # Barcodes from path, incase one wishes to use MethylumIDATepic manually.
    # bar <- dir(path, recursive=T)[grepl('R0[12345678]C0[12]_(Red|Grn).idat',
    # dir(path, recursive=T))]
    bar <- dir(path, recursive = TRUE)[grepl("_(Red|Grn).idat", dir(path, recursive = TRUE))]
    bar <- gsub("_(Red|Grn).idat", "", bar, ignore.case = TRUE)
    # Automatically sift out duplicated filenames!
    bar <- bar[!duplicated(basename(bar))]
    return(bar)
}#}}}

generateManifest <- function(anno = c("450k", "EPIC", "EPICv2")) { #{{{
    anno <- match.arg(anno)
    # Possible to add more manifests here!
    anno <- switch(anno, 
        `450k` = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        EPIC   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", 
        EPICv2 = "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
    )
    man <- getAnnotationObject(anno)
    x <- getAnnotation(man)[, c("Name", "AddressB", "AddressA", "Type", "Color")]
    # Name = Name AddressB = M AddressA = U Type = DESIGN Color = COLOR_CHANNEL
    # Generate manifest.
    snpI  <- getProbeInfo(man, type = "SnpI" )[, c(1, 2, 3, 4)]
    snpII <- getProbeInfo(man, type = "SnpII")[, c(1, 2)]
    snpI  <- cbind(snpI[, c("Name", "AddressB", "AddressA")], rep("I", nrow(snpI)),
        snpI[, "Color"])
    snpII <- cbind(snpII[, "Name"], rep("", nrow(snpII)), snpII[, "AddressA"], rep("II",
        nrow(snpII)), rep("", nrow(snpII)))
    colnames(snpI) <- colnames(snpII) <- colnames(x)
    x1 <- rbind(data.frame(x, stringsAsFactors = F), data.frame(snpI, stringsAsFactors = F),
        data.frame(snpII, stringsAsFactors = F))
    x1$col <- x1$Color
    x1$Color[x1$Color == ""   ] <- "Both"
    x1$col  [x1$Color == "Red"] <- "R"
    x1$col  [x1$Color == "Grn"] <- "G"
    is.na(x1$col) <- x1$Color == "Both"
    colnames(x1) <- c("Probe_ID", "M", "U", "DESIGN", "COLOR_CHANNEL", "col")
    x1$COLOR_CHANNEL <- factor(x1$COLOR_CHANNEL)
    x1$col <- factor(x1$col)
    return(x1)
} #}}}

# memory hole {{{
## for HM27k, all probes are of "design I": single-channel, two-address pairings
## for HM450k, probes are either "design I" or "design II" as noted in manifest!
## for epic, probes are either "design I" or "design II" as noted in manifest
## IDATs from GoldenGate methylation arrays are not supported at this time.
#cy3 <- function(object) { #
#  if(is.element('Color_Channel', fvarLabels(object)) &&
#     !is.element('COLOR_CHANNEL', fvarLabels(object))) {
#    fvarLabels(object)<-gsub('Color_Channel','COLOR_CHANNEL',fvarLabels(object))
#  }
#  if(!is.element('COLOR_CHANNEL', fvarLabels(object))) {
#    annotChip <- paste(annotation(object),'db',sep='.')
#    object <- addColorChannelInfo(object, annotChip)
#  }
#  return(which(fData(object)$COLOR_CHANNEL=='Grn'))
#} #


#cy5 <- function(object) { #
#  if(is.element('Color_Channel', fvarLabels(object)) &&
#     !is.element('COLOR_CHANNEL', fvarLabels(object))) {
#    fvarLabels(object)<-gsub('Color_Channel','COLOR_CHANNEL',fvarLabels(object))
#  }
#  if(!is.element('COLOR_CHANNEL', fvarLabels(object))) {
#    annotChip <- paste(annotation(object),'db',sep='.')
#    object <- addColorChannelInfo(object, annotChip)
#  }
#  return(which(fData(object)$COLOR_CHANNEL=='Red'))
#} # }}}


