function (Proj = NULL, groupBy = "Sample", chromSizes = getChromSizes(ArchRProj), 
    TSS = getTSS(ArchRProj), flank = 2000, norm = 100, smooth = 11, 
    pal = NULL, returnDF = FALSE, threads = getArchRThreads(), 
    logFile = createLogFile("plotTSSEnrichment")) 
{

    tstart <- Sys.time()
    chr <- paste0(seqnames(chromSizes))
    chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
    TSS <- sort(sortSeqlevels(TSS))
    splitTSS <- split(resize(TSS, 1, "start"), seqnames(TSS))[chr]
    window <- 2 * flank + 1
    groups <- Proj@meta.data[,groupBy]
    uniqGroups <- gtools::mixedsort(unique(groups))
    frag <- as.data.frame(fread(fragment,header=F))
    grFrag <- GRanges(seqnames = frag$V1,ranges = IRanges(start = frag$V2, end = frag$V3),strand = "*")
    dfTSS <- ArchR:::.safelapply(seq_along(uniqGroups), function(x) {
        cellx <- rownames(groups)[which(paste0(groups[, 1]) == 
            uniqGroups[x])]
        for (i in seq_along(chr)) {
            TSSi <- splitTSS[[chr[i]]]
            covi <- unlist(suppressMessages(getFragmentsFromProject(ArchRProj = proj, 
                subsetBy = chromSizes[paste0(seqnames(chromSizes)) %in% 
                  chr[i]], cellNames = cellx)), 
                use.names = FALSE) %>% sort %>% {
                coverage(IRanges(c(start(.), end(.)), width = 1))
            }

            if (i == 1) {
                sumTSS <- ArchR:::rleSumsStranded(list(chr1 = covi), 
                  list(chr1 = TSSi), window, as.integer)
            }
            else {
                sumTSS <- sumTSS + ArchR:::rleSumsStranded(list(chr1 = covi), 
                  list(chr1 = TSSi), window, as.integer)
            }
        }
        normBy <- mean(sumTSS[c(1:norm, (flank * 2 - norm + 1):(flank * 
            2 + 1))])
        df <- DataFrame(group = uniqGroups[x], x = seq_along(sumTSS) - 
            flank - 1, value = sumTSS, normValue = sumTSS/normBy, 
            smoothValue = ArchR:::.centerRollMean(sumTSS/normBy, 11))
        df
    }) %>% Reduce("rbind", .)

    if (returnDF) {
        return(dfTSS)
    }
    else {
        plotDF <- data.frame(x = dfTSS$x, v = dfTSS$smoothValue, 
            group = dfTSS$group)
        plotDF <- plotDF[sort(unique(c(1, seq(1, nrow(plotDF), 
            11), nrow(plotDF)))), , drop = FALSE]
        if (is.null(pal)) {
            pal <- paletteDiscrete(values = unique(plotDF$group))
        }
        p <- ggplot(plotDF, aes(x, v, color = group)) + geom_line(size = 1) + 
            theme_ArchR() + xlab("Distance From Center (bp)") + 
            ylab("Normalized Insertion Profile") + scale_color_manual(values = pal) + 
            scale_y_continuous(limits = c(0, max(plotDF$v) * 
                1.05), expand = c(0, 0)) + scale_x_continuous(limits = c(min(plotDF$x), 
            max(plotDF$x)), expand = c(0, 0))
        p
    }
}