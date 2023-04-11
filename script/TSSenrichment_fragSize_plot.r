parser = argparse::ArgumentParser(description=cat("Script to plot TSSEnrichment\n\n# Authors:   Maven\n# Contact information:   mawen3@genomics.cn\n# Date: 2022-04-7\n# R package version: cowplot_1.1.1  data.table_1.14.2  RColorBrewer_1.1-2  ggplot2_3.3.5  dplyr_1.0.7\n\n"))
parser$add_argument('-T', '--tssbed', help='input tss bed file')
parser$add_argument('-F', '--fragment', help='input fragment file')
parser$add_argument('-G', '--group',default='sample', help='sample name [default = \"%(default)s\"')
parser$add_argument('-Fl', '--flank',default=2000,type="integer", help='flank [default = %(default)s]')
parser$add_argument('-S', '--smooth',default=11,type="integer", help='smooth [default = %(default)s]')
parser$add_argument('-M', '--maxSize',default=750,type="integer", help='max fragment size [default = %(default)s]')
parser$add_argument('-O', '--out', default='.', help='Out directory for results [default = .]')

args = parser$parse_args()

library(ArchR)
library(cowplot)
## TSS 
tmp <- as.data.frame(fread(args$tssbed))
TSS <- GRanges(seqnames = tmp$V1, ranges=IRanges(start = tmp$V2,width=1),strand = tmp$V6)
values(tmp_TSS) <- DataFrame(tx_id = c(1:length(tmp$V1)),tx_name = tmp$V4)
TSS <- sort(sortSeqlevels(TSS))

chr <- paste0(seqnames(TSS))
chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
splitTSS <- split(resize(TSS, 1, "start"), seqnames(TSS))[chr]
window <- 2 * flank + 1

## get fragment
frag <- as.data.frame(fread(args$fragment,header=F))
grFrag <- GRanges(seqnames = frag$V1,ranges = IRanges(start = frag$V2, end = frag$V3),strand = "*")
values(grFrag) <- DataFrame(RG = frag$V4)

for (i in seq_along(chr)) {
    TSSi <- splitTSS[[chr[i]]]
    covi <-  grFrag[which(seqnames(grFrag) == chr[i]),] %>% sort %>% {
        coverage(IRanges(c(start(.), end(.)), width = 1))
    }
    if (i == 1) {
        sumTSS <- ArchR:::rleSumsStranded(list(chr1 = covi), 
          list(chr1 = TSSi), window, as.integer)
        fsi <- grFrag[which(seqnames(grFrag) == chr[i]),] %>% width %>% tabulate(nbins = argas$maxSize)

    }
    else {
        sumTSS <- sumTSS + ArchR:::rleSumsStranded(list(chr1 = covi), 
          list(chr1 = TSSi), window, as.integer)
        fsi <- fsi + grFrag[which(seqnames(grFrag) == chr[i]),] %>% width %>% tabulate(nbins = argas$maxSize)
    }
}
normBy <- mean(sumTSS[c(1:norm, (flank * 2 - norm + 1):(flank * 
    2 + 1))])
dfTSS <- DataFrame(group = args$group, x = seq_along(sumTSS) - 
    flank - 1, value = sumTSS, normValue = sumTSS/normBy, 
    smoothValue = ArchR:::.centerRollMean(sumTSS/normBy, 11))
dfFS <- DataFrame(group = args$group, fragmentSize = seq_along(fsi),
    fragmentPercent = round(100 * fsi/sum(fsi), 4))


plotTSSDF <- data.frame(x = dfTSS$x, v = dfTSS$smoothValue, 
    group = dfTSS$group)
plotTSSDF <- plotDF[sort(unique(c(1, seq(1, nrow(plotDF), 
    11), nrow(plotDF)))), , drop = FALSE]

plotFSDF <- data.frame(dfFS)

pal <- paletteDiscrete(values = unique(plotTSSDF$group))

p <- ggplot(plotTSSDF, aes(x, v, color = group)) + geom_line(size = 1) + 
    theme_ArchR()+ 
    guides(color= "none")+
    xlab("Distance From Center (bp)") + 
    ylab("Normalized Insertion Profile") + scale_color_manual(values = pal) + 
    scale_y_continuous(limits = c(0, max(plotDF$v) * 1.05), expand = c(0, 0)) + 
    scale_x_continuous(limits = c(min(plotDF$x), max(plotDF$x)), expand = c(0, 0))

ggsave(paste0(args$out,"/report/plot6_TSS.png"),p,width=5,height=5)

svg(paste0(args$out,"/report/plot6_TSS.png"),width=5,height=5)
print(p)
dev.off()

pal <- paletteDiscrete(values = unique(plotFSDF$group))

p1 <- ggplot(plotFSDF, aes(fragmentSize, fragmentPercent,
    color = group)) + geom_line(size = 1) + theme_ArchR() +
    guides(color= "none")+
    xlab("ATAC-seq Fragment Size (bp)") + ylab("Percentage of Fragments") +
    scale_color_manual(values = pal) + scale_y_continuous(limits = c(0,
    max(plotFSDF$fragmentPercent) * 1.05), expand = c(0,
    0)) + scale_x_continuous(limits = c(min(plotFSDF$fragmentSize),
    max(plotFSDF$fragmentSize)), expand = c(0, 0))

ggsave(paste0(args$out,"/report/plot5_InterSize.png"),p1,width=5,height=5)

svg(paste0(args$out,"/report/plot5_InterSize.svg"),width=5,height=5)
print(p1)
dev.off()