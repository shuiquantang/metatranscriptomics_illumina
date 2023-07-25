#karyoploteR
library(getopt)
library(karyoploteR)

spec = matrix(c(
 'ref_label', 'l', 1, "character",
 'read_length', 'r', 1, "numeric",
 'ref_size','s', 1, "numeric",
 'bam_file', 'b', 1, "character",
 'genome_name', 'n', 1, "character",
 'output_file', 'o', 1, "character"
 ), byrow=TRUE, ncol=4
)

opt = getopt(spec)

ref <- toGRanges(data.frame(chr=c(opt$ref_label), start=c(1), end=c(opt$ref_size)))

pdf(opt$output_file, width=10, height=6.4)

kp <- plotKaryotype(genome = ref)

tickDist = round(opt$ref_size/10000)*1000

kpAddBaseNumbers(kp, tick.dist = tickDist, add.units = TRUE)

kpPlotMarkers(kp, chr=opt$ref_label, x=(opt$ref_size/2), labels=(opt$genome_name),
              y=0.1, marker.parts = c(0, 0, 0),r0=-0.3, r1=0.3,
              line.color = "#FFAA22", label.color = "black", text.orientation="horizontal")

window_size <- 1e3

kp <- kpPlotBAMDensity(kp, data=opt$bam_file, window.size=window_size, col="red", data.panel=1, r0=-0.09)

kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density*opt$read_length/window_size, r0=-0.09)

kpAddLabels(kp, labels = "Read Depth", srt=90, pos=1, label.margin = 0.08)

print(kp)

dev.off()
