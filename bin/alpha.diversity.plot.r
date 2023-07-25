#!/usr/bin/Rscript

library(ggplot2)
library(getopt)

spec = matrix(c(
 'input_file', 'i', 1, "character",
 'cat_pdf', 'a', 1, "character",
 'cat_png', 'b', 1, "character",
 'bar_pdf', 'c', 1, "character",
 'bar_png', 'd', 1, "character",
 'title', 't', 1, "character"
 ), byrow=TRUE, ncol=4
)

opt = getopt(spec)

# read the table
#print (opt$input_file)


rawdata <- read.delim(opt$input_file, sep='\t', header=TRUE)

rawdata$category <- as.character(rawdata$category)

rawdata$category <- factor(rawdata$category, levels=unique(rawdata$category))

rawdata$sample_id <- as.character(rawdata$sample_id)

rawdata$sample_id <- factor(rawdata$sample_id, levels=unique(rawdata$sample_id))

p<-ggplot(rawdata, aes(x=sample_id, y=shannon_index)) + geom_bar(stat="identity") + coord_flip() + xlab(opt$title) + theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), axis.text=element_text(size=14)) 

sample_size <- length(rawdata$category)

pdf_width <- 2+0.25*sample_size

pdf(opt$bar_pdf, width=7, height=pdf_width)

print(p)

dev.off()


png_width <- 4+0.5*sample_size

png(opt$bar_png, width=14, height=png_width, units = "cm", pointsize = 10, bg = "white", res=1200)

print(p)

dev.off()


category_size <- length(unique(rawdata$category))

if (category_size > 1){


pdf_width <- 4+category_size

p<-ggplot(rawdata, aes(x=category, y=shannon_index)) + geom_boxplot(color="red") + geom_dotplot(binaxis='y', stackdir='center')  + ylab(opt$title) + theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text=element_text(size=14)) + expand_limits(y=0)

pdf(opt$cat_pdf, width=pdf_width, height=5)

print(p)

dev.off()

png_width <- 6+2*category_size

png(opt$cat_png, width = png_width, height = 8, units = "cm", pointsize = 10, bg = "white", res=1200)

print(p)

dev.off()

}
