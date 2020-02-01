#! /usr/bin/env Rscript

rm(list=ls())
library(ggplot2)
library(RColorBrewer) 
library(argparse)
library(reshape2)
library(stringr)
library(dplyr)
library(ggsci)


parser <- ArgumentParser(description="generate abundance plots")
parser$add_argument("-d", "--directory", help="location of count files")
parser$add_argument("-f","--feature", default="organism", help="colname to be sorted by. Used for annotating the plots")
parser$add_argument("-o","--outdir", default=getwd(), help="Output directory, default is working directory.")
parser$add_argument("-t","--top", default=10, help="number of features to display for each sample. Default = 10")
args <- parser$parse_args()

#To run on the command line comment out the this section
############################################# Interactive options ###########################################################
args$feature <- "Organisms"                                                   
args$directory <- "/home/drewx/Documents/DevOps/assembly-freeX/organism"
args$outdir <- "/home/drewx/Documents/DevOps/assembly-freeX/analysis.organism"

args$feature <- "Function"
args$directory <- "/home/drewx/Documents/DevOps/assembly-freeX/functions"
args$outdir <- "/home/drewx/Documents/DevOps/assembly-freeX/analysis.functions"
args$top = 50
############################################################################################################################

#https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

directory <- args$directory
feature <- args$feature
top <- args$top
outdir <- args$outdir
top <- args$top


dir.create(outdir)
setwd(outdir)

source("/home/drewx/Documents/misc-scripts/skunkworks.R")

merged_csv <- merge_data(directory, val_col = 2)

#merged_df  <- merged_csv %>% arrange(desc(A2a.tab_function.tsv)) 
merged_df  <- merged_csv %>% arrange(desc(A2a.tab_organism.tsv)) 

filename <- paste0(feature, "_top_perc", top, ".tsv")
write.table(merged_df, file = filename,  sep = "\t", row.names = FALSE, quote = F)  
colnames(merged_df)[1] <- c(feature)
colnames(merged_df)[c(2:6)] <- c("D1rep1","D2rep2","D3rep3","D2","D3")
top_v  = 50
for (col in 2:ncol(merged_df)){
    data <- head(merged_df[feature][order(-merged_df[col]),], top_v)
  if (col == 2)
  {
    top <- data
  }
  else{
    top <- union(top, data)
  }
}

topx  <- data.frame(top)

subset_merged <- merged_df[merged_df[[feature]] %in% topx$top,]
nr <- nrow(subset_merged)
subset_merged[nr+1,1] <- "Other"

for (col in 2:ncol(subset_merged)){
    subset_merged[nr+1,col] <- 100 - sum(subset_merged[c(1:nr), col])

}

#subset_merged <- subset_merged[-c(1,2,6, nrow(subset_merged)),] %>%  arrange(`D1-rep1`)
#subset_merged <- subset_merged[-c(nrow(subset_merged)),] %>%  arrange(`D1-rep1`)

subset_merged <- subset_merged[-c(nrow(subset_merged)),] %>%  
                                    arrange(desc(D1rep1)) %>%
                                    rowwise() %>%
                                    mutate(D1_avg = mean(c(D1rep1, D1rep2, D1-rep3))) 


%>%  
                                    select(D1_avg,  D2, D3)

subset_merged
subset_merged[,feature] <- factor(subset_merged[,feature], levels = subset_merged[,feature])
#subset_merged <- subset_merged[-c(nrow(subset_merged)),] %>%  arrange(`D1-rep1`)
#subset_merged$Function <- factor(subset_merged$Function,  levels = subset_merged$Function[order(subset_merged$Function)])  
#df$derma <- factor(df$derma, levels = df$derma[order(df$prevalence)])

melt_merged_df <- melt(data = subset_merged, id.vars = feature )
Colors =c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075')
colourCount = length(unique(melt_merged_df[,feature]))
#getPalette = colorRampPalette(brewer.pal(8, Dark))
getPalette = colorRampPalette(Colors)
mycols <- getPalette(colourCount)
names(mycols) <- unique(melt_merged_df[,feature])
mycols["Other"] <- "#ffffff"

p <- ggplot(data = melt_merged_df, aes_string(x = "variable" , y = "value", fill=feature)) + 
        geom_bar(stat = 'identity',  colour="black") + 
        scale_fill_manual(values = mycols, name = Hmisc::capitalize(feature)) +
        xlab("Sample") + 
        ylab("Relative abundance (%)" ) +
        theme_bw() +
        scale_fill_simpsons(name = Hmisc::capitalize(feature)) +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=16,  vjust = 1, face="bold"),
              legend.text=element_text(size=16, face = "italic"),
              legend.title=element_text(size=18, face="bold"),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              panel.grid = element_blank())
#
# 
fname <- paste0("top_",feature,"_1x.svg")

setwd("/home/drewx/Dropbox/PhDX/Chapter4/Figures/samsa_function")

ggsave(fname,
       width = 200,
       units = "mm")

#
#labels = function(x) str_wrap(x, width = 40
p  <- ggplot(data = melt_merged_df, aes_string(x = "variable" , y = "value", fill=feature)) + 
        scale_fill_manual(values = mycols, name = Hmisc::capitalize(feature)) +
        geom_bar(stat = 'identity', colour="black") + 
        xlab("Sample") + 
        ylab("Relative abundance (%)" ) +
        theme_bw() +
        theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16,  vjust = 1, face="bold"),
            legend.title=element_text(size=14, face="bold"),
            legend.text = element_text(size=12),
            panel.grid = element_blank(),
            legend.direction="vertical",
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        guides(linetype=guide_legend(nrow=2)) 
  
p

fname <- paste0("top_",feature,"_2x.tiff")

setwd("/home/drewx/Dropbox/PhDX/Chapter4/Figures/samsa_function")

ggsave(fname,
       units = "mm",
       dpi = 1000, 
       compression = "lzw")




# i = 1
# col_names <- list()
# for (count_file in feature_files){
#   count_df <- read.table(count_file, sep="\t", stringsAsFactors = FALSE, quote = "")
#   colnames(count_df) <- c("percent","count", feature )
#   top_feature <- head(count_df, top)
#   #print(count_file)
#   #print(top_feature)
#   top_feature[nrow(top_feature)+1,] <- list(100 - sum(top_feature$percent), 0, "Other")
#   if (i > 1) {
#     merged_df <- merge(merged_df, temp_df, by=feature, all = TRUE)
#   }
#   else{
#     merged_df = top_feature[,c(1,3)]
#   }
#   temp_df <-  top_feature[,c(1,3)]
#   basefname <- tools::file_path_sans_ext(basename(count_file))
#   col_names[i] <- unlist(strsplit(basefname,"[.]"))[1]
#   i = i + 1
#   print(head(merged_df))  
# }
# 
# 
# merged_df[is.na(merged_df)] <- 0
# colnames(merged_df)[c(2:ncol(merged_df))] <- col_names 
# merged_df <- merged_df[order(-merged_df[,2]),]
# other_index <- which(merged_df[feature]=="Other")
# merged_df <- rbind(merged_df[-other_index,], merged_df[other_index,] )
# merged_csv <- merged_df
# merged_csv[c(2:ncol(merged_df))] <- round(merged_df[c(2:ncol(merged_df))], 2)

