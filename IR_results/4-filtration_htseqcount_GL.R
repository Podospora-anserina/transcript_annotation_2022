# Detection and filtration of putative intron retention
###################################################################
# The aim of this script is to select and filterate putative intron retention


# Calls the needed functions 
source("functions.R")

# Defines the average and standard deviation thresholds
#thresold.intron.average <- 10
#thresold.intron.sd <- 3

# Reads average, standard deviation and reads per kilobase for CDS matrixes
df.average.counts <- read.table(file = "average_htseqcounts.tab",
                                      sep = "\t",
                                      header = T)
colnames(df.average.counts) <- unlist(strsplit(colnames(df.average.counts), ".txt"))

df.sd.counts <- as.matrix(read.table(file = "standard-deviation_htseqcounts.tab",
                                    sep = "\t",
                                    header = T))
colnames(df.sd.counts) <- unlist(strsplit(colnames(df.sd.counts), ".txt"))

df.rpk.CDS <- as.matrix(read.table(file = "RPK_CDS.tab",
                                  sep = "\t",
                                  header = T))

# Collapses average and standard deviation matrix
df.average.counts <- collapse_matrix(df.average.counts)
df.sd.counts <- collapse_matrix(df.sd.counts)


# Creates a data frame that will store intron features
df.select.RI <- data.frame(intron = df.average.counts[,1],
                           average = as.numeric(df.average.counts[,2]),
                           sd = as.numeric(df.sd.counts[,2]),
                           rpk = rep(0,times = nrow(df.average.counts)),
                           RNAseq = df.average.counts[,3])

# Search for IR in data, based on two thresholds, the number of
# detected IR is store in a table :
res_table <- NULL

for(T1 in c(5, 10, 20, 40)){
  for(T2 in c(1, 2, 3, 4, 5)){
    
    thresold.intron.average <- T1 
    thresold.intron.sd      <- T2    
    
    temp <- df.select.RI[(df.select.RI$average >= thresold.intron.average) &
                           (df.select.RI$sd <= thresold.intron.sd),]
    num_IR <- nrow(temp)
    
    # Search for a positive control : Pa_6_990.G_intron_1
    control = "Pa_6_990.G_intron_1" %in% temp[,1]
    
    res_table <- rbind(res_table, c(thresold.intron.average,
                                  thresold.intron.sd,
                                  num_IR, control))
  }
 
}

colnames(res_table) <- c("average", "sd", "#selected", "positive control")

write.table(res_table, file = "param_eval.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

thresold.intron.average <- 30
thresold.intron.sd <- 20

# Filters intron that doesn't match the average and standard deviation threshold
df.select.RI <- df.select.RI[(df.select.RI$average >= thresold.intron.average) &
                             (df.select.RI$sd <= thresold.intron.sd),]


# Loops for each intron
for(intron in 1:nrow(df.select.RI))
  
{
  # Gets the gene associated with the intron
  name.gene <- paste0(sub(pattern = ".G_intron.*",
                   df.select.RI[intron,1],
                   replacement = "\\1"), ".1")
  
  # Stores the read per kilobase (CDS) in the intron column of df.select.RI
  if(name.gene %in% row.names(df.rpk.CDS)){
    df.select.RI$rpk[intron] <- as.numeric(df.rpk.CDS[name.gene,df.select.RI$RNAseq[intron]])    
  }else{
    df.select.RI$rpk[intron] <- NA
  }

}

# Filters intron that has a gene with a rpk below 100
df.select.RI <- df.select.RI[df.select.RI$rpk >= 200,]
# Another filtering condition....
df.select.RI <- df.select.RI[(df.select.RI$average / df.select.RI$rpk) >= 0.1,]

# Search for the positive control
df.select.RI[df.select.RI[,1] == "Pa_6_990.G_intron_1",]

# Writes the selected intron retention events in a table
write.table(file = "selected_RI.tab",
            x = df.select.RI,
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)

# Lists gene names that show intron retention given the thresholds used
list.genes.RI <- unique(paste0(sub(pattern = ".G_intron.*",
                            x = unique(df.select.RI[,1]),
                            replacement = "\\1"), ".1"))

# Writes the intron retention associated gene names in a list
write.table(file = "list_genes_RI.txt",
            x = list.genes.RI,
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)

# Comparison with previous results
damien <- as.vector(as.matrix(read.table("list_genes_RI_Damien")))
sum(list.genes.RI %in% damien)/length(list.genes.RI)
sum(damien %in% list.genes.RI)/length(damien)
length(damien)


