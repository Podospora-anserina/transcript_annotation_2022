# Processing of compiled htseqcount data
###################################################################
# The aim of this script is to calculate the average and square-deviation of read
# counts respectively for each annotated intron.


# Calls the needed functions 
source("functions.R")

# Reads the compiled counts table
df.introns.counts <- as.matrix(read.table(file ="compiled_htseqcounts.tab",
                                        header = T,
                                        sep = "\t"))

# Sorts the compiled counts table per intron names
df.introns.counts <- df.introns.counts[order(rownames(df.introns.counts)),]

# Creates the matrixes that will store respectively average and standard deviation values
df.average.counts <- matrix(data = NA,
                        nrow = nrow(df.introns.counts),
                        ncol = ncol(df.introns.counts),
                        dimnames = list(rownames(df.introns.counts),
                                        colnames(df.introns.counts)))

df.sd.counts <- matrix(data = NA,
                      nrow = nrow(df.introns.counts),
                      ncol = ncol(df.introns.counts),
                      dimnames = list(rownames(df.introns.counts),
                                      colnames(df.introns.counts)))

# Loops for each RNA-seq data
for(name.RNAseq in colnames(df.introns.counts))
  
{

  # Computes and stores the average read counts value respectively for each intron
  df.average.counts[,name.RNAseq] <- by(df.introns.counts[,name.RNAseq],
                                        rownames(df.introns.counts),
                                        mean_values)
  
  # Computes and stores the standard deviation of read counts value respectively for each intron
 df.sd.counts[,name.RNAseq] <- by(df.introns.counts[,name.RNAseq],
                                  rownames(df.introns.counts),
                                  sd_values)
 
  
}

# Writes the average counts for each intron in a table
write.table(x = df.average.counts,
            file = "average_htseqcounts.tab",
            quote = F,
            row.names = T,
            col.names = T,
            sep = "\t")

# Writes the standard deviation counts for each intron in a table
write.table(x = df.sd.counts,
            file = "standard-deviation_htseqcounts.tab",
            quote = F,
            row.names = T,
            col.names = T,
            sep = "\t")


