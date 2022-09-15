# The aim of this script is to compile all tophat results in one data frame that
# will be next filtered from low quality splicing events

# Calls the needed functions
source("07-TopHat/splicing_condition/scripts/functions.R")

# Gets the tophat bed file names
list.tophat.files <- list.files("07-TopHat/splicing_condition/tophat")

# Initialization of the empty dataframe that will contain the tophat junctions bed files compilation
tophat.compilation <- NULL

# Loops for each tophat bed files
for(tophat.file in list.tophat.files)
{
  # Reads the tophat junctions bed file
  temp.bed.junction <- read.delim(paste0("07-TopHat/splicing_condition/tophat/",tophat.file),
                                  header = F,
                                  sep = "\t")[-1,]
  
  # Adds the sample which the tophat junctions file comes from.
  temp.bed.junction <- cbind(temp.bed.junction,
                             sub(pattern = "_junctions.bed",
                                 replacement = "\\1",
                                 x = tophat.file))
  
  # Compiles the current tophat junctions data frame to the main data frame tophat.compilation
  tophat.compilation <- rbind(tophat.compilation,
                              temp.bed.junction)
  
}

# Sorts the compiled tophat junctions data frame by chromosome and start position
tophat.compilation <- tophat.compilation[order(tophat.compilation[,1],
                                               tophat.compilation[,2]),]

# Extracts the left and right overhang
values.overhang <- as.numeric(unlist(strsplit(tophat.compilation[,11],split = ",")))

# Computes the start position of tophat junctions
tophat.compilation[,2] <- tophat.compilation[,2] + values.overhang[seq(1,length(values.overhang)-1,2)] + 1

# Computes the end position of tophat junctions
tophat.compilation[,3] <- tophat.compilation[,3] - values.overhang[seq(2,length(values.overhang),2)]

# Adds a new NA column for next junctions identifiers
tophat.compilation <- cbind(tophat.compilation,NA)

# Adds junctions identifiers in the 14th column
tophat.compilation[,14] <- paste(tophat.compilation[,1],
                                 tophat.compilation[,2],
                                 tophat.compilation[,3],
                                 sep = "_")

# Filters uninformative columns
tophat.compilation <- tophat.compilation[,c(1,2,3,5,6,13,14)]

# Defines columns names to tophat.compilation
colnames(tophat.compilation) <- c("Chromosome","Start","End","Coverage","Strand","Sample","Identifier")

# Writes tophat.compilation table into a file needed for the next script
write.table(x = tophat.compilation,
            file = "07-TopHat/splicing_condition/tables/compiled_tophat_junctions.tab",
            sep = "\t",
            row.names = F,
            col.names = T)

