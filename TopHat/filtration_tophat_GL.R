# The aim of this script is to filtrate low quality tophat junctions in the 
# compiled junctions data frame

# Calls the needed functions
source("splicing_condition/scripts/functions.R")

# Annotation file of Podospora anserina introns made from the 2016 version of 
# Podospora anserina genome annotation
gff.introns <- read.table("annotations/introns_annotation_PODANS_v2016.gff")

# Adds a new NA column for next introns identifiers
gff.introns <- cbind(gff.introns,NA)

# Defines the 14th column name
colnames(gff.introns)[10] <- "Identifier"

# Adds introns identifiers in the 10th column
gff.introns[,10] <- paste(gff.introns[,1],
                          gff.introns[,4],
                          gff.introns[,5],
                          sep = "_")


# Reads the table with compilated tophat junctions made from compilation_tophat.R
tophat.compilation <- read.table(file = "splicing_condition/tables/compiled_tophat_junctions.tab",
                                 sep = "\t",
                                 header = T)

# Lists all junctions that are present in at least two RNA-seq for a first filtering
id.junctions.selected <- matrix(unique(tophat.compilation[duplicated(tophat.compilation$Identifier),"Identifier"]),
                                     ncol = 1)

# Filters all junctions present in only one RNA-seq
tophat.filtration <- tophat.compilation[tophat.compilation$Identifier %in% id.junctions.selected[,1],]

# Compiles tophat junctions
#tophat.filtration <- as.data.frame(t(apply(id.junctions.selected,1,compile_junctions,tophat.filtration)))

# GL :: Same calculation, but this time we know how long we have to wait 
# for final results. 
temp = NULL
for(i in 1:nrow(id.junctions.selected)){
  print(paste(i, "/", nrow(id.junctions.selected)))
  temp <- rbind(temp, compile_junctions(id.junctions.selected[i,], tophat.filtration))
}

tophat.filtration <- as.data.frame(temp)

# Gets the identifier of previous annotated introns present in tophat junctions
id.introns.junctions <- gff.introns[gff.introns[,10] %in% tophat.filtration[,7],10]

# Filters the annotated introns from the tophat junctions matrix
tophat.filtration <- tophat.filtration[!tophat.filtration[,7] %in% id.introns.junctions,]

# Adds two columns that will store the coverage and dataset results of filtering tests
tophat.filtration <- cbind(tophat.filtration,NA,NA)

# Uses the functions to filtrate the low quality tophat junctions that have 
# coverage below 5 and less than 2 independents RNA-seq (from 2 different dataset)
tophat.filtration[,8]<- apply(tophat.filtration,1,max_values)

# Filters tophat junctions having max coverage value below 5 and with only one
# dataset where it was detected
tophat.filtration <- tophat.filtration[tophat.filtration[,8] >= 5 ,] 

# Writes the matrix of tophat junctions that passed through the filtration step.
write.table(x = tophat.filtration,
            file = "splicing_condition/tables/selected_tophat_junctions.tab",
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)
