### Sequential Coverage Values (SCV) methodology script  ###
###################################################################
# The aim of this script is to extract and select StringTie transcript predictions 
# according to the 2016 version of Podospora anserina genome annotation

# Calls the needed functions 
source("functions.R")

# Lists the column names of gff format annotation
list.gff.colnames <- c("chromosome","source","feature","start","end","score","strand","phase","attributes")

# Gets the stringtie gff file names
list.stringtie.files <- list.files("stringtie_dataset-B")

# Gets the RNA-seq names
list.RNAseq <- sub(pattern = ".gtf",
                   x = list.stringtie.files,
                   replacement = "\\1")

# Reads the Podospora anserina genome annotation (v.2016)
gff.PODANS <- read.table("genome_annotation_PODANS_v2016.gff",
                         sep = "\t",
                         header = F,
                         quote = "\"" ,
                         comment.char = "",
                         col.names = list.gff.colnames)

# Filters CDS & transcript annotations that have redundant information with gene annotations
gff.PODANS <- gff.PODANS[gff.PODANS$feature != "CDS" ,]

# Lists all gene names
list.genes <- extract_names(gff.PODANS[gff.PODANS[,3] == "gene",9],"ID=gene:"," ;.*","\\1")

# Lists the column names of the next matrix
list.colnames <- c("pos.TSS","pos.TES","nb.TSS","nb.TES","RNA-seq","coverage")

# Creates two empty vectors that will store coding transcripts and non-coding transcripts annotations
gff.NAT <- NULL

# Loops in each stringtie files to intersect transcripts with annotated genes and
# extract their genomic coordinates
for(stringtie.file in list.stringtie.files){
  
  # Displays which stringtie files is currently processed
  print(stringtie.file)
  
  # Selects the RNA-seq name from which the stringtie files comes from
  sample.name <- list.RNAseq[stringtie.file == list.stringtie.files]
  
  # Reads the gff stringtie file
  gff.stringtie <- read.table(file = paste0("stringtie_dataset-B/",stringtie.file),
                              sep = "\t",
                              quote = "\"" ,
                              header = F,
                              col.names = list.gff.colnames)
  
  # Selects only transcripts stringtie annotation in the gff file
  gff.stringtie <- gff.stringtie[gff.stringtie$feature == "transcript",]
  
  # Adds the RNA-seq name associated with the processed stringtie file
  gff.stringtie$attributes <- paste(gff.stringtie$attributes," note ",sample.name,sep = "")
  
  # Loops in each chromosome to intersect stringtie transcripts with annotated genes
  for(chromosome in unique(gff.stringtie$chromosome)){
    
    # Selects the annotations from gff.stringtie and gff.PODANS corresponding to the
    # current processed chromosome in the loop
    gff.sub.stringtie <- gff.stringtie[gff.stringtie$chromosome == chromosome,]
    gff.sub.PODANS <- gff.PODANS[gff.PODANS$chromosome == chromosome,]
    
    # Loops in each annotated transcripts in the current sub stringtie file processed
    for(transcript in 1:nrow(gff.sub.stringtie)){
      
      # Gets the genomic coordinates of the annotated stringtie transcript
      tmp.start <- gff.sub.stringtie$start[transcript]
      tmp.end <- gff.sub.stringtie$end[transcript]
      
      # Intersects the stringtie transcript cordinates with the annotated elements cordinates
      tmp.elements <- gff.sub.PODANS[gff.sub.PODANS$start > tmp.start &
                                       gff.sub.PODANS$start < tmp.end |
                                       gff.sub.PODANS$end > tmp.start &
                                       gff.sub.PODANS$end < tmp.end,]
      
      # If there is only one annotated element intersecting with the stringtie transcript
      if(nrow(tmp.elements) == 1){
        
        # If the annotated element is an annotated gene and the predicted transcript fully covers it
        if(tmp.elements[1,3] == "gene" & tmp.elements$start[1] > tmp.start & tmp.elements$end[1] < tmp.end){
            ## GL :: Write code here
            if(gff.sub.stringtie[transcript,7] != as.character(tmp.elements[7])){
              gff.NAT <- rbind(gff.NAT, c(gff.sub.stringtie[transcript,],
                                          tmp.elements, "over"))          
            }
        }  
      # Else if there is no element intersecting with the stringtie transcript
      }else if(nrow(tmp.elements) == 0){
        # Intersects annotated elements coordinates with the stringtie transcript cordinates
        tmp.elements <- gff.sub.PODANS[gff.sub.PODANS$start < tmp.start &
                                         gff.sub.PODANS$end > tmp.start |
                                         gff.sub.PODANS$start < tmp.end &
                                         gff.sub.PODANS$end > tmp.end,]
        
        # If there is one genomic element, and only one...
        if(nrow(tmp.elements) == 1){
          ## GL : Write code here
          if(gff.sub.stringtie[transcript,7] != as.character(tmp.elements[7])){
            gff.NAT <- rbind(gff.NAT, c(gff.sub.stringtie[transcript,],
                                        tmp.elements, "inside"))          
          }
        }
      }
    } # End of transcript loop
  } # End of chromosome loop
} # End of stringtie file loop

# Sorting of the results
gff.NAT <- gff.NAT[order(unlist(gff.NAT[,1]), unlist(gff.NAT[,4])),]

# Write the potential NAT in a file 
write.table(x = gff.NAT,
            file = "potential_NAT-V2.tab",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)
