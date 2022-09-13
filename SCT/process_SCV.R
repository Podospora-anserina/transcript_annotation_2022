### Sequential Coverage Values (SCV) methodology script  ###
###################################################################
# The aim of this script is to extract and select StringTie transcript predictions 
# according to the 2016 version of Podospora anserina genome annotation

# Calls the needed functions 
source("UTR_annotation/scripts/functions.R")

# Lists the column names of gff format annotation
list.gff.colnames <- c("chromosome","source","feature","start","end","score","strand","phase","attributes")

# Gets the stringtie gff file names
list.stringtie.files <- list.files("UTR_annotation/stringtie")

# Gets the RNA-seq names
list.RNAseq <- sub(pattern = ".gtf",
                   x = list.stringtie.files,
                   replacement = "\\1")

# Reads the Podospora anserina genome annotation (v.2016)
gff.PODANS <- read.table("annotations/genome_annotation_PODANS_v2016.gff",
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

# Creates an NA matrix that will stock UTR information for each gene
tab.SCV <- matrix(data = NA,
                  nrow = length(list.genes),
                  ncol = length(list.colnames),
                  dimnames = list(list.genes,list.colnames))


# Creates an NA matrix for TSS, TES and coverage data that will be extracted from stringtie files
tab.TSS <- matrix(data = NA,
                  nrow = length(list.genes),
                  ncol = length(list.RNAseq),
                  dimnames = list(list.genes,list.RNAseq))

tab.TES <- matrix(data = NA,
                  nrow = length(list.genes),
                  ncol = length(list.RNAseq),
                  dimnames = list(list.genes,list.RNAseq))

tab.cov <- matrix(data = NA,
                  nrow = length(list.genes),
                  ncol = length(list.RNAseq),
                  dimnames = list(list.genes,list.RNAseq))

# Creates two empty vectors that will store coding transcripts and non-coding transcripts annotations
gff.mRNA <- NULL
gff.nTARs <- NULL


# Loops in each stringtie files to intersect transcripts with annotated genes and
# extract their genomic coordinates
for(stringtie.file in list.stringtie.files)
  
{
  # Displays which stringtie files is currently processed
  print(stringtie.file)
  
  # Selects the RNA-seq name from which the stringtie files comes from
  sample.name <- list.RNAseq[stringtie.file == list.stringtie.files]
  
  # Reads the gff stringtie file
  gff.stringtie <- read.table(file = paste0("UTR_annotation/stringtie/",stringtie.file),
                              sep = "\t",
                              quote = "\"" ,
                              header = F,
                              col.names = list.gff.colnames)
  
  # Selects only transcripts stringtie annotation in the gff file
  gff.stringtie <- gff.stringtie[gff.stringtie$feature == "transcript",]
  
  # Adds the RNA-seq name associated with the processed stringtie file
  gff.stringtie$attributes <- paste(gff.stringtie$attributes," note ",sample.name,sep = "")
  
  # Loops in each chromosome to intersect stringtie transcripts with annotated genes
  for(chromosome in unique(gff.stringtie$chromosome))
    
  {
    # Selects the annotations from gff.stringtie and gff.PODANS corresponding to the
    # current processed chromosome in the loop
    gff.sub.stringtie <- gff.stringtie[gff.stringtie$chromosome == chromosome,]
    gff.sub.PODANS <- gff.PODANS[gff.PODANS$chromosome == chromosome,]
    
    
    # Loops in each annotated transcripts in the current sub stringtie file processed
    for(transcript in 1:nrow(gff.sub.stringtie))
    {
      
      # Gets the genomic coordinates of the annotated stringtie transcript
      tmp.start <- gff.sub.stringtie$start[transcript]
      tmp.end <- gff.sub.stringtie$end[transcript]
      
      # Intersects the stringtie transcript cordinates with the annotated elements cordinates
      tmp.elements <- gff.sub.PODANS[gff.sub.PODANS$start > tmp.start &
                                       gff.sub.PODANS$start < tmp.end |
                                       gff.sub.PODANS$end > tmp.start &
                                       gff.sub.PODANS$end < tmp.end,]
      
      # If there is only one annotated element intersecting with the stringtie transcript
      if(nrow(tmp.elements) == 1)
      {
        # If the annotated element is an annotated gene and the predicted transcript fully covers it
        if(tmp.elements[1,3] == "gene" & tmp.elements$start[1] > tmp.start & tmp.elements$end[1] < tmp.end)
          
        {
          # Gets the gene name intersecting with the stringtie transcript
          tmp.gene <- extract_names(tmp.elements[1,9],"ID=gene:"," ;.*","\\1")
          
          # Gets the coverage of the stringtie transcript intersecting with the gene
          tmp.cov <- as.numeric(extract_names(gff.sub.stringtie$attributes[transcript],".*cov ","; FPKM.*","\\1"))
          
          # Stocks the coverage value associated with the identified gene in tab.cov
          tab.cov[tmp.gene,sample.name] <- tmp.cov
          
          # If the strand of the identified gene is forward
          if(tmp.elements$strand[1] == "+")
          {
            # Stocks tmp.start as the TSS of the identified gene in tab.TSS
            tab.TSS[tmp.gene,sample.name] <- tmp.start
            
            # Stocks tmp.end as the TES of the identified gene in tab.TES
            tab.TES[tmp.gene,sample.name] <- tmp.end
            
          # Else if the strand of the identified gene is reverse  
          }else if(tmp.elements$strand[1] == "-")
          {
            # Stocks tmp.end as the TSS of the identified gene in tab.TSS
            tab.TSS[tmp.gene,sample.name] <- tmp.end
            
            # Stocks tmp.start as the TES of the identified gene in tab.TES
            tab.TES[tmp.gene,sample.name] <- tmp.start 
            
          }
        }
      # Else if there is no element intersecting with the stringtie transcript
      }else if(nrow(tmp.elements) == 0) 
      {
        # Intersects annotated elements cordinates with the stringtie transcript cordinates
        tmp.elements <- gff.sub.PODANS[gff.sub.PODANS$start < tmp.start &
                                         gff.sub.PODANS$end > tmp.start |
                                         gff.sub.PODANS$start < tmp.end &
                                         gff.sub.PODANS$end > tmp.end,]
        
        # If there is no genomic element where the stringtie transcript is nested
        if(nrow(tmp.elements) == 0)
        {
          
          # Adds the stringtie transcript annotation to the nTARs transcripts gff
          gff.nTARs <- rbind(gff.nTARs,
                             gff.sub.stringtie[transcript,])
        }
      }
    } # End of transcript loop
  } # End of chromosome loop
} # End of stringtie file loop


# Sorts gff.nTARs by chromosome and genomic cordinates
gff.nTARs <- gff.nTARs[order(gff.nTARs$chromosome,
                             gff.nTARs$start),]


# Writes the nTARs gff and bed file from gff.nTARs
write.table(x = gff.nTARs,
            file = "annotations/nTARs_raw.gff",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

# Converts gff.nTARs into a bed format
bed.nTARs <- cbind(gff.nTARs$chromosome,
                   as.numeric(gff.nTARs$start),
                   as.numeric(gff.nTARs$end),
                   ".",
                   500,
                   gff.nTARs[,7],
                   sub(pattern = ".*note ",
                       replacement = "\\1",
                       gff.nTARs[,9]))

# Writes bed.nTARs table in a bed file
write.table(x = bed.nTARs,
            file = "annotations/nTARs_raw.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)


# Initializes the threshold coverage values
threshold <- c(1,10,25,50,75,
               100,250,500,750,
               1000,2500,5000,7500,
               10000,15000,20000)


# Loops in each annotated genes
for(gene in rownames(tab.SCV))
  
{
  # Tests if there is at least one transcript predicted for the current gene
  if(sum(!is.na(tab.cov[gene,])) != 0)
    
  {
    # Gets the max coverage value of stringtie transcripts.
    coverage.max <- max(as.numeric(tab.cov[gene,!is.na(tab.cov[gene,])]))
    
    # Computes the differences between the max coverage value and the threshold values
    threshold.calc <- coverage.max - threshold
    
    # Selects the highest threshold chosen to select the stringtie transcripts
    threshold.max <- threshold[threshold.calc == min(threshold.calc[threshold.calc > 0])]
    
    # Selects the RNA-seq from where the stringtie transcript is/are above the threshold.
    RNAseq <- colnames(tab.cov)[tab.cov[gene,] >= threshold.max & !is.na(tab.cov[gene,])]
    
    # Gets the TSS positions and TES of selected stringtie transcripts
    TSS.position <- tab.TSS[gene,RNAseq]
    TES.position <- tab.TES[gene,RNAseq]
    
    # Computes the number of unique TSS and TES
    TSS.nb <- length(unique(TSS.position))
    TES.nb <- length(unique(TES.position))
    
    # Collapses the TSS and TES positions, as for the RNA-seq samples, respectively in one respective string
    TSS.position <- paste(TSS.position,collapse = ";")
    TES.position <- paste(TES.position,collapse = ";")
    RNAseq <- paste(RNAseq,collapse = ";")
    
    # Adds the values in the SCV matrix
    tab.SCV[gene,] <- c(TSS.position,TES.position,TSS.nb,TES.nb,RNAseq,threshold.max)
    
    
  }
  
  
}

# Writes the final SCV matrix
write.table(x = tab.SCV,
            file = "UTR_annotation/tables/SCV_UTR_PODANS.tab",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = T)


