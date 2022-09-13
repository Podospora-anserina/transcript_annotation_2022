### Analysis of SCV methodology results ###
###################################################################
# The aim of this script is to analyze SCV results from stringTie transcript 
# predictions obtained through the process_SCV.R script

# Calls the needed functions 
source("./06-SCT_results/UTR_annotation/scripts/functions.R")

# Lists the column names of gff format annotation
list.gff.colnames <- c("chromosome","source","feature","start","end","score","strand","phase","attributes")

# Gets the Podospora anserina genome annotation (v.2016)
gff.PODANS <- read.table("./06-SCT_results/annotations/genome_annotation_PODANS_v2016.gff",
                         sep = "\t",
                         header = F,
                         quote = "\"" ,
                         comment.char = "",
                         col.names = list.gff.colnames)

# Gets the SCV results matrix 
tab.SCV <- read.table(file = "./06-SCT_results/UTR_annotation/tables/SCV_UTR_PODANS.tab",
                         header = T,
                         sep = "\t")

# Filters genes that have no predicted transcript
tab.SCV <- tab.SCV[!is.na(tab.SCV[,1]),]

# Creates an empty vector that will contains a matrix
tab.SCV.RNAseq <- NULL

# Creates an empty vector that will contains 5'UTR coordinates
coord_5UTR <- NULL

# Creates an empty object that will contains the gff annotation of TSS, TES and UTRs
gff.stringtie <- NULL

# Loops in each gene with at least a predicted transcript
for(gene in 1:nrow(tab.SCV)){
  # Gets the gene name
  gene.name <- rownames(tab.SCV)[gene]
  
  # Gets the pattern of the gene name
  element.name <- sub(pattern = "\\.G",
                      replacement = "\\1",
                      x = gene.name)
  
  # Gets the gene and CDS annotations from gff.PODANS
  gff.gene <- gff.PODANS[grepl(pattern = gene.name,x = gff.PODANS$attributes,fixed = T)&
                           gff.PODANS$feature == "gene",]
  
  gff.CDS <- gff.PODANS[grepl(pattern = gene.name,x = gff.PODANS$attributes,fixed = T)&
                          gff.PODANS$feature == "CDS",]
  
  # Gets the chromosome, strand, TSS, TES, Start and stop positions,RNA-seq 
  # sources and the number of predicted transcripts
  chromosome <- gff.gene$chromosome[1]
  TSS.position<- as.numeric(extract_values(tab.SCV[gene,"pos.TSS"],";"))
  TES.position<- as.numeric(extract_values(tab.SCV[gene,"pos.TES"],";"))
  RNAseq <- extract_values(tab.SCV[gene,"RNA.seq"],";")
  Nb.transcripts <- extract_length(tab.SCV[gene,"RNA.seq"],";")
  Strand <- gff.gene$strand[1]
  Nb.exons <- nrow(gff.CDS)
  seq.exons <- 1:Nb.exons
  
  # Selects the position of the TSS and TES according to the gene strand
  if(Strand == "+") # Tests if the gene is in the forward strand
  {
    start.position <- as.numeric(gff.gene$start[1])
    stop.position <- as.numeric(gff.gene$end[1])
    
    CDS.position.1 <- as.numeric(gff.CDS$end[1])
    CDS.position.2 <- as.numeric(gff.CDS$start[nrow(gff.CDS)])

  }else # Else if the gene is in the reverse strand
  {
    start.position <- as.numeric(gff.gene$end[1])
    stop.position <- as.numeric(gff.gene$start[1])
    
    CDS.position.1 <- as.numeric(gff.CDS$start[nrow(gff.CDS)])
    CDS.position.2 <- as.numeric(gff.CDS$end[1])

  }
  
  # Creates a gff annotation line for each TSS and TES, also one 5'UTR and 3'UTR 
  # with the furthest TSS/TES from the CDS
  gff.TSS <- cbind(chromosome,"Stringtie","TSS",TSS.position,TSS.position,".",Strand,".",
                   paste0("Parent=gene:",gene.name,"; note ",RNAseq,"; color #004c00;"))
  
  gff.TES <- cbind(chromosome,"Stringtie","TES",TES.position,TES.position,".",Strand,".",
                   paste0("Parent=gene:",gene.name,"; note ",RNAseq,"; color #004c00;"))
  
  colnames(gff.TSS) <- list.gff.colnames
  colnames(gff.TES) <- list.gff.colnames
  
  # Computes UTRs size 
  size.5UTR <- abs(TSS.position - start.position) + 1
  size.3UTR <- abs(TES.position - stop.position) + 1
  
  #get 5UTR coordinates
  fiveUTR_coord <- data.frame(gene,TSS.position,start.position,size.5UTR,gene.name)
  coord_5UTR <- rbind(coord_5UTR, fiveUTR_coord)
  
  
  # Selects the furthest TSS from the start codon and the furthest TES from the stop codon
  start.5UTR <- unique(TSS.position[size.5UTR == max(size.5UTR)])
  end.3UTR <- unique(TES.position[size.3UTR == max(size.3UTR)])
  
  # Creates a gff annotation for CDS features of the gene
  gff.CDS$attributes <- paste0("ID=CDS:",element.name,".P1"," ;Parent=mRNA:",element.name,".T1",
                               sub(x = gff.CDS$attributes,
                                   pattern = ".*;function",
                                   replacement = " ;function"))
  
  # Writes attributes column of the gene and its elements 
  if(Strand == "+")
  {
    gff.mRNA <- cbind(chromosome,"Stringtie","mRNA",start.5UTR,end.3UTR,".",Strand,".",
                      paste0("ID=mRNA:",element.name,".T1 ;Parent=gene:",gene.name))
    
    gff.gene <- cbind(chromosome,"Stringtie","gene",start.5UTR,end.3UTR,".",Strand,".",gff.gene$attributes[1])
    
    colnames(gff.mRNA) <- list.gff.colnames
    colnames(gff.gene) <- list.gff.colnames
    
    if(Nb.exons == 1)
    {
      gff.exon <- cbind(chromosome,"Stringtie","exon",start.5UTR,end.3UTR,".",Strand,".",
                        paste0("ID=exon:",element.name,".T1-1 ;Parent=mRNA:",element.name,".T1"))

    }else if(Nb.exons > 1)
    {
      gff.exon <- rbind(cbind(chromosome,"Stringtie","exon",start.5UTR,CDS.position.1,".",Strand,".",
                              paste0("ID=exon:",element.name,".T1-1 ;Parent=mRNA:",element.name,".T1")),
                        cbind(chromosome,"Stringtie","exon",CDS.position.2,end.3UTR,".",Strand,".",
                              paste0("ID=exon:",element.name,".T1-",Nb.exons," ;Parent=mRNA:",element.name,".T1")))
      
      if(Nb.exons > 2)
      {
        colnames(gff.exon) <- list.gff.colnames
        sub.seq.exons <- seq.exons[c(-1,-Nb.exons)]
        sub.exons <- gff.CDS[sub.seq.exons,]
        sub.exons$feature <- "exon"
        sub.exons[,9] <- paste0("ID=exon:",element.name,".T1-",sub.seq.exons," ;Parent=mRNA:",element.name,".T1")
        gff.exon <- rbind(gff.exon,sub.exons)
      }
      
    }

  }else
  {
    gff.mRNA <- cbind(chromosome,"Stringtie","mRNA",end.3UTR,start.5UTR,".",Strand,".",
                      paste0("ID=mRNA:",element.name,".T1 ;Parent=gene:",gene.name))
    
    gff.gene <- cbind(chromosome,"Stringtie","gene",end.3UTR,start.5UTR,".",Strand,".",gff.gene[1,9])
    
    
    if(Nb.exons == 1)
    {

    gff.exon <- cbind(chromosome,"Stringtie","exon",end.3UTR,start.5UTR,".",Strand,".",
                      paste0("ID=exon:",element.name,".T1-1 ;Parent=mRNA:",element.name,".T1"))

    }else if(Nb.exons > 1)
    {

      gff.exon <- rbind(cbind(chromosome,"Stringtie","exon",CDS.position.1,start.5UTR,".",Strand,".",
                              paste0("ID=exon:",element.name,".T1-1 ;Parent=mRNA:",element.name,".T1")),
                        cbind(chromosome,"Stringtie","exon",end.3UTR,CDS.position.2,".",Strand,".",
                              paste0("ID=exon:",element.name,".T1-",Nb.exons," ;Parent=mRNA:",element.name,".T1")))
    
      if(Nb.exons > 2)
      {
        colnames(gff.exon) <- list.gff.colnames
        sub.seq.exons <- seq.exons[c(-1,-Nb.exons)]
        sub.exons <- gff.CDS[sub.seq.exons,]
        sub.seq.exons <- rev(sub.seq.exons)
        sub.exons$feature <- "exon"
        sub.exons[,9] <- paste0("ID=exon:",element.name,".T1-",sub.seq.exons," ;Parent=mRNA:",element.name,".T1")
        gff.exon <- rbind(gff.exon,sub.exons)
      }
    }
  }
  
  colnames(gff.gene) <- list.gff.colnames
  colnames(gff.mRNA) <- list.gff.colnames
  colnames(gff.exon) <- list.gff.colnames
  
  # Binds all gff annotations in gff.stringtie
  gff.stringtie <- rbind(gff.stringtie,
                         gff.gene,
                         gff.mRNA,
                         gff.CDS,
                         gff.exon,
                         gff.TSS,
                         gff.TES)
  
  # Binds all values in one row and then in the final matrix tab.SCV.RNAseq
  sub_SCV.RNAseq <- cbind(gene,TSS.position,TES.position,RNAseq,size.5UTR,size.3UTR,Nb.transcripts,Strand)
  tab.SCV.RNAseq <- rbind(tab.SCV.RNAseq,sub_SCV.RNAseq)
  
  
}

# Order the matrix containing gff annotation
gff.stringtie <- gff.stringtie[order(gff.stringtie$chromosome,
                                     as.numeric(gff.stringtie$start)),]

# Write the stringtie gff annotation in a file
write.table(x = gff.stringtie,
            file = "./06-SCT_results/annotations/genome_annotation_PODANS_SCV.gff",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

# Gets all gene names of P. anserina
list.genes <- extract_names(gff.PODANS[gff.PODANS$feature == "gene",9],
                            "ID=gene:"," ;function.*",replacement = "\\1")


# Converts tab.SCV.RNAseq into a dataframe
df.SCV.RNAseq <-data.frame(Gene = tab.SCV.RNAseq[,1],
                           TSS.pos = as.numeric(tab.SCV.RNAseq[,2]),
                           TES.pos = as.numeric(tab.SCV.RNAseq[,3]),
                           RNAseq = tab.SCV.RNAseq[,4],
                           Size.5UTR = as.numeric(tab.SCV.RNAseq[,5]),
                           Size.3UTR = as.numeric(tab.SCV.RNAseq[,6]),
                           Transcripts.Nb = as.numeric(tab.SCV.RNAseq[,7]),
                           Categories = NA,
                           Strand = tab.SCV.RNAseq[,8],
                           Dataset = NA)

# Lists the genes with at least one transcript predicted by stringtie
list.genes.stringtie <- unique(df.SCV.RNAseq[,1])

# Lists the genes without any transcript prediction
list.genes.no.prediction <- list.genes[!list.genes %in% list.genes.stringtie]

# Defines tags for each categories based on the number of predicted transcripts
SCV.categories <- c(paste("n = ",1:9,sep = ""),
                    paste("n ",sprintf("\u2265 "),10,sep = ""))


# Defines dataset names
SCV.dataset <- c("A","B","C")

# Defines tags for each dataset
SCV.RNAseq <- c("ERR","SRR6","SRR3")

# Creates a dataframe that will contain the number of predictions according to their dataset and category
df.predictions <- data.frame(Categories = rep(SCV.categories,3),
                             Number = 0,
                             Sample = rep(SCV.dataset,
                                          rep(10,3)))

# Creates a dataframe that will contain distance between closest and furthest UTR
df.distance.UTR <- data.frame(gene = rep(list.genes.stringtie,
                                         2),
                              UTR = rep(c("TSS","TES"),
                                        c(length(list.genes.stringtie),
                                          length(list.genes.stringtie))),
                              distance = 0)

# Loops in each gene name
for(gene in unique(df.SCV.RNAseq[,1]))
{
  # Extracts the rows from df.SCV.RNAseq associated with gene
  sub.df <- df.SCV.RNAseq[df.SCV.RNAseq[,1] == gene,]
  
  # Gets the number of predictions
  predictions.size <- nrow(sub.df)
  
  # Computes the distance between the furthest and the closest TSS from the Start codon
  df.distance.UTR[df.distance.UTR[,1] == gene &
                    df.distance.UTR[,2] == "TSS",3] <- max(sub.df$Size.5UTR) - min(sub.df$Size.5UTR)
  
  # Computes the distance between the furthest and the closest TES from the STOP
  df.distance.UTR[df.distance.UTR[,1] == gene &
                    df.distance.UTR[,2] == "TES",3] <- max(sub.df$Size.3UTR) - min(sub.df$Size.3UTR)

  # Loops in each RNAseq pattern
  for(RNAseq in SCV.RNAseq)
  {
    # Computes the number of predictions associated with the RNAseq pattern
    sample.size <- sum(grepl(pattern = RNAseq,sub.df[,4],fixed = T))
    
    # Computes the total of predictions associated with the dataset/RNAseq pattern following its category
    
    # If the number of predictions is between 1 and 9 
    if(predictions.size > 0 & predictions.size < 10)
    {
      df.predictions[predictions.size+10*(which(RNAseq == SCV.RNAseq)-1),"Number"] <- df.predictions[predictions.size+10*(which(RNAseq == SCV.RNAseq)-1),"Number"] + sample.size
      
    # Else if the number of predictions is higher or equal to 10
    }else if(predictions.size >= 10)
      
    {
      df.predictions[10+10*(which(RNAseq == SCV.RNAseq)-1),"Number"] <- df.predictions[10+10*(which(RNAseq == SCV.RNAseq)-1),"Number"] + sample.size
    }
  }
}

# Sorts the factor levels of df.predictions
df.predictions$Categories <- factor(df.predictions$Categories,
                                    levels = SCV.categories)

# Gets the number of predicted transcripts for each gene
nb.transcripts <- by(tab.SCV[,"RNA.seq"],rownames(tab.SCV),extract_length,";")

# Creates a dataframe that will contain total number of predictions according to their category
df.total <- data.frame(Categories = SCV.categories,
                       Number = rep(0,10))

# Loops in each category
for(category in 1:10)
  
{
  # Computes the total number of predictions per category
  if(category < 10){ df.total[category,2] <- sum(nb.transcripts == category)
  }else if(category == 10)
  {
    df.total[category,2] <- sum(nb.transcripts >= category)
    
  }
}

# Creates a tibble grouped by UTR type
tibble(df.distance.UTR) %>% group_by(UTR) -> tibble.distance.UTR

# Computes mean distance between furthest and closest TSS/TES from CDS
tibble.distance.UTR %>% summarize(distance.mean = mean(distance)) -> df.mean.distance

# Filters values higher than the 99th centile
tibble.distance.UTR %>%  filter(distance <= quantile(distance,0.99)) -> tibble.distance.UTR

# Sorts the factor levels in tibble.distance.UTR
tibble.distance.UTR$UTR <- factor(tibble.distance.UTR$UTR,
                                  levels = c("TSS","TES"))


# Loops in each RNAseq pattern
for(RNAseq in SCV.RNAseq)
{
  # Stocks the dataset name associated with the RNAseq pattern
  df.SCV.RNAseq[grepl(pattern = RNAseq,
                      x = df.SCV.RNAseq$RNAseq,
                      fixed = T),"Dataset"] <- SCV.dataset[RNAseq == SCV.RNAseq]
  
}


# Defines the categories of transcript with only one or multiple predictions
df.SCV.RNAseq[df.SCV.RNAseq[,7] == 1,8] <- "N = 1"
df.SCV.RNAseq[df.SCV.RNAseq[,7] > 1,8] <- "N > 1"

# Creates a data frame used for Figure 4.A (UTR sizes considering the number of predictions)
df.UTR <- data.frame(UTR.size = c(df.SCV.RNAseq$Size.5UTR,
                                  df.SCV.RNAseq$Size.3UTR),
                     UTR.end = c(rep("5'UTR",nrow(df.SCV.RNAseq)),
                                 rep("3'UTR",nrow(df.SCV.RNAseq))),
                     UTR.categories = rep(df.SCV.RNAseq$Categories,2))


# Creates a tibble grouped by UTRs and categories
tibble(UTR.size = df.UTR$UTR.size,
       UTR.end = as.character(df.UTR$UTR.end),
       UTR.categories = df.UTR$UTR.categories) %>% group_by(UTR.end,UTR.categories) -> tibble.UTR


# Computes average for each categories of UTR
tibble.UTR %>% group_by(UTR.end,UTR.categories) %>% summarize(UTR.mean = mean(UTR.size)) -> tibble.UTR.mean

# Filters values higher than the 99th centile 
tibble.UTR %>% filter(UTR.size <= quantile(UTR.size,0.99)) -> tibble.UTR

# Sorts the factor levels of tibble.UTR
tibble.UTR$UTR.end <- factor(tibble.UTR$UTR.end,
                             levels = c("5'UTR","3'UTR"))

# For figure 3.B

# Lists the categories of transcripts with same of different TSS/TES
list.categories <- c("TSS_0_TES_0","TSS.1_TES_0","TSS_0_TES.1")

# Selects only genes with multiple predictions
sub.df.SCV.RNAseq <- df.SCV.RNAseq[df.SCV.RNAseq$Transcripts.Nb > 1,]

# Lists the gene names in sub.df.SCV.RNAseq
list.gene.query <- unique(sub.df.SCV.RNAseq[,1])

# Creates a data frame used for figure 3.B
df.dup.positions <- data.frame(Gene = unique(sub.df.SCV.RNAseq$Gene),
                                   cond.TSS = F,
                                   cond.TES = F)

# For each gene with multiple predictions
for(gene in list.gene.query)
{
  # Gets the TSS and TES of the gene
  gene.TSS <- sub.df.SCV.RNAseq[sub.df.SCV.RNAseq[,1] == gene,"TSS.pos"]
  gene.TES <- sub.df.SCV.RNAseq[sub.df.SCV.RNAseq[,1] == gene,"TES.pos"]
  
  # Tests if there are TSS or TES are common between predictions
  df.dup.positions[df.dup.positions$Gene == gene,"cond.TSS"] <- sum(duplicated(gene.TSS)) > 0
  df.dup.positions[df.dup.positions$Gene == gene,"cond.TES"] <- sum(duplicated(gene.TES)) > 0 

  
}

# Creates a data frame used for figure 3B
df.categories <- data.frame(category = list.categories,
                         cond.TSS =c (F,T,F),
                         cond.TES = c(F,F,T),
                         size = 0)

# Loops for each category in df.categories
#for(category in df.categories$category)
#  {
    # Computes the number of genes associated with the category
#    Nb.genes <-  sum(df.dup.positions.all$cond.TSS == df.categories[df.categories$category == category,2] &
#                     df.dup.positions.all$cond.TES == df.categories[df.categories$category == category,3])
    
    
    # Stocks in df.categories the number of gene computed.
#    df.categories$size[df.categories$category == category] <- Nb.genes 
    
#  }
  
# Sorts the factor levels of df.categories
df.categories$category <- factor(df.categories$category,
                                 levels = list.categories[c(1,3,2)])



# Creates a data.frame used for figure 3.C
df.genes.venn <- list(A = unique(df.SCV.RNAseq[df.SCV.RNAseq$Dataset == "A",1]),
                      B = unique(df.SCV.RNAseq[df.SCV.RNAseq$Dataset == "B",1]),
                      C = unique(df.SCV.RNAseq[df.SCV.RNAseq$Dataset == "C",1]))


#PIERRE
#Damien only numbered the genes in his output file, the easiest way I found is to use a dataframe I made to get 5'UTR coordinates which as both the numbering and gene names
#not the best way to do it but the fastest way at that time
#add gene name in df.SCV.RNAseq
coord_5UTR2 <- coord_5UTR[,c(1,5)]
colnames(coord_5UTR2) <- c("Gene","gene.name")
coord_5UTR2$Gene <- as.character(coord_5UTR2$Gene)
coord_5UTR2 <- distinct(coord_5UTR2)
df.SCV.RNAseq <- left_join(df.SCV.RNAseq,coord_5UTR2)



# Writes the SCV transcripts dataframe in a file
write.table(x = df.SCV.RNAseq,
            file = "./06-SCT_results/UTR_annotation/tables/SCV_mRNA_PODANS.tab",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

