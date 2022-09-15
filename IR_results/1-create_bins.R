## Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>

###
# Script to divide introns into bins.
###

# Data reading
intron_pos <- read.table("introns_annotation_PODANS_v2016.gff", sep = "\t")

# Choose bin size
bin_size <- 8

# table to store intron
intron_bins <- NULL

for(i in 1:nrow(intron_pos)){
  print(paste(i, "/", nrow(intron_pos)))
  # from ...
  pos1 <- intron_pos[i,4]
  # ... to
  pos2 <- intron_pos[i,5]
  
  # calculate number of bin
  bin_num <- floor((pos2 - pos1)/bin_size)
  
  # calculate position of bin
  for(j in 1:bin_num){
    
    if(j == 1){
      bin_pos1 <- pos1
    }else{
      bin_pos1 <- bin_pos2 + 1
    }
    
    bin_pos2 <- bin_pos1 + bin_size
    temp <- c(intron_pos[i,1], bin_pos1, bin_pos2, "intron.bin",
              intron_pos[i,6],  intron_pos[i,7], intron_pos[i,8],
              paste0("intron_id=", intron_pos[i,9], ".bin", j))    
    
    intron_bins <- rbind(intron_bins, temp)
  }
  
}

row.names(intron_bins) <- intron_bins[,4]

write.table(cbind(intron_bins[,1], rep("r-script", nrow(intron_bins)),
                  intron_bins[,4], intron_bins[,2], intron_bins[,3],
                  intron_bins[,5:8]),
                  file = "bins_annotated_introns_v2.gff", 
            quote = F, row.names = F, col.names = F, sep = "\t")

write.table(cbind(intron_bins[,1], rep("r-script", nrow(intron_bins)),
                  intron_bins[,4], intron_bins[,2], intron_bins[,3],
                  intron_bins[,5:8])[1:20,],
            file = "bins_annotated_introns_v2_20lines.gff", 
            quote = F, row.names = F, col.names = F, sep = "\t")
