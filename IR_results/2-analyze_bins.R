# Gaëlle Lelandais <gaelle.lelandais@universite-paris-sacay.fr>

# All results files are in the same folder
all_files <- list.files("./htseqcount")

combined_results <- NULL
sample_names <- NULL

for(f in all_files){

  print(f)
  # sample name
  sample_names <- c(sample_names, unlist(strsplit(f, "-"))[2])
  
  # Data reading
  res_file <- read.table(paste0("./htseqcount/", f))
  res_file <- res_file[c(-nrow(res_file), 
                         -(nrow(res_file)-1),
                         -(nrow(res_file)-2),
                         -(nrow(res_file)-3),
                         -(nrow(res_file)-4)),]
  
  # Extract value for each bin, for each intron...
  bin_val <- NULL
  
  all_bins <- cbind(matrix(unlist(strsplit(res_file[,1], ".bin")), 
                           ncol = 2, byrow = T), res_file[,2])
  all_introns <- NULL
  
  # List of different introns
  uniq_intron <- unique(all_bins[,1])
  top <- length(uniq_intron)
  count <- 1
  
  for(i in uniq_intron){
    
 #   print(paste(count, "/", top))
    count <- count + 1
    
    # results for bins are merged for each intron
    all_introns <- c(all_introns, 
                     paste(all_bins[which(all_bins[,1] == i),3], collapse = ";"))
  }
  
  # Intron names are used 
  names(all_introns) <- uniq_intron
  
  combined_results <- cbind(combined_results, all_introns)
# end of for()  
}

# Names of samples are used
colnames(combined_results) <- sample_names

# Writing of the results
write.table(combined_results, file = "compiled_htseqcounts.tab",
            quote = F, row.names = T, col.names = T, sep = "\t")
