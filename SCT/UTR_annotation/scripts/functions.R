# Functions called by the scripts used in this study

# Calls for ggplot library
library(ggplot2)
library(ggsci)
library(ggvenn)
library(UpSetR)

# Deletes patterns before and after the characters chain of interest
extract_names <- function(x, pattern_1, pattern_2, replacement)
  
{
  tmp <- sub(pattern = pattern_1,
             x = x,
             replacement = replacement)
  
  tmp <- sub(pattern = pattern_2,
             x = tmp,
             replacement = replacement)
  
  return(tmp)

}

# Extracts values from a collapsed string 
extract_values <- function(data,split_character = ";")
  
{
  extracted_values <- unlist(strsplit(as.character(data),split = split_character))

  return(extracted_values)
}


# Extracts the number of values from a collapsed string
extract_length <- function(data,split_character = ";")
  
{
  extracted_length <- length(unlist(strsplit(as.character(data),split = split_character)))
  
  return(extracted_length)
}



# Creates a custom data frame according to the respective size of variables
create_custom_dataframe <- function(variable_1 = 1,
                                    variable_2 = 1,
                                    variable_3 = 1,
                                    col_names = NULL)
{
  
  size_var_1 <- length(variable_1)
  
  size_var_2 <- length(variable_2)
  
  size_var_3 <- length(variable_3)
  
  
  df_custom <- data.frame(V1 = rep(variable_1,size_var_2*size_var_3),
                          V2 = rep(rep(variable_2,rep(size_var_1,size_var_2)),
                                   size_var_3),
                          V3 = rep(variable_3,rep(size_var_1*size_var_2,size_var_3)))
  
  colnames(df_custom) <- col_names
  
  return(df_custom)
}

# Calculates the average from compiled values
mean_values <- function(data,split_character = ";")
  
{mean(as.numeric(unlist(strsplit(as.character(data),split = split_character))))}


# Calculates the square deviation from compiled values
sd_values <- function(data,split_character = ";")
  
{sd(as.numeric(unlist(strsplit(as.character(data),split = split_character))))}

# Computes the max value of the given compiled values
max_values <- function(data,split_character = ";",n_col = 4)
  
{
  # Extracts each value in compiled data and calculates the maximum value within
  v_max <- max(as.numeric(unlist(strsplit(data[4],split = split_character))))
  
  return(v_max)
  
}

# Collapses a matrix 
collapse_matrix <- function(matrix.ini)
  
{
  matrix.final <- NULL
  
  for(i in 1:length(matrix.ini[1,]))
    
  {
    matrix.tmp <- matrix(c(rownames(matrix.ini),
                           as.numeric(matrix.ini[,i]),
                           rep(colnames(matrix.ini)[i],times = length(matrix.ini[,i]))),
                         ncol = 3)
    
    matrix.final <- rbind(matrix.final,
                      matrix.tmp)
  }
  
  return(matrix.final)
}


# Compiles data for each respective junction in one row
compile_junctions <- function(identifier,matrix_junctions)
{
  
  # Extracts the information at the row corresponding to the junction identifier
  temp_matrix <- matrix_junctions[matrix_junctions[,7] %in% identifier,]
  
  # Extracts in seperate vectors the chromosome, positions, compiled coverage 
  # and compiled samples of the tophat junction
  temp_chromosome <- temp_matrix[1,1]
  temp_start <- temp_matrix[1,2]
  temp_end <- temp_matrix[1,3]
  temp_cov <- paste(temp_matrix[,4],collapse = ";")
  temp_strand <- temp_matrix[1,5]
  temp_samples <- paste(temp_matrix[,6],collapse = ";")
  
  # Builds a matrix with the compiled information of the tophat junction
  matrix_junction <- matrix(data = c(temp_chromosome,
                                     temp_start,
                                     temp_end,
                                     temp_cov,
                                     temp_strand,
                                     temp_samples,
                                     identifier),
                            nrow = 1)
  
  # Returns the matrix which has tophat junction information
  return(matrix_junction)
}


#
theme_custom_A <- function(police = 16,v_angle = 0,legend.position = "right")
  
{
  theme(axis.title.x = element_text(size = police),
        axis.title.y = element_text(size = police),
        axis.title = element_text(size = police + 4),
        axis.text.x = element_text(size = police - 4,angle = v_angle,vjust = 0.75),
        axis.text.y = element_text(size = police - 4),
        legend.title = element_text(size = police - 4),
        legend.text = element_text(size = police - 8),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "darkgrey"),
        legend.position = legend.position)
  
  
  
}

theme_custom_B <- function(police = 16,v_angle = 0)
  
{
  theme(plot.title = element_text(size = police + 8,
                                  hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = police),
        legend.text = element_text(size = police -2),
        panel.background = element_rect(fill = "white"),
        axis.line = element_blank())
  
  
  
}
