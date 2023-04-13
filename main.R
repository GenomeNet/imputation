Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3) # suppress TF messages
suppressPackageStartupMessages(suppressWarnings(library(deepG)))
library(magrittr)
a = Sys.time()
# Load file from stdin and save to disk
message("Loading file")
file.stdin <- file("stdin")
fileConn <- file("file.fasta")
writeLines(readLines(file.stdin), fileConn)
close(fileConn)
close(file.stdin)

args = commandArgs(trailingOnly=TRUE)
model <- keras::load_model_hdf5("VirusNetGenus.hdf5", compile = FALSE)
genus_labels <- readRDS("labels.rds")
message("Processing file")

pred <- predict_model(vocabulary=c("A", "C", "G", "T", "N"),output_format = "one_seq",
                      model = model1,
                      layer_name = "flatten", 
                      sequence = NULL,
                      path_input = "~/file.fasta",
                      round_digits = 4,
                      filename = NULL,
                      step = 100,
                      batch_size = 2,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard", 
                      mode = "label",
                      format = "fasta")

fasta_file <- "~/file.fasta"

# Read in the fasta file as a character vector
fasta_lines <- readLines(fasta_file)

# Remove any empty lines
fasta_lines <- fasta_lines[fasta_lines != ""]

# Initialize empty variables for the header and sequence
header <- ""
sequence <- ""

# Initialize an empty list to store the fasta entries
fasta_entries <- list()

# Loop through each line in the fasta file
for (line in fasta_lines) {
  
  # If the line is a header line, save the current sequence and reset the variables
  if (substr(line, 1, 1) == ">") {
    if (header != "") {
      fasta_entries[[header]] <- sequence
      sequence <- ""
    }
    header <- line
  }
  
  # If the line is a sequence line, append it to the current sequence
  else {
    sequence <- paste0(sequence, line)
  }
  
}

# Add the last sequence to the fasta entries list
fasta_entries[[header]] <- sequence

# Print the fasta entries
for (header in names(fasta_entries)) {
  sequence <- fasta_entries[[header]]
  #cat(header, "\n", sequence, "\n")
  n_positions <- which(strsplit(sequence, "")[[1]] == "N")
  if (length(n_positions)==0){
    cat(sequence,"\n")}
  i=1 
  while (i<=length(n_positions)){
    j <- n_positions[i] 
    if  (j<1001){
      a = which.max(pred$state[1,j*4-3:j*4])
      }
      if (a==1){
        sequence <- paste0(substr(sequence, 1, j-1), "A", substr(sequence, j+1, nchar(sequence)))}
      else if (a==2){
        sequence <- paste0(substr(sequence, 1, j-1), "C", substr(sequence, j+1, nchar(sequence)))}
      else if (a==3){
        sequence <- paste0(substr(sequence, 1, j-1), "G", substr(sequence, j+1, nchar(sequence)))}
      else if (a==4){
        sequence <- paste0(substr(sequence, 1, j-1), "T", substr(sequence, j+1, nchar(sequence)))}
    i <- i+1
  }
}



write.csv(data.frame(sequence), file = stdout(), row.names = FALSE)

b = Sys.time()
message(paste0("Prediction took ", round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 3), " seconds"))

agg <- colMeans(df)
agg_o <- agg[order(agg, decreasing = T)]
for (i in 1:5){
   message(paste0(names(agg_o[i]), ": " , round(agg_o[i]* 100, digits = 1), "%"))
}
