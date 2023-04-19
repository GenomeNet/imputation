library(deepG)
library(magrittr)
library(ggplot2)
library(keras)

sequence <- microseq::readFasta("example.fasta")$Sequence

#for testing a sequence without Ns, this mask 1000 nucleotides N:
#indices <- sample(1:nchar(sequence), 1000)
#for (i in indices) {
#  substring(sequence, i, 20000000L) <- "N"
#}

model <- keras::load_model_hdf5("model_imputation.hdf5", compile = FALSE)
len <- nchar(sequence)

pred <- predict_model(vocabulary=c("a", "c", "g", "t", "n"),
                      output_format = "one_seq",
                      model = model,
                      layer_name = "flatten",
                      path_input = "example.fasta",
                      sequence = sequence,
                      round_digits = 4,
                      filename = NULL,
                      step = 1000,
                      batch_size = 64,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard",
                      mode = "label",
                      format = "fasta")

#this is for the last prediction, for example if the sequence is 2700 nt long, the first 2000 is classified using the function above
#the last 700 is using this.
if (len>1000){
  pred1 <- predict_model(vocabulary=c("a", "c", "g", "t", "n"),
                      output_format = "one_seq",
                      model = model,
                      layer_name = "flatten",
                      path_input = "example.fasta",
                      sequence = substr(sequence, len-999, len),
                      round_digits = 4,
                      filename = NULL,
                      step = 1000,
                      batch_size = 4,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard",
                      mode = "label",
                      format = "fasta")
  k <- 0
} else{
  k <- (1000-len)*4
}
num_of_pred <- len %/% 1000
last_pred <- num_of_pred * 1000 
last_piece <- len %% 1000
last_addition <- (1000 - last_piece)*4


n_positions <- which(strsplit(sequence, "")[[1]] == "N")
i=1 
while (i<=length(n_positions)){
  j <- n_positions[i] 
  j_rem_1000 <- ((j-1) %% 1000) +1
  if  (j<=last_pred){
    a = which.max(pred$state[((j-1)%/%1000)+1,(k+j_rem_1000*4-3):(k+j_rem_1000*4)])
    if (a==1){
      substring(sequence, j) <- "A"
    }else if (a==2){
      substring(sequence, j) <- "C"
    }else if (a==3){
      substring(sequence, j) <- "G"
    }else if (a==4){
      substring(sequence, j) <- "T"
    } 
  } else {
    a = which.max(pred1$state[1,(last_addition+j_rem_1000*4-3):(last_addition+j_rem_1000*4)])
    if (a==1){
      substring(sequence, j) <- "A"
    }else if (a==2){
      substring(sequence, j) <- "C"
    }else if (a==3){
      substring(sequence, j) <- "G"
    }else if (a==4){
      substring(sequence, j) <- "T"
    }
  } 
  i <- i+1
}

#here is to test how many of mismatches there are out of 1000. there should be no Ns in your fasta file to test this
#sequence1 <- microseq::readFasta("/content/example.fasta")$Sequence
#num_mismatches <- sum(strsplit(sequence, "")[[1]] != strsplit(sequence1, "")[[1]])
#cat("Number of mismatches: ", num_mismatches)

seq <- microseq::readFasta("example.fasta")
seq$Sequence <- sequence
writeFasta(seq,"imputed_file.fasta")

