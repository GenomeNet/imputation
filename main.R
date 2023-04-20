library(deepG)

sequence <- microseq::readFasta("example.fasta")$Sequence

#for testing a sequence without Ns, this mask 5000 nucleotides N:
#indices <- sample(1:nchar(sequence), 5000)
#for (i in indices) {
#  substring(sequence, i, 20000000L) <- "N"
#}

model <- keras::load_model_hdf5("/content/model_imputation_maxlen100.hdf5", compile = FALSE)

len <- nchar(sequence)

maxlen <- 100
pred <- predict_model(vocabulary=c("a", "c", "g", "t", "n"),
                      output_format = "one_seq",
                      model = model1,
                      layer_name = "flatten_4",
                      path_input = "example.fasta",
                      sequence = sequence,
                      round_digits = 4,
                      filename = NULL,
                      step = maxlen,
                      batch_size = 512,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard",
                      mode = "label",
                      format = "fasta")

#this is for the last prediction, for example if the sequence is 2770 nt long, the first 2700 is classified using the function above
#the last 70 is using this.
if (len>maxlen){
  pred1 <- predict_model(vocabulary=c("a", "c", "g", "t", "n"),
                      output_format = "one_seq",
                      model = model1,
                      layer_name = "flatten_4",
                      path_input = "example.fasta",
                      sequence = substr(sequence, len-maxlen+1, len),
                      round_digits = 4,
                      filename = NULL,
                      step = maxlen,
                      batch_size = 4,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard",
                      mode = "label",
                      format = "fasta")
} else{
  k <- (maxlen-len)*4
}
num_of_pred <- len %/% maxlen
last_pred <- num_of_pred * maxlen 
last_piece <- len %% maxlen
last_addition <- (maxlen - last_piece)*4


n_positions <- which(strsplit(sequence, "")[[1]] == "N")
i=1 

if (len>maxlen){
  while (i<=length(n_positions)){
    j <- n_positions[i] 
    j_rem_maxlen <- ((j-1) %% maxlen) +1
    if  (j<=last_pred){
      a = which.max(pred$state[((j-1)%/%maxlen)+1,(j_rem_maxlen*4-3):(j_rem_maxlen*4)])
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
      a = which.max(pred1$state[1,(last_addition+j_rem_maxlen*4-3):(last_addition+j_rem_maxlen*4)])
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
} else {
  while (i<=length(n_positions)){
    j <- n_positions[i] 
    j_rem_maxlen <- ((j-1) %% maxlen) +1
    a = which.max(pred$state[1,(k+j_rem_maxlen*4-3):(k+j_rem_maxlen*4)])
    if (a==1){
      substring(sequence, j) <- "A"
    }else if (a==2){
      substring(sequence, j) <- "C"
    }else if (a==3){
      substring(sequence, j) <- "G"
    }else if (a==4){
      substring(sequence, j) <- "T"
    }
  i <- i+1
  } 
}

#here is to test how many of mismatches there are out of 5000. there should be no Ns in your fasta file to test this
#sequence1 <- microseq::readFasta("/content/example.fasta")$Sequence
#num_mismatches <- sum(strsplit(sequence, "")[[1]] != strsplit(sequence1, "")[[1]])
#cat("Number of mismatches: ", num_mismatches)

seq <- microseq::readFasta("example.fasta")
seq$Sequence <- sequence
writeFasta(seq,"imputed_file.fasta")

