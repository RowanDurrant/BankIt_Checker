library(Biostrings)

Percent_Ns = function(filename){
  x = readDNAStringSet(filename)

  for(i in 1:length(x)){
    freqs = alphabetFrequency(x[i], baseOnly = F, as.prob = F)
    N_percent = (sum(freqs[colnames(freqs) %in% c("n", "N")])/sum(freqs)) * 100
      message("\n", names(x[i]), " contains ", N_percent,"% Ns.")

  }
  message("\nEnd of warnings.")
  
}
