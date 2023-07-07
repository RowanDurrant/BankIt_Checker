library(Biostrings)

Percent_Ns = function(filename, output.as.df = F){
  x = readDNAStringSet(filename)
  
  if(output.as.df == F){
    for(i in 1:length(x)){
      x[i] = gsub("-","",x[i])
      freqs = alphabetFrequency(x[i], baseOnly = F, as.prob = F)
      N_percent = (sum(freqs[colnames(freqs) %in% c("n", "N")])/sum(freqs)) * 100
      message("\n", names(x[i]), " contains ", N_percent,"% Ns.")

  }
  message("\nEnd of warnings.")
    
  }
  if(output.as.df == T){
    ID = names(x)
    N_percent = rep(NA, length(x))
    for(i in 1:length(x)){
      x[i] = gsub("-","",x[i])
      freqs = alphabetFrequency(x[i], baseOnly = F, as.prob = F)
      N_percent[i] = (sum(freqs[colnames(freqs) %in% c("n", "N")])/sum(freqs)) * 100
      
    }
    df = as.data.frame(cbind(ID,N_percent))
    return(df)
  }
  
}
