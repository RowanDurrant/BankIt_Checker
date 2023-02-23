`%!in%` = Negate(`%in%`)
library(Biostrings)

BankItCheck = function(filename){
  x = readDNAStringSet(filename)

  
  if(length(names(x)[duplicated(names(x))]) > 0){
    message("The following sequences have duplicate IDs: \n", paste0(names(x)[duplicated(names(x))], collapse="\n"))
    
    
  }else{message("No duplicate sequence IDs detected")}
  
  
  for(i in 1:length(x)){
    if(nchar(names(x[i])) > 25){
      message("\nSequence ID ", names(x[i]), " is too long: must be 25 characters or less.")
    }
    
    if(grepl("^[a-zA-Z0-9_]*$", x=names(x[i])) == FALSE){
      message("\nSequence ID ", names(x[i]), " contains special characters.\n
              Please check that these are accepted by GenBank (- _ . : * #).")
      
    }
    
    freqs = alphabetFrequency(x[i], baseOnly = F, as.prob = F)
    
    if(sum(freqs) < 200){
      message("\n",names(x[i]), " is too short (fewer than 200 bases).")
    }
    unacceptable = sum(freqs[colnames(freqs) %!in% c("a", "c", "g", "t", "n",
                                                  "A", "C", "G", "T", "N", "-",
                                                  "U", "u", "R", "r", "Y", "y",
                                                  "S", "s", "W", "w", "K", "k",
                                                  "M", "m", "B", "b", "D", "d",
                                                  "H", "h", "V", "v")])
    if(unacceptable > 0){
      message("\n",names(x[i]), " contains these non-IUPAC characters: \n", 
              paste0(colnames(freqs)[which(colnames(freqs) %!in% c("a", "c", "g", "t", "n",
                                                                    "A", "C", "G", "T", "N", "-",
                                                                    "U", "u", "R", "r", "Y", "y",
                                                                    "S", "s", "W", "w", "K", "k",
                                                                    "M", "m", "B", "b", "D", "d",
                                                                    "H", "h", "V", "v") & 
                                                                        freqs > 0) ], collapse="\n"))
      
      
      }
    else if(sum(freqs[colnames(freqs) %!in% 
                   c("a", "c", "g", "t", "n","A", "C", "G", "T", "N", "-")]) > 0){
        message("Ambiguous nucleotides present in sequence '", names(x[i]),"': ", 
                paste0(colnames(freqs)[which(colnames(freqs) %!in% c("a", "c", "g", "t", "n",
                                                                     "A", "C", "G", "T", "N", "-") & 
                                               freqs > 0) ], collapse="\n"), ". Please check if these are allowed.")
      }
    if(sum(freqs[colnames(freqs) %in% c("n", "N")]) > sum(freqs)/2){
      message("\n", names(x[i]), " contains over 50% Ns, and will not be accepted as it is too low quality.")
      
    }
  }
  message("\nEnd of warnings.")
  
}

