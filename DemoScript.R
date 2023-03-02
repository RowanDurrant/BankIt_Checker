filename = "example.fasta"

source("BankItCheck.R")
BankItCheck(filename) #full diagnostics

source("Percent_Ns.R")
Percent_Ns(filename) #prints the % Ns for each sequence
