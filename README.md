# BankIt_Checker
R function that checks if .fasta files are suitable for GenBank submission, based on the errors and warnings I've found so far using BankIt. The checks are non-exhaustive, and sequence length limits may vary between organsisms and submssion types. If you are aware of any rules that GenBank submission has that this script doesn't currently check for, feel free to open an issue!

BTW: this works best if you're using a *non-aligned* fasta file. Weird alignments can make your sequences appear higher quality than they actually are!

02 March 2023 Update: added a separate function just to display the %Ns for each sequence.
