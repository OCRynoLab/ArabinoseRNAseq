# ArabinoseRNAseq
RNA-seq workflow for 2025 manuscript on L-arabinose impact on E. coli gene expression

Ryno Lab 2025

Authors: Lisa Ryno, Jason Kuchtey, Isabella Moppel

Process:
1. Download all files from GitHub
2. Make sure the following programs are installed and up-to-date on your computer: gzip, sickle, salmon, R, RStudio (some of these are big and might take a long time to install). Note: sickle and salmon MUST be installed in a Linux environment, can be on a Linux VM if necessary – check that you can run them from the command line with “sickle -h” and “salmon -h”.
3. IN A LINUX ENVIRONMENT: If the paired read files are gzipped, unzip them with “gzip -d <filename>”
4. IN A LINUX ENVIRONMENT: Run cleanQuant.sh on each set of paired reads: move eColi_transcriptome_fasta.fa and cleanQuant.sh into the parent directory containing the paired reads, and run "bash cleanQuant.sh" from the command line (from that directory). Will take 15-30mins for each set of paired reads. A mapping rate of 60-70% is pretty common and perfectly okay. The output that you care about is the directory named “quant” that will appear. After running the “bash cleanQuant.sh” from the terminal, you can go ahead and close that terminal – all further analysis will be done from the GUI and R.
5. Open either R workflow in RStudio. Set the working directory, and make sure that BOTH EColi_k12.gff3.gz and EColi_k12.gff3 are in this working directory. Run the code chunk by chunk (by selecting and *control* *enter*). PLEASE read the comments above each chunk of code, sometimes there are things you'll have to do before you run it.
Analyses are optional, take them or leave them. Read the comments!
