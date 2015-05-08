# Orthology-with-UBLAST
Workflow for Orthology analysis using UBLAST (faster, more sensitive than BLAST) and Reciprocal Best Hits

* Find __genome-wide orthologs__ with UBLAST (Linux bash script "ublast.sh")
      * the UBLAST algorithm implemented in USEARCH is faster and more sensitive than BLAST 
      * [USEARCH >> UBLAST documentation](http://www.drive5.com/usearch/)
      
* Select the __Best Reciprocal Hits__ from the UBLAST output (R script "RBH_ublast.R")
      
