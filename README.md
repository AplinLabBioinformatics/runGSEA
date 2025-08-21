# runGSEA
Run GSEA via Command Prompt

If you have not already downloaded GSEA, go to the GSEA website (https://www.gsea-msigdb.org/gsea) and download the windows version (GSEA_Win_4.X.X-installer.exe). Currently, runGSEA only works for windows.

Download the complete gene set collection (.gmt files) for the current release. The filename will look like this:
msigdb_v202X.X.Hs_files_to_download_locally.zip

Download the appropriate chip file (i.e. Human_Ensembl_Gene_ID_MSigDB.v20XX.X.Hs.chip) for the dataset and the corresponding msigdb release version. Note that this is not needed for GEO2R data.

It is good to familiarize yourself with the GSEA wiki page on how to run RNASeq (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA).

### To install, run
 
    devtools::install_github("AplinLabBioinformatics/runGSEA")
  
### GEO2R example 1 - RNAseq data
For the GEO2R_example1.Rmd file, you will need to change the paths (i.e. "C:/GSEA_R_Demo/"). 

The GEO2R_example1.Rmd has GEO2R code for GSE181194 (RNAseq dataset). 

To replace it with code for a different dataset:
1. follow the 'Analyze with GEO2R' link on a GEO Accession webpage (with RNAseq data). 
2. assigning samples into groups. 
3. click on the R script tab.
