# MxTOP
Sample codes to reproduce the analysis results and figures.

Requirements
R (at least v4.0.3)
Bioconductor (v3.12)

Getting Started
These codes demonstrate how to reproduce most findings from Skardziute et al. paper. There are three Rmd files in the code directory. prep_data.Rmd is used to prepare main data matrix (gene level). Analysis.Rmd and AnalysisSup.Rmd shows basic principles of data acquisition. Note that example files contain only two samples and data from a single chromosome. Full data sets were deposited to NCBI GEO (GSE185551).

Input data
The main input is the coverage matrix, which represents target coverages matrix. Data is read in as an RDS file, but it is possible to change this to any other common table format. Besides the coverage matrix, there are several other input files describing different genomic features - see example files for the exact format.

Final figures
Codes and data to generate figures from the article are present. There are individual files for each main figure and one file for all supplemental figures. Figures not generated with R are not included.
