# One-pot trimodal mapping of unmethylated, hydroxymethylated and open chromatin sites unveils distinctive 5hmC roles at dynamic chromatin loci  
Kotryna Skardžiūtė, Kotryna Kvederavičiūtė, Inga Pečiulienė, Milda Narmontė, Povilas Gibas, Janina Ličytė, Saulius Klimašauskas, Edita Kriukienė

Sample codes to reproduce the analysis results.

## Requirements
R (at least v4.0.3)
Bioconductor (v3.12)

## Getting Started
These codes demonstrate how to reproduce most findings from Skardziute et al. paper. There are multiple R files in the code directory that correspond to individual figures.

## Input data for regions
For regions creation, the main input is the coverage matrix, which represents target coverages matrix and density matrix. Data is read in as an RDS and bw files, but it is possible to change this to any other common table format. 

## Input data for final figures
Codes and data to generate figures from the article are present. There are individual files for each main figure and one file for all supplemental figures. Figures not generated with R are not included.
