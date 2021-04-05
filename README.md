# ERK-inhibitors
Contains Python, R code and CellProfiler pipelines.

Python code:
Diff Mutations CancerRX ERKi- code to identify mutations associated with increased or decreased sensitivity to ERK inhibitors.
Elastic Net regression CancerRX ERKi- code to calculate gene expression impact on IC50 values for ERK inhibitors.
Genomic, transcriptomic and IC50 data was downloaded from Genomics of Drug Sensitivity in Cancer https://www.cancerrxgene.org/

CellProfiler pipelines:
H1299 ERK-KTR, SH-SY5Y ERK-KTR, and HCT-116 ERK-KTR- pipelines used to calculate ERK-KTR cytoplasm to nucleus ratios.
Illumination Correction- pipeline used to correct illumination for mClover.

R:
ERKi Heatmap- code used for unsupervised clustering and heatmap creation.
Example file for H1299 is provided ('Combined H1299 bins.xlsx').

Dependencies:
Python 3.8
CellProfiler 4.0.5
RStudio 1.4.1106

Python libraries:
numpy
pandas
scipy
statsmodels
math
sklearn
os- optional

R libraries:
ComplexHeatmap
openxlsx

Other:
PANCANCER_Genetic_feature file used for mutational analysis can be downloaded from https://www.cancerrxgene.org/downloads/genetic_features
