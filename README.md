# ERK-inhibitors
Contains Python code and CellProfiler pipelines.

Python code:
Diff Mutations CancerRX ERKi- code to identify mutations associated with increased or decreased sensitivity to ERK inhibitors.
Elastic Net regression CancerRX ERKi- code to calculate gene expression impact on IC50 values for ERK inhibitors.
Genomic, transcriptomic and IC50 data was downloaded from Genomics of Drug Sensitivity in Cancer https://www.cancerrxgene.org/

CellProfiler pipelines:
H1299 ERK-KTR, SH-SY5Y ERK-KTR, and HCT-116 ERK-KTR- pipelines used to calculate ERK-KTR cytoplasm to nucleus ratios.
Illumination Correction- pipeline used to correct illumination for mClover.

Dependencies:
Python 3.8
CellProfiler 4.0.5

Python libraries:
numpy
pandas
scipy
statsmodels
math
sklearn
os- optional

Other:
PANCANCER_Genetic_feature file used for mutational analysis can be downloaded from https://www.cancerrxgene.org/downloads/genetic_features
