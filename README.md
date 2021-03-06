# OperonIdentifier
OperonIdentifier is a function for predicting operon structures. The algorithm is based on calculations of the correlation of gene expression information and is implemented in Pyrhon.

OperonIdentifier both refines already predicted operons with gene expression information (method I), but also predicts operons based on gene expression information only (method II). The input to the function is the gene expression information in xlsx format and the prediction of operon structures using another available operon prediction tool also in xlsx format.

Different parameters need to be set in the OperonIdentifier function for different organisms. Here, 2 code modifications are available for the organisms - Escherichia coli BW25113 and Clostridium beijerinckii NRRL B-598. The E_coli and C_beijerinckii folders contain all the necessary function input files. 
