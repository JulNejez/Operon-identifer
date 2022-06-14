# import library
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics
from tabulate import tabulate
from pandas import DataFrame

# load data
gene_expression = pd.read_excel("E.coli.xlsx")
operons = pd.read_excel("E_coli_OM.xlsx")

# Extraction names of gene
names = gene_expression["#name"]
gene_expression.pop("#name")

# Save  ID gene
id_gene = operons["IdGene"]

GE = np.zeros((4481,18))
for i in range(len(names)):
    gene = names[i]

    j = 0
    while j < len(id_gene):
        if id_gene[j] == gene:
            break
        else:
            j = j+1

    if j < len(id_gene):
        GE[j, :] = data_coli.loc[i]

# Save to excel
excel = pd.DataFrame(GE)
excel.to_excel("GE_OM_Ecoli.xlsx")