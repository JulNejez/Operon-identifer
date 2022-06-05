# Import packages
import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
import statistics
from tabulate import tabulate
from pandas import DataFrame

# Function for operon prediction for C. beijerinckii NRRL B-598
def OperonIdentifier(operon_prediction, gene_expression, method):
    # Load data and operons
    data = pd.read_excel(gene_expression)
    operons = pd.read_excel(operon_prediction)

    # To be changed according to the tool used
    # Column indexing, must be modified according to the files used

    # Operon-mapper as default operon prediction file
    operons = operons['OM']

    # Output file naming
    name = "operon_prediction_I_OM.xlsx"

    # FGENESB as default operon prediction file
    #operons = operons['FGENESB']

    # Output file naming
    #name = "operon_prediction_I_FGENESB.xlsx"


    gene_id = data['Gene_id']
    data.pop('Gene_id')

    ## Method 1 - Refinement of operon structures
    if method == 1:

        # Initiating a list where newly predicted operons will be stored
        operons_new = []
        limit_all = 0.83
        i = 0

        # While cyclus, which gradually passes through all rows of the dataset ("length" of data only the same as "length" of the predicted sheet with operons)
        while i <= (len(operons) - 1):

            # Initiation of auxiliary variables and variables where the currently analyzed operons will be stored
            provisional_operon = []
            start = i

            if i == (len(operons) - 1):
                operons_new.append(operons_new[-1] + 1)
                i = i + 1

            if i > (len(operons) - 1):
                break

            # If the operons do not equal in the rows below each other, calculate the correlation of 2 consecutive lines to see if they happen to correlate
            if (operons[i] != operons[i + 1]):
                provisional_operon.append(data.loc[i])
                provisional_operon.append(data.loc[i + 1])
                provisional_operon = pd.DataFrame(provisional_operon).T
                correlation = provisional_operon.corr()

                # Calculation of average correlation
                mat = np.matrix(correlation)
                mean_matrix = mat.mean()

                # If the average correlation is less than limit_all, add another operon to the operons_new
                if mean_matrix < limit_all:
                    limit = 0.58

                    if operons_new == []:
                        operons_new.append(1)
                        i += 1

                    else:
                        operons_new.append(operons_new[-1] + 1)
                        i += 1

                else:
                    limit = mean_matrix * 0.9
                    if operons_new == []:
                        operons_new.append(1)
                        operons_new.append(operons_new[-1])
                    else:
                        operons_new.append(operons_new[-1] + 1)
                        operons_new.append(operons_new[-1])

                    i = i + 2
                    provisional_operon = []
                    mean_last_row = 1
                    while (mean_last_row >= limit) and (i <= (len(operons) - 1)):
                        provisional_operon = (data.loc[start:i + 1])
                        provisional_operon = pd.DataFrame(provisional_operon).T
                        correlation = provisional_operon.corr()
                        last_row_corr = correlation[i]
                        mean_last_row = statistics.mean(last_row_corr)

                        if mean_last_row >= limit:
                            operons_new.append(operons_new[-1])
                            i += 1
                        else:
                            i = i
                    if i > (len(operons) - 1):
                        break

            else:
                while (operons[i] == operons[i + 1]) and (i < (len(operons) - 2)):
                    provisional_operon.append(data.loc[i])
                    i += 1

                # Adding the last gene that belongs to the operon to the operon
                provisional_operon.append(data.loc[i])

                ## Correlation in operon
                # First I swap the rows and columns, then I perform the correlation calculation

                provisional_operon = pd.DataFrame(provisional_operon).T
                correlation = provisional_operon.corr()
                matrix_size = len(correlation)

                # Calculation of the average value of correlation in operon, determination of the minimum value of correlation (95% of the average value)
                mat = np.matrix(correlation)
                mean_matrix = mat.mean()
                limit = mean_matrix

                if limit < limit_all:
                    limit = 0.58
                else:
                    limit = limit * 0.9

                # If the list with new operons is empty, the operon number is 1
                if operons_new == []:
                    operons_new.append(1)

                    for j in range(1, matrix_size):
                        mean_row = statistics.mean(correlation[j])

                        if mean_row >= limit:
                            operons_new.append(operons_new[-1])
                        else:
                            operons_new.append(operons_new[-1] + 1)

                # Otherwise, start numbering the operon from operons_new[-1] + 1
                else:
                    operons_new.append(operons_new[-1] + 1)

                    for j in range(start + 1, i + 1):
                        mean_row = statistics.mean(correlation[i])

                        if mean_row >= limit:
                            operons_new.append(operons_new[-1])
                        else:
                            operons_new.append(operons_new[-1] + 1)
                            start = start - 1

                            # Update of the provisional operon and recalculation of correlation
                            provisional_operon = data.loc[start:i]
                            provisional_operon = pd.DataFrame(provisional_operon).T
                            correlation = provisional_operon.corr()
                            mat = np.matrix(correlation)
                            mean_matrix = mat.mean()
                            # limit = 0.95 * mean_matrix
                            limit = mean_matrix * 0.9

                            if limit < limit_all:
                                limit = 0.58

                ## I'll see if the following gene belongs to the operon
                # I have to increase the index by one, because I already have a provisional operon
                if i <= (len(operons) - 1):
                    i = i + 1
                else:
                    break

                mean_last_row = 1

                while (mean_last_row >= limit) and (i <= (len(operons) - 1)):
                    provisional_operon = (data.loc[start:i])
                    provisional_operon = pd.DataFrame(provisional_operon).T
                    correlation = provisional_operon.corr()
                    last_row_corr = correlation[i]
                    mean_last_row = statistics.mean(last_row_corr)

                    if mean_last_row >= limit:
                        operons_new.append(operons_new[-1])
                        i += 1
                    else:
                        i = i

                if i > (len(operons) - 1):
                    break

    ## Method 2 - Prediction of operon structures
    elif method == 2:

        # output file naming
        name = "operon_prediction_II.xlsx"

        # Initiating a worksheet where newly predicted operons will be stored
        operons_new = []
        length = len(data)
        limit_all = 0.85
        i = 0

        while i < length:
            start = i
            provisional_operon = []

            # I take the first two lines, I calculate the correlation, if the correlation is greater than or equal to the limit, then the genes form one operon, if not, the genes form a different operon
            provisional_operon = data.loc[i:(i + 1)]
            provisional_operon = pd.DataFrame(provisional_operon).T
            correlation = provisional_operon.corr()
            mat = np.matrix(correlation)
            mean_matrix = mat.mean()
            limit = mean_matrix * 0.9

            if limit < limit_all:
                limit = limit_all

            if mean_matrix < limit_all:
                if i == 0:
                    operons_new.append(1)
                    i += 1
                else:
                    operons_new.append(operons_new[-1] + 1)
                    i += 1

            else:
                if i == 0:
                    operons_new.append(1)
                    operons_new.append(1)

                else:
                    operons_new.append(operons_new[-1] + 1)
                    operons_new.append(operons_new[-1])

                i += 2

                mean_last_row = 1
                while (mean_last_row >= limit) and (i < length):
                    provisional_operon = data.loc[start:i]
                    provisional_operon = pd.DataFrame(provisional_operon).T
                    correlation = provisional_operon.corr()
                    last_row = correlation[i]
                    mean_last_row = statistics.mean(last_row)

                    if mean_last_row > limit_all:
                        operons_new.append(operons_new[-1])
                        i += 1

    # Add information transcript unit (TU) or operon
    TU_or_operon = []
    for i in range(len(operons_new)):
        if i == 0:
            if operons_new[i] == operons_new[i + 1]:
                TU_or_operon.append("Operon")
            else:
                TU_or_operon.append("TU")
        elif i == (len(operons_new) - 1):
            if operons_new[i] == operons_new[i - 1]:
                TU_or_operon.append("Operon")
            else:
                TU_or_operon.append("TU")
        else:
            if (operons_new[i] == operons_new[i - 1]) or (operons_new[i] == operons_new[i + 1]):
                TU_or_operon.append("Operon")
            else:
                TU_or_operon.append("TU")

    # Save gene_id, operons_new and TU_or_operon in one table and save as xlsx
    final_table = pd.DataFrame({'ID_gene': gene_id,
                                'Operon_prediction': operons_new,
                                'Operon/TU': TU_or_operon})
    final_table.to_excel(name)

    return (final_table)

# Call function
if __name__ == "__main__":
    # Method 1, input Operon-mapper
    OperonIdentifier("operon_prediction_OM.xlsx", "gene_expression.xlsx", 1)

    # Method 1, input FGENESB
    #OperonIdentifier("operon_prediction_FGENESB.xlsx", "gene_expression.xlsx", 1)

    #Method 2
    #OperonIdentifier("operon_prediction_FGENESB.xlsx", "gene_expression.xlsx", 2)
