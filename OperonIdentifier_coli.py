# Import packages
import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
import statistics
from tabulate import tabulate
from pandas import DataFrame

# Function for operon prediction for E. coli BW25113
def OperonIdentifier(operon_prediction, gene_expression, method):

    # Load data and operons
    data = pd.read_excel(gene_expression)
    operons = pd.read_excel(operon_prediction)

    #Column indexing, must be modified according to the files used
    operons = operons['OM']
    gene_id = data['Gene_id']
    data.pop('Gene_id')

    ## Method 1 - Refinement of operon structures
    if method == 1:

        # Output file naming
        name = "operon_prediction_I.xlsx"

        mean_matrix_save = []

        # Initiating a list where newly predicted operons will be stored
        operons_new = []
        limit_all = 0.75
        i = 0

        # While cyclus, which gradually passes through all rows of the dataset ("length" of data only the same as "length" of the predicted sheet with operons)
        while i <= (len(operons ) -1):

            # Initiation of auxiliary variables and variables where the currently analyzed operons will be stored
            provisional_operon = []
            start = i

            if i == (len(operons ) -1):
                operons_new.append(operons_new[-1] + 1)
                i = i+ 1

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

                if mean_matrix > 0.46:
                    limit = mean_matrix * 0.9
                else:
                    limit = limit_all

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
                        mean_row = statistics.mean(correlation[j])
                        if j == 305:
                            print(mean_row)

                        if mean_row >= limit:
                            operons_new.append(operons_new[-1])
                        else:
                            operons_new.append(operons_new[-1] + 1)
                            start = start - 1

                            # Update of the provisional operon and recalculation of correlation
                            provisional_operon = data.loc[start:i]
                            provisional_operon = pd.DataFrame(provisional_operon).T
                            correlation = provisional_operon.corr()
                            matrix_size = len(correlation)
                            mat = np.matrix(correlation)
                            mean_matrix = mat.mean()

                            if mean_matrix < 0.46:
                                limit = limit_all

                ## I'll see if the following gene belongs to the operon
                # I have to increase the index by one, because I already have a provisional operon
                if i <= (len(operons) - 1):
                    i = i + 1
                else:
                    break

                mean_matrix_save.append(mean_matrix)
                mean_last_row = 1

                while (mean_last_row >= limit) and (i <= (len(operons) - 1)):
                    provisional_operon = (data.loc[start:i])
                    provisional_operon = pd.DataFrame(provisional_operon).T
                    correlation = provisional_operon.corr()
                    last_row_corr = correlation[i]
                    mean_last_row = statistics.mean(last_row_corr)

                    if (mean_matrix * 1.1 < mean_last_row):
                        mean_last_row = 0
                    if mean_last_row >= limit:
                        operons_new.append(operons_new[-1])
                        i += 1
                    else:
                        i = i

                if i > (len(operons) - 1):
                    break

    ## Method 2 - Prediction of operon structures
    elif method == 2:

        # Output file naming
        name = "operon_prediction_II.xlsx"

        mean_matrix_save = []

        # Iniciace novýc proměnných
        operons_new = []
        limit_all = 0.70
        length = len(data)
        i = 0

        while (i < (length)):
            start = i

            # Take first and second rows
            provisional_operon = data.loc[i:i + 1]
            provisional_operon = pd.DataFrame(provisional_operon).T
            correlation = provisional_operon.corr()
            mat = np.matrix(correlation)
            mean_matrix = mat.mean()

            if mean_matrix > limit_all:
                if mean_matrix > 0.75:
                    limit = mean_matrix * 0.9
                else:
                    limit = 0.55

                # If operons_new is empty
                if operons_new == []:
                    operons_new.append(1)
                    operons_new.append(1)
                    i = i + 2

                else:
                    operons_new.append(operons_new[-1] + 1)
                    operons_new.append(operons_new[-1])
                    i = i + 2

                mean_last_row = 1
                mean_matrix_save.append(mean_matrix)

                while (mean_last_row >= limit) and (i < length):
                    provisional_operon = data.loc[start:i]
                    provisional_operon = pd.DataFrame(provisional_operon).T
                    correlation = provisional_operon.corr()
                    mat = np.matrix(correlation)
                    mean_matrix = mat.mean()
                    last_row = correlation[i]
                    mean_last_row = statistics.mean(last_row)
                    mean_matrix_save.append(mean_last_row)

                    if ((mean_matrix_save[-2] * 1.1) < mean_last_row):
                        mean_last_row = 0

                    if mean_last_row >= limit:
                        operons_new.append(operons_new[-1])
                        i = i + 1

                    elif (i == length - 1) and mean_last_row < limit:
                        operons_new.append(operons_new[-1] + 1)
                        i = i + 1

            else:

                if operons_new == []:
                    operons_new.append(1)
                    i = i + 1
                elif (i) == (length - 2):
                    operons_new.append(operons_new[-1] + 1)
                    operons_new.append(operons_new[-1] + 1)
                else:
                    operons_new.append(operons_new[-1] + 1)
                    i = i + 1

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
# Import packages
if __name__ == "__main__":
    # Method 1
    OperonIdentifier("operon_prediction_OM_coli.xlsx", "gene_expression_coli.xlsx", 1)

    # Method 2
    #OperonIdentifier("operon_prediction_OM_coli.xlsx", "gene_expression_coli.xlsx", 2)
