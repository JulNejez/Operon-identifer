# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics
from tabulate import tabulate
from pandas import DataFrame


def OperonIdentifier(operons, data, method):
    # Load data and operons
    data = pd.read_excel(data)
    operons = pd.read_excel(operons)

    operons = operons['FGENESB']

    ## Method 1
    if method == 1:

        # Initiating a list where newly predicted operons will be stored
        operons_new = []
        limit_all = 0.77
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

                    if operons_new == []:
                        operons_new.append(1)
                        i += 1

                    else:
                        operons_new.append(operons_new[-1] + 1)
                        i += 1

                else:
                    if operons_new == []:
                        operons_new.append(1)
                        operons_new.append(operons_new[-1])
                    else:
                        operons_new.append(operons_new[-1] + 1)
                        operons_new.append(operons_new[-1])

                    i = i + 2
                    provisional_operon = []
                    limit = mean_matrix
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
                            matrix_size = len(correlation)
                            mat = np.matrix(correlation)
                            mean_matrix = mat.mean()
                            # limit = 0.95 * mean_matrix
                            limit = mean_matrix

                            if limit < limit_all:
                                limit = limit_all

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

    ## Method 2
    elif method == 2:

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

    return (operons_new)