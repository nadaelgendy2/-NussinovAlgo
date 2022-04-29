import numpy
import numpy as np

Pairsdict = {'G,C': 1, 'A,U': 1, 'C,G': 1, 'U,A': 1, 'G,U': 1, 'U,G': 1}


def ZeroMatrix(seq):
    matrix = np.zeros((len(seq), len(seq)), int)
    np.fill_diagonal(matrix, 0)

    return matrix


def Check(seq):
    if seq in Pairsdict:
        return Pairsdict[seq]
    else:
        return -1


def MatrixFilling(matrix, RnaSeq):
    DirectionList = {}

    for dindex in range(1, len(RnaSeq)):
        for row in range(len(RnaSeq) - dindex):
            column = row + dindex

            if column - row >= 0:
                Leftindex = matrix[row][column - 1]
                belowindex = matrix[row + 1][column]
                seq = RnaSeq[row] + ',' + RnaSeq[column]
                Diagnoalindex = matrix[row + 1][column - 1] + Check(seq)

                bir = []
                birmax = 0

                for bindex in range(row, column):
                    bir.append(matrix[row][bindex] + matrix[bindex + 1][column])

                for bindex in range(row, column):
                    birmax = max(bir)

                matrix[row][column] = max(Leftindex, belowindex, Diagnoalindex, birmax)
                if Leftindex == matrix[row][column]:
                    DirectionList[str(row), str(column)] = 'Left'

                elif belowindex == matrix[row][column]:
                    DirectionList[str(row), str(column)] = 'Below'

                elif Diagnoalindex == matrix[row][column]:
                    DirectionList[str(row), str(column)] = 'Diagonal'

                else:
                    DirectionList[str(row), str(column)] = 'Burification'
            else:
                matrix[row][column] = -1

    tracebackList = []
    outputDots = []

    tracebackList = traceback(RnaSeq, matrix, tracebackList, 0, len(RnaSeq) - 1, DirectionList)

    for element in range(len(RnaSeq)):
        outputDots.append('.')
    for element in tracebackList:
        outputDots[min(element)] = "("
        outputDots[max(element)] = ")"

    Output = "".join(outputDots)
    print(Output)

    return matrix


def traceback(RnaSeq, matrix, TracebackList, i, j, DirectionList):
    if i < j:
        if DirectionList[str(i), str(j)] == 'Below':

            traceback(RnaSeq, matrix, TracebackList, i + 1, j, DirectionList)

        elif DirectionList[str(i), str(j)] == 'Left':

            traceback(RnaSeq, matrix, TracebackList, i, j - 1, DirectionList)

        elif DirectionList[str(i), str(j)] == 'Diagonal':

            TracebackList.append((i, j))

            traceback(RnaSeq, matrix, TracebackList, i + 1, j - 1, DirectionList)

        elif DirectionList[str(i), str(j)] == 'Burification':
            maximumK = i + 1

            while maximumK < j:

                if matrix[i][j] == matrix[i][maximumK] + matrix[maximumK + 1][j]:
                    traceback(RnaSeq, matrix, TracebackList, i, maximumK, DirectionList)

                    traceback(RnaSeq, matrix, TracebackList, maximumK + 1, j, DirectionList)

                    break
                maximumK = maximumK + 1

    elif i > j:
        return TracebackList

    return TracebackList


# "GGGAAAUCC"
# "CGGACCCAGACUUUC"

rnaseq = input("Enter RNA Sequence : ").upper()
matrix = ZeroMatrix(rnaseq)
matrix = MatrixFilling(matrix, rnaseq)
print(matrix)
