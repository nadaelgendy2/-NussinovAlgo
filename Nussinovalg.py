import numpy
import numpy as np
#MA7DSH YL3AB F AY HAGA 3ASHAN HAYDRBKO ALLAHUMA BALGHT!
Pairsdict={  'G,C' : 1,  'A,U':1,  'G,U':1, 'C,G' : 1   ,'U,A':1,  'U,G':1}

def ZeroMatrix(seq):

    #print( [[0 for row in range(row)] for column in range(column)])
    matrix=np.zeros((len(seq),len(seq)),int)
    np.fill_diagonal(matrix,0)
   # print(matrix)
    return matrix


def Check(seq):
    if seq in Pairsdict:
        return Pairsdict[seq]
    else:
        return -1


def MatrixFilling(matrix,RnaSeq):
    for index in range(1, len(RnaSeq)):
        for row in range(len(RnaSeq) - index):
            column = row + index
            if column - row >= 0:
                Leftindex = matrix[row][column - 1]
                belowindex = matrix[row + 1][column]
                seq=RnaSeq[row]+','+ RnaSeq[column]
                Diagnoalindex = matrix[row + 1][column - 1] +Check(seq)
                bir=[]
                birmax=0
                for index in range(row,column):
                  bir.append(matrix[row][index]+matrix[index + 1][column])
                for index in range(row,column):
                    birmax=max(bir)

                matrix[row][column] = max(Leftindex,belowindex,Diagnoalindex,birmax)  # max of all

                #saving the indices of the max each time
                leftind={}
                belowind = {}
                diagind = {}
                birind = {}
                if matrix[row][column] == Leftindex:
                    leftind={row,column}

                elif matrix[row][column] == belowindex:
                    belowind = {row, column}

                elif matrix[row][column] == Diagnoalindex:
                    diagind = {row,column}

                else:
                    for k in range(row + 1, column - 1):

                        if matrix[row][column] == birmax:
                            birind = {row,column}

                #print(birind,diagind,belowind,leftind)

            else:
                matrix[row][column] = -1
    return matrix


#"GGGAAAUCC"
#"CGGACCCAGACUUUC"

rnaseq=input("Enter RNA Sequence : ").upper()
matrix=ZeroMatrix(rnaseq)
matrix=MatrixFilling(matrix,rnaseq)
print(matrix)

