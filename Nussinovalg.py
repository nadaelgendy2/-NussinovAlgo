import numpy as np

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
                for t in range(row,column):
                  bir.append(matrix[row][t]+matrix[t + 1][column])
                for t in range(row,column):
                    birmax=max(bir)
                matrix[row][column] = max(Leftindex,belowindex,Diagnoalindex,birmax)  # max of all
            else:
                matrix[row][column] = -1
    return matrix


#"GGGAAAUCC"
#"CGGACCCAGACUUUC"

rnaseq=input("Enter RNA Sequence : ").upper()
matrix=ZeroMatrix(rnaseq)
matrix=MatrixFilling(matrix,rnaseq)
print(matrix)

