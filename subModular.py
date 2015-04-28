import numpy as np
import random


def checkForSubMudularity(mat):
    (nRows,nCols) = mat.shape
    notSubMod = False
    for row1 in range (0, nRows - 1):
        for row2 in range (row1 + 1, nRows):
            for col1 in range (0, nCols - 1):
                for col2 in range (col1 + 1, nCols):
                    if ( mat.item(row1,col1) + mat.item(row2,col2) \
                          >  mat.item(row1,col2) + mat.item(row2,col1) ):
                        notSubMod = True
                    if notSubMod:
                        break
                if notSubMod:
                    break
            if notSubMod:
                break
        if notSubMod:
            break
    
    return not notSubMod
                
def creatRandMatrix(nRows,nCols):
    
    rows = []
    for _ in range(0,nRows):
        row = []
        for _ in range(0,nCols):
            row.append(random.random() * 10 )
        rows.append(row)
    
    return np.matrix(rows)

def swapVal(mat,r1,c1,r2,c2):
    v1 = mat.item(r1,c1)
    v2 = mat.item(r2,c2)
    mat.itemset(r1,c1,v2)
    mat.itemset(r2,c2,v1)

# def sortMatrix(mat):
#     allVals = []
#     for row in mat.tolist():
#         allVals += row
#     allVals.sort()
#     iter1 = 0
#     iter2 = len(allVals) - 1
#     for 
#     
#     