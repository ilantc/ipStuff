import numpy as np
import random
import time
import math
import sys

def unique_permutations(seq):
    """
    Yield only unique permutations of seq in an efficient way.

    A python implementation of Knuth's "Algorithm L", also known from the 
    std::next_permutation function of C++, and as the permutation algorithm 
    of Narayana Pandita.
    
    code taken from:
    http://stackoverflow.com/questions/12836385/how
    -can-i-interleave-or-create-unique-permutations
    -of-two-stings-without-recurs/12837695#12837695
    """

    # Precalculate the indices we'll be iterating over for speed
    i_indices = range(len(seq) - 1, -1, -1)
    k_indices = i_indices[1:]

    # The algorithm specifies to start with a sorted version
    seq = sorted(seq)

    while True:
        yield seq

        # Working backwards from the last-but-one index,           k
        # we find the index of the first decrease in value.  0 0 1 0 1 1 1 0
        for k in k_indices:
            if seq[k] < seq[k + 1]:
                break
        else:
            # Introducing the slightly unknown python for-else syntax:
            # else is executed only if the break statement was never reached.
            # If this is the case, seq is weakly decreasing, and we're done.
            return

        # Get item from sequence only once, for speed
        k_val = seq[k]

        # Working backwards starting with the last item,           k     i
        # find the first one greater than the one at k       0 0 1 0 1 1 1 0
        for i in i_indices:
            if k_val < seq[i]:
                break

        # Swap them in the most efficient way               k     i
        (seq[k], seq[i]) = (seq[i], seq[k])         # 0 0 1 1 1 1 0 0       

        # Reverse the part after but not                           k
        # including k, also efficiently.                     0 0 1 1 0 0 1 1
        seq[k + 1:] = seq[-1:k:-1]

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

def createMatrixFromList(vals,nRows,nCols):
    rows = []
    valIter = 0
    for _ in range(0,nRows):
        row = []
        for _ in range(0,nCols):
            row.append(vals[valIter])
            valIter += 1
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

if __name__ == '__main__':
    
    numbers = [1,1,1,1,8,8,16,2,2,2,2,3,3,15,3,3]
    numbers = [random.random() * 10 for _ in range(16)]
    numbers[0] = numbers[1] = numbers[2] = numbers[3]
    numbers[4] = numbers[5] = numbers[6] = numbers[7]
    numbers[8] = numbers[9]
    numbers[10] = numbers[11]
    numbers[12] = numbers[13]
    numbers[14] = numbers[15] 
#     numbers[4] = numbers[3] = numbers[2] = numbers[1] = numbers[0]
#     numbers[10] = numbers[9] = numbers[8] = numbers[7] = numbers[6] = numbers[5] 
    
    vals = set(numbers)
    
    countDenominator = 1
    for val in vals:
        countDenominator *= math.factorial(numbers.count(val))
    totalIterNum = math.factorial(len(numbers)) / countDenominator 
    
    permsGenerator = unique_permutations(numbers)
    
    iterNum = 0
    startTime = time.clock()
    while True:
        try:
            currPerm = permsGenerator.next()
            mat = createMatrixFromList(currPerm,4,4)
            if (checkForSubMudularity(mat)):
                print "found sub modular perm, iter num", iterNum, "perm is:\n", currPerm
                break
        except StopIteration:
            print "finished all perms!\n none found for numbers:\n", numbers
            break
        iterNum += 1
        if iterNum % 10000 == 0:
            totalTime = time.clock() - startTime
            iterTime = totalTime / iterNum
            print "iter", iterNum, "out of", totalIterNum, "\naverage iter time =", iterTime, \
                "elapsed time =", totalTime, \
                "estimated time to finish =", iterTime * (totalIterNum - iterNum)
         
        