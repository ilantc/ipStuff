import numpy as np
import random
import time
import math
import gurobipy as gp

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

def createILP(vals,nRows,nCols):
    n = len(vals)
    model = gp.Model()
    # variables - x[i,j] = 1 iff the ith item in the permutation is vals[j]
    x = {}
    for i in range(n):
        for j in range(n):
            x[i,j] = model.addVar(vtype=gp.GRB.BINARY,     name=('A_%s = vals_%s' % (i,j)))
    model.update()
    for r1 in range(nRows - 1):
        for r2 in range(r1 + 1,nRows):
            for c1 in range(nCols - 1):
                A1 = (r1*nCols) + (c1)
                A3 = (r2*nCols) + (c1) 
                for c2 in range(c1 + 1,nCols):
                    A2 = (r1*nCols) + (c2)
                    A4 = (r2*nCols) + (c2)
                    model.addConstr(gp.quicksum(x[A1,j]*vals[j] for j in range(n)) + \
                                    gp.quicksum(x[A4,j]*vals[j] for j in range(n)),gp.GRB.LESS_EQUAL, \
                                    gp.quicksum(x[A2,j]*vals[j] for j in range(n)) + \
                                    gp.quicksum(x[A3,j]*vals[j] for j in range(n)),'r_%s_c_%s_r_%s_c_%s' % (r1,c1,r2,c2))
    
    for i in range(n):
        model.addConstr(gp.quicksum(x[i,j] for j in range(n)),gp.GRB.EQUAL,1)
        model.addConstr(gp.quicksum(x[j,i] for j in range(n)),gp.GRB.EQUAL,1)
    
    firstColIndices = filter(lambda num: num % nCols == 0, range(n))
    minIndex = vals.index(min(vals))
    maxIndex = vals.index(max(vals))
    
#     # min and max elements in the first col
#     model.addConstr(gp.quicksum(x[i,minIndex] for i in firstColIndices),gp.GRB.GREATER_EQUAL,1)
#     model.addConstr(gp.quicksum(x[i,maxIndex] for i in firstColIndices),gp.GRB.GREATER_EQUAL,1)

#     # min element top left corner, max element bottom left corner
#     model.addConstr(x[n - nCols,maxIndex],gp.GRB.GREATER_EQUAL,1)
#     model.addConstr(x[0,minIndex],gp.GRB.GREATER_EQUAL,1)

#     # every col is non decreasing
#     for col in range(nCols):
#         for row in range(nRows - 1):
#             index1 = (row * nRows) + col
#             index2 = ((row + 1) * nRows) + col
#             model.addConstr(gp.quicksum(x[index1,j] * vals[j] for j in range(n)),gp.GRB.LESS_EQUAL,\
#                             gp.quicksum(x[index2,k] * vals[k] for k in range(n)))
    
    model.setObjective(1,gp.GRB.MAXIMIZE)
    model.setParam( 'OutputFlag', False )
    model.update()
    model.optimize()
    if model.getAttr('status') != gp.GRB.OPTIMAL:
        print "\t", vals
    
    return (model, x)

# def solveAlg1(vals,nRows,nCols):
#     for i in range(nCols):
#         # find the nRow numbers that has 
#     

def vars2matrix(vals,lpvars,nRows,nCols):
    newPermutation = []
    for i in range(nRows*nCols):
        for j in range(nRows*nCols):
            if lpvars[i,j].x > 0:
                newPermutation.append(vals[j])
    return createMatrixFromList(newPermutation,nRows,nCols)

if __name__ == '__main__':
    
    def expDistribution():
        lam = 0.01
        return -(math.log(random.random()))/lam
    
    def uniDistribution():
        maxVal = 10
        minVal = -10
        return (random.random() * (maxVal - minVal)) + minVal
    
    nRows = 4
    nCols = 5
    niter = 50
    dist = expDistribution
    ncorrect = 0
    for _ in range(niter):
#         numbers = [-(math.log(random.random()))/0.5  for _ in range(nRows * nCols)]
        numbers = [dist() for _ in range(nRows * nCols)]
#     numbers[0] = numbers[1] = numbers[2] = numbers[3]
#     numbers[4] = numbers[5] = numbers[6] = numbers[7]
#     numbers[8] = numbers[9]
#     numbers[10] = numbers[11]
#     numbers[12] = numbers[13]
#     numbers[14] = numbers[15] 
#     numbers[4] = numbers[3] = numbers[2] = numbers[1] = numbers[0]
#     numbers[10] = numbers[9] = numbers[8] = numbers[7] = numbers[6] = numbers[5] 
    
        (model,lpVars) = createILP(numbers, nRows,nCols)
        mat = vars2matrix(numbers,lpVars, nRows,nCols)
        numbers.sort()
        print ["%.3f" % t for t in numbers]
        print mat
        isSubModular = checkForSubMudularity(mat)
        print "is sub modular =",isSubModular
        if isSubModular:
            ncorrect += 1
        else:
            print numbers
    print ncorrect, "out of",niter,"were subModular"
    exit()
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
         
        