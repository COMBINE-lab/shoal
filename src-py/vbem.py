from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
from scipy.special import digamma
#from math import exp
import numpy as np

minIter = 100
maxIter = 1000
minAlpha = 1e-8
alphaCheckCutoff = 10**-2       #https://github.com/COMBINE-lab/salmon/blob/94fddc82fc694baee0c12160cc176baa1bd2468c/src/CollapsedEMOptimizer.cpp#L894
relDiffTolerance =  0.01
digammaMin = 10**-10            #Not sure why this value https://github.com/COMBINE-lab/salmon/blob/develop/src/CollapsedEMOptimizer.cpp#L39
minEQClassWeight = 2**-1074     #http://en.cppreference.com/w/cpp/types/numeric_limits/denorm_min

def dumpSFfile(outfile,
               tnames,
               txpLengthDict,
               effLengthDict,
               alphasList):
    '''
    Dumps the values into quant.sf format
    input:
        outfile: name of the file to be dumped.
        tnames: List of txpNames
        txpLengthDict: Dictionary of txp Length
        effLengthDict: Dictionary of effective txp length
        TPMlist: List of TPM values
        alphasList: List of alphas
    '''
    print(len(tnames))
    print(len(alphasList))
    print(len(effLengthDict))
    norm = 0.0
    for index, txp in enumerate(tnames):
        if index < len(alphasList) and index < len(effLengthDict):
            norm += old_div(alphasList[index], effLengthDict[txp])
    norm /= 10**6
    with open(outfile, 'w') as oFile:
        oFile.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for index, txp in enumerate(tnames):
            if index < len(alphasList) and index < len(effLengthDict):
                tpm = old_div((old_div(alphasList[index], effLengthDict[txp])), norm)
                oFile.write("{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\n".format(txp, txpLengthDict[txp], effLengthDict[txp], tpm, alphasList[index]))

def VBEMUpdate(alphaIn,
               txpGroupCombinedWeights,
               txpGroupCounts,
               priorAlphas):
    '''
    input:
        txpGroupLabels: Represents eqClass i.e. tuples of txp Ids
        alphas: estimated counts after runing online phase
        txpGroupCombinedWeights: Dictionary with eqClass label to w_i^j weight for each txp in eqClass.
        txpGroupCounts: Dictionary with eqClass label to d^j counts of reads for eqClass
        priorAlphas: prior values for estimated counts
    output:
        alphas: estimated counts after one round of VBEM optimization.
    Variable names are consistent with https://github.com/COMBINE-lab/salmon/blob/develop/src/CollapsedEMOptimizer.cpp#L166
    for simplicity.
    '''

    alphaSum = np.sum(alphaIn)
    logNorm = digamma(alphaSum)

    expTheta = np.exp(digamma(alphaIn) - logNorm)
    alphasOut = np.copy(priorAlphas)
    expTheta[np.isnan(expTheta)] = 0
    expTheta[np.isinf(expTheta)] = 0

    for txps, count in txpGroupCounts.items():
        denom = 0.0
        auxs = txpGroupCombinedWeights[txps]

        if len(txps) > 1:
            for index, txp in enumerate(txps):
                aux = auxs[index]
                if expTheta[txp] > 0.0:
                    denom += expTheta[txp] * aux
            if denom > minEQClassWeight:
                invDenom = old_div(count, denom)
                for index, txp in enumerate(txps):
                    aux = auxs[index]
                    if expTheta[txp] > 0.0:
                        alphasOut[txp] += expTheta[txp] * aux * invDenom
        else:
            txp = txps[0]
            alphasOut[txp] += count

    return alphasOut

def runVBEM(alphas,
            effLens,
            txpGroupCombinedWeights,
            txpGroupCounts,
            priorAlphas):
    '''
    Applying variational Bayesian EM method on rich equivalence classes.
        input:
            txpGroupLabels: Represents eqClass i.e. tuples of txp Ids
            alphas: estimated counts after runing online phase
            txpGroupCombinedWeights: Dictionary with eqClass label to w_i^j weight for each txp in eqClass.
            txpGroupCounts: Dictionary with eqClass label to d^j counts of reads for eqClass
            priorAlphas: prior values for estimated counts

        output:
            alphas: estimated counts after VBEM optimization.
    '''
    for k, v in txpGroupCombinedWeights.items():
        count = txpGroupCounts[k]
        newWeights = []
        w = 0
        for t in k:
            w += old_div(count, effLens[t])
        norm = old_div(1.0, w)
        for t in k:
            newWeights.append((old_div(count, effLens[t])) * norm)
        txpGroupCombinedWeights[k] = np.array(newWeights)
 
    converged = False
    itNum = 0
    alphasDPrime = np.zeros(len(priorAlphas))
    while itNum < minIter or (itNum < maxIter and not converged):
        itNum += 1
        alphasPrime = VBEMUpdate(alphas, txpGroupCombinedWeights, txpGroupCounts, priorAlphas)
        if np.array_equal(alphasDPrime, alphasPrime):
            break

        converged = True
        maxRelDiff = -np.inf
        for txpId in range(len(alphas)):
            if alphasPrime[txpId] > alphaCheckCutoff:
                relDiff = old_div(abs(alphas[txpId] - alphasPrime[txpId]), alphasPrime[txpId])
                maxRelDiff = max(relDiff, maxRelDiff)
                if relDiff > relDiffTolerance:
                    converged = False
                    break
        alphasDPrime = np.copy(alphas)
        alphas = np.copy(alphasPrime)

        if maxRelDiff < relDiffTolerance:
            break
            #alphasPrime[txpId] = 0.0

        if itNum %10 == 0:
            print("iteration = {} | Max rel diff. = {}".format(itNum, maxRelDiff))

    alphas[alphas <= (priorAlphas + minAlpha)] = 0.0
    print ("Total {} iterations done".format(itNum))
    return alphas
