import os
import sys
import pandas as pd
import numpy as np
import logging
#import equiv_class_tools
from readBootstraps import getBootstraps

sep = os.path.sep

def computeFromBoot(sampleDirs, args):
    logger = logging.getLogger('shoal')
    bootstrapMeans = []
    bootstrapVars = []
    quant = None
    for sd in sampleDirs: #Adds the values from different samples into one data frame
        sffile = sep.join([sd, "quant.sf"])
        print("Importing file: {};".format(sffile))
        if quant is None:
            quant = pd.read_table(sffile).set_index('Name')
        print("importing bootstraps")
        bootstraps = np.stack(getBootstraps(sd), axis=-1)
        bootstrapMeans.append(bootstraps.mean(axis=1))
        bootstrapVars.append(bootstraps.var(axis=1))
    bootstrapMeans = np.stack(bootstrapMeans, axis=-1)
    bootstrapVars = np.stack(bootstrapVars, axis=-1)
    logger.info(bootstrapMeans)
    logger.info(bootstrapVars)
    means = bootstrapMeans.mean(axis=1)

    #import scipy.stats
    #for sd in sampleDirs: #Adds the values from different samples into one data frame
    #    sffile = sep.join([sd, "quant.sf"])
    #    print("Importing file: {};".format(sffile))
    #    print("importing bootstraps")
    #    bootstraps = np.stack(getBootstraps(sd), axis=-1)
    #    c = []
    #    for bs in xrange(bootstraps.shape[1]):
    #        c.append(np.linalg.norm(bootstraps[:,bs] - means))
    #    print(sorted(c))

    stddev = np.sqrt( np.sum(np.square(means - bootstrapMeans.T) + bootstrapVars.T, axis=0)  / len(sampleDirs) )
    factors = np.array([np.exp(-stddev[i]/(means[i] + 1e-5)) for i in xrange(len(means))])
    names = quant.index
    d = pd.DataFrame({'Name' : names, 'mean' : means, 'stddev' : stddev, 'factor' : factors}).set_index('Name')
    return d


def normFactors(readVals):
    import scipy.stats.mstats
    gmeans = scipy.stats.mstats.gmean(readVals, axis=-1)
    print(gmeans)
    sjs = []
    for samp in xrange(readVals.shape[1]):
        sj = readVals[:,samp] / gmeans
        sj = sj[~np.isnan(sj)]
        sjs.append(np.median(sj))
    sjs = np.array(sjs)
    return sjs

def computeFromPointEst(sampleDirs, args):
    logger = logging.getLogger('shoal')
    readVals = []
    quant = None
    for sd in sampleDirs: #Adds the values from different samples into one data frame
        sffile = sep.join([sd, "quant.sf"])
        print("Importing file: {};".format(sffile))
        if quant is None:
            quant = pd.read_table(sffile).set_index('Name')
            readVals.append(quant['NumReads'].values)
        else:
            quant2 = pd.read_table(sffile)
            quant2.set_index('Name', inplace=True)
            quant += quant2
            readVals.append(quant2['NumReads'].values)
    readVals = np.stack(readVals, axis=-1)

    ##sjs = normFactors(readVals)
    ##logger.info("Normalization factors = {}".format(sjs))
    ##readVals /= sjs
    #logger.info(readVals)
    #logger.info(readVals.shape)
    #m = sklearn.decomposition.NMF(2)
    #m.fit(readVals.T)
    #logger.info(m.components_.shape)
    #logger.info(m.components_)
    #logger.info(m.reconstruction_err_)
    #import scipy.stats
    #import scipy.optimize
    #for i in xrange(20):
    #    ma = 0
    #    mj = -1
    #    c = []
    #    for j in xrange(2):
    #       s = scipy.optimize.lsq_linear(m.components_.T, readVals[:,i])
    #       #print(s)
    #       #c.append(scipy.stats.spearmanr(readVals[:,i], m.components_[j,:]).correlation)
    #       y = s.x[0] * m.components_[0,:] + s.x[1] * m.components_[1,:]
    #       c.append(scipy.stats.spearmanr(readVals[:,i], y).correlation)
    #    logger.info("rhos = {}".format(c[0]))
    means = readVals.mean(axis=1)
    #logger.info("simple prior")
    #c = []
    #for i in xrange(20):
    #    c.append(scipy.stats.spearmanr(readVals[:,i], means).correlation)
    #for e in c:
    #    logger.info("rhos = {}".format(e))
    stddev = np.sqrt(readVals.var(axis=1))
    N = readVals.shape[1]
    cvs = np.array([ (1 + 1/(4*N)) * (stddev[i]/(means[i] + 1e-5)) for i in xrange(len(means))])
    factors = np.exp(-cvs)
    names = quant.index
    d = pd.DataFrame({'Name' : names, 'mean' : means, 'stddev' : stddev, 'factor' : factors}).set_index('Name')
    return d

def createPrior(args):
    '''
    quantification output file and equivalence class file parsing code, and calling vbem module
    some part from RapClust https://github.com/COMBINE-lab/RapClust/blob/master/bin/RapClust
    input:
        sampleDirs: path to the directory of the samples
    output:
        NONE
    side effect:
        writes prior file to outdir
    '''
    # Create a logger object.
    logger = logging.getLogger('shoal')

    # Initialize coloredlogs.
    import coloredlogs
    coloredlogs.install(level='DEBUG')

    sep = os.path.sep
    # A list of all sample directories
    snames = args.samples.strip().split(',')
    sampleDirs = [sep.join([args.basepath, sample]) for sample in snames]
    outdir = args.outdir

    useBootstraps = False if args.noboot else True
    for sd in sampleDirs:
        sfFile = sep.join([sd, "quant.sf"])
        bsFile = sep.join([sd, "aux_info", "bootstrap", "bootstraps.gz"])
        if not os.path.exists(sfFile):
            logger.critical("Could not find quant file {}. Aborting".format(sfFile))
            sys.exit(1)
        if useBootstraps and not os.path.exists(bsFile):
            logger.info("Sample {} did not contain bootstrap information, disabling use of bootstrap variance".format(bsFile))
            useBootstraps = False

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            print("The output directory already exists, and is not a directory!")
            sys.exit(1)
    else:
        # create it
        os.makedirs(outdir)

    logger.info("Reading in the quant files")
#    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampleDirs]
    d = None
    if useBootstraps:
        d = computeFromBoot(sampleDirs, args)
    else:
        d = computeFromPointEst(sampleDirs, args)
    d.to_csv(outdir + '/prior.tsv', sep='\t')


#def VBEMcaller(args):
#    '''
#    quantification output file and equivalence class file parsing code, and calling vbem module
#    some part from RapClust https://github.com/COMBINE-lab/RapClust/blob/master/bin/RapClust
#    input:
#        sampleDirs: path to the directory of the samples
#    output:
#        alphas: estimated counts after VBEM optimization
#    '''
#    # Create a logger object.
#    import logging
#    logger = logging.getLogger('shoal')
#
#    # Initialize coloredlogs.
#    import coloredlogs
#    coloredlogs.install(level='DEBUG')
#
#    weight = float(args.weight)
#    sep = os.path.sep
#
#    # A list of all sample directories
#    sd = sep.join([args.basepath, 'salmonData', 'quant', args.sample])
#    outdir = args.outdir
#
#    if os.path.exists(outdir):
#        if not os.path.isdir(outdir):
#            print("The output directory already exists, and is not a directory!")
#            sys.exit(1)
#    else:
#        # create it
#        os.makedirs(outdir)
#
#    eqfile = sep.join([sd, 'aux_info', 'eq_classes.txt'])
#    quantfile = sep.join([sd, 'quant.sf'])
#
#    eqCollection = eqtools.EquivCollection(True)
#    eqCollection.fromFile(eqfile)
#    logger.info("Imported file: {}; # eq = {}".format(eqfile, len(eqCollection.eqClasses)))
#
#    priorTable = pd.read_table(args.prior).set_index('Name')
#    inputTable = pd.read_table(quantfile).set_index('Name')
#
#    prior = []
#    alphasIn = []
#    tnames = eqCollection.tnames
#    effLens = np.array(inputTable.loc[tnames,'EffectiveLength'].values)
#    flatPrior = 1e-3 * np.array(inputTable.loc[tnames, 'EffectiveLength'].values)
#    alphas = np.array(inputTable.loc[tnames, 'NumReads'].values)
#    factors = np.array(priorTable.loc[tnames, 'factor'].values)
#    means = np.array(priorTable.loc[tnames, 'mean'].values)
#    prior = flatPrior #+ ((1.0 - factors) * (weight * means))
#    alphasIn = alphas
#    logger.info(prior)
#    #for txp in tqdm(eqCollection.tnames):
#    #    flatPrior = 1e-3 * inputTable.loc[txp, 'EffectiveLength']
#    #    factor = priorTable.loc[txp, 'factor']
#    #    mean = priorTable.loc[txp, 'mean']
#    #    alpha = inputTable.loc[txp, 'NumReads']
#    #    prior.append((factor * flatPrior) + (1.0 - factor) * (weight * mean))
#    #    alphasIn.append(alpha)
#
#    print ("All files imported; Starting VBEM")
#    alphasOut = vbem.runVBEM(np.array(alphasIn), effLens, eqCollection.auxDict, eqCollection.eqClasses, np.array(prior))
#    print ("Optimization Complete; Dumping output")
#    vbem.dumpSFfile(sep.join([outdir, sd.split(sep)[-1]+'_{}.sf'.format(weight)]),
#                    eqCollection.tnames,
#                    inputTable['Length'].to_dict(),
#                    inputTable['EffectiveLength'].to_dict(),#txpLength,
#                    alphasOut)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="shoal")

    subparsers = parser.add_subparsers(help="sub-command help")
    parserCreatePrior = subparsers.add_parser('create-prior', help='build prior')
    parserCreatePrior.add_argument('-n', '--noboot', action='store_true', help="Don't use bootstraps even if they are available")
    parserCreatePrior.add_argument('--samples', help="Comma separated list of sample names")
    parserCreatePrior.add_argument('--basepath', help="Path where the yaml file resides")
    parserCreatePrior.add_argument('--outdir', help="Location where the prior should be written")
    parserCreatePrior.set_defaults(func=createPrior)

    #parserEstimate= subparsers.add_parser('estimate', help='estimate abundance')
    #parserEstimate.add_argument('--sample', help="A single sample on which to run the VBEM")
    #parserEstimate.add_argument('--basepath', help="Path where the yaml file resides")
    #parserEstimate.add_argument('--outdir', help="Path to output abundance estimate")
    #parserEstimate.add_argument('--weight', help="weight of the prior")
    #parserEstimate.add_argument('--prior', help="Path to prior file")
    #parserEstimate.set_defaults(func=VBEMcaller)

    args = parser.parse_args()
    args.func(args)
