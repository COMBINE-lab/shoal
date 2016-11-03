import gzip
import struct
import os
import logging
import logging.handlers
import sys
import json

def getBootstraps(quantDir):
    logging.basicConfig(level=logging.INFO)

    bootstrapFile = os.path.sep.join([quantDir, "aux_info", "bootstrap", "bootstraps.gz"])
    nameFile = os.path.sep.join([quantDir, "aux_info", "bootstrap", "names.tsv.gz"])
    if not os.path.isfile(bootstrapFile):
       logging.error("The required bootstrap file {} doesn't appear to exist".format(bootstrapFile))
       sys.exit(1)
    if not os.path.isfile(nameFile):
       logging.error("The required transcript name file {} doesn't appear to exist".format(nameFile))
       sys.exit(1)

    with gzip.open(nameFile) as nf:
        txpNames = nf.read().strip().split('\t')

    ntxp = len(txpNames)
    logging.info("Expecting bootstrap info for {} transcripts".format(ntxp))

    with open(os.path.sep.join([quantDir, "aux_info", "meta_info.json"])) as fh:
        meta_info = json.load(fh)

    if meta_info['samp_type'] == 'gibbs':
        s = struct.Struct('<' + 'i' * ntxp)
    elif meta_info['samp_type'] == 'bootstrap':
        s = struct.Struct('@' + 'd' * ntxp)
    else:
        logging.error("Unknown sampling method: {}".format(meta_info['samp_type']))
        sys.exit(1)

    numBoot = 0

    listBootstrap = []
    # Now, iterate over the bootstrap samples and write each
    with gzip.open(bootstrapFile) as bf:
        while True:
            try:
                x = s.unpack_from(bf.read(s.size))
                listBootstrap.append(x)
                numBoot += 1
            except:
                logging.info("read all bootstrap values")
                break

    logging.info("read bootstraps successfully.")
    return listBootstrap
