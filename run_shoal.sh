#!/bin/bash
# Using getopt


if [ "$#" -ne 4 ]
then   # Script needs exactly two command-line argument.
  echo "Usage $0 -[options -q <quant directory>, -o <shoal output path>]"
  exit 1
fi

while getopts ":q:o:" opt; do
  case $opt in
    q)
      echo "-q<Input quant files Path>: $OPTARG"
	  BASEPATH=$OPTARG
      ;;
	o)
	  echo "-o<Output shoal files Path>: $OPTARG"
	  OUT=$OPTARG
	  ;;
    \?)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
      exit 1
      ;;
  esac
done

###############################################
#make c++ code
cd src-cpp
make
cd ..
###############################################

SHOALDIR=$(pwd)/

#souce files
SHOAL_CPP=$SHOALDIR"src-cpp/shoal"
SHOAL_PY=$SHOALDIR"src-py/shoal.py"

#salmon quant files path
#BASEPATH=$SHOALDIR"quant/"

#salmon quants path
QUANTFILES=$BASEPATH"*"

#output path
#OUT=$SHOALDIR"output/"

#prior file path
PRIOR=$SHOALDIR"prior"

SAMPLE=`ls $BASEPATH | paste -s -d ,`

python $SHOAL_PY create-prior \
	--samples $SAMPLE \
	--basepath $BASEPATH \
	--outdir $PRIOR

parallel -j 4 "$SHOAL_CPP -s {} -o $OUT{/}_noprior.sf" ::: $QUANTFILES
parallel -j 4 "$SHOAL_CPP -p $PRIOR/prior.tsv -s {} -o $OUT{/}_adapt.sf -t adapt-prior" ::: $QUANTFILES
