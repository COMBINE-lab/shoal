#!/usr/bin/env bash
# Using getopt


if [ "$#" -ne 4 ]
then   # Script needs exactly two command-line argument.
  echo "Usage $0 -[options -q <quant directory>, -o <shoal output path>]"
  exit 1
fi

while getopts ":q:o:" opt; do
  case $opt in
    q)
      echo "-q<Input quant files Path> = $OPTARG"
	  BASEPATH=$OPTARG
      ;;
	o)
	  echo "-o<Output shoal files Path> = $OPTARG"
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
# make c++ executable if it doesn't yet exist
if [ ! -f src-cpp/shoal ]
then
    cd src-cpp
    make
    cd ..
fi
###############################################

SHOALDIR=$(pwd)/

#souce files
SHOAL_CPP=$SHOALDIR"src-cpp/shoal"
SHOAL_PY=$SHOALDIR"src-py/shoal.py"


#salmon quant files path
#BASEPATH=$SHOALDIR"quant/"

#salmon quants path
QUANTFILES=`find $BASEPATH -name quant.sf`
QUANTNAMES=()
QUANTDIRS=()
{
for qf in $QUANTFILES
do
    qd=$(dirname $qf)
QUANTDIRS+=($qd)
QUANTNAMES+=($(basename $qd))
done
}

# This is slick, but python is so much les cryptic
# http://superuser.com/questions/461981/how-do-i-convert-a-bash-array-variable-to-a-string-delimited-with-newlines
SAMPLE=$(IFS=$','; echo "${QUANTNAMES[*]}")

#
# Check for quant files
#
for qf in $QUANTFILES
do
    if [ ! -f $qf ] 
    then
        echo "Error: The file ${qf} should exist but does not."
        exit 1
    fi
done

for qd in "${QUANTDIRS[@]}"
do
    if [ ! -f "$qd/aux_info/eq_classes.txt" ]
    then
        echo "The quantification folder $qd is missing an equivalence class file."
        echo "Please make sure you run salmon with --dumpEqWeights."
        exit 1
    fi
done

#
# Check that the output directory is a directory
#
if [ -f $OUT ]
then
    echo "the provided output directory, $OUT, already exists, and is a file."
    exit 1
fi

# If it doesn't exist, then make it
if [ ! -e $OUT ]
then
    mkdir -p $OUT
fi

# Canonicalize the output directory
OUT=$(readlink -m $OUT)"/"

#prior file path
PRIOR=$OUT"/prior/"

python $SHOAL_PY create-prior \
  -n \
	--samples $SAMPLE \
	--basepath $BASEPATH \
	--outdir $PRIOR

parallel -j 4 "$SHOAL_CPP -p $PRIOR/prior.tsv -s {} -o $OUT{/}_adapt.sf -t adapt-prior" ::: ${QUANTDIRS[*]}
