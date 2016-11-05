# shoal 
Improved multi-sample transcript abundance estimates using adaptive priors
--------------------------------------------------------------------------

![A shoal](https://upload.wikimedia.org/wikipedia/commons/b/b1/School_jacks_klein.JPG)<sup id="a1">[1](#f1)</sup>


### What is shoal?
shoal is a tool which jointly quantify transcript abundances across multiple samples.
Specifically, shoal learns an empirical prior on transcript-level abundances across all of the 
samples in an experiment, and subsequently applies a variant of the variational Bayesian 
expectation maximization algorithm to apply this prior adaptively across multi-mapping groups 
of reads.  
shoal can increase quantification accuracy, inter- sample consistency, and reduce false positives in 
downstream differential analysis when applied to multi-condition RNA-seq experiments. 
Moreover, shoal, runs downstream of Salmon and requires less 
than a minute per-sample to re-estimate transcript abundances while accounting for the learned empirical prior.

### Using shoal
Shoal requires to have salmon output of all the samples in the experiment
separately using the latest version of Salmon (from the [develop branch of the Salmon repo](https://github.com/COMBINE-lab/salmon/tree/develop); a tagged release and pre-compiled binary will be made available soon).  Please run Salmon with the `--dumpEqWeights` option, which will produce output suitable for shoal.

* clone shoal into your local machine:
```bash
git clone https://github.com/COMBINE-lab/shoal.git
```

* run shoal<sup id="a2">[2](#f2)</sup>:
```bash
./run_shoal.sh -q <salmon_quant_directory_path> -o <output_directory_path>
```
This script assumes that all of the Salmon quantification directories are _subdirectories_ of the path that you provide via the `-q` option.  So, e.g., if you have an experiment with six samples across 2 conditions (say, `A{1,2,3}` and `B{1,2,3}`), then the shoal script would expect a layout like:

```
exp_quants
  |
  |--- A1
     |
     |--- quant.sf
  |--- A2
    |
    |--- quant.sf
  |--- A3
    |
    |--- quant.sf
  |--- B1
    |
    |--- quant.sf
  |--- B2
    |
    |--- quant.sf
  |--- B3
    |
    |--- quant.sf
```

the script would then be invoked by passing `-q exp_quants` to provide the top-level quantification directory for the entire experiment.
Specifically, a command like `./run_shoal -q exp_quants -o exp_shoal_quants` would produce a modified (Salmon-format) quantification file for each of the samples (`{A,B}{1,2,3}`) in the directory `exp_shoal_quants` as described below (the script will create the output directory if it does not already exist). 

* shoal output:  
-- shoal generates `.sf` files for each sample in the experiment with naming convention as follows:
```bash
<output_directoty>/<sample_name>_adapt.sf
```

Footnotes:
----------
<b id="f1">1</b> This image is from [the wikipedia artical on shoaling](https://en.wikipedia.org/wiki/Shoaling_and_schooling#/media/File:School_jacks_klein.JPG). It is licensed under CC-BY-SA.[↩](#a1)

<b id="f2">2</b> shell script can be given executable permission with command: `chmod +x run_shoal.sh`
