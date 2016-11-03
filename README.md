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
Shoal requires to have salmon output from version --: xx

* clone shoal into your local machine:
```bash
git clone https://github.com/COMBINE-lab/shoal.git
```

* change current working directory to shoal and make folder with name `quant` inside it
```bash
cd shoal; mkdir quant;
```

* move the salmon output ***directory*** of all the samples in the experiment into the `quant` directory
```bash
mv ~/ConA* quant/
mv ~/ConB* quant/
```

* run shell script
```bash
./run_shoal.sh
```
NOTE: shell script can be given executable permission with command: ```chmod +x run_shoal.sh```

Footnotes:
----------
<b id="f1">1</b> This image is from [the wikipedia artical on shoaling](https://en.wikipedia.org/wiki/Shoaling_and_schooling#/media/File:School_jacks_klein.JPG). It is licensed under CC-BY-SA.[↩](#a1)
