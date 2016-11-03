# shoal 
Improved multi-sample transcript abundance estimates using adaptive priors
--------------------------------------------------------------------------

![A shoal](https://upload.wikimedia.org/wikipedia/commons/b/b1/School_jacks_klein.JPG)<sup id="a1">[1](#f1)</sup>


# What is shoal?


# Using shoal
Shoal requires to have salmon output from version --: xx

* clone shoal into your local machine:
```
git clone https://github.com/COMBINE-lab/shoal.git
```

* change current working directory to shoal and make folder with name `quant` inside it
```
cd shoal; mkdir quant;
```

* move the salmon output *directory* of all the samples in the experiment into the `quant` directory
```
mv ~/ConA* quant/
mc ~/ConB* quant/
```

* run shell script
```
./run_shoal.sh
```
NOTE: shell script can be given executable permission with command: ```chmod +x run_shoal.sh```

Footnotes:
----------
<b id="f1">1</b> This image is from [the wikipedia artical on shoaling](https://en.wikipedia.org/wiki/Shoaling_and_schooling#/media/File:School_jacks_klein.JPG). It is licensed under CC-BY-SA.[↩](#a1)
