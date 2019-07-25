# mDNACN
Association of mitochondrial DNA CN with Alzheimer's disease


Ricardos AMPAD files

```/sc/orga/projects/LOAD/Data/AMPAD```

Location of AMPAD bim files.

```
/sc/orga/projects/AMPADWGS
/sc/orga/projects/AMPADWGS/RawData
/sc/orga/projects/AMPADWGS/RawDataSinai
/sc/orga/projects/AMPADWGS/RawDataMayo
```

## fastMitoCalc
https://lgsun.irp.nia.nih.gov/hsgu/software/mitoAnalyzer/index.html

Qian, Yong, Thomas J. Butler, Krista Opsahl-Ong, Nicholas S. Giroux, Carlo Sidore, Ramaiah Nagaraja, Francesco Cucca, et al. 2017. “fastMitoCalc: An Ultra-Fast Program to Estimate Mitochondrial DNA Copy Number from Whole-Genome Sequences.” Bioinformatics  33 (9): 1399–1401.


##
Useage

```
bsub -J "mDNACN" -P acc_load -q premium -n 2 -R span[hosts=1] -R rusage[mem=16000] -W 01:00 -L /bin/bash -o data/mDNACN.stdout -eo data/mDNACN.stderr \
"perl src/fastMitoCalc/fastMitoCalc.pl -f raw/MAP15387421.final.bam -w data -p src/fastMitoCalc/BaseCoverage"
```

```
src/fastMitoCalc/fastMitoCalc.pl -f raw/MAP15387421.final.bam -w data -p src/fastMitoCalc/BaseCoverage
```
