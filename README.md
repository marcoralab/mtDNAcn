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

```
/sc/orga/projects/LOAD/shea/Projects/mDNACN
```

## ROS/MAP

cogdx

| Value | Coding | n |
| ----- | ------ | - |
| 1     | NCI: No cognitive impairment (No impaired domains) | 374 |
| 2     | MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI | 273 |
| 3     | MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI | 20 |
| 4     | 	AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD) | 433 |
| 5     | 	AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD) | 57 |
| 6     | 	Other dementia: Other primary cause of dementia | 21 |

msex

| Value | n |
| ----- | - |
| F | 779 |
| M | 399 |

race

| Value | Coding | n |
| ----- | ------ | - |
| 1 | White | 1177 |
| 2 | Black | 1 |

apoe_genotype

| Value | n |
| ----- | - |
| 22 | 7 |
| 23 | 146 |
| 24 | 22 |
| 33 | 709 |
| 34 | 263 |
| 44 | 18 |
| NA | 13 |

## MSBB

SEX

| Value | n |
| ----- | - |
| F | 216 |
| M | 117 |

apoe

| Value | n |
| ----- | - |
| 22 | 2 |
| 23 | 19 |
| 24 | 1 |
| 33 | 96 |
| 34 | 53 |
| 44 | 3 |
| NA | 159 |

RACE

| Value | Coding | n |
| ----- | ------ | - |
| A | Asian | 1 |
| B | Black | 31 |
| H | Hispanic | 23 |
| U | Uknown | 1 |
| W | White | 277 |

dx

| Value | n |
| ----- | - |
| AD | 156 |
| CONTROL | 55 |
| OTHER | 122 |

## Mayo

diagnosis

| Value | n |
| ----- | - |
| AD | 92 |
| Control | 100 |
| Pathologic Aging | 74 |
| PSP | 83 |

Sex

| Value | n |
| ----- | - |
| F | 182 |
| M | 167 |

ApoE

| Value | n |
| ----- | - |
| 22 | 1 |
| 23 | 36 |
| 24 | 2 |
| 33 | 176 |
| 34 | 69 |
| 44 | 9 |
| NA | 56 |

WGS_Source_Tissue_Type

| Value | n |
| ----- | - |
| Cerebellar Cortex | 8 |
| Temporal Cortex | 341 |

## fastMitoCalc
https://lgsun.irp.nia.nih.gov/hsgu/software/mitoAnalyzer/index.html

Qian, Yong, Thomas J. Butler, Krista Opsahl-Ong, Nicholas S. Giroux, Carlo Sidore, Ramaiah Nagaraja, Francesco Cucca, et al. 2017. “fastMitoCalc: An Ultra-Fast Program to Estimate Mitochondrial DNA Copy Number from Whole-Genome Sequences.” Bioinformatics  33 (9): 1399–1401.


## Useage
### Whole Genome Sequencing
```
bsub -J "mDNACN" -P acc_load -q premium -n 2 -R span[hosts=1] -R rusage[mem=16000] -W 01:00 -L /bin/bash -o data/mDNACN.stdout -eo data/mDNACN.stderr \
"perl src/fastMitoCalc/fastMitoCalc.pl -f raw/MAP15387421.final.bam -w data -p src/fastMitoCalc/BaseCoverage"
```

```
src/fastMitoCalc/fastMitoCalc.pl -f raw/MAP15387421.final.bam -w data -p src/fastMitoCalc/BaseCoverage
```

### Whole Exome Sequencing

If sequence reads from the mtDNA are not targeted, need to use off target reads only. The requires .bed file of the targeted captures regions for the scpefic exome kit. For example the [IDT xGen Exome Research Panel from the UKB](http://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=3801).

Extract all off target reads
```
bedtools intersect -v -abam the.bam -b the.bed -bed | awk '{print $1, $2, $3}' OFS="\t" > test.bed
```

Extract all off target reads that 50kb away from target reads
- extract off target reads
- write out chr, start, end
- remove values of -1
- remove off target reads within 50kb

```
bedtools intersect -v -abam WES.bam -b Targets.bed -bed | \
awk '{print $1, $2, $3}' OFS="\t" | awk '$2!=-1' | \
bedtools window -w 50000 -v -a stdin -b Targets.bed > OffTargets.bed
```

Subset the number of reads from the .bed file to reduce computation time
- Select reads maping to autosomes only
- calculat length of reads and remove any with less then 100bp
- randomly sample 2% of reads from each chromosome
- Write out

```
read_table2('OffTargets.bed', col_names = F, col_types = list(X1 = col_character())) %>%
  filter(X1 %in% c(1:22)) %>%
  mutate(length = X3 - X2) %>%
  filter(length >= 100) %>%
  group_by(X1) %>%
  sample_frac(0.02, replace = F) %>%
  ungroup() %>%
  write_tsv('OffTargetsClean.bed', col_names = F)
```

fastMitoCalc limited to the regions in the bed file only
```
perl scr/fastMitoCalc.pl -f WES.bam -b OffTargetsClean.bed -w data -p src/BaseCoverage
```
