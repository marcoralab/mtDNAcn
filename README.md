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

NP.1 Neuropathology Category as measured by CERAD

| Value | Coding | n |
| ----- | ------ | - |
| 1     | Normal | 91 |
| 2     | Definite AD | 161 |
| 3     | probable AD | 41 |
| 4     | possible AD | 40 |

CDR Clinical Dementia Rating

| Value | Coding | n |
| ----- | ------ | - |
| 0     | no dementia | 39 |
| 0.5     | Questionable dementia | 43 |
| 1     | Mild dementia | 36 |
| 2     | Moderate Dementia | 41 |
| 3     | Severe dementia | 74 |
| 4     | Profound dementia | 45 |
| 5     | Terminal Dementia | 55 |

bbscore - Braak Score

| Value | Coding | n |
| ----- | ------ | - |
| 0     |  | 13 |
| 1     | I - transentorhinal | 27 |
| 2     | II - transentorhinal | 43 |
| 3     | III - limbic | 49 |
| 4     | IV - limbic | 32 |
| 5     | V - neocortical | 37 |
| 6     | VI  - neocortical | 100 |
| NA     |  | 22 |

NIA-Reagan (Kurt)

| CERAD / Braak | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
| ---------- | - | - | - | - | - | - | - |
| Normal | Normal | Other | Other | Other | Other | Other | Other |
| Possible | Low | Low | Low | Low | Low | Low | Low |
| Probable | Low | Low | Low | Intermediate | Intermediate | Intermediate | Intermediate |
| Definite | Low | Low | Low | Intermediate | Intermediate | High | High |

NIA-Reagan

| CERAD / Braak | 1 | 2 | 3 | 4 | 5 | 6 |
| ---------- | - | - | - | - | - | - |
| Normal | No AD | No AD | - | - | - | - |
| Possible | Low | Low | - | - | - | - |
| Probable | - | - | Intermediate | Intermediate | - | - |
| Definite | - | - | - | - | High | High |

| CERAD / Braak | 0 | 1 | 2 | 3 | 4 | 5 | 6 | NA |
| ---------- | -- | -- | -- | -- | -- | -- | -- | -- |
| Normal | 12 | 21 | 30 | 20 | 2 | NA | NA | 6 |
| Possible | 1  | 6  | 4 | 16 | 10 | 2 | 1 | NA |
| Probable | NA | NA | 4 | 12 | 6  | 8  | 8  | 3  |
| Definite | NA | NA | 5 | 11 | 14 | 27 | 91 | 13 |


| Value | Coding | CERD | Braak | n |
| ----- | ------ | ---- | ----- | - |
| 0     | No AD | Normal | I/II |
| 1     | Low | Possible AD | I/II | 27 |
| 2     | Intermediate | Probable AD | III/IV | 43 |
| 3     | High | Definite AD | V/V1 | 49 |

NIA-AA

| Thal | CERAD               | Braak: | None or I/II | III/IV | V/VI | 
| -    | ----------          | ------ | - | - | - |
| 0    | None                |        | Other<sup>§</sup> | Other<sup>§</sup>  | Other<sup>§</sup>  |
| 1/2  | None - Sparse       |        | Low | Low | Low<sup>¶</sup>  |
|      | Modearte - Frequent |        | Low<sup>†</sup>  | Intermediate | Intermediate<sup>¶</sup>  |
| 3    | Any                 |        | Low<sup>†</sup>  | Intermediate | Intermediate<sup>¶</sup>  |
| 4/5  | None - Sparse       |        | Low<sup>†</sup>  | Intermediate | Intermediate<sup>¶</sup>  |
|      | Modearte - Frequent |        | Low<sup>†</sup>  | Intermediate | High |

<sup>†</sup>
§
¶
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

| rs429358	| rs7412	| Name |
| -------- | ------ | ------ |
| 19:45411941:T:C | 19:45412079:C:T | NA |
| C	| T	| ε1 |
| T	| T	| ε2 |
| T	| C	| ε3 |
| C	| C	| ε4 |

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
