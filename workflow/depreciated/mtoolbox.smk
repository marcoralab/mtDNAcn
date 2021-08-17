import os
from itertools import product
import pandas as pd
RWD = os.getcwd()

shell.prefix('module load R/3.6.3 samtools singularity freebayes/1.3.2 vcflib bcftools tabix bedtools; ')

ID = ["HG00119", "HG00124", "SM-CTENH", "1000137_23163_0_0", "1000198_23163_0_0", "1000421_23163_0_0", "1000516_23163_0_0"]

rule all:
    input:
        expand("data/mtoolbox/{id}/mtoolbox_config_{id}.sh", id = ID),
        expand(RWD + "/data/mtoolbox/{id}/logassemble.txt", id = ID),
        expand(RWD + "/data/mtoolbox/{id}/mt_classification_best_results.csv", id = ID),
        expand(RWD + "/data/mtoolbox/{id}/prioritized_variants.txt", id = ID),
        expand(RWD + "/data/mtoolbox/{id}/sample.vcf", id = ID),

rule mtoolbox_config:
    input: "sandbox/data/{id}.bam"
    output: "data/mtoolbox/{id}/mtoolbox_config_{id}.sh"
    params:
        input_path = RWD + '/sandbox/data',
        list = '{id}',
        output_name = RWD + '/data/mtoolbox/{id}'
    shell:
        """
        printf '%s\n' \
            'input_type=bam' \
            'ref=RCRS' \
            'input_path={params.input_path}' \
            'list={params.list}' \
            'output_name={params.output_name}' > \
            {output}
        """

rule mtoolbox:
    input:
        bam = "sandbox/data/{id}.bam",
        config = "data/mtoolbox/{id}/mtoolbox_config_{id}.sh"
    output:
        RWD + "/data/mtoolbox/{id}/logassemble.txt",
        RWD + "/data/mtoolbox/{id}/mt_classification_best_results.csv",
        RWD + "/data/mtoolbox/{id}/prioritized_variants.txt",
        RWD + "/data/mtoolbox/{id}/sample.vcf",
        directory(RWD + "/data/mtoolbox/{id}/OUT_{id}"),
        temp(RWD + "/data/mtoolbox/{id}/processed_fastq.tar.gz"),
        temp(RWD + "/data/mtoolbox/{id}/VCF_dict_tmp"),
    conda: "envs/environment.yaml"
    shell: "MToolBox.sh -i {input.config}"
