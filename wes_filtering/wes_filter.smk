import os

if config == {}:
    configfile: "config.json"

DATA_DIR = config["data_dir"]
REFERENCE = config["reference"]
FAMILY_FILE=config['famfile']
OUT_DIR=config["out_dir"]
CRAM_LOCATIONS=config["cramloc"]
LOWER=config["lower"]
UPPER=config["upper"]
CASE_SPECIFIC=config["case_specific"]
INTERVAL=config['interval']
CHANGE=config['change']
CHANGE=int(CHANGE)

SAMPLES = []
with open(FAMILY_FILE) as input:
    for line in input:
        data=line.strip().split(',')
        SAMPLES.append(data[-1])

print(SAMPLES)

rule all: 
    input: expand("%s/{sample}_mpileup"%OUT_DIR, sample=SAMPLES), expand("%s/{sample}_mpileup_filtered_checkme.txt"%OUT_DIR,sample=SAMPLES), "mod_bedfile/changed_bedfile_pm%s.bed"%CHANGE

rule modify_bedfile:
    input: "%s"%INTERVAL
    output: "mod_bedfile/changed_bedfile_pm%s.bed"%CHANGE
    shell: """
    export PATH=/opt/conda/envs/wes_filter/bin:$PATH
     
    mkdir -p mod_bedfile

    python modify_bed.py -i {input} -c {CHANGE}
    """

rule mpileup:
    input: "%s/{sample}.glnexus.family.combined_intersection_filtered_gq_20_depth_10_position.vcf"%DATA_DIR, "mod_bedfile/changed_bedfile_pm%s.bed"%CHANGE
    output: "%s/{sample}_mpileup"%OUT_DIR, "%s/{sample}_mpileup_filtered_checkme.txt"%OUT_DIR
    params: prefix="{sample}"
    shell: """

    export PATH=/opt/conda/envs/wes_filter/bin:$PATH

    fam=$( grep -w {params.prefix} {FAMILY_FILE} )


    if [ {CASE_SPECIFIC} == "False" ]
    then

    python get_pileup_updated.py -v {input[0]} -p {CRAM_LOCATIONS} -r {REFERENCE} -f $fam -o {OUT_DIR} -l {LOWER} -u {UPPER} -i {input[1]}
    else

    python get_pileup_updated.py -v {input[0]} -p {CRAM_LOCATIONS} -r {REFERENCE} -f $fam -o {OUT_DIR} -l {LOWER} -u {UPPER} -c -i {input[1]}
    fi
    rm -r {OUT_DIR}/{params.prefix}_bed_coor.txt
    """
