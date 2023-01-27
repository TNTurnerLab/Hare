import os
import sys

configfile: "/dnv_wf_cpu/config.json"

REFERENCE = config["reference"]
FAMILY_FILE=config["family_file"]
GATK_OUT=config["gatk_out"]
DV_OUT=config["dv_out"]
GLNEXUS_FILE_DV_BCF=config["glnexus_file_dv_bcf"]
GLNEXUS_FILE_HC_BCF=config["glnexus_file_hc_bcf"]
GLNEXUS_FILE_DV_VCF=config["glnexus_file_dv_vcf"]
GLNEXUS_FILE_HC_VCF=config["glnexus_file_hc_vcf"]
CHROM_LENGTH=config["chrom_length"]
REGIONS=config["regions"]
GQ_VALUE=config["gq_value"]
DEPTH_VALUE=config["depth_value"]
OUT_DIR=config["out_dir"]
SUFFIX_DV=config["suffix_dv"]
SUFFIX_HC=config["suffix_hc"]
GLNEXUS_DV_MODEL=config["glnexus_dv_model"]

if not os.path.isfile(REFERENCE):
    sys.exit("Please update REFERENCE entry in the config file, the file cannot be seen or doesn't exist")
if not os.path.isfile(FAMILY_FILE):
    sys.exit("Please update FAMILY_FILE entry in the config file, the file cannot be seen or doesn't exist")
if GATK_OUT is None:
    sys.exit("Please update gatk_out entry in the config file")
if DV_OUT is None:
   sys.exit("Please update deep_out entry in the config file")
if GLNEXUS_FILE_DV_BCF is None:
    sys.exit("Please update glnexus_file_dv_bcf entry in the config file")
if GLNEXUS_FILE_HC_BCF is None:
   sys.exit("Please update glnexus_file_hc_bcf entry in the config file")
if GLNEXUS_FILE_DV_VCF is None:
    sys.exit("Please update glnexus_file_dv_vcf entry in the config file")
if GLNEXUS_FILE_HC_VCF is None:
    sys.exit("Please update glnexus_file_hc_vcf entry in the config file")
if REGIONS is None:
    REGIONS=1
if GQ_VALUE is None:
    sys.exit("Please update gq_value entry in the config file")
if DEPTH_VALUE is None:
    sys.exit("Please update depth_value entry in the config file")
if OUT_DIR is None:
    sys.exit("Please update out_dir entry in the config file")
if not os.path.isfile('/dnv_wf_cpu/filter_glnexuscombined_updated.py'):
    sys.exit("filter_glnexuscombined.py file is missing, please make sure you downloaded the git directory correctly")
if not os.path.isfile('/dnv_wf_cpu/naive_inheritance_trio_py2.py'):
    sys.exit("naive_inheritance_trio_py2.py is missing, please make sure you downloaded the git directory correctly")
if DV_OUT.endswith('/') and DV_OUT!='':
    DV_OUT=DV_OUT[:-1]
if GLNEXUS_FILE_DV_BCF.endswith('/') and GLNEXUS_FILE_DV_BCF!='':
    GLNEXUS_FILE_DV_BCF=GLNEXUS_FILE_DV_BCF[:-1]
if GLNEXUS_FILE_HC_BCF.endswith('/') and GLNEXUS_FILE_HC_BCF!='':
    GLNEXUS_FILE_HC_BCF=GLNEXUS_FILE_HC_BCF[:-1]
if GLNEXUS_FILE_DV_VCF.endswith('/') and GLNEXUS_FILE_DV_VCF!='':
    GLNEXUS_FILE_DV_VCF=GLNEXUS_FILE_DV_VCF[:-1]
if GLNEXUS_FILE_HC_VCF.endswith('/') and GLNEXUS_FILE_HC_VCF!='':
    GLNEXUS_FILE_HC_VCF=GLNEXUS_FILE_HC_VCF[:-1]
if REGIONS!=1:
    if REGIONS.endswith('/'):
        REGIONS=REGIONS[:-1]
if OUT_DIR.endswith('/') and OUT_DIR!='':
    OUT_DIR=OUT_DIR[:-1]
if GATK_OUT.endswith('/') and GATK_OUT!='':
    GATK_OUT=GATK_OUT[:-1]
if CHROM_LENGTH is None:
    CHROM_LENGTH==1
if GLNEXUS_DV_MODEL is None:
    GLNEXUS_DV_MODEL="DeepVariant"

FULL={}
FAMILIES = []
INDIVIDUALS=[]
with open(FAMILY_FILE) as f:
    for line in f:
        data=line.strip().split(',')
       # print(data)
        if len(data) ==3:
            FAMILIES.append(line.strip().split(',')[2])
            FULL[data[-1]]=data
            i=0
            while i< len(data):
                INDIVIDUALS.append(data[i])
                i+=1
#FAMILIES=["NA12878"]
GQ_VALUES=[GQ_VALUE]
DEPTH_VALUES=[DEPTH_VALUE]
print(FAMILIES)
def glnexus_intro_hc(wildcards):
    temp=[]
    FAMILY=wildcards.family
    print(FAMILY)
    print(FULL)
    caller='hc'
    for thing in FULL[FAMILY]:
        if caller=='hc':
            temp.append(GATK_OUT+'/'+thing+SUFFIX_HC)

    return (temp)
def glnexus_intro_dv(wildcards):
    temp=[]
    FAMILY=wildcards.family
    print(FAMILY)
    print(FULL)
    caller='dv'
    for thing in FULL[FAMILY]:

        if caller=='dv':
            temp.append(DV_OUT+'/'+thing+SUFFIX_DV)
    return (temp)

rule all:
    input: expand("%s/{family}.glnexus.trio.dv.bcf" % GLNEXUS_FILE_DV_BCF, family = FAMILIES),expand("%s/{family}.glnexus.trio.hc.bcf" % GLNEXUS_FILE_HC_BCF, family = FAMILIES),expand("%s/{family}.trio.hc.vcf.gz" % GLNEXUS_FILE_HC_VCF, family=FAMILIES), expand("%s/{family}.trio.hc.vcf.gz.tbi" % GLNEXUS_FILE_HC_VCF, family=FAMILIES), expand("%s/{family}.trio.dv.vcf.gz" % GLNEXUS_FILE_DV_VCF, family=FAMILIES), expand("%s/{family}.trio.dv.vcf.gz.tbi" % GLNEXUS_FILE_DV_VCF, family=FAMILIES), expand("%s/{family}/{family}.glnexus_denovo_actual.dv.vcf.gz" % OUT_DIR, family = FAMILIES),expand("%s/{family}/{family}.glnexus_denovo_actual.hc.vcf.gz" % OUT_DIR, family = FAMILIES), expand("%s/{family}/{family}.glnexus.family.combined_intersection_filtered_gq_{gq_val}_depth_{depth_val}_position.vcf" %OUT_DIR,family = FAMILIES,gq_val=GQ_VALUES,depth_val=DEPTH_VALUES)


rule glnexus_hc:
    input: glnexus_intro_hc
    output: "%s/{family}.glnexus.trio.hc.bcf" % GLNEXUS_FILE_HC_BCF, "%s/{family}.trio.hc.vcf.gz" % GLNEXUS_FILE_HC_VCF, "%s/{family}.trio.hc.vcf.gz.tbi" % GLNEXUS_FILE_HC_VCF
    params: prefix="{family}"
    shell: """
    mkdir -p {GLNEXUS_FILE_HC_BCF}
    export PATH=/opt/conda/bin:$PATH
    /glnexus_cli --config gatk --mem-gbytes 64 --dir /dnv_wf_cpu/{params.prefix}_hc {input} > {output[0]}
    bcftools view {output[0]} | bgzip -c > {output[1]}
    tabix {output[1]}
    rm -r /dnv_wf_cpu/{params.prefix}_hc
    """
rule glnexus_dv:
    input: glnexus_intro_dv
    output: "%s/{family}.glnexus.trio.dv.bcf" % GLNEXUS_FILE_DV_BCF,"%s/{family}.trio.dv.vcf.gz" % GLNEXUS_FILE_DV_VCF, "%s/{family}.trio.dv.vcf.gz.tbi" % GLNEXUS_FILE_DV_VCF
    params: prefix="{family}"
    shell: """
    mkdir -p {GLNEXUS_FILE_DV_BCF}
    export PATH=/opt/conda/bin:$PATH
    /glnexus_cli --config {GLNEXUS_DV_MODEL} --mem-gbytes 64 --dir /dnv_wf_cpu/{params.prefix}_dv {input} > {output[0]}
    bcftools view {output[0]} | bgzip -c > {output[1]}
    tabix {output[1]}
    rm -r /dnv_wf_cpu/{params.prefix}_dv
    """

rule grabFamilies_DV:
    input:  "%s/{family}.trio.dv.vcf.gz" % GLNEXUS_FILE_DV_VCF
    output:  "%s/{family}/{family}.glnexus_denovo_actual.dv.vcf.gz" % OUT_DIR
    params:  prefix="{family}"
    shell: """
    export PATH=/opt/conda/envs/py2/bin:$PATH
    export PATH=/opt/conda/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib/openblas-base/
    mkdir -p {OUT_DIR}
    trio=$( grep {params.prefix} {FAMILY_FILE} )
    echo $trio
    echo "Keep Samples DV"
    mkdir -p {OUT_DIR}/{params.prefix}
    #Grabs the trio from the GLnexus DeepVariant (DV) file, then compresses and creates an index
    bcftools view -Oz --threads 4 -s $trio {input}  > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.dv.vcf.gz

    tabix {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.dv.vcf.gz


    echo "Naive search DV"
    #Searches family trio .vcf file for possible denovo varints, adding the label TRANSMITTED=no;INH=denovo_pro if the variant is a candidate.
    #It then pulls out the varints with that label into a separate .vcf file.

    python2 /dnv_wf_cpu/naive_inheritance_trio_py2.py  {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.dv.vcf.gz > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo.dv.vcf
    grep '#' {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo.dv.vcf > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.dv.vcf
    grep -w 'TRANSMITTED=no;INH=denovo_pro' {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo.dv.vcf >> {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.dv.vcf
    bgzip {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.dv.vcf
    tabix {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.dv.vcf.gz




    ls {output}
    """

rule grabFamilies_HC:
    input: "%s/{family}.trio.hc.vcf.gz" % GLNEXUS_FILE_HC_VCF
    output:  "%s/{family}/{family}.glnexus_denovo_actual.hc.vcf.gz" % OUT_DIR
    params:  prefix="{family}"
    shell: """
    mkdir -p {OUT_DIR}
    export PATH=/opt/conda/envs/py2/bin:$PATH
    export PATH=/opt/conda/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib/openblas-base/
    mkdir -p {OUT_DIR}/{params.prefix}
    trio=$( grep {params.prefix} {FAMILY_FILE} )
    echo $trio
    echo "Keep Samples HC"
    #Same steps as above, only with the HaplotypeCaller (HC) GLnexus output

    bcftools view -Oz --threads 4 -s $trio {input}  > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.hc.vcf.gz

    tabix {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.hc.vcf.gz


    echo "Naive search HC"
    python2 /dnv_wf_cpu/naive_inheritance_trio_py2.py  {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.hc.vcf.gz > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo.hc.vcf

    grep '#' {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo.hc.vcf > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.hc.vcf
    grep -w 'TRANSMITTED=no;INH=denovo_pro' {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo.hc.vcf >> {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.hc.vcf
    bgzip {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.hc.vcf
    tabix {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.hc.vcf.gz

    """

rule combinedAndFilter:
    input: "%s/{family}/{family}.glnexus_denovo_actual.hc.vcf.gz" % OUT_DIR,"%s/{family}/{family}.glnexus_denovo_actual.dv.vcf.gz" % OUT_DIR
    output: "%s/{family}/{family}.glnexus.family.combined_intersection_filtered_gq_{gq_val}_depth_{depth_val}_position.vcf" % OUT_DIR
    params: prefix="{family}", gq="{gq_val}", depth="{depth_val}"
    shell:"""
    #Filters out -L chromosomes that are not -L chr1-22, grabbing variants with an allele count of 1, variants that were found from both DV and HC (defined as set=Intersection), and removes variants that are either 10 A's or T's in a row.
    echo "Combining files"
    
    ref=$( head -n 1 {REFERENCE} | cut -d' ' -f 1)
    echo $ref
    /opt/conda/bin/python /dnv_wf_cpu/test_intersect.py -g {input[0]} -d {input[1]} -r $ref -c {CHROM_LENGTH}
    cat {OUT_DIR}/{params.prefix}/{params.prefix}_combined_out.vcf |  awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.combined.vcf
    zcat {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.dv.vcf.gz  | grep '#' > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection.vcf
    grep -v 'chrUn' {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus_denovo_actual.combined.vcf | grep -v '_random' | grep -v '_alt'  | grep -v 'chrY' | grep -v 'chrM' | grep  'AC=1' |  egrep -v 'AAAAAAAAAA|TTTTTTTTTT' >> {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection.vcf


    #Finds the order of the family position found within the combined .vcf file and what the order is in the family file (which should be in order of father, mother, child)
    echo 'Set up filter script'
    actual=$( cat {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection.vcf | grep '#' | tail -n 1 | cut -f 10-13 )
    echo $actual
    temp=$( grep {params.prefix} {FAMILY_FILE} )
    temp=$( echo $temp | cut -d',' -f 1-3 )
    trio=$( echo $temp | tr ',' ' ')
    echo $trio
    echo "Run filter script"
    #Python script that filters for parents with no alt allele, depth of set value, GQ of set value, and allele balance of .25
    /opt/conda/bin/python /dnv_wf_cpu/filter_glnexuscombined_updated.py {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection.vcf $trio $actual {params.gq} {params.depth}

    echo "Make position file"
    cat {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_20_depth_10.vcf | grep '#' > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_{params.gq}_depth_{params.depth}_position.vcf
    echo "Make CpG file"
    cat {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_20_depth_10.vcf | grep '#' > {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_{params.gq}_depth_{params.depth}_position_all.vcf
    echo "Filter by position"

    if [ "{REGIONS}" != "1" ]
    then
        /opt/conda/bin/bedtools intersect -v -a {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_{params.gq}_depth_{params.depth}.vcf -b {REGIONS}/LCR-hs38-5bp-buffer.bed.gz    {REGIONS}/hg38_centromeres_09252018.bed.gz  {REGIONS}/recent_repeat_b38-5bp-buffer.bed.gz  >> {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_{params.gq}_depth_{params.depth}_position.vcf
        /opt/conda/bin/bedtools intersect -a {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_{params.gq}_depth_{params.depth}_position.vcf -b {REGIONS}/CpG_sites_sorted_b38.bed.gz >> {OUT_DIR}/{params.prefix}/{params.prefix}.glnexus.family.combined_intersection_filtered_gq_{params.gq}_depth_{params.depth}_position_all.vcf
    fi
    ls {output}
    """
