# Hare
#### Developed and maintained by Jeffrey K. Ng
#### Maintainer Email:  jeffrey.ng@wustl.edu
#### Washington University in St. Louis Medical School
#### Tychele N. Turner, Ph.D., Lab


Hare is a *de novo* variant caller leveraging the power of Parabricks GPU accelerated variant calling.  This is an updated version of the [workflow described in Ng et al. 2022.](https://doi.org/10.1002/humu.24455)  The original code of the workflow can be [found here.](https://github.com/TNTurnerLab/GPU_accelerated_de_novo_workflow)  This version has been tested with the output from Parabricks v.3.0.0, as well as, the now free version of [Parabricks v4.0.0.0-1](https://docs.nvidia.com/clara/parabricks/4.0.0/index.html).  You can find dependancies and instructions on how to run Parabricks [here.](https://docs.nvidia.com/clara/parabricks/4.0.0/GettingStarted.html)

# How to Run
## Input

Three main inputs:
1)  .bam or .cram files for the trio(s)
2)  A comma delimited text file, with one trio per line, with sample IDs formatted in the following way:  Father,Mother,Child
3)  The reference genome .fasta used when running GATK and DeepVariant

You'll first need to run the your crams through Parabricks GATK Haplotypecaller and DeepVariant, the instructions of which can be found above.  When doing so, please ensure you are using the ```--gvcf``` flag.You'll use that output for the Snakemake file found in this repo.  

**NOTE** This pipeline has been specifically tested using output from Parabricks GATK Haplotypecaller and DeepVariant.  Other g.vcf output may not work properly

This workflow also makes use of specific RepeatMasker files, links to which can be found below.
 
 
#### Centromeres
```
wget -q https://de.cyverse.org/dl/d/B42A0F3D-C402-4D5F-BBD5-F0E61BE2F4AC/hg38_centromeres_09252018.bed.gz
wget -q https://de.cyverse.org/dl/d/37B13DB5-0478-4C4B-B18D-33AFB742E782/hg38_centromeres_09252018.bed.gz.tbi
```
 
#### LCR plus 5 bp buffer
```
wget -q https://de.cyverse.org/dl/d/870755FF-CD04-4010-A1EC-658D7E1151EF/LCR-hs38-5bp-buffer.bed.gz
wget -q https://de.cyverse.org/dl/d/01D038EA-51CC-4750-9814-0BB3784E808E/LCR-hs38-5bp-buffer.bed.gz.tbi
```
 
#### Recent repeats plus 5 bp buffer
```
wget -q https://de.cyverse.org/dl/d/185DA9BC-E13D-429B-94EA-632BDAB4F8ED/recent_repeat_b38-5bp-buffer.bed.gz
wget -q https://de.cyverse.org/dl/d/4A6AF6EF-D3F0-4339-9B8E-3E9E83638F00/recent_repeat_b38-5bp-buffer.bed.gz.tbi
```
 
#### CpG Locations
```
wget -q https://de.cyverse.org/dl/d/786D1640-3A26-4A1C-B96F-425065FBC6B7/CpG_sites_sorted_b38.bed.gz
wget -q https://de.cyverse.org/dl/d/713F020E-246B-4C47-BBC3-D4BB86BFB6E9/CpG_sites_sorted_b38.bed.gz.tbi
```



## Running
 
While the use of Docker is highly recommended, the workflow is able to run outside of a docker envionment, please see the software dependancies below.  If you would like to run Parabricks 4.0.0, having a system able to run Docker **is** a requirement.

The workflow Docker image can be pulled from here:
```
tnturnerlab/hare:v1.1
```
 
### Snakemake

#### Setting up the config.json file
Before running, please make any necessary changes to these options below in the config.json. 
 ```
* regions:  "/region" *If you don't have the RepeatMasker files, please make this entry blank*
* gq_value:  20 *Default gq value filter*
* depth_value: 10 *Default depth value filter*
* suffix_dv:  *Suffix of the DeepVariant data files.  Assumes input files are \<sample\_name\>\<suffix\>* 
* suffix_hc:  *Suffix of the GATK Haplotypecaller data files.  Assumes input files are \<sample\_name\>\<suffix\>* 
* family_file: "/dnv_wf_cpu/<your_family_file>"
* chrom_length:  *Optional chromosome length file, use if you are not using human reference GRCh38.  Can leave blank if using GRCh38.  Please make this a two column, tab delimited file, with the first chromosome and the second column the length of the chromosome*
```
Below is an example Docker run command:

```
docker run -v "/path/to/hare/code:/dnv_wf_cpu" -v "/path/to/reference:/reference" -v "/path/to/deepvariant/output:/dv" -v "/path/to/gatk/output:/gatk"  -v "/path/to/RepeatMasker/region/files:/region" tnturnerlab/hare:v1.1 /opt/conda/envs/snake/bin/snakemake -s /dnv_wf_cpu/hare_1.1.smk -j 6 --cores -k --rerun-incomplete -w 120 
```


### Cromwell 

We also provide this workflow in a .wdl format.  Unlike the Snakemake, you will be able to run Parabricks directly from this workflow, instead of separately.  You can also run this workflow in the cloud.  To run this, you'll need to download the [Cromwell .jar found here](https://github.com/broadinstitute/cromwell/releases).  This wdl was specificlly tested on cromwell-83.  

The basic config file looks like this:
```
{
  "jumping_hare.num_ram_hc": "Int (optional, default = 120)",
  "jumping_hare.extra_mem_hc": "Int (optional, default = 65)",
  "jumping_hare.maxPreemptAttempts": "Int (optional, default = 3)",
  "jumping_hare.cpu_hc": "Int (optional, default = 24)",
  "jumping_hare.glnexus_deep_model": "String (optional, default = \"DeepVariant\")",
  "jumping_hare.test_intersect": "File", #pathway to the test_intersect.py file
  "jumping_hare.deep_model": "String (optional, default = \"shortread\")",
  "jumping_hare.gpuDriverVersion_DV": "String (optional, default = \"460.73.01\")",
  "jumping_hare.sample_suffix": "String", #suffix of the input cram file.  If your sample was NA12878.final.cram, you would put ".final.cram" here
  "jumping_hare.typeOfGPU_HC": "String (optional, default = \"nvidia-tesla-t4\")",
  "jumping_hare.gq": "Int (optional, default = 20)",
  "jumping_hare.glnexus_ram_dv": "Int (optional, default = 100)",
  "jumping_hare.filter_glnexuscombined_updated": "File",  #pathway to filter_glnexuscombined_updated.py
  "jumping_hare.num_gpu_HC": "Int (optional, default = 2)",
  "jumping_hare.num_ram_dv": "Int (optional, default = 120)",
  "jumping_hare.naive_inheritance_trio_py2": "File", #pathway to naive_inheritance_trio_py2.py
  "jumping_hare.num_gpu_dv": "Int (optional, default = 4)",
  "jumping_hare.glnexus_DV.extramem_GLDV": "Int? (optional)",
  "jumping_hare.extra_mem_dv": "Int (optional, default = 65)",
  "jumping_hare.typeOfGPU_DV": "String (optional, default = \"nvidia-tesla-t4\")",
  "jumping_hare.combinedAndFilter.extramem_GLDV": "Int? (optional)",
  "jumping_hare.glnexus_ram_hc": "Int (optional, default = 100)",
  "jumping_hare.deep_docker": "String (optional, default = \"nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1\")",
  "jumping_hare.pathToReference": "File",  #pathway to tarball of reference information
  "jumping_hare.wes": "Boolean (optional, default = false)", #Please set this to true if you are analyzing WES data
  "jumping_hare.glnexus_cpu": "Int (optional, default = 32)",
  "jumping_hare.gpuDriverVersion_HC": "String (optional, default = \"460.73.01\")",
  "jumping_hare.cram_files": "Array[Array[WomCompositeType {\n cram -> File\ncrai -> File \n}]]", #cram/bam file input, please see example for formating
  "jumping_hare.cpu_dv": "Int (optional, default = 24)",
  "jumping_hare.glnexus_HC.extramem_GLDV": "Int? (optional)",
  "jumping_hare.interval_file": "String (optional, default = \"None\")",
  "jumping_hare.depth": "Int (optional, default = 10)",
  "jumping_hare.reference": "String", #name of reference fasta
  "jumping_hare.regions": "File? (optional)",
  "jumping_hare.hare_docker": "String (optional, default = \"tnturnerlab/hare:v1.1\")",
  "jumping_hare.trios": "Array[WomCompositeType {\n father -> String\nmother -> String\nchild -> String \n}]", #trios, MUST be in same order as trios in cram_file
  "jumping_hare.chrom_length": "File? (optional)" #Optional chromosome length file if you are not using Human build GRCh38
}
```
Required arguments are highlighted in comments above.  We have provided an example config to help with formatting. Please modify the computational requirements to fit your HPC.  If you are running it on Google CLoud Platform, you may keep the computation settings. Requirements are based on [NVIDIA's own workflows found here.](https://github.com/clara-parabricks-workflows/parabricks-wdl)  If you are going to use this wdl, please tarball your reference files.  If you are running WES data, please include your capture region in this tarball.
```
tar -jcf reference.tar.bz2 reference.fa reference.fa.fai reference.dict
```

We also provide the Dockerfile if you would like to make modifications. 

 # Output
 Below is a brief descirption of the main output folders from Hare:
 * dv_bcf: Output folder for GLnexus .bcf files for DeepVariant output.
 * hc_bcf: Output folder for GLnexus .bcf files for GATK HaplotypeCaller output.
 * dv_vcf: Output folder for converted GLnexus .bcf files to .vcf.gz files for DeepVariant output.  Also includes tabix index file.
 * hc_vcf: Output folder for converted GLnexus .bcf files to .vcf.gz files for GATK HaplotypeCaller output.  Also includes tabix index file.
 * out_hare:  Output folder where the *de novo* variant files can be found.  If you are running multiple trios, each trio will have an individual folder, identified by the child ID.
 
The main output files are:

* out_hare/<child_name>/<child_name>.glnexus.family.combined_intersection_filtered_gq_<gq_value>_depth_<depth_value>_position.vcf
 
  * This file holds the *de novo* variants

* out_hare/<child_name>/<child_name>.glnexus.family.combined_intersection_filtered_gq_<gq_value>_depth_<depth_value>_position_all.vcf

  * This file holds the *de novo* variants specifically within CpG regions.
  
This current version has able to find almost all of the same de novo variants found from the original pipeline.  The NA12878 trio from the 1000 Genomes Project 30x WGS data is used as an example:

![NA12878](https://github.com/TNTurnerLab/Hare/blob/main/docs/compare_old_pipeline_to_1.1.png)
 
 ### Software Used:
* bcftools v1.11 
* python v3.9.7
* tabix v1.11 
* vcflib v1.0.0-rc0 
* bedtools v2.29.2 
* samtools v1.11 
* snakemake v7.15.2-0
* python v2.7
* GLnexus v1.4.1
* pytabix v0.1

