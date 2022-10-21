# Hare
#### Developed and maintained by Jeffrey K. Ng
#### Maintainer Email:  jeffrey.ng@wustl.edu
#### Washington University in St. Louis Medical School
#### Tychele N. Turner, Ph.D., Lab


A de novo variant caller leveraging Parabricks GPU accelerated variant calling.  This is an updated version of the [workflow described in Ng et al. 2022.](https://doi.org/10.1002/humu.24455).  The original code of the workflow can be [found here.](https://github.com/TNTurnerLab/GPU_accelerated_de_novo_workflow)  This version has been tested with the output from Parabricks v.3.0.0 as well as the now free version of [Parabricks v4.0.0.0-1](https://docs.nvidia.com/clara/parabricks/4.0.0/index.html).

# How to Run
## Input

Two main inputs:
1) The two input for this workflow are the output g.vcf files from Parabricks accelerated from GATK and DeepVariant. 
2)  A comma delimited text file, with one trio per line, with sample IDs formatted in the following way:  Father,Mother,Child

If you would like to download the RepeatMasker files, please use the following links:
 
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

#### Setting up the config.json file
Before running, please make any necessary changes to these options below in the config.json. 
 
* regions:  "/region" *If you don't have the RepeatMasker files, please make this entry blank*
* gq_value:  20 *gq value filter*
* depth_value: 10 *Depth value filter*
* suffix_dv:  *Suffix of the DeepVariant data files.  Assumes input files are \<sample\_name\>\<suffix\>* 
* suffix_hc:  *Suffix of the GATK Haplotypecaller data files.  Assumes input files are \<sample\_name\>\<suffix\>* 

## Running
 
While the use of docker is highly recommended, but the workflow is able to run outside of a docker envionment, please see the software dependancies below.  If you like to run Parabricks 4.0.0 having a system able to run Docker **is** a requirement.

The workflow Docker image can be pulled from here:
```
tnturnerlab/hare:v1.1
```
We also provide the Dockerfile if you would like to make modifications.  

```
docker run -v "/path/to/hare/code:/dnv_wf_cpu" -v "/path/to/reference:/reference" -v "/path/to/deepvariant/output:/dv" -v "/path/to/gatk/output:/gatk"  -v "/path/to/RepeatMasker/region/files:region" tnturnerlab/hare:v1.1 /opt/conda/envs/snake/bin/snakemake -s /dnv_wf_cpu/hare_1.1.smk -j 6 --cores -k --rerun-incomplete -w 120 
```
 # Output
 Below is a brief descirption of the main output folders from Hare:
 * dv_vcf:  Output folder for converted GLnexus .bcf files to .vcf.gz files for DeepVariant output.  Also includes tabix index file.
 * hc_vcf:  Output folder for converted GLnexus .bcf files to .vcf.gz files for GATK HaplotypeCaller output.  Also includes tabix index file.
 * dv_bcf: Output folder for GLnexus .bcf files for DeepVariant output.
 * hc_bcf: Output folder for GLnexus .bcf files for GATK HaplotypeCaller output.
 * out_hare:  Output folder where the *de novo* variant files can be found.  If you are running multiple trios, each trio will have an individual folder, identified by the child ID.
 
The main output files are:

* out_hare/<child_name>/<child_name>.glnexus.family.combined_intersection_filtered_gq_<gq_value>_depth_<depth_value>_position.vcf
 
  * This file holds the *de novo* variants

* out_hare/<child_name>/<child_name>.glnexus.family.combined_intersection_filtered_gq_<gq_value>_depth_<depth_value>_position_all.vcf

  * This file holds the *de novo* variants specifically within CpG regions.
  
This current version has able to find almost all of the same de novo variants found from the original pipeline.  The NA12878 trio from the 1000 Genomes Collection 30x WGS data is used as an example:

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

