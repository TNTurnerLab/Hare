version 1.0

struct cram_pairs {
    File cram
    File crai
}
struct family {
    String father
    String mother
    String child
}

workflow jumping_hare {
    input {
        Array[Array[cram_pairs]] cram_files
        Array[family] trios
        File pathToReference
        String reference
        String basename_split=".cram"
        String deep_docker="nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1"
        String typeOfGPU_DV="nvidia-tesla-t4"
        String gpuDriverVersion_DV = "460.73.01"
        String typeOfGPU_HC="nvidia-tesla-t4"
        String gpuDriverVersion_HC = "460.73.01"
        String hare_docker="tnturnerlab/hare:v1.1"
        File naive_inheritance_trio_py2
        File test_intersect
        File filter_glnexuscombined_updated
        String deep_model='shortread'
        Boolean wes=false
        String glnexus_deep_model='DeepVariant'
        
        Int cpu_dv=24
        Int num_ram_dv=120
        Int extra_mem_dv=65
        Int glnexus_cpu=32
        Int glnexus_ram_dv=100
        Int num_gpu_dv=4

        Int num_gpu_HC=2
        Int cpu_hc=24
        Int num_ram_hc=120
        Int extra_mem_hc=65
        Int glnexus_ram_hc=100
        String sample_suffix
        Int maxPreemptAttempts=3
        File? chrom_length
        File? regions
        String interval_file="None"
        Int gq=20
        Int depth=10
    }

    Int region_check = if !defined(regions) then 1 else 0
    scatter (i in range (length(cram_files)))
    {
        Array[cram_pairs] pairs = cram_files[i]
        
        scatter(c_pair in pairs) {
            call deepVariant {
                input:
                    deep_docker=deep_docker,
                    member=c_pair,
                    reference=reference,
                    basename_split=basename_split,
                    deep_model=deep_model,
                    pathToReference=pathToReference,
                    nThreads = cpu_dv,
                    gbRAM = num_ram_dv,
                    gpuDriverVersion=gpuDriverVersion_DV,
                    typeOfGPU=typeOfGPU_DV,
                    num_gpu_dv=num_gpu_dv,
                    extra_mem_dv = extra_mem_dv,
                    maxPreemptAttempts = maxPreemptAttempts,
                    suffix=sample_suffix,
                    interval_file=interval_file,
                    wes=wes
        }
            call gatk {
                input:
                    gatk_docker=deep_docker,
                    member=c_pair,
                    reference=reference,
                    basename_split=basename_split,
                    pathToReference=pathToReference,
                    nThreads = cpu_hc,
                    gbRAM = num_ram_hc,
                    extra_mem_dv = extra_mem_hc,
                    maxPreemptAttempts = maxPreemptAttempts,
                    suffix=sample_suffix,
                    gpuDriverVersion=gpuDriverVersion_HC,
                    typeOfGPU=typeOfGPU_HC,
                    num_gpu_HC=num_gpu_HC,
                    interval_file=interval_file,
                    wes=wes
            }

        }
        call glnexus_DV {
            input:
                fam=trios[i],
                glnexus_cpu=glnexus_cpu,
                glnexus_ram=glnexus_ram_dv,
                hare_docker=hare_docker,
                maxPreemptAttempts=maxPreemptAttempts,
                deepout=deepVariant.deep_out,
                glnexus_deep_model=glnexus_deep_model,
                naive_inheritance_trio_py2=naive_inheritance_trio_py2
        }   
        call glnexus_HC {
            input:
                fam=trios[i],
                glnexus_cpu=glnexus_cpu,
                glnexus_ram=glnexus_ram_hc,
                hare_docker=hare_docker,
                maxPreemptAttempts=maxPreemptAttempts,
                gatkout=gatk.gatk_out,
                naive_inheritance_trio_py2=naive_inheritance_trio_py2
        }   
        call combinedAndFilter{
            input:
            fam=trios[i],
            chrom_length=chrom_length,
            pathToReference=pathToReference,
            reference=reference,
            hare_docker=hare_docker,
            maxPreemptAttempts=maxPreemptAttempts,
            test_intersect=test_intersect,
            filter_glnexuscombined_updated=filter_glnexuscombined_updated,
            region_check=region_check,
            regions=regions,
            depth=depth,
            gq=gq,
            dv_vcf_out=glnexus_DV.dnv_vcf,
            dv_vcf_out_tbi=glnexus_DV.dnv_vcf_tbi,
            hc_vcf_out=glnexus_HC.dnv_vcf,
            hc_vcf_out_tbi=glnexus_HC.dnv_vcf_tbi,
            dv_vcf_out_full=glnexus_DV.glenexus_vcf_dv,
            hc_vcf_out_full=glnexus_HC.glenexus_vcf_hc,
            dv_vcf_out_full_tbi=glnexus_DV.out_vcf_dv_tbi,
            hc_vcf_out_full_tbi=glnexus_HC.out_vcf_hc_tbi


        }
    }
}
task deepVariant {
    input{
    cram_pairs member
    String deep_docker
    File pathToReference
    Int nThreads 
    String reference
    Int gbRAM 
    String gpuDriverVersion
    String deep_model
    String typeOfGPU
    Int? extra_mem_dv 
    Int maxPreemptAttempts
    String basename_split
    String suffix
    Int num_gpu_dv
    String interval_file
    
    Boolean wes
    }
    String wes_mode= if wes then "--use-wes-model" else ""
    String region= if wes then "--interval-file " +interval_file else ""
    Int tmpsize=select_first([extra_mem_dv,65])
    Int disk_size = ceil(size(member.cram,"GB")+size(pathToReference,"GB")+tmpsize)
    String file=basename(member.cram,basename_split)
    String ref_base=basename(pathToReference)
    String name=basename(file,suffix)
    command <<<
    
    set -euo pipefail
    tar -xf ~{pathToReference}
    pbrun deepvariant --num-gpus ~{num_gpu_dv} --ref ~{reference}  \
    --in-bam ~{member.cram} \
    --out-variants ~{name}.dv.g.vcf.gz \
    --mode ~{deep_model} \
    --num-gpus ~{num_gpu_dv} \
    --gvcf  ~{wes_mode} ~{region} 


    >>>

    runtime {
        docker : "~{deep_docker}"
        disks : "local-disk ~{disk_size} SSD"
        cpu : nThreads
        gpuType : "~{typeOfGPU}"
        gpuCount : "~{num_gpu_dv}"
        gpu: "num=~{num_gpu_dv}:gmodel=~{typeOfGPU}:j_exclusive=yes"
        memory : "~{gbRAM}GB"
        resource: "gpuhost span[hosts=1]"
        nvidiaDriverVersion : "~{gpuDriverVersion}"
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }

    output {
        File deep_out="~{name}.dv.g.vcf.gz"
    }
}


task gatk {
    input{
    cram_pairs member
    String gatk_docker
    File pathToReference
    Int nThreads 
    String reference
    Int gbRAM 
    Int? extra_mem_dv 
    Int maxPreemptAttempts
    String basename_split
    String typeOfGPU
    String suffix
    String gpuDriverVersion
    Int num_gpu_HC
    String interval_file
    Boolean wes
    }
    
    String region= if wes then "--interval-file " +interval_file else ""
    Int tmpsize=select_first([extra_mem_dv,65])
    Int disk_size = ceil(size(member.cram,"GB")+size(pathToReference,"GB")+tmpsize)
    String file=basename(member.cram,basename_split)
    String name=basename(file,suffix)
    command <<<
    
    set -euo pipefail
    tar -xf ~{pathToReference}
    pbrun haplotypecaller --num-gpus ~{num_gpu_HC} --ref ~{reference} --in-bam ~{member.cram} --out-variants ~{name}.hc.g.vcf.gz --haplotypecaller-options "-standard-min-confidence-threshold-for-calling 30" --gvcf ~{region} 

    >>>

    runtime {
        docker : "~{gatk_docker}"
        disks : "local-disk ~{disk_size} SSD"
        cpu : nThreads
        memory : "~{gbRAM}GB"
        gpuType : "~{typeOfGPU}"
        gpuCount : "~{num_gpu_HC}"
        gpu: "num=~{num_gpu_HC}:gmodel=~{typeOfGPU}:j_exclusive=yes"
        resource: "gpuhost span[hosts=1]"
        nvidiaDriverVersion : "~{gpuDriverVersion}"
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }

    output {
        File gatk_out="~{name}.hc.g.vcf.gz"
    }
}

task glnexus_DV {
    input{
        family fam
        String hare_docker
        Array[File] deepout
        Int? extramem_GLDV
        Int glnexus_cpu
        Int glnexus_ram
        Int maxPreemptAttempts
        File naive_inheritance_trio_py2
        String glnexus_deep_model
    }

    
    Int tmpsize=select_first([extramem_GLDV,50])
   
    Int disk_size = ceil(size(deepout,"GB")+tmpsize)
    command <<<
    set -euo pipefail
    export PATH=/opt/conda/bin:$PATH
    for f in ~{sep=' ' deepout}
    do
    echo $f
    /opt/conda/bin/tabix $f
    done
    
    /glnexus_cli --config ~{glnexus_deep_model} --mem-gbytes ~{glnexus_ram} --dir ~{fam.child}_glnexus_dv ~{sep=' ' deepout} > ~{fam.child}.trio.glnexus.dv.bcf
    rm -r  ~{fam.child}_glnexus_dv 
    /opt/conda/bin/bcftools view ~{fam.child}.trio.glnexus.dv.bcf | bgzip -c > ~{fam.child}.trio.dv.vcf.gz
    /opt/conda/bin/tabix ~{fam.child}.trio.dv.vcf.gz
    export PATH=/opt/conda/envs/py2/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib/openblas-base/
    trio="~{fam.father},~{fam.mother},~{fam.child}"
    bcftools view -Oz --threads 4 -s $trio ~{fam.child}.trio.dv.vcf.gz  > ~{fam.child}.glnexus.dv.vcf.gz

    tabix ~{fam.child}.glnexus.dv.vcf.gz


    echo "Naive search DV"
    #Searches family trio .vcf file for possible denovo varints, adding the label TRANSMITTED=no;INH=denovo_pro if the variant is a candidate.
    #It then pulls out the varints with that label into a separate .vcf file.

    python2 ~{naive_inheritance_trio_py2}  ~{fam.child}.glnexus.dv.vcf.gz > ~{fam.child}.glnexus_denovo.dv.vcf
    grep '#'  ~{fam.child}.glnexus_denovo.dv.vcf > ~{fam.child}.glnexus_denovo_actual.dv.vcf
    grep -w 'TRANSMITTED=no;INH=denovo_pro' ~{fam.child}.glnexus_denovo.dv.vcf >> ~{fam.child}.glnexus_denovo_actual.dv.vcf
    bgzip ~{fam.child}.glnexus_denovo_actual.dv.vcf
    tabix ~{fam.child}.glnexus_denovo_actual.dv.vcf.gz
    >>>
    runtime {
        docker : "~{hare_docker}"
        disks : "local-disk ~{disk_size} SSD"
        cpu : glnexus_cpu
        memory : "~{glnexus_ram}GB"
        preemptible : maxPreemptAttempts
    }

    output{
        File glenexus_vcf_dv="~{fam.child}.trio.dv.vcf.gz"
        File out_vcf_dv_tbi="~{fam.child}.trio.dv.vcf.gz.tbi"
        File dnv_vcf="~{fam.child}.glnexus_denovo_actual.dv.vcf.gz"
        File dnv_vcf_tbi="~{fam.child}.glnexus_denovo_actual.dv.vcf.gz.tbi"

    }
}


task glnexus_HC {
    input{
        family fam
        String hare_docker
        Array[File] gatkout
        Int? extramem_GLDV
        Int glnexus_cpu
        Int glnexus_ram
        Int maxPreemptAttempts
        File naive_inheritance_trio_py2
    }

    
    Int tmpsize=select_first([extramem_GLDV,50])
   
    Int disk_size = ceil(size(gatkout,"GB")+tmpsize)
    command <<<
    set -euo pipefail
    export PATH=/opt/conda/bin:$PATH
    for f in ~{sep=' ' gatkout}
    do
    echo $f
    /opt/conda/bin/tabix $f
    done
    
    /glnexus_cli --config gatk --mem-gbytes ~{glnexus_ram} --dir ~{fam.child}_glnexus_hc ~{sep=' ' gatkout} > ~{fam.child}.trio.glnexus.hc.bcf
    rm -r  ~{fam.child}_glnexus_hc
    /opt/conda/bin/bcftools view ~{fam.child}.trio.glnexus.hc.bcf | bgzip -c > ~{fam.child}.trio.hc.vcf.gz
    /opt/conda/bin/tabix ~{fam.child}.trio.hc.vcf.gz
    export PATH=/opt/conda/envs/py2/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib/openblas-base/
    trio="~{fam.father},~{fam.mother},~{fam.child}"
    bcftools view -Oz --threads 4 -s $trio ~{fam.child}.trio.hc.vcf.gz  > ~{fam.child}.glnexus.hc.vcf.gz

    tabix ~{fam.child}.glnexus.hc.vcf.gz


    echo "Naive search hc"
    #Searches family trio .vcf file for possible denovo varints, adding the label TRANSMITTED=no;INH=denovo_pro if the variant is a candidate.
    #It then pulls out the varints with that label into a separate .vcf file.

    python2 ~{naive_inheritance_trio_py2}  ~{fam.child}.glnexus.hc.vcf.gz > ~{fam.child}.glnexus_denovo.hc.vcf
    grep '#'  ~{fam.child}.glnexus_denovo.hc.vcf > ~{fam.child}.glnexus_denovo_actual.hc.vcf
    grep -w 'TRANSMITTED=no;INH=denovo_pro' ~{fam.child}.glnexus_denovo.hc.vcf >> ~{fam.child}.glnexus_denovo_actual.hc.vcf
    bgzip ~{fam.child}.glnexus_denovo_actual.hc.vcf
    tabix ~{fam.child}.glnexus_denovo_actual.hc.vcf.gz
    >>>
    runtime {
        docker : "~{hare_docker}"
        disks : "local-disk ~{disk_size} SSD"
        cpu : glnexus_cpu
        memory : "~{glnexus_ram}GB"
        preemptible : maxPreemptAttempts
    }

    output{
        File glenexus_vcf_hc="~{fam.child}.trio.hc.vcf.gz"
        File out_vcf_hc_tbi="~{fam.child}.trio.hc.vcf.gz.tbi"
        File dnv_vcf="~{fam.child}.glnexus_denovo_actual.hc.vcf.gz"
        File dnv_vcf_tbi="~{fam.child}.glnexus_denovo_actual.hc.vcf.gz.tbi"

    }
}

task combinedAndFilter {
    input{
        family fam
        File dv_vcf_out_full
        File hc_vcf_out_full
        File dv_vcf_out_full_tbi
        File hc_vcf_out_full_tbi
        File dv_vcf_out
        File hc_vcf_out
        File dv_vcf_out_tbi
        File hc_vcf_out_tbi
        File test_intersect
        File filter_glnexuscombined_updated
        Int? extramem_GLDV
        File pathToReference
        Int maxPreemptAttempts
        String hare_docker
        File? regions
        String reference
        Int region_check
        Int gq
        Int depth
        File? chrom_length
        

    }

    Int chrom_check=if !defined(chrom_length) then 0 else 1

    
    Int tmpsize=select_first([extramem_GLDV,50])
   
    Int disk_size = ceil(size(dv_vcf_out,"GB")+size(hc_vcf_out,"GB")+tmpsize+size(dv_vcf_out_full,"GB")+size(hc_vcf_out_full,"GB"))
    command <<<
    set -euo pipefail
    #Filters out -L chromosomes that are not -L chr1-22, grabbing variants with an allele count of 1, variants that were found from both DV and HC (defined as set=Intersection), and removes variants that are either 10 A's or T's in a row.
    echo "Combining files"
    tar -xf ~{pathToReference}
    mv ~{hc_vcf_out}* .
    mv ~{dv_vcf_out}* .
    ref=$( head -n 1 ~{reference} | cut -d' ' -f 1)
    if [ "~{chrom_check}" == "1" ]
    then
    /opt/conda/bin/python ~{test_intersect} -g ~{fam.child}.glnexus_denovo_actual.hc.vcf.gz -d ~{fam.child}.glnexus_denovo_actual.dv.vcf.gz -r $ref -c ~{chrom_length}
    else
    /opt/conda/bin/python ~{test_intersect} -g ~{fam.child}.glnexus_denovo_actual.hc.vcf.gz -d ~{fam.child}.glnexus_denovo_actual.dv.vcf.gz -r $ref -c None
    fi
    cat ~{fam.child}_combined_out.vcf |  awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ~{fam.child}.glnexus_denovo_actual.combined.vcf
    zcat ~{fam.child}.glnexus_denovo_actual.dv.vcf.gz  | grep '#' > ~{fam.child}.glnexus.family.combined_intersection.vcf
    grep -v 'chrUn' ~{fam.child}.glnexus_denovo_actual.combined.vcf | grep -v '_random' | grep -v '_alt'  | grep -v 'chrY' | grep -v 'chrM' | grep  'AC=1' |  egrep -v 'AAAAAAAAAA|TTTTTTTTTT' >> ~{fam.child}.glnexus.family.combined_intersection.vcf


    #Finds the order of the family position found within the combined .vcf file and what the order is in the family file (which should be in order of father, mother, child)
    echo 'Set up filter script'
    actual=$( cat ~{fam.child}.glnexus.family.combined_intersection.vcf | grep '#' | tail -n 1 | cut -f 10-13 )
    echo $actual
    trio="~{fam.father} ~{fam.mother} ~{fam.child}"
    echo $trio
    echo "Run filter script"
    #Python script that filters for parents with no alt allele, depth of set value, GQ of set value, and allele balance of .25
    /opt/conda/bin/python ~{filter_glnexuscombined_updated} ~{fam.child}.glnexus.family.combined_intersection.vcf $trio $actual ~{gq} ~{depth}

    echo "Make position file"
    cat ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_20_depth_10.vcf | grep '#' > ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}_position.vcf
    echo "Make CpG file"
    cat ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_20_depth_10.vcf | grep '#' > ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}_position_all.vcf
    echo "Filter by position"

    if [ "~{region_check}" != "1" ]
    then
        
        tar -xf ~{regions}
      
        /opt/conda/bin/bedtools intersect -v -a ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}.vcf -b regions/LCR-hs38-5bp-buffer.bed.gz    regions/hg38_centromeres_09252018.bed.gz  regions/recent_repeat_b38-5bp-buffer.bed.gz  >> ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}_position.vcf
        /opt/conda/bin/bedtools intersect -a ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}_position.vcf -b regions/CpG_sites_sorted_b38.bed.gz >> ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}_position_all.vcf
    fi
    tar -cjf ~{fam.child}.hat.output.tar.bz2 ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}_position.vcf ~{fam.child}.glnexus.family.combined_intersection_filtered_gq_~{gq}_depth_~{depth}_position_all.vcf ~{dv_vcf_out_full} ~{hc_vcf_out_full} ~{dv_vcf_out_full_tbi} ~{hc_vcf_out_full_tbi}
    >>>
    runtime {
        docker : "~{hare_docker}"
        disks : "local-disk ~{disk_size} SSD"
        cpu : 1
        memory : "16GB"
        preemptible : maxPreemptAttempts
    }
    
    output{
        File dnv_all="~{fam.child}.hat.output.tar.bz2"
    }
}
