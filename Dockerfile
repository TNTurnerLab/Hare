FROM continuumio/miniconda3:4.12.0 as conda_setup 
RUN conda install -c bioconda bcftools=1.11 \
    && conda install -c bioconda tabix=1.11 \
    && conda install -c bioconda vcflib=1.0.0-rc0 \
    && conda install -c bioconda bedtools=2.29.2 \
    && conda install -c bioconda samtools=1.11 
RUN conda create -n snake -c bioconda -c conda-forge snakemake=7.15.2-0
RUN conda create -n py2 python=2.7
WORKDIR /
RUN wget https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli && chmod +x glnexus_cli
RUN apt-get update && apt-get install -y --no-install-recommends gcc libz-dev libopenblas-base && pip3 install pytabix

CMD ["/bin/bash"]
