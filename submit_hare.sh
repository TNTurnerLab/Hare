#!/bin/bash


export LSF_DOCKER_VOLUMES="/storage1/fs1/tychele/Active/testing/ris_testing/parabricks/HAT_1.1:/dnv_wf_cpu /storage1/fs1/tychele/Active/testing/reference:/reference /storage1/fs1/tychele/Active/testing/ris_testing/parabricks/test_4.0.0/dv:/dv /storage1/fs1/tychele/Active/testing/ris_testing/parabricks/test_4.0.0/gatk:/gatk /storage1/fs1/tychele/Active/testing/ris_testing/parabricks/glnexus_1kg_combinedtriotest/regions:/region  /scratch1/fs1/tychele:/scratch1/fs1/tychele $HOME:$HOME"

bsub -oo %J.log.main.snake.txt -m compute1-exec-356.ris.wustl.edu -q tychele  -g /jeffrey.ng/200test -R 'span[hosts=1] rusage[mem=5G]' -n 1 -M 5G  -G compute-tychele -a 'docker(jng2/testrepo2_actual:deep1.4)' /opt/conda/envs/snake/bin/snakemake -s /dnv_wf_cpu/hare_1.1.smk -j 3 --cluster-config cluster.json --cluster "bsub -q tychele -G compute-tychele -n {cluster.n} -oo %J.lsf.oo -R 'span[hosts=1] rusage[mem={cluster.mem}]' -g /jeffrey.ng/200test  -a 'docker(jng2/testrepo2_actual:deep1.4)' -M {cluster.mem}" -k --rerun-incomplete -w 120


#bsub -oo %J.log.main.snake.txt -q tychele -m compute1-exec-356.ris.wustl.edu  -g /jeffrey.ng/200test -R 'span[hosts=1] rusage[mem=5G]' -n 1 -M 5G  -G compute-tychele -a 'docker(jng2/testrepo2_actual:deep1.4)' /opt/conda/envs/snake/bin/snakemake -s /dnv_wf_cpu/tortoise_1.1.smk -j 10 --cluster-config cluster.json --cluster "bsub -q tychele -G compute-tychele -n {cluster.n} -oo %J.lsf.oo -R 'span[hosts=1] rusage[mem={cluster.mem}]' -g /jeffrey.ng/200test  -a 'docker(jng2/testrepo2_actual:deep1.4)' -M {cluster.mem}" -k --rerun-incomplete -w 120 --unlock
