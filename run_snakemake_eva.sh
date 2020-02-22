mkdir -p autosnake/
snakemake -w 100 --jobname "_{rule}_{jobid}" -j 1000 \
    --cluster-config config/cluster.yaml \
    --resources io=100 \
    --cluster "qsub -r yes -q all.q -l h_vmem={cluster.mem} -o autosnake/ -e autosnake/ -cwd -V -S /bin/bash" $*  

