
## Invoking Snakemake run

/mnt/storage/apps/miniconda3/bin/snakemake -s path/to/mutect2.snakefile --cluster-config path/to/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete
