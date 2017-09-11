# Bagnoli_2017
Code for analysis of single-cell RNA-seq data of Bagnoli et al., 2017

## Preprocessing
All scRNA-seq data was preprocessed with [zUMIs](https://github.com/sdparekh/zUMIs/) (Parekh et al., 2017).

The command was run as follows:

```
bash zUMIs-master.sh -f JM8.read1.fastq.gz -r JM8.read2.fastq.gz -n JM8 -g /data/ngs/genomes/Mouse/mm10/STAR5idx_ERCC_noGTF/ -a /data/ngs/genomes/Mouse/mm10/Mus_musculus.GRCm38.75.clean.spike.gtf -c 1-14 -m 15-24 -l 50 -z 2 -u 3 -p 16 -R yes -d 10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000
```

This resulted in the [JM8.rds object found in this repository](Data/JM8.rds).
