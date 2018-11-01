# TelescopeEncode

## Setup

### Data / Metadata

ENCODE data is obtained from the Sequence Read Archive.
The data analyzed comes from two projects: 

| Name           | GEO      | SRA       | Samples | Runs |
| -------------  | -------- | --------- | :-----: | :--: |
| ENCODE/CSHL    | GSE30567 | SRP007461 | 99      | [205](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP007461) |
| ENCODE/Caltech | GSE33480 | SRP014320 | 24      | [104](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014320) |


The RunInfo table for all SRA runs can be downloaded using the 
[SRA run selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP007461,SRP014320)
and is saved as [`metadata/SraRunTable.txt`](metadata/SraRunTable.txt). There are a total of 309 runs found that correspond to 123 samples. A [python script](scripts/parse_SraRunTable.py) is used to parse this run table and create three analysis sets. 

+ [**set1**](metadata/set1.tsv) ENCODE Tiers 1 and 2 common cell types, whole-cell bulk RNA-seq
+ [**set2**](metadata/set2.tsv) ENCODE Tiers 1 and 2 common cell types, subcellular RNA fractions
+ [**set3**](metadata/set3.tsv) Other cell types not in tiers 1 or 2


**ENCODE Tier 1 Cell Types:**

+ H1-hESC
+ GM12878
+ K562

**ENCODE Tier 2 Cell Types:**

+ HeLa-S3
+ HepG2
+ HUVEC
+ SK-N-SH
+ IMR90
+ A549
+ MCF-7
+ CD20+
+ Monocytes-CD14+
+ NHEK





```bash
curl -O -L 'https://www.encodeproject.org/metadata/type=Experiment&assay_slims=Transcription&assay_title=polyA+RNA-seq&assay_title=total+RNA-seq&award.project=ENCODE&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=fastq/metadata.tsv'
```




```bash
mkdir refs
mkdir 
```

polyA ENCODE paired end
Download metadata corresponding to polyA RNA-seq fastq files.
Only select runs with paired-end data and human

```
curl $(head -n1 metadata/files.paired_polyA.txt) > metadata/metadata.paired_polyA.txt
```

```python
lines = 
```




# refs

```
cd refs
ln -s /lustre/groups/cbi/shared/References/Homo_sapiens/UCSC/hg38
cd ..
```


GRCh38 XY reference genome (ENCODE3 used only one reference genome for analysis)

```
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
wget https://www.encodeproject.org/files/ENCFF908UQN/@@download/ENCFF908UQN.fasta.gz

```




GENCODE transcripts:

```
https://www.encodeproject.org/files/ENCFF824ZKD/@@download/ENCFF824ZKD.gtf.gz
```

Spike-in
```


```

```
http://labshare.cshl.edu/shares/mhammelllab/www-data/TEToolkit/TE_GTF/hg38_rmsk_TE.gtf.gz
```

```
ccmd='sbatch {cluster.args} -N {cluster.N} -p {cluster.p} -t {cluster.t} -o {cluster.out} -J {cluster.name}'
snakemake \
  -j 10000 \
  --local-cores 12 \
  -c "$ccmd" \
  -u cluster.yaml \
  --cluster-status 'scripts/slurm_status.py' \
  -T -k --ri \
  all


```


```
conda create -n TETools \
  python=3.6 future cython snakemake appdirs=1.4.3 requests readline \
  biopython numpy scipy sra-tools hisat2 bowtie2 bowtie samtools picard bedtools \
  stringtie HTSeq htslib=1.5 pysam=0.12.0.1 parallel-fastq-dump \
  fastqc multiqc docopt pandas

source activate TETools

pip install git+ssh://git@github.com/mlbendall/telescope.git@master



conda create -n testrepenrich python=2.7.3 biopython bowtie bedtools=2.20 samtools=0.1.19
conda update -n base conda
```





```
cp metadata/sampletable_long.tsv results
cp metadata/set1.txt results
cp refs/HERV_rmsk.hg38.v2.gtf

cat cur_samples.txt | while read s; do
    cp samples/$s/inform-telescope_report.tsv results/$s.inform-telescope_report.tsv
    cp samples/$s/bt2_multi.summary.txt results/$s.bt2_multi.summary.txt
done
cp cur_samples.txt results/samples.txt
```
    
    
    
for f in samples/*/inform-telescope_report.tsv; do
    cp $f results/$(basename $(dirname $f)).inform-telescope_report.tsv
    
    cp $f results/$(basename $(dirname $f)).inform-telescope_report.tsv
    
done

```