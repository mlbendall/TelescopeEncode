#! /usr/bin/env python
from future import standard_library
standard_library.install_aliases()

from os import path

import snakemake
from snakemake.exceptions import RuleException
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

""" Functions """
def dict_args(d):
    """ Convert a dictionary to command line arguments"""
    ret = ''
    for k,v in d.items():
        if len(k) == 1:
            ret += '-{} {} '.format(k,v)
        else:
            ret += '--{} {} '.format(k,v)
    return ret

def read_value(fn):
    """ Read a value from text file.
        
        The file should contain a single word indicating the value. For example, the 
        library layout of a run could be stored in a file named 'LAYOUT.txt', and the 
        value should be either 'PAIRED' or 'SINGLE'.
    """
    with open(fn, 'rU') as fh:
        return fh.read().strip()

configfile: "config.yaml"

wildcard_constraints:
    sampid="GSM\d+",
    sraid="[SED]RR\d+",
    sraproj="[SED]RP\d+",



''' Constants '''
SAMP_TABLE = [l.strip('\n').split('\t') for l in open(config['sample_table'], 'rU')][1:]
SAMPLES = [l.strip() for l in open('cur_samples.txt', 'rU')] #  [r[0] for r in SAMP_TABLE]
RUNS = {r[0]:r[1].split(',') for r in SAMP_TABLE}

# print(RUNS)

localrules: all, complete_sample, make_sampdir, make_rundir

rule all:
    input:
        expand("samples/{s}/completed.txt", s=SAMPLES)

rule complete_sample:
    input:
        "samples/{sampid}/bt2_multi.unsorted.bam",
        "samples/{sampid}/bt2_multi.summary.txt",       
        "samples/{sampid}/inform-telescope_report.tsv"
    output:
        touch("samples/{sampid}/completed.txt")
        
rule make_sampdir:
    output:
        "samples/{sampid}/"
    shell:
        "mkdir -p {output}"

# rule make_rundir:
#     output:
#         "runs/{sraid}/"
#     shell:
#         "mkdir -p {output}"

''' References '''
localrules: download_references, download_telescope_annotation

rule download_references:
    input:
        ancient("refs/")
    output:
        "refs/download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
        "refs/download/ENCFF908UQN.fasta.gz",
        "refs/download/ENCFF824ZKD.gtf.gz"
    shell:
        "mkdir -p refs/download && "
        "wget -O {output[0]} https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz && "
        "wget -O {output[1]} https://www.encodeproject.org/files/ENCFF908UQN/@@download/ENCFF908UQN.fasta.gz && "
        "wget -O {output[2]} https://www.encodeproject.org/files/ENCFF824ZKD/@@download/ENCFF824ZKD.gtf.gz"


rule download_telescope_annotation:
    output:
        config['herv_annotation']
    shell:
        "wget -O {output} https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf "


rule build_bowtie2_index:
    input:
        genome = "refs/download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
        ctrls = "refs/download/ENCFF908UQN.fasta.gz"
    output:
        expand(config['bt2idx'] + ".{i}.bt2", i = range(1, 5)),
        expand(config['bt2idx'] + ".rev.{i}.bt2", i = range(1, 3))
    threads: snakemake.utils.available_cpu_count()        
    shell:
        "tmpidx=$(mktemp -d) && tmpfasta=$(mktemp) && "
        "cat {input.genome} {input.ctrls} | gunzip > $tmpfasta && "
        "bowtie2-build --threads {threads} -f $tmpfasta $tmpidx/$(basename {config[bt2idx]}) && "
        "mkdir -p $(dirname {config[bt2idx]}) && "
        "mv $tmpidx/* $(dirname {config[bt2idx]})"

rule build_bowtie_index:
    input:
        genome = "refs/download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
        ctrls = "refs/download/ENCFF908UQN.fasta.gz"
    output:
        expand(config['btidx'] + ".{i}.ebwt", i = range(1, 5)),
        expand(config['btidx'] + ".rev.{i}.ebwt", i = range(1, 3))
    threads: snakemake.utils.available_cpu_count()        
    shell:
        "tmpidx=$(mktemp -d) && tmpfasta=$(mktemp) && "
        "cat {input.genome} {input.ctrls} | gunzip > $tmpfasta && "
        "bowtie-build -f $tmpfasta $tmpidx/$(basename {config[btidx]}) && "
        "mkdir -p $(dirname {config[btidx]}) && "
        "mv $tmpidx/* $(dirname {config[btidx]})"


''' Run-level download and preprocessing '''

def sra_url(sraid):
    base = "ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra"
    return '/'.join([base, sraid[:3], sraid[:6], sraid, '%s.sra' % sraid])


rule rundata:
    output:
        sra  = "runs/{sraid}/{sraid}.sra",
        r1   = "runs/{sraid}/reads_1.fastq.gz",
        r2   = "runs/{sraid}/reads_2.fastq.gz",
        lout = "runs/{sraid}/LAYOUT.txt"
    params:
        url = lambda x: sra_url(x.sraid),
    threads: snakemake.utils.available_cpu_count()
    run:
        # Create run directory
        # Download SRA file
        # Dump SRA to FASTQ
        shell(
            "mkdir -p runs/{wildcards.sraid} && "
            "curl -o {output.sra} {params.url} && "
            "parallel-fastq-dump "
            "-t {threads} "
            "--tmpdir $TMPDIR "
            "-O runs/{wildcards.sraid} "
            "--gzip "
            "--split-files --skip-technical "
            " -F -B -Q 33 --defline-qual '+' "
            "-s {output.sra} "
        )
        
        # Rename FASTQ to reads_1.fastq.gz and reads_2.fastq.gz (if paired)
        prefix = 'runs/%s/%s' % (wildcards.sraid, wildcards.sraid)
        f1, f2, f3 = [path.exists('%s_%d.fastq.gz' % (prefix, i)) for i in [1,2,3]]
        if not f1:
            raise RuleException("Expected %s_1.fastq.gz." % (prefix))
        if f2 and f3:    
            raise RuleException("Too many reads dumped." % (prefix))
        elif f2:
            shell(
                "mv runs/{wildcards.sraid}/{wildcards.sraid}_1.fastq.gz {output.r1} && "
                "mv runs/{wildcards.sraid}/{wildcards.sraid}_2.fastq.gz {output.r2} && "
                "echo 'PAIRED' > {output.lout} "
            )
        elif f3:
            shell(
                "mv runs/{wildcards.sraid}/{wildcards.sraid}_1.fastq.gz {output.r1} && "
                "mv runs/{wildcards.sraid}/{wildcards.sraid}_3.fastq.gz {output.r2} && "
                "echo 'PAIRED' > {output.lout} "
            )
        else:
            shell(
                "mv runs/{wildcards.sraid}/{wildcards.sraid}_1.fastq.gz {output.r1} && "
                "touch {output.r2} && "
                "echo 'SINGLE' > {output.lout} "
            )


rule flexbar:
    input:
        "runs/{sraid}/reads_1.fastq.gz",
        "runs/{sraid}/reads_2.fastq.gz",
        "runs/{sraid}/LAYOUT.txt"
    output:
        "runs/{sraid}/clean_1.fastq.gz",
        "runs/{sraid}/clean_2.fastq.gz",
        "runs/{sraid}/clean.log"
    threads: snakemake.utils.available_cpu_count()
    run:
        # Parse flexbar arguments from config.yaml
        flexbar_args = dict_args(config['flexbar'])
        
        # Setup adapters and inputs depending on library layout
        layout = read_value(input[2])
        if layout == 'SINGLE':
            adapters = config['adapters']['SE']
            rargs = '--reads {}'.format(input[0])
        else:
            adapters = config['adapters']['PE']
            rargs = '--reads {} --reads2 {}'.format(input[0], input[1])
            
        # Run flexbar
        shell(
            "module load flexbar/3.0.3 && "
            "flexbar "
            "--threads {threads} "
            "{flexbar_args} "
            "--adapters {adapters} "
            "{rargs} "
            "--target runs/{wildcards.sraid}/clean "
        )
        # Adjust outputs for single-end reads
        if layout == 'SINGLE':
            shell(
                "mv runs/{wildcards.sraid}/clean.fastq.gz {output[0]} && "
                "touch {output[1]} "
            )


''' Run-level alignment '''
rule bowtie2:
    input:
        "runs/{sraid}/clean_1.fastq.gz",
        "runs/{sraid}/clean_2.fastq.gz",
        "runs/{sraid}/LAYOUT.txt",          
        expand(config['bt2idx'] + ".{i}.bt2", i = range(1, 5)),
        expand(config['bt2idx'] + ".rev.{i}.bt2", i = range(1, 3))
    output:
        "runs/{sraid}/bt2_{preset, \w+}.unsorted.bam",
        "runs/{sraid}/bt2_{preset, \w+}.summary.txt"
    threads: snakemake.utils.available_cpu_count()
    run:
        # Parse bowtie2 arguments from config.yaml based on "preset"
        if wildcards.preset in config['bowtie2']:
            bowtie2_args = dict_args(config['bowtie2'][wildcards.preset])
        else:
            bowtie2_args = ""

        # Setup inputs depending on library layout
        layout = read_value(input[2])
        if layout == 'SINGLE':
            rargs = '-U {}'.format(input[0])
        else:
            rargs = '-1 {} -2 {}'.format(input[0], input[1])
        
        # Run bowtie2
        shell(
            "(bowtie2 "
            "-p {threads} "
            "{bowtie2_args} "
            "--rg-id {wildcards.sraid} "
            "-x {config[bt2idx]} "
            "{rargs} | "
            "samtools view -b > {output[0]}"
            ") 3>&1 1>&2 2>&3 | tee {output[1]} "
        )


''' Sample-level processing '''
localrules: concat_alignments, set_sample_layout, combine_bt2_summary, concat_fastq

rule concat_alignments:
    input:
        bams = lambda wildcards: expand("runs/{r}/bt2_{p}.unsorted.bam", 
                                         p=wildcards.preset, r=RUNS[wildcards.sampid]),
        logs = lambda wildcards: expand("runs/{r}/bt2_{p}.summary.txt", 
                                         p=wildcards.preset, r=RUNS[wildcards.sampid])                                     
    output:
        "samples/{sampid}/bt2_{preset, \w+}.unsorted.bam"
    shell:
        "samtools cat -o {output} {input.bams}"


rule set_sample_layout:
    input:
        lout = lambda wildcards: expand("runs/{r}/LAYOUT.txt", r=RUNS[wildcards.sampid])
    output:
        "samples/{sampid}/LAYOUT.txt"
    run:
        # Check library layout
        layouts = [read_value(f) for f in input.lout]
        if layouts[0] == 'SINGLE':
            if not all(l == 'SINGLE' for l in layouts):
                raise RuleException("Layouts do not match.")
            shell("echo 'SINGLE' > {output}")
        else:
            if not all(l == 'PAIRED' for l in layouts):
                raise RuleException("Layouts do not match.")        
            shell("echo 'PAIRED' > {output}")


rule combine_bt2_summary:
    input:
        logs = lambda wildcards: expand("runs/{r}/bt2_{p}.summary.txt", 
                                         p=wildcards.preset, r=RUNS[wildcards.sampid]),
        lout = "samples/{sampid}/LAYOUT.txt"                                     
    output:
        "samples/{sampid}/bt2_{preset, \w+}.summary.txt"
    run:
        layout = read_value(input.lout)
        shell(
            "python scripts/combine_bowtie2_logs.py --layout {layout} {input.logs} > {output}"
        )

rule concat_fastq:
    input:
        r1 =   lambda wildcards: expand("runs/{r}/clean_1.fastq.gz", 
                                         r=RUNS[wildcards.sampid]),
        r2 =   lambda wildcards: expand("runs/{r}/clean_2.fastq.gz", 
                                         r=RUNS[wildcards.sampid]),
        lout = "samples/{sampid}/LAYOUT.txt"
    output:
        "samples/{sampid}/clean_1.fastq.gz",
        "samples/{sampid}/clean_2.fastq.gz",
    run:
        layout = read_value(input.lout)
        if layout == 'SINGLE':
            shell(
                "cat {input.r1} > {output[0]} && "
                "touch {output[1]}"
            )
        else:
            shell(
                "cat {input.r1} > {output[0]} && "
                "cat {input.r2} > {output[1]}"
            )



''' Telescope '''
rule telescope:
    input:
        "samples/{sampid}/bt2_multi.unsorted.bam",
        config['herv_annotation']
    output:
        "samples/{sampid}/{preset,\w+}-telescope_report.tsv",
        "samples/{sampid}/{preset,\w+}.log"
    run:
        if wildcards.preset in config['telescope']:
            telescope_args = dict_args(config['telescope'][wildcards.preset])
        else:
            telescope_args = ""
        shell(
            "telescope "
            "assign "
            "{telescope_args} "            
            "--outdir $(dirname {output[0]}) "
            "{input[0]} "
            "{input[1]} "
            "2>&1 | tee {output[1]}"
        )

''' Unique counts '''


''' Best counts '''

''' RepEnrich '''
# rule bowtie:
#     input:
#         "runs/{sraid}/clean_1.fastq.gz",
#         "runs/{sraid}/clean_2.fastq.gz",
#         "runs/{sraid}/LAYOUT.txt",          
#         expand(config['bt2idx'] + ".{i}.bt2", i = range(1, 5)),
#         expand(config['bt2idx'] + ".rev.{i}.bt2", i = range(1, 3))
#     output:
#         "runs/{sraid}/bt_{preset, \w+}.unsorted.sam",
#         "runs/{sraid}/bt_{preset, \w+}.multimap.fastq",
#         "runs/{sraid}/bt_{preset, \w+}.summary.txt"
#     threads: snakemake.utils.available_cpu_count()
#     run:
#         # Parse bowtie arguments from config.yaml based on "preset"
#         if wildcards.preset in config['bowtie']:
#             bowtie_args = dict_args(config['bowtie'][wildcards.preset])
#         else:
#             bowtie_args = ""
# 
#         # Setup inputs depending on library layout
#         layout = read_value(input[2])
#         if layout == 'SINGLE':
#             rargs = '{}'.format(input[0])
#         else:
#             rargs = '-1 {} -2 {}'.format(input[0], input[1])
#         
#         # Run bowtie
#         shell(
#             "(bowtie "
#             "--threads {threads} "
#             "{bowtie_args} "
#             "--max {output[1]} "
#             "{config[btidx]} "
#             "{rargs} "
#             "{output[0]} "
#             ") 2> {output[2]}"
#             "samtools view -b > {output[0]}"
#             ") 3>&1 1>&2 2>&3 | tee {output[1]} "
#         )

''' TEtranscripts '''
