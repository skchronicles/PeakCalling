##########################################################################################################
# Benchmarking of ChIP-Seq peak Callers
# Author: Skyler Kuhn (NIH/NCI) [C]
# CCR Collaborative Bioinformatics Resource
# Version 1.0.0
# See README.MD for more information
# USAGE:
#   sbatch --cpus-per-task=8 --mem=16g snakemake.sh
##########################################################################################################

import sys

shell.prefix("set -eo pipefail; ")

configfile: "config.yaml"
localrules: all

controls = config["controls"]
if controls is None:
    sys.exit("Controls are needed")

samples_narrow = config["samples_narrow"]
if samples_narrow is None:
    samples_narrow = []

samples_broad = config["samples_broad"]
if samples_broad is None:
    samples_broad = []

ALL_SAMPLES = samples_narrow + samples_broad + controls

ALL_BAM = expand("bam/{sample}.sorted.Q5DD.bam", sample = ALL_SAMPLES)
ALL_BAM.extend(expand("bam/{sample}.sorted.Q5DD.bam.bai", sample = ALL_SAMPLES))
ALL_BAM.extend(expand("bam/{sample}.sorted.Q5DD.bam.flagstat", sample = ALL_SAMPLES))
ALL_BAM.extend(("bam/control.sorted.Q5DD.bam", "bam/control.sorted.Q5DD.bam.bai", "bam/control.sorted.Q5DD.bam.flagstat"))


ALL_PEAKS = expand("peaks/mac2/narrow/{sample}_peaks.narrowPeak", sample = samples_narrow) + \
expand("peaks/mac2/broad/{sample}_peaks.broadPeak", sample = samples_broad) + \
expand("peaks/ranger/{sample}_summit.bed", sample = samples_narrow) + \
expand("peaks/bcp/{sample}_region.bed", sample = samples_broad) + \
expand("peaks/ccat/{sample}_summit.bed", sample = samples_broad) + \
expand("peaks/sicer/narrow/{sample}.sorted.Q5DD-W100-normalized.wig", sample = samples_narrow) + \
expand("peaks/sicer/broad/{sample}.sorted.Q5DD-W200-normalized.wig", sample = samples_broad)


ALL_BED = expand("bed/{sample}.sorted.Q5DD.bed", sample = ALL_SAMPLES)
ALL_BED.extend(expand("bed/{sample}.sorted.bed", sample = ALL_SAMPLES))

ALL_NGSQC = expand("ngsqc/{sample}/NGSQC_report.txt", sample = ALL_SAMPLES)
ALL_NGSQC.extend(expand("ngsqc/{sample}Q5DD/NGSQC_report.txt", sample = ALL_SAMPLES))

rule all:
    input: ALL_PEAKS + ALL_BAM + ALL_BED + ALL_NGSQC

rule merge_controls:
    input:   bam = expand("bam/{sample}.sorted.Q5DD.bam", sample = controls),
             bai = expand("bam/{sample}.sorted.Q5DD.bam.bai", sample = controls),
             flagstat = expand("bam/{sample}.sorted.Q5DD.bam.flagstat", sample = controls) 
    output:  bam = "bam/control.sorted.Q5DD.bam",
             bai = "bam/control.sorted.Q5DD.bam.bai",
             flagstat = "bam/control.sorted.Q5DD.bam.flagstat"
    log:     "log/merge_controls"
    threads: 2
    shell:
        '''
        inbam=( {input.bam} )
        if [[ ${{#inbam[@]}} -eq 1 ]]; then
            ln -s $(cd $(dirname {input.bam}) && pwd)/$(basename {input.bam}) {output.bam}
            ln -s $(cd $(dirname {input.bai}) && pwd)/$(basename {input.bai}) {output.bai}
            ln -s $(cd $(dirname {input.flagstat}) && pwd)/$(basename {input.flagstat}) {output.flagstat}
        else
            module load samtools/1.2
            samtools merge -r -@{threads} {output.bam} {input.bam}
        fi
        '''


rule MAC2_narrowPeaks:
    input:  "bam/{sample}.sorted.Q5DD.bam", "bam/control.sorted.Q5DD.bam", "bam/{sample}.sorted.Q5DD.ppqt"
    output: "peaks/mac2/narrow/{sample}_model.r", "peaks/mac2/narrow/{sample}_peaks.narrowPeak",
            "peaks/mac2/narrow/{sample}_peaks.xls", "peaks/mac2/narrow/{sample}_summits.bed",
            "peaks/mac2/narrow/{sample}_model.pdf"
    log:    "log/mac2/{sample}.find_narrow_peaks"
    threads: 2
    run:
        fh = open(input[2])
        fragsize = 150
        for line in fh:
            try:
                linelist = line.strip().split("\t")
                fragsize = int(linelist[2].split(",")[0])
            except IndexError:
                pass
        fh.close()
        # took out -B
        shell('''
        module load macs/2.1.0.20150420 R
        macs2 callpeak -t {input[0]} \
        -c {input[1]} -f BAM -g {config[macs_g]} \
        --outdir peaks/mac2/narrow -n {wildcards.sample} \
        --nomodel --extsize {fragsize} -q 0.01 &> {log}
        cd peaks/mac2/narrow && Rscript {wildcards.sample}_model.r
        ''')

rule MAC2_broadPeaks:
    input:  "bam/{sample}.sorted.Q5DD.bam", "bam/control.sorted.Q5DD.bam", "bam/{sample}.sorted.Q5DD.ppqt"
    output: "peaks/mac2/broad/{sample}_peaks.xls", "peaks/mac2/broad/{sample}_peaks.broadPeak"
    log:    "log/mac2/{sample}.find_broad_peaks"
    threads: 2
    run:
        fh = open(input[2])
        fragsize = 150
        for line in fh:
            try:
                linelist = line.strip().split("\t")
                fragsize = int(linelist[2].split(",")[0])
            except IndexError:
                pass
        fh.close()

        shell('''
        module load macs/2.1.0.20150420
        macs2 callpeak -t {input[0]} \
        -c {input[1]} -f BAM -g {config[macs_g]} \
        --broad --broad-cutoff 0.1 --nomodel --extsize {fragsize} \
        --outdir peaks/mac2/broad -n {wildcards.sample} -q 0.001 &> {log}
        ''')
    

rule bams2beds:
    input:  "bam/{sample}.sorted.bam", "bam/{sample}.sorted.Q5DD.bam"
    output: "bed/{sample}.sorted.bed", "bed/{sample}.sorted.Q5DD.bed"
    log:    "log/bed2bam/"
    threads: 2
    shell:
        '''
        module load bedtools
        bedtools bamtobed -i {input[0]} > {output[0]} 2> {log}{wildcards.sample}.bam2bed
        bedtools bamtobed -i {input[1]} > {output[1]} 2> {log}{wildcards.sample}.Q5DD.bam2bed
        '''

rule ngsqc:
    input:  "bed/{sample}.sorted.bed", "bed/{sample}.sorted.Q5DD.bed"
    output: "ngsqc/{sample}/NGSQC_report.txt", "ngsqc/{sample}Q5DD/NGSQC_report.txt"
    log:    "log/ngsqc/"
    threads: 2
    shell:
        '''
        module load bedtools
        /scratch/ChIPSeqBenchmarking/NGSQC_linux_x86_64 -v -o ngsqc/{wildcards.sample} {input[0]} /scratch/ChIPSeqBenchmarking/genomes/hg19.genome 2> {log}{wildcards.sample}.ngsqc
        /scratch/ChIPSeqBenchmarking/NGSQC_linux_x86_64 -v -o ngsqc/{wildcards.sample}Q5DD {input[1]} /scratch/ChIPSeqBenchmarking/genomes/hg19.genome 2> {log}{wildcards.sample}Q5DD.ngsqc
        '''

rule ranger:
    input:  "bam/{sample}.sorted.Q5DD.bam","bam/control.sorted.Q5DD.bam"
    output: "peaks/ranger/{sample}_summit.bed", "peaks/ranger/{sample}_region.bed"
    log:    "log/ranger/"
    threads: 4
    shell:
        '''
        module load peakranger
        mkdir --p peaks/ranger
        peakranger ranger --format bam --data {input[0]} --control {input[1]} --output peaks/ranger/{wildcards.sample}
        '''

rule bcp:
    input:  "bam/{sample}.sorted.Q5DD.bam","bam/control.sorted.Q5DD.bam"
    output: "peaks/bcp/{sample}_region.bed"
    log:    "log/bcp/"
    threads: 4
    shell:
        '''
        module load peakranger
        mkdir --p peaks/bcp
        peakranger bcp --format bam --data {input[0]} --control {input[1]} --output peaks/bcp/{wildcards.sample}
        '''

rule ccat:
    input:  "bam/{sample}.sorted.Q5DD.bam","bam/control.sorted.Q5DD.bam"
    output: "peaks/ccat/{sample}_summit.bed", "peaks/ccat/{sample}_region.bed"
    log:    "log/ccat/"
    threads: 4
    shell:
        '''
        module load peakranger
        mkdir --p peaks/ccat
        peakranger ccat --format bam --data {input[0]} --control {input[1]} --output peaks/ccat/{wildcards.sample}
        '''

rule sicer_narrow:
    input:  "bed/{sample}.sorted.Q5DD.bed", "bam/{sample}.sorted.Q5DD.ppqt"
    output: "peaks/sicer/narrow/{sample}.sorted.Q5DD-W100-normalized.wig"
    log:    "log/sicer/{sample}.find_narrow_peaks"
    threads: 2
    run:
        fh = open(input[1])
        fragsize = 150
        for line in fh:
            try:
                linelist = line.strip().split("\t")
                fragsize = int(linelist[2].split(",")[0])
            except IndexError:
                pass
        fh.close()
        
        shell('''
        module load sicer bedtools
        if [ ! -f bed/control.sorted.Q5DD.bed ]; then
            bedtools bamtobed -i bam/control.sorted.Q5DD.bam > bed/control.sorted.Q5DD.bed
        fi
        bash SICER.sh ./bed {wildcards.sample}.sorted.Q5DD.bed control.sorted.Q5DD.bed ./peaks/sicer/narrow hg19 1 100 {fragsize} 0.79 200 0.01
        ''')

rule sicer_broad:
    input:  "bed/{sample}.sorted.Q5DD.bed", "bam/{sample}.sorted.Q5DD.ppqt"
    output: "peaks/sicer/broad/{sample}.sorted.Q5DD-W200-normalized.wig"
    log:    "log/sicer/{sample}.find_broad_peaks"
    threads: 2
    run:
        fh = open(input[1])
        fragsize = 150
        for line in fh:
            try:
                linelist = line.strip().split("\t")
                fragsize = int(linelist[2].split(",")[0])
            except IndexError:
                pass
        fh.close()
        
        shell('''
        module load sicer bedtools
        if [ ! -f bed/control.sorted.Q5DD.bed ]; then
            bedtools bamtobed -i bam/control.sorted.Q5DD.bam > bed/control.sorted.Q5DD.bed
        fi
        bash SICER.sh ./bed {wildcards.sample}.sorted.Q5DD.bed control.sorted.Q5DD.bed ./peaks/sicer/broad hg19 1 200 {fragsize} 0.79 400 0.01
        ''')
