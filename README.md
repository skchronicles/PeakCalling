# ChIP-Seq Peakcalling Benchmarking
   

## 1. PeakRanger   
##### Description:  
PeakRanger is a multi-purporse software suite for analyzing next-generation sequencing (NGS) data. 
It contains the following tools:
1. `nr`: a noise ratio estimator useful for QC statistics. Estimates signal to noise ratio which is an indicator for ChIP 
enrichment.
2. `lc`: library complexity calculator useful for QC statistics. Calculates the ratio of unique reads over total reads. 
Only accepts bam files.
3. `ranger`: ChIP-Seq peak caller. Ranger servers better as a narrow-peak caller. It behaves in a conservative but 
sensitive way compared to similar algorithms. It is able to identify enriched genomic regions while at the same time 
discover summits within these regions. 
Ranger supports HTML-based annotation reports.
4. `bcp`: ChIP-Seq peak caller. Tuned for the discovery of broad peaks. BCP supports HTML-based annotation reports.  
5. `ccat`: ChIP-Seq peak caller. Tuned for the discovery of broad peaks. CCAT supports HTML-based annotation reports.
  
Peakranger is installed on [Biowulf.](https://hpc.nih.gov/apps/peakranger.html) 

##### Loading PeakRanger on Biowulf:  

    module load peakranger

##### Running NR (Noise Ratio Estimator):  

    peakranger nr \
    --format bam \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    
##### Running LC (Library Complexity Calculator):  

    peakranger lc \
    --format bam \
    {*.bam} \
    --output bcp_results  

##### Running Ranger (Narrow Peak Caller):  

    peakranger ranger \
    --format bam \
    --report \
    --plot_region 10000 \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    -t 4
 
##### Running BCP (Broad Peak Caller):  

    peakranger bcp \
    --format bam \
    --report \
    --plot_region 10000 \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    -t 4

##### Running CCAT (Broad Peak Caller):  

    peakranger ccat \
    --format bam \
    --report \
    --plot_region 10000 \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    -t 4
    

## 2. MAC2 
##### Description:
MACS empirically models the length of the sequenced ChIP fragments, which tends to be shorter than sonication or library 
construction size estimates, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a 
dynamic Poisson distribution to effectively capture local biases in the genome sequence, allowing for more sensitive and 
robust prediction. MACS compares favorably to existing ChIP-Seq peak-finding algorithms and can be used for ChIP-Seq with 
or without control samples.MAC2 is installed on [Biowulf.](https://hpc.nih.gov/apps/macs.html)  

##### Loading MACs on Biowulf:

    module load macs


##### Running MAC2 (Narrow Peak Mode):  

    module load macs/2.1.0.20150420 R
    macs2 callpeak -t {input[0]} \
    -c {input[1]} -f BAM -g {config[macs_g]} \
    --outdir peaks/mac2/narrow -n {wildcards.sample} \
    --nomodel --extsize {usePhantomPeaks.Rscript} -B -q 0.01 &> {log}
	cd peaks/mac2/narrow && Rscript {wildcards.sample}_model.r
    
##### Running MAC2 (Broad Peak Mode):  
 
    module load macs/2.1.0.20150420
        macs2 callpeak -t {input[0]} \
        -c {input[1]} -f BAM -g {config[macs_g]} \
        --broad --broad-cutoff 0.1 --nomodel --extsize {usePhantomPeaks.Rscript} \
        --outdir peaks/mac2/broad -n {wildcards.sample} -q 0.001 &> {log}

<!--- *need to use the --nomodel --extsize $extsize. get extsize by running phantompeakqualtools ...something like this > Rscript /data/CCBR_Pipeliner/3.0/Pipeliner/Results-template/Scripts/phantompeakqualtools/run_spp_nodups.R -c=$bamfile -p=${NTHREADS} -savp=${CC_PLOT} -out=${CC_SCORES}... the CC_Scores files will have the peaks, pick the first one for extsize* --->


## 3. SICER
##### Description: 
Sicer is a clustering approach for identification of enriched domains from histone modification ChIP-Seq data.

##### Loading SICER on Biowulf:

    module load sicer

##### Running SICER with controls:  

    bash {params.SICERDIR}/SICER.sh ./ {wildcards.name}.bed {params.ctrl}.bed ./ hg18 1 300 300 0.75 600 1E-2
Example:  sh DIR/SICER.sh ["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy
threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"]
["FDR"]   

Meanings of the parameters that are not self-explanatory: 
   * Species: allowed species and genome versions are listed in GenomeData.py. You can add your own species and/or genome versions and relevant data there. Redundancy Threshold: The number of copies of identical reads allowed in a
library.
   * Window size: resolution of SICER algorithm. For histone modifications, one can use
200 bp
   * Fragment size: is for determination of the amount of shift from the beginning of a
read to the center of the DNA fragment represented by the read.
FRAGMENT_SIZE=150 means the shift is 75.
   * Effective genome fraction: Effective Genome as fraction of the genome size. It
depends on read length.
   * Gap size: needs to be multiples of window size. Namely if the window size is 200,
the gap size should be 0, 200, 400, 600, ….

##### Running SICER without controls:
    bash {params.SICERDIR}/SICER-rb.sh ./ {wildcards.name}.bed ./ hg18 1 300 300 0.75 600 100

    

## 4. GEM 
##### Description: 
GEM is a high-resolution peak calling and motif discovery tool for ChIP-seq and ChIP-exo data. GEM only supports BED and SAM 
alignment file formats. 
GEM is installed on [Biowulf.](https://hpc.nih.gov/apps/gem.html)

##### Loading GEM on Biowulf:

    module load gem

##### Running GEM:  

    java -Xmx10g -jar $GEMJAR --t 24 \
    --d ./Read_Distribution_default.txt \
    --g ./mm10.chrom.sizes 
    --genome /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/ \
    --s 2000000000 
    --expt SRX000540_mES_CTCF.bed \
    --ctrl SRX000543_mES_GFP.bed \
    --f BED \
    --out mouseCTCF --k_min 6 --k_max 13

## 5. MUSIC
##### Description: 
MUSIC is a tool for identification of enriched regions at multiple scales in the read depth signals from ChIP-Seq experiments. 
MUSIC is installed on [Biowulf.](https://hpc.nih.gov/apps/music.html)

##### Loading Music on Biowulf:

    module load samtools  # needed to convert to sam format
    module load music

##### Running Music:
    mkdir chip; mkdir input
    samtools view chip.bam | MUSIC -preprocess SAM stdin chip/ 
    samtools view input.bam | MUSIC -preprocess SAM stdin input/
    samtools view /directory/to/chip.bam | MUSIC -preprocess SAM stdin chip/ 
    samtools view /directory/to/input.bam | MUSIC -preprocess SAM stdin input/
    mkdir chip/sorted;mkdir chip/dedup;mkdir input/sorted;mkdir input/dedup
    MUSIC -sort_reads chip chip/sorted 
    MUSIC -sort_reads input input/sorted 
    MUSIC -remove_duplicates chip/sorted 2 chip/dedup 
    MUSIC -remove_duplicates input/sorted 2 input/dedup
      
    MUSIC -get_multiscale_broad_ERs \
    -chip chip/dedup \
    -control input/dedup \
    -mapp Mappability_36bp \
    -l_mapp 36 \
    -begin_l 1000 \
    -end_l 16000 \
    -step 1.5
   * This code tells MUSIC to identify the enriched regions starting from 1kb smoothing window length upto 16kb with 
   multiplicative factor of 1.5 using the default parameters for the remaining parameters. The ERs for each scale are dumped.
   
* You will have to generate appropriate mappability files (https://github.com/gersteinlab/MUSIC#multi-mappability-profile-generation) Our read lengths are 100 or 125, so generate for those. Also, make sure you have bowtie2 indices for the genomes. Look at https://github.com/gersteinlab/MUSIC#running-music-with-default-parameters-and-automatic-selection-of-l_p-parameter- to automate parameter selection.*
    

## 6. PePr
##### Description:
PePr is a ChIP-Seq Peak-calling and Prioritization pipeline that uses a sliding window approach and models read counts 
across replicates and between groups with a negative binomial distribution. PePr empirically estimates the optimal 
shift/fragment size and sliding window width, and estimates dispersion from the local genomic area. Regions with less 
variability across replicates are ranked more favorably than regions with greater variability. Optional post-processing 
steps are also made available to filter out peaks not exhibiting the expected shift size and/or to narrow the width of peaks. 
PePr is installed on [Biowulf.](https://hpc.nih.gov/apps/PePr.html)

##### Loading PePr on Biowulf:  
    module load PePr

##### Running Pepr:
    PePr -c chip_rep1.bam,chip_rep2.bam \
    -i input_rep1.bam,input_rep2.bam \
    -f bam \
    -n {expname}
    
 * --shiftsize Half the fragment size.. again comes for ppqt... this should be half ext size that we used for macs
 -f needs to be bampe for PE data...not to worry about this now --> this seems to get *
 

## 6. DFilter
##### Description:
DFilter has been made to detect regulatory regions and enriched sites using tag count data. It has been made using 
a generalized approach so that data from multiple kinds of assays can be analyzed. The raw tags files can be in 6-column 
bed file, bedgraph, bam or sam format. For more information, read through
DFilter's [documentation](http://collaborations.gis.a-star.edu.sg/~cmb6/kumarv1/dfilter/tutorial.html).

##### Location of DFilter:  
    /data/CCBR_Pipeliner/db/PipeDB/bin/DFilter1.6

##### Running DFilter:
    enter command
