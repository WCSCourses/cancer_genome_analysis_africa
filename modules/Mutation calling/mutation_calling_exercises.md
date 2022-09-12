Somatic Mutation Calling
------------------------

# Overview
In somatic mutation calling, we aim to locate
and genotype mutations that have occurred in a collection of
somatic cells (such as a tumor). This involves calling the mutations
present in those cells and then removing the _germline background_ that
is present in all cells within the organism.


In this practical, we will:

1. See how to generate a BAM file for a tumor and normal sample
2. Prepare this file for variant calling
3. Call variants with GATK HaplotypeCaller
4. Look at our variants in IGV to assess their quality.
5. Build an intuition for variant quality control.

# Input data

In the VM, we have a directory titled `TODO REPLACE FOLDER PATH`.
This directory contains data from the [Texas Cancer Research Biobank Open Access project.](http://stegg.hgsc.bcm.edu/open.html).
We're using whole-genome sequenced data from Case 006, a woman in her 60s who presented with neuroendocrine carcinoma
of the pancreas and received no prior treatment.

```bash

```
Files present in the `TODO REPLACE FOLDER PATH` directory.


The `TODO REPLACE FOLDER PATH` directory contains the FASTQs, BAMs, and example VCFs for this patient's tumor and normal sample.
As we go through this tutorial, you'll run the commands to generate each of these files, but you can also skip over long-running
computational tasks since the data needed is already present. As this is standard whole-genome data, you can also use this for 
testing and learning other software so long as you follow the data access agreement rules.

# Practical
Below, we ask some concept and knowledge questions before moving on to analyzing a match tumor-normal pair.

## Concepts and knowledge questions

1. What is variant calling?
```



```

2. Approximately how many variants do we expect for a given normal sample if we know:  
- the average genome differs from the reference approximately every 1000bp
- the human genome is roughly 3.2 Billion basepairs in length?
```


```

3. Let's say our tumor has a mutation burden of 1.5 mutations per _megabase_. Approximately how many mutations does
this tumor have in total?

```



```

4. To call somatic variants, we do the following:
- Call all mutations in the tumor
- Call all mutations in the normal
- Subtract out the germline background to generate somatic alls.
Knowing this, which of the following Venn Diagrams best represents our
data and the expected number of variants in the germline and somatic VCFs?
![](images/somatic_venns.png)

## Preprocessing (read alignment, duplicate marking, sorting, indexing)
Our raw sequencing reads are random pieces of DNA derived from our tumor and normal tissues.
To make use of them, we need to align our sequences to a reference genome to produce _alignments_.
This will help us make sense of _where_ each read comes from in the reference and what sites in our
sample are different, or _variant_, from the reference genome.

We will also perform duplicate marking, base quality score recalibration, sorting, and indexing after alignment.
These processes will do the following:
1. Remove reads that are optical duplicates (i.e., technical noise) that could bias our variant calls.
2. Normalize our base quality scores, which helps improve the accuracy of variant calls.
3. Sort and index our file so we a) save disk space storing it and 2) programs can make use of efficient random access into it.

### Staying organized
We'll run our tutorial in a new directory so we can keep our work organized. We'll use some basic linux commands to set up a directory.

```bash
cd ~

mkdir mutation-calling-exercise
cd mutation-calling-exercise

```

### Alignment
We use a program called [Burrows-Wheeler Alginer (BWA)](https://github.com/lh3/bwa) for aligning reads. While there are many, many programs for sequence alignment,
BWA is generally considered to be both fast and accurate and is widely used within the community. It is also free and open-source.


BWA contains several algorithms for read alignment. In this tutorial, we'll use the `mem` algorithm, which leverages Maximal Exact Matches
and clever heuristics to be both faster and more accurate than the `aln/samse/sampe` and `bwasw` algorithms.


BWA requires an input reference genome as well as several indexes it uses to efficiently align reads. Today, we'll be using the `Homo_sapiens_assembly38.fasta` reference from the Broad Institute. These indexes are provided for you
in the `~/references/` directory and are also downloadable from the [Broad Resource Bundle site](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false). If you needed to generate them for a new reference genome, you could use the `bwa index` command like so:

```bash
bwa index ~/references/Homo_sapiens_assembly38.fasta
```


Once we have our indices, we are ready to align reads. Remember, we'll use the mem algorithm. The basic input form of BWA is like so:

```bash
bwa mem <reference> <inputFASTQ_1> <inputFASTQ_2>
```

to see a full list of options:
```bash
bwa mem
```

We'll also tune the performance of BWA by using the following parameters, such as the number of processing threads and the 
number of reads we keep in memory at any given time. We'll also pass the `-Y ` flag to enable softclipping on supplementary
alignments.

We will also add a Read Group to our data. This step is essential for making our BAM useful in downstream calling - it annotates
each read with the Read Group so that callers know how to use reads in their statistical models. The Read Group takes the form
of a string so we must enclose in single quotes.

To align our reads, we run the following command:
```bash
bwa mem \
    -Y \
    -t 2 \
    -K 100000 \
    -R '@RG\tID:TCRBOA6-Normal-RG1\tLB:lib1\tPL:Illumina\tSM:TCRBOA6-Normal\tPU:TCRBOA6-Normal-RG1' \
    ~/references/Homo_sapiens_assembly38.fasta \
    TODO REPLACE FOLDER PATH /chr22.TCRBOA6-Normal_1.fastq.gz TODO REPLACE FOLDER PATH /chr22.TCRBOA6-Normal_2.fastq.gz \
    | samtools sort \
    -@ 2 \
    -o chr22.TCRBOA6-Normal.bam -
```
Expected runtime: 30-50 minutes.

The `samtools sort` command will write a sorted BAM file (since the default output of BWA is unsorted SAM).
Sorting means that alignments appear in the file in their order along the genome (i.e., the 1st position of chromosome 1 is
at the start of the file and the last alignment in the file will be on the last chromosome). This makes it possible to 
index our BAM later.

Let's check if our BAM is valid - we can do so with `samtools quickcheck`:

```bash
samtools quickcheck chr22.TCRBOA6-Normal.bam
```

Expected runtime: 1 second.

This command should return nothing and report nothing if everything is in order.

Are there other commands we might use for BAM quality control?

```


```


Lastly, we need to do the same process for our tumor. This will take longer to align - expect roughly an hour (but
remember, these files have already been provided for you).

```bash
bwa mem \
    -Y \
    -t 2 \
    -K 100000 \
    -R '@RG\tID:TCRBOA6-Tumor-RG1\tLB:lib1\tPL:Illumina\tSM:TCRBOA6-Tumor\tPU:TCRBOA6-Tumor-RG1' \
    ~/references/Homo_sapiens_assembly38.fasta \
    TODO REPLACE FOLDER PATH /chr22.TCRBOA6-Tumor_1.fastq.gz TODO REPLACE FOLDER PATH /chr22.TCRBOA6-Tumor_2.fastq.gz \
    | samtools sort \
    -@ 2 \
    -o chr22.TCRBOA6-Tumor.bam -
```

#### Duplicate Marking
Once we've aligned our reads we need to perform some additional steps to normalize our data. The first of these is duplicate
marking, where we remove reads that are likely optical duplicates generated during the sequencing process. These duplicates
can distort our variant calls downstream if not removed.

We'll use the Picard MarkDuplicates tool to mark duplicates. Conveniently, this is now integrated into the Genome Analysis Toolkit.

```bash
gatk MarkDuplicates \
    --java-options -Xmx4g \
    -I chr22.TCRBOA6-Normal.bam\
    -O chr22.TCRBOA6-Normal.markdups.bam \
    -M chr22.TCRBOA6-Normal.metrics.txt
```

We'll need to do the same for our tumor sample:

```bash
gatk MarkDuplicates \
    --java-options -Xmx4g \
    -I chr22.TCRBOA6-Tumor.bam\
    -O chr22.TCRBOA6-Tumor.markdups.bam \
    -M chr22.TCRBOA6-Tumor.metrics.txt
```

Expected runtime: 1 minute per sample

The MarkDuplicates tool prints some statistics to stdout / stderr when it's finished. How long did the run take?

```


```

How many unmapped reads were present? (hint: check the metrics.txt file)

```


```

What percent of reads were optical duplicates?

```


```

#### Base Quality Score Recalibration

The next step in the process is to generate a Base Quality Score Recalibration (BQSR) report.
The BQSR process will normalize the quality scores within a BAM file based on a set of known variants
(especially indels), resulting in more accurate variant calls downstream.

```bash
gatk BaseRecalibrator \
    --java-options -Xmx4g \
    --input chr22.TCRBOA6-Normal.markdups.bam \
    --output chr22.TCRBOA6-Normal.markdups.BQSR-REPORT.txt \
    --known-sites ~/references/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --reference ~/references/Homo_sapiens_assembly38.fasta
```


```bash
gatk BaseRecalibrator \
    --java-options -Xmx4g \
    --input chr22.TCRBOA6-Tumor.markdups.bam \
    --output chr22.TCRBOA6-Tumor.markdups.BQSR-REPORT.txt \
    --known-sites ~/references/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --reference ~/references/Homo_sapiens_assembly38.fasta
```



Expected runtime: 5 minutes for normal, 10 minutes for tumor.



What is the name of the read group in the BQSR Report for the normal sample? Note that
it should be the same as we generated for alignment.

```



```

## Indexing BAM files

To call variants, we'll need to index our BAM files. We can use the `samtools index`
command to do so:

```bash
samtools index chr22.TCRBOA6-Tumor.markdups.bam

samtools index chr22.TCRBOA6-Normal.markdups.bam
```

Our BAM index is kind of like the index of a book. What do you think its coordinate system is?

```



```


## Variant calling

We can now run MuTect2 to call somatic variants in our matched tumor and normal samples.
Since we're only interested in chromosome 22 (at least for this analysis), we can 
tell MuTect2 to only call that interval using `-L chr22`. 

```bash
gatk Mutect2 \
    -R ~/references/Homo_sapiens_assembly38.fasta \
    --input chr22.TCRBOA6-Tumor.markdups.bam \
    --tumor-sample TCRBOA6-Tumor \
    --input chr22.TCRBOA6-Normal.markdups.bam \
    --normal-sample chr22.TCRBOA6-Normal \
    -L chr22 \
    --output chr22.TCRBOA6-Tumor.TCRBOA6-Normal.vcf
```

Expected runtime: 30-60 minutes


We now have a somatic VCF file from MuTect2. But before we start looking for driver mutations or 
mutational signatures, we should run some basic quality control and assess some of our variants in IGV.
The following section will give a brief overview of quality control and assessment techniques.



## Variant assessment and quality control

There are many ways to do variant quality control. Below, we demonstrate how to use Integrated Genomics Viewer (IGV),
which is perhaps the most commonly-used program for analyzing variants in BAMs/VCFs.

We also introduce some packages for basic quality control, but we enourage you to focus on _what's_ being plotted
or analyzed rather than how to do it. Programs change throughout time; focusing on the basic parameters to think about when
checking your data - depth, quality, filter fields, etc. - will help you build an intuition for how to do quality control regardless
of what pipeline or software your use.

How many variants are in our VCF file?  
(Hint: we can use `grep -c "<pattern>" <file>` to count the number of lines that match a pattern in a file.
If we want the number of lines that _don't_ match a pattern, we can use `grep -c -v "<pattern>" <file>`).

```


```


## Variant annotation with Funcotator

To make use of our variants, we'll want to annotate them. This process uses databases 

## Manual Review

Manual review of variants is an essential step of the variant calling process. While variant callers use complex statistical
models and clever heuristics, they still very often make erroneous calls (poor specificity) or miss true variants (poor sensitivity).
We must filter our variants to generate the most specific, sensitive set of variants we can.


We can view our variants using the Unix program `less`. 

```bash
less -S chr22.TCRBOA6-Tumor.TCRBOA6-Normal.vcf
```

The VCF header (lines beginning with `#`) has information about the fields within the file.
Each variant caller will use its own fields, though there are some common ones across callers
that are useful for filtering:

- `GT`: Genotype. In tumors, we generally expect true somatic variants to have a `0/1` genotype and a `0/0` genotype in the normal sample.
- `AD`: Allele Depth. This field provides one value per allele and contains the number of reads supporting that allele. The two AD values will generally add up to the DP value at a given position. AD is a very useful field - we can use it to determine our CCF, 
tumor purity, and true positive variants.
- `DP`: Depth. The number of reads covering a given position. We should expect this value to be close to the estimated depth of our sequencing.
- `QUAL`: the VCF quality field is filled by every caller, though it doesn't have a strict definition. It can't be compared across variant callers and a filter that works for one callset may not work well for another.
- For mutect2, `F1R2` and `F2R1`: these describe the number of reads in proper versus improper orientation that supports the variant. Too many reads in the improper orientation can indicate that the variant we're examining is unlikely to be real.


Let's look at an example line from our VCF (and the last line in our VCF header):
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCRBOA6-Normal  TCRBOA6-Tumor
chr22   10573224        .       C       CT      .       .       AS_SB_TABLE=6,19|3,2;DP=30;ECNT=1;MBQ=35,35;MFRL=546,354;MMQ=40,40;MPOS=38;NALOD=0.997;NLOD=2.66;POPAF=6.00;RPA=2,3;RU=T;STR;TLOD=10.46 GT:AD:AF:DP:F1R2:F2R1:FAD:SB     0/0:9,0:0.092:9:3,0:5,0:9,0:2,7,0,0     0/1:16,5:0.263:21:6,2:10,3:16,5:4,12,3,2
```

This variant is an indel - specificially, an insertion of a single T after a C at position 10573224 on chromosome 22. In our normal sample, mutect2 called this variant homozygous reference, with nine reads supporting the reference allele. In the tumor,
mutect2 called the variant heterozygous, with 16 reference-supporting reads and 5 alternate-supporting reads. The Allele Fraction (AF) field in the tumor is 0.263 (5 / (5 + 16)); this is a pretty low allele fraction for a real variant. The variant has a depth of 9 in
the normal, which is relatively low for this high-depth sequencing. While this variant seems slightly suspicious, there's nothing
here immediately flags it for removal outside of the low allele fraction.


To better assay the variant, let's look at it in IGV. IGV review is performed in nearly every study. There's a lot of nuance
and intuition used in determining whether a variant is real or not in IGV. As you view many true variants, take note of what
a real variant looks like and what characteristics are shared among false variants.


We can open up IGV, change the genome to human HG38, and load our files for the tumor and normal. Next, we'll use the location
bar to navigate to the variant we reviewed above. We can then click the coverage bars at the location to open up the detailed
allele fraction view. The resulting IGV screen will look roughly like the screenshot below:

![](images/igv_3.png)