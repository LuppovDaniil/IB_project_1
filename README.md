# Lab Notebook

Here you can find the detailed report including raw data, commands, intermediate and final results. Moreover, you can use this report as a manual and repeat steps written below to gain your own results and compare them to ours.


*10/19/2022*

Firstly, we studied our model organism *E.Coli* making notes in our Lab report. The thing of our interest is the certain k-12 strain, which is non-resistant to antibiotics.

The raw data is stored in *rawdata/*. It includes the sequence of the strain (GCF_000005845.2_ASM584v2_genomic.fna.gz) and the annotation (GCF_000005845.2_ASM584v2_genomic.gff.gz), both were taken from [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/). Also you can find the folder *rawdata/reads/* with paired raw Illumina reads of an resistant to ampicillin *E. coli* strain, whose were taken from [here](https://doi.org/10.6084/m9.figshare.10006541.v3).

We downloaded required data and extracted all files using tar:

```bash
tar -xf archive.tar.gz
```

Then we looked at the first 20 lines in each file to make sure that we received correct types of data.

```bash
head -20 filename.format
```

Then we checked the number of reads in forward and reverse files.

```bash
wc -l filename.fastq
```

We got 1823504 reads for the forward file and 1823504 for the reverse file. Thus, we were sure that numbers of reads match each other.

Then we installed [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). This program calculates the quality of sequencing, so we wanted to make sure that our data is qualitative enough to work with.

```bash
apt-get install fastqc
fastqc -h
```

The second command was run to make sure that the program was installed successfully and look at the short manual.

Then we ran the program with our forward and reverse files. The flag `-o` was used just for outputting results to the same directory.

```bash
fastqc -o . /pathtofile1/file1.fastq /pathtofile2/file2.fastq 
```

Then we opened reports made with fastqc. You can find them in the folder *results/*. Some interesting artifacts were found. First, the quality of reads which are placed at the end of each file was not as fine as at the start. But it is the usual thing and these "obstacles" are not significant. Then we discovered another "problem" at the page "Per tile sequence quality". We consider it was caused by bubbles in a flow cell. Also these artifacts are not critical. All pictures are presented in our lab report.

Then we decided to use [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) program to improve the quality. We set the following parameters: cutted bases off the start of a read if their quality was below 20 (LEADING:20), cutted bases off the end of a read if their quality was below 20 (TRAILING:20), trimmed reads with the window size 10 and the average quality within the window 20 (SLIDINGWINDOW:10:20), dropped reads if they were below length 20 (MINLEN:20).

```bash
conda install trimmomatic
trimmomatic PE -phred33 amp_res_1.fastq amp_res_2.fastq _1P.fq _1U.fq _2P.fq _2U.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20
```
Trimmomatic report: Input Read Pairs: 455876 Both Surviving: 446259 (97,89%) Forward Only Surviving: 9216 (2,02%) Reverse Only Surviving: 273 (0,06%) Dropped: 128 (0,03%)

Then we checked the quality with fastqc again, results were better than previous ones. Pictures are shown in our lab report too.

Then we wanted to align our reads using [BWA-MEM](https://github.com/lh3/bwa) program, which is optimized for long NGS reads of 100bp or more. It looks for the maximum exact matches of seeds, and then it extends the seed using fitting or local alignment to map the entire read. First, we needed to index our WGS.

```bash
apt-get install bwa
bwa index GCF_000005845.2_ASM584v2_genomic.fna
```

The next command was run to make alighnment.

```bash
bwa mem reference_file forward_reads reverse_reads > alignment.sam
```

BWA-MEM created the sam file, which we wanted to compress. Then we used [samtools](http://www.htslib.org/) program. 

```bash
apt-get install samtools
samtools view -S -b alignment.sam > alignment.bam
samtools flagstat alignment.bam
```
The last command gave us statistics. 891649 reads were mapped.


*10/23/2022*

*10/26/2022*
