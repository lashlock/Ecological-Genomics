

# 2017 Ecological Genomics Course

### Author: Lauren Ashlock     


## Overall Description of notebook      
Notes from class material,and class project will populate this notebook. 

## Date started: (2017-02-01)   
## Date end:   ongoing    


### Table of contents for 60 entries (Format is *Page: Date(with year-month-day). Title*)        
* [Page 1: 2017-02-01](#id-section1). Sequencing strategies applied to biological questions
* [Page 2: 2017-02-06](#id-section2).  RNAseq   
* [Page 3: 2017-02-08](#id-section3). Transcriptomics
* [Page 4: 2017-02-13](#id-section4) . RNAseq mapping
* [Page 5: 2017-02-15](#id-section5). SNPs and Population genetics
* [Page 6: 2017-02-22](#id-section6). DESEQ 2 Tutorial
* [Page 7: 2017-02-27](#id-section7). Scott Edwards and Differential Expression Analysis
* [Page 8: 2017-03-01](#id-section8). Differential expression - Catch up Day
* [Page 9:2017-03-06](#id-section9). Population Genomics
* [Page 10:2017-03-08](#id-section10).Effective population size
* [Page 11:2017-03-08](#id-section11). R script for homework 2
* [Page 12:2017-03-20](#id-section12). Population genetic structure
* [Page 13:2017-22-03](#id-section13). Population Genetic Structure
* [Page 14 :2017-27-03](#id-section14). Species Divergence
* [Page 15:2017-29-03](#id-section15).Identifying local adaptation
* [Page 16:2017-03-04](#id-section16).Fst
* [Page 17:2017-04-04](#id-section17).HW #3 notes
* [Page 18:2017-05-04](#id-section18). Assignment 3 R code
* [Page 19:2017-10-04](#id-section19).Metagenomics
* [Page 20:2017-12-04](#id-section20).Microbiome
* [Page 21:2017-12-04](#id-section21).Review of Pipeline
* [Page 22:2017-17-04](#id-section22).Microbiome
* [Page 23:2017-19-04](#id-section23).Microbiome
* [Page 24:2017-02-05](#id-section24).Final Project Script
* [Page 25:](#id-section25).
* [Page 26:](#id-section26).
* [Page 27:](#id-section27).
* [Page 28:](#id-section28).
* [Page 29:](#id-section29).
* [Page 30:](#id-section30).
* [Page 31:](#id-section31).
* [Page 32:](#id-section32).
* [Page 33:](#id-section33).
* [Page 34:](#id-section34).
* [Page 35:](#id-section35).
* [Page 36:](#id-section36).
* [Page 37:](#id-section37).
* [Page 38:](#id-section38).
* [Page 39:](#id-section39).
* [Page 40:](#id-section40).
* [Page 41:](#id-section41).
* [Page 42:](#id-section42).
* [Page 43:](#id-section43).
* [Page 44:](#id-section44).
* [Page 45:](#id-section45).
* [Page 46:](#id-section46).
* [Page 47:](#id-section47).
* [Page 48:](#id-section48).
* [Page 49:](#id-section49).
* [Page 50:](#id-section50).
* [Page 51:](#id-section51).
* [Page 52:](#id-section52).
* [Page 53:](#id-section53).
* [Page 54:](#id-section54).
* [Page 55:](#id-section55).
* [Page 56:](#id-section56).
* [Page 57:](#id-section57).
* [Page 58:](#id-section58).
* [Page 59:](#id-section59).
* [Page 60:](#id-section60).

------
<div id='id-section1'/>
### Page 1: 

**Info update: Whole Genome Sequencing**

> Applications   
>
> Prior considerations
>
> Methods

- High power and resolution 

Reference Genome necessary?

- No? -> De Novo assembly
  - Adaptations
  - Genome expression
- Yes? -> Looking for important variation
  - Epigenetic modifications
  - Analyze DNA --> protein process
  - Gene expression studies

Do you have the necessary resources?

- Money
- Computational resources
  - Core facility
  - Command line
  - Programming language (Python/Perl)

Limitations of your data:

- Polymorphic genes or Paralogs (core, highly conserved genes)
- Rapidly evolving genes  and Large gene families (poor representation)
- Sequencing from one individual  or pooled samples
- Impossible to sequence all aspects of the genome
  - heterochromatic regions  (DNA too tied up in chromatin to be sequenced)
  - highly repetitive regions (hard to tell where it should go in your assembly)

Prior Considerations:

- Sequencing platform
  - Short reads (Illumina Seq 150bp, or SOLID)
  - Long reads (PacBio 5kb, Iion Torrent 500bp, Illumina Moleculo up to 10 kb)
- Knowledge of organism
  - Genome size (K-mer approach)
  - Repeat content
  - Error rate of sequencing 
  - Degree of genome duplications (polyploidy)
- Methods
  - Wet lab procedures
    - Sample (high quality DNA), avoid energetically active tissue (lots of mitochondrial DNA), avoid gut and skin (DNA from other orgs)
    - Quantity: 1mg -> 6ug
  - Library Prep
    - single or paired end sequences
    - Mate pair (ends far apart you put together)
    - Assemble paired reads into a contig
    - Assemble scaffold
    - N approximation of gaps using long reads
  - Before Library Prep
    - GC content
    - Repeat abundane
    - Duplicate reads
  - Quality control
    - Reaper
    - N50 contig in the middle (50% longer 50% shorter) - the number of bp in that contig)
  - Annotation
    - Automated or Manual

**Info Update: RNA seq**

> Advantage   
>
> Limitation
>
> Workflow
>
> > set-up
> >
> > wet lab
> >
> > seq strategy
> >
> > bio info
> >
> > statistical measures

Advantages:

- Differential gene expression
- Allele specific expression (environmental response/ adaptation)
- Functionally relevant subset of the genome 
- Wide range of expression value
- Info re splicing events

Limitations: 

- Does not necessarily correlate with protein abundance
- Tissue specific

Work Flow:

- Purpose:
  - coding or regulatory non coding? 
  - Reference genome?
  - Alternative splicing?
  - Technology?
  - Population level or treatment level comparison
- Stats:
  - Biological replication
  - Tissue choice
  - Coordinated time of sampling
  - Small organism (pooled samples... )
- Wet Lab:
  - RNAse free environment
  - Treat with DNAse
  - eliminate rRNA 
  - Enrich RNA with poly-A tail
  - cDNA synthesis: oligodT primer (complementary to poly-A tail). Also add reverse transcriptase, to make DNA from the mRNA
  - Library Prep: single vs paired end. 
- Sequencing platform
  - pyrosequencing (Roche), Ion Torrent, Illumina HiSeq
    - Roche - incorrect homopolymer
    - GC content template 
- Sequencing coverage
  - @ least 100 million bp (100bp)
  - 10 million
- Programs: 
  - unix
  - Python
  - R

**Info Update: Amplicon Seq**

Library Prep

- Extract genomic DNA
- 1st PCR
  - amplify targeted gene using specified primers
    - (using 16s as an example)
    - variable and conserved regions (200-600bp)
- Clean up
  - using gel extraction or column
- PCR
  - adding bar codes and adapters for sequencing
- Clean up
- Pool samples
- Sequencing

Sequencing: 

- 454 (outdated)
- miseq: paired end sequencing, 300bp

Data Analysis:

- Trim adapters
- Align reads using conserved sequences

**Info Update: GBS/RAD Seq**



------
<div id='id-section2'/>
### Page 2:   

**RNAseq workflow**

1. Clean reads *fastq*

   1. Adapters
   2. Nucleotide quality
   3. Length

2. Evaluate quality

3. Assemble De Novo transcriptome *fasta*

   1. Evaluate assembly
      1. compare to closely related species or core set of genes

4. Annotate reference

   1. nr (gene annotation)
   2. uniprot (protein) database
   3. Gene Ontology (GO) 
   4. ​

5. Map reads to de novo transcriptome

   1. Generates a lot of alignment files *sam*

6. Extract read count info (# of reads that map to each contig or each sample) and Identify SNPs 

7. Differential gene expression analyses 

   1. co expression network analyses

8. Population genomics 

   1. genetic differentiation
   2. population structure
   3. demographic history
   4. signatures of selection

   -----------------------------------

   Logging onto the server: Host name: pbio381.uvm.edu

   UN and PW same as netID

   - Each file has the individual number, disease level (0-5), sample date, paired files (left and right paired ends)
   - my files to clean 08__5-14__S_1_R1 and R2
   - seeing the top of a zipped file 
     - zcat filename | head
   - First line is unique identifier (sequencer, adaptors)
   - Starting at second line is our data
   - Line after plus sign indicates quality of that nucleotide (in one character per quality score so it alligns perfectly with your sequence)
     - phred score
   - Input Read Pairs: 32170422 Both Surviving: 27242302 (84.68%) Forward Only Surviving: 3624951 (11.27%) Reverse Only Surviving: 601581 (1.87%) Dropped: 701588 (2.18%)
   - Put clean paired fastqc files in fastqc_out directory
   - also put them on my zoo.uvm.edu~/lashlock

   ​

   ​

------
<div id='id-section3'/>
### Page 3:   

**Info update: transcriptomics**

- Useful when studying non-model organisms  or non traditional model organisms
- Helpful for finding genes responding to multiple stimuli
- Novel transcripts w/o homologs in closely related model organisms

Methods

- Microarrays
  - Older technology
  - Easy to use, and simple analyses
- RNAseq
  - Newer technology
  - genome wide ecological transcripts
  - Data processing is more intensive and time consuming (because more data?)

Questions:

1) How much variation is there in gene expression and how is it structured?

- Evolutionary Processes - Gene expression is heritable and natural selection can act on it
    - can also be facilitated by epigenetics
    - Qst and Fst comparisons
    - eQTL
- Macroevolution - 
    - Drift, selection, bottlenecks

2) How do environmental stimuli affect gene expression?

- Abiotic stress
- Environmental heterogeneity in time or space
- host-parasite interactions 
- Selective biotic and abiotic interactions
  - Digging deeper to identify the molecular basis of response
  - how the environment affects the genome
  - how the environment affects phenotypic plasticity
- Considerations: need to be flash frozen... and must consider that it is only a snapshot of expression
  - can use time course analyses to get around this

3) How does gene expression affect phenotype?

- Alternate phenotypes
- Moving from correlation to causation
  - transgenics, RNAi, Crisper/CAS9

Future directions

- RNAseq
- Database for proposed ecological annotations

Primary problems

- Bias in signal
- Heterologous arrays
- Polyploidy
- RNA pooling
- Statistical analyses
- Unannotated genes


------
<div id='id-section4'/>
### Page 4:    

Missed class this day

Tutorials and scripts can be found [here](https://adnguyen.github.io/2017_Ecological_Genomics/Tutorial/2017-02-13_RNAseq_Mapping.html).



------
<div id='id-section5'/>
### Page 5:   

**Info Update: SNPs and Population Genomics**

- SNP data - expressed sequences

1.    Tissue
      1. breadth of tissue from different developmental stages
         1. This controls for exon skipping
2. Pool samples and create your sequence libraries (30-100Million p.e. long reads)
3. Process raw sequence data
      1. Important for SNP detection
4. Digital normalization
      1. Remove high coverage reads and associated errors
         1. Loss of quantitative info
5. Assemble clean p.e. reads
6. Prune assembled transcripts
      1. Reduce DNA contaminations, noncoding DNA, and gene fragments
7. Assembly evaluation using either a reference genome or conserved genes in other eukaryotic organisms
8. SNP detection
      1. Software: constant patterns of sequence variation
         1. sequence errors - hopefully your software will eliminate reads of low frequency
         2. Errors can also be found in homogeneous regions... amplified by PCR... can filter for these
         3. Artifacts caused by InDels - filter SNP clusters near indels, quality scores
9. SNP validation - primers
      1. Use Sanger sequencing or mass spec to quality control a portion of your sequence data
10. Applications
11. Differences in population structure
12. How natural selection is acting on particular loci
13. Methods for Applications
14. Outlier - for a given locus, whats the level of differentiation compared to differences across the genome? Using Fst
15. Non outlier - Tests high Fst loci for other features associated with selection
            1. Fitness advantage
            2. Functional enrichment

Command Line notes:

Getting up to date on data processing

opened sam file in my scripts directory

```
cd scripts 
ll
```

Save the tail of the SAM file 

```
tail -n 100 filename > tail.sam
vim tail.sam
:set nowrap
```

Reading a SAM file

Col1: the read/query/name

Col2: Flag (corresponds to information regarding the mapping success)

- 77 all bad
- 141 second in pair and all bad
- 113 this is the first read of a pair, both reads in pair were flipped and both mapped
- 177 this is the second read of a pair, both reads in pair were flipped and both mapped
- 65 this is first read in pair and both reads aligned the forward strand
- 129 This is second read of pair and both reads aligned the forward strand.

Col3: The reference sequence name to which the read mapped

Col4: The left most position in the reference where the read mapped

Col5: Phred score

Col6: CIGAR string that gives alignment information

- how many bases match (M) where there's an insertion (I) or deletion (D)

Col7: an '=', mate position, inferred insert size (col7,8,9)

Col8: The query sequence and Phred quality score from the FASTQ file (cols 10,11)

Col9: information found in the tags at the end, if the read mapped, including whether it is a complete read (XT:A:U), the number of best hits (X0:i:1), the number of suboptimal hits (X1:i:0)



Use the grep command to tell you how many reads mapped uniquely

```
grep -c XT:A:U 08_5-14_S_1_bwaaln.sam
```

5991502 uniquely mapped reads

```
grep -c X0:i:1 08_5-14_S_1_bwaaln.sam

```

6014626 best hits

The below renaming of SAM file needs to be done in the directory of your SAM file (I think it is in my scripts file) 

Going to use regex to rename SAM files

```
sed -i 's/::/\_/g' filename.sam
		#'search/for ::/replace with_ (/ to take literally)/g is an option to do it globally
		
```

Find and copy python script into your scripts directory. Run python script.

```
cd /data/scripts
ll
cp countxpression_PE.py ~/scripts
python countxpression_PE.py 20 35 countstatssummary.txt YOURFILENAME.sam
```



------
<div id='id-section6'/>
### Page 6:

Notes in r script file

------
<div id='id-section7'/>
### Page 7:

Differential Expression Analysis using DESEQ

- Models run today will only run through a random 10% sample of our data
- Notes in annotated script DESEQ2_SSW_round2.R
- ​

------
<div id='id-section8'/>
### Page 8:

**Differential Expression Analysis HW due next Wednesday**

Info Update: WGCNA (weighted gene correlation network analysis)

- An R package applies correlation network methods to describe correlation (co-expression) patterns among genes in  gene expression data
  - Network construction 
    - Each node is a gene
    - Edges are showing how strongly correlated they are by expression
      - Signed/unsigned networks excluding or including positive and negative correlates
    - Input - rows are genes (read counts for these genes) and columns are individuals
    - This matrix is used to evaluate correlations between each set of genes 
    - This leads you to your next step of module identification (genes that are highly correlated to each other and are clustering together)
    - This program uses a weighted analysis
      - has a soft (continuous) threshold for determining whether or not genes are linked. Does this by value, not a dichotomous pass/fail that an unweighted analysis would have
  - Module Identification
    - Unsupervised clustering (not a biased set of genes, doesn't choose input)
    - removes weak edges (correlations)
    - Summarize profiles of modules - Eigengenes which collapses module network into a principal component
  - Relationship of modules to external information
    - Gene significance assigns a positive number with each gene
  - Relationship between/within modules
  - Finding key drivers in modules of interest


------
<div id='id-section9'/>
### Page 9:

Population genomics: 

- SNPs and lots of them 
- Processes: population structure, diversity within populations, selection (positive, negative)

Pipeline: 

Raw reads -> Clean -> Assemble "draft transcriptome" -> Map reads -> Transcriptomics ->

Counts # reads/transcript -> Differential gene expression ... now we are looking at population genomics instead

Take the same mapped reads -> Call SNPs and genotypes -> Allele frequencies/SFS/nucleotide diversity

Complications with calling SNPs:

- Sequencing Error (Illumina 1:100)
  - Filters: minor allele frequencies (how many individuals is it found in), eliminate low frequency reads
    - Depth of coverage

Complications with calling genotypes: Assigning probabilities to each genotype and paralogs 

Diversity: 

- pi: pairwise nucleotide diversity
  - equivalent (at the sequence level) to the expected level of heterozygosity

Working with our data:

- Using VCF tools to filter SNPs by depth/quality/etc
  - We are going to treat repeated samples of each individual as technical replicates and compare them with each other

Looked at the summary of our vcf file

```
vcftools --vcf SSW_bamlist.txt.vcf
```

- After filtering, kept 24 out of 24 Individuals
- After filtering, kept all 7.47M  

Looking at how many sites were unresolved

```
 grep "unres" SSW_bamlist.txt.vcf | wc

```

5.63 million

Looking at how many sites were paralogs

```
grep "para" SSW_bamlist.txt.vcf | wc

```

4,354

Now, we are going to filter for biallelic sites

- highly unlikely sites will have more than two alleles, likely due to sequencing error

```
vcftools --vcf SSW_bamlist.txt.vcf --min-alleles 2 --max-alleles 2

```

After filtering, kept 24 out of 24 Individuals
After filtering, kept 20319 out of a possible 7472775 Sites



Now, we will filter for a minor allele frequency of 0.02 (eliminating sites that only show up once in one individual...?)

```
vcftools --vcf SSW_bamlist.txt.vcf --maf 0.02
```

After filtering, kept 24 out of 24 Individuals
After filtering, kept 5656584 out of a possible 7472775 Sites

Now, filter for individuals that have less than 20% missing data

```
vcftools --vcf SSW_bamlist.txt.vcf --max-missing 0.8

```

After filtering, kept 24 out of 24 Individuals
After filtering, kept 100219 out of a possible 7472775 Sites



Now, output the resulting filtered data as a new vcf file in my data directory

```
vcftools --vcf SSW_bamlist.txt.vcf --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/mydata/biallelic.MAF0.2.Miss0.8

```

Parameters as interpreted:

        --vcf SSW_bamlist.txt.vcf
        --maf 0.02
        --max-alleles 2
        --min-alleles 2
        --max-missing 0.8
        --out /users/l/a/lashlock/mydata/biallelic.MAF0.2.Miss0.8
        --recode

After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 1180 out of a possible 7472775 Sites

Then check it out by going to your mydata directory and vim the file

```
cd ~
cd mydata
ll
vim biallelic.MAF0.2.Miss0.8.recode.vcf
:set nowrap

```



Alright, so now conduct Hardy Weinberg on your filtered file, to see if your loci are in HWE

```
--vcf biallelic.MAF0.2.Miss0.8.recode.vcf --hardy

```

After filtering, kept 24 out of 24 Individuals
Outputting HWE statistics (but only for biallelic loci)

        HWE: Only using fully diploid SNPs.
After filtering, kept 1180 out of a possible 1180 Sites

All in HWE

Head the hwe output file

```
 head out.hwe
```

Switch to R mode (to exit R mode type 'q()')

```
R
```

Read in our output table and look at the structure of the file

```
> hardy <- read.table("out.hwe", header=T)
> str(hardy)
```

'data.frame':   442 obs. of  8 variables:
 $ CHR               : Factor w/ 111 levels "TRINITY_DN35598_c0_g1_TRINITY_DN35598_c0_g1_i1_g.5802_m.5802",..: 65 65 100 100 100 100 100 100 88 88 ...
 $ POS               : int  4566 4665 978 1404 1722 3426 3729 3912 115 141 ...
 $ OBS.HOM1.HET.HOM2.: Factor w/ 27 levels "10/11/3","11/0/13",..: 27 22 27 27 20 27 22 18 18 27 ...
 $ E.HOM1.HET.HOM2.  : Factor w/ 16 levels "10.01/10.98/3.01",..: 14 12 14 14 11 14 12 10 10 14 ...
 $ ChiSq_HWE         : num  0.0109 0.1067 0.0109 0.0109 0.1983 ...
 $ P_HWE             : num  1 1 1 1 1 1 1 1 1 1 ...
 $ P_HET_DEFICIT     : num  1 1 1 1 1 1 1 1 1 1 ...
 $ P_HET_EXCESS      : num  1 0.936 1 1 0.874 ...

Now let's look at the rows in this file that are not in HWE

```
 hardy[which(hardy$P_HET_EXCESS<0.001),]
```

[1] CHR                POS                OBS.HOM1.HET.HOM2. E.HOM1.HET.HOM2.
[5] ChiSq_HWE          P_HWE              P_HET_DEFICIT      P_HET_EXCESS
<0 rows> (or 0-length row.names)

```
hardy[which(hardy$P_HET_DEFICIT<0.001),]
```

                                                                 CHR POS
291 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742  99
293 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742 138
401     TRINITY_DN39079_c3_g1_TRINITY_DN39079_c3_g1_i1_g.8354_m.8354 244
406     TRINITY_DN39696_c4_g1_TRINITY_DN39696_c4_g1_i1_g.8926_m.8926 283
    OBS.HOM1.HET.HOM2. E.HOM1.HET.HOM2. ChiSq_HWE        P_HWE P_HET_DEFICIT
291            11/0/13  5.04/11.92/7.04        24 9.114786e-08  9.114786e-08
293             19/0/5  15.04/7.92/1.04        24 6.498371e-06  6.498371e-06
401            13/0/11  7.04/11.92/5.04        24 9.114786e-08  9.114786e-08
406            13/0/11  7.04/11.92/5.04        24 9.114786e-08  9.114786e-08
    P_HET_EXCESS
291            1
293            1
401            1
406            1













------
<div id='id-section10'/>
### Page 10:

Rate of evolution due to the relationship of effective population size and the substitution rate



4 methods of measuring effective population size

- From species life history
- From variance in allele frequency between populations
- From genetic polymorphism data
- Correlated to body saize

Effective population size varies across species

Signatures in the genome:

- genetic hitchhiking = selective sweep
- Background selection (opposite of hitchhiking elimination of adaptive or neutral alleles due to linkage with deleterious alleles)
- Smaller number of sex chromosomes and larger amount of autosomes

Mutations

- across gene or chromosome
  - Duplication, inversion, deletion, insertion, translocation
- base level (point mutation)
  - Substitutions
    - Transitions : purine for purine and pyrimidine for pyrimidine
    - Transversion: purine for pyrimidine or vice versa
  - Synonymous (silent) mutations
    - Natural selection does not act on these mutations, the only reason they increase in frequency in a population is due to drift
  - Nonsynonymous mutations (missense/replacement)
    - Natural selection acts on these: purifying selection (removal of deleterious alleles) or positive selection (increase in frequency of adaptive alleles)

Classes of Mutations

- Neutral
  - no natural selection, drift acts on these, more drift in larger genomes w=0, <1/Ne
- slightly deleterious
  - Small effect of w
  - Natural selection and drift act here
  - haploid = 1/Ne, diploid = 2Ne
- slightly advantageous
- Deleterious
  - big effect of selection >1/Ne
  - negative (NeRR)
- Advantageous
  - big effect of selection >1/Ne
  - positive (NeRR)

Variation in mutation rate 

- generation time
  - shorter generation time = higher mutation rate
- Selection 
  - lowers mutation rate (controversial)

NeRR and Linkage

- selective sweep
- clonal interference (competition for fixation between adaptive alleles)

Fitness landscapes

- Populations further away from an adaptive peak will (theoretically) have a higher rate of mutation, because they are working towards reaching the adaptive peak
- Populations close to an adaptive peak will (theoretically) have a lower mutation rate because mutations will drift them away from that peak

**Popgen Pt 2**

Path to data:

```
cd /data/project_data/snps/reads2snps/
ll
```

Looking into our zipped vcf file

```
vcftools --gzvcf SSW_byind.txt.vcf.gz

```



After filtering, kept 22 out of 22 Individuals
After filtering, kept 7485987 out of a possible 7485987 Sites



Filtering Data

Filtering for biallelic loci

```
vcftools --gzvcf SSW_byind.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.8 
```

After filtering, kept 22 out of 22 Individuals
Outputting VCF file...
After filtering, kept 5565 out of a possible 7485987 Sites

This output went to my home directory

cd to home directory

```
cd ~/
ll
```

Zip the output file

```
gzip SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf
```

This next term calls vcf tools and reads in new files, and looking for files that match HWE 

```
vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --hardy
```

After filtering, kept 22 out of 22 Individuals
Outputting HWE statistics (but only for biallelic loci)

        HWE: Only using fully diploid SNPs.
After filtering, kept 5565 out of a possible 5565 Sites



Now open R on the terminal

```
R
```

Read file into R

```
hwe <- read.table("out.hwe", header=T)
```

Check out the structure and content of the file

```
str(hwe)
summary(hwe)
```

Now, let's search for snps that deviate from HWE 

```
which(hwe$P_HET_DEFICIT<0.01)
```

[1] 1001 1021 1023 1300 1302 1320 1407 1409

Now, look at all columns

```
hwe[which(hwe$P_HET_DEFICIT<0.01),]

```

                                                                  CHR POS
1001 TRINITY_DN45155_c27_g2_TRINITY_DN45155_c27_g2_i2_g.18743_m.18743 216
1021 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742  99
1023 TRINITY_DN45155_c27_g1_TRINITY_DN45155_c27_g1_i1_g.18742_m.18742 138
1300     TRINITY_DN39079_c3_g1_TRINITY_DN39079_c3_g1_i1_g.8354_m.8354 244
1302     TRINITY_DN39079_c3_g1_TRINITY_DN39079_c3_g1_i1_g.8354_m.8354 279
1320     TRINITY_DN39696_c4_g1_TRINITY_DN39696_c4_g1_i1_g.8926_m.8926 283
1407   TRINITY_DN42225_c1_g1_TRINITY_DN42225_c1_g1_i1_g.12458_m.12458 220
1409   TRINITY_DN42225_c1_g1_TRINITY_DN42225_c1_g1_i1_g.12458_m.12458 255
     OBS.HOM1.HET.HOM2. E.HOM1.HET.HOM2. ChiSq_HWE        P_HWE P_HET_DEFICIT
1001             20/0/2  18.18/3.64/0.18        22 1.701645e-03  1.701645e-03
1021            10/0/12  4.55/10.91/6.55        22 3.671957e-07  3.671957e-07
1023             17/0/5  13.14/7.73/1.14        22 1.061317e-05  1.061317e-05
1300            12/0/10  6.55/10.91/4.55        22 3.671957e-07  3.671957e-07
1302             20/0/2  18.18/3.64/0.18        22 1.701645e-03  1.701645e-03
1320            12/0/10  6.55/10.91/4.55        22 3.671957e-07  3.671957e-07
1407            10/0/12  4.55/10.91/6.55        22 3.671957e-07  3.671957e-07
1409             20/0/2  18.18/3.64/0.18        22 1.701645e-03  1.701645e-03
     P_HET_EXCESS
1001            1
1021            1
1023            1
1300            1
1302            1
1320            1
1407            1
1409            1

Running a comparison between samples in vcf tools

You need to provide the way these samples are differentiated

Navigate back to where you got the vcf file

```
cd /data/project_data/snps/reads2snps/
ll
```

enter into vim for text file

```
vim ssw_healthloc.txt
```

Or you can copy to screen using cat

```
cat ssw_healthloc.txt
```

Grab the healthy individuals and output a txt file to home directory

```
grep "HH" ssw_healthloc.txt > ~/HOneSampPerInd.txt
```

Grab the sick individuals

```
grep "SS" ssw_healthloc.txt > ~/SOneSampPerInd.txt
```

Append (using >>) HS individuals to SS text file 

```
grep "HS" ssw_healthloc.txt >> ~/SOneSampPerInd.txt
```

cd to home dir

```
cd ~/
ll
```

The files are there

VCF tools wants only the sample ID so we will select the column we want

```
cut -f 1 HOneSampPerInd.txt >HOneSampPerInd2.txt
```

cat to see if it worked

```
 cat HOneSampPerInd2.txt
```

Now do it again for the sick file

```
cut -f 1 SOneSampPerInd.txt >SOneSampPerInd2.txt
```

Now use these files to run the vcf allele frequency command

```
 
vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --freq2 --keep HOneSampPerInd2.txt --out H_AlleleFreqs
```

After filtering, kept 6 out of 22 Individuals

```
vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --freq2 --keep SOneSampPerInd2.txt --out S_AlleleFreqs
```

After filtering, kept 14 out of 22 Individuals

Now hop into your output files and edit the headers

- take out freq
- replace with MAJOR [TAB] MINOR
- Do for both files








------
<div id='id-section11'/>
### Page 11:

R script for Homework #2

```
##Differential Expression Assignment
##Ecological Genomics
##Lauren Ashlock
##3_06_2017


library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)
ncol(countData)
dim(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)
dim(colData)


#comparing healthy vs sick individuals within the intertidal


colInt<-subset(colData,location=="int")
colSub <- subset(colData, location=="sub")

nrow(colInt)
#output
#48

countDataInt<-countData[, which(colnames(countData) %in% row.names(colInt))]
countDataSub<-countData[, -which(colnames(countData) %in% row.names(colInt))]
dim(countDataInt)
dim(countDataSub)

ddsInt <- DESeqDataSetFromMatrix(countData = countDataInt, colData = colInt ,design = ~ health)

ddsInt <- ddsInt[ rowSums(counts(ddsInt)) > 100, ]

colData(ddsInt)$health <- factor(colData(ddsInt)$health, levels=c("H","S"))

ddsInt <- DESeq(ddsInt) 

resInt <- results(ddsInt)

resInt <- resInt[order(resInt$padj),]

head(resInt)

summary(resInt)

#output
# out of 12399 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 205, 1.7% 
# LFC < 0 (down)   : 37, 0.3% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 9381, 76% 
# (mean count < 41)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


#Comparing healthy vs sick individuals in subtidal

ddsSub <- DESeqDataSetFromMatrix(countData = countDataSub, colData = colSub ,design = ~ health)

ddsSub <- ddsSub[ rowSums(counts(ddsSub)) > 100, ]

colData(ddsSub)$health <- factor(colData(ddsSub)$health, levels=c("H","S"))

ddsSub <- DESeq(ddsSub) 

resSub <- results(ddsSub)

resSub <- resSub[order(resSub$padj),]

head(resSub)

summary(resSub)


#output
# out of 12392 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 20, 0.16% 
# LFC < 0 (down)   : 113, 0.91% 
# outliers [1]     : 647, 5.2% 
# low counts [2]   : 4289, 35% 
# (mean count < 13)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

#Comparing healthy vs sick controlling for location 

ddsLoc <- DESeqDataSetFromMatrix(countData = countData, colData = colData ,design = ~ location + health)

ddsLoc <- ddsLoc[ rowSums(counts(ddsLoc)) > 100, ]

colData(ddsLoc)$health <- factor(colData(ddsLoc)$health, levels=c("H","S"))

ddsLoc <- DESeq(ddsLoc) 

resLoc <- results(ddsLoc)

resLoc <- resLoc[order(resLoc$padj),]

head(resLoc)

summary(resLoc)

#output
# out of 12947 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 209, 1.6% 
# LFC < 0 (down)   : 65, 0.5% 
# outliers [1]     : 400, 3.1% 
# low counts [2]   : 7679, 59% 
# (mean count < 23)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

#Table of up/down regulated genes in three models
vec<- c(205, 37, 20, 113, 209, 65)

summaryMatrix<- matrix(data=vec, nrow=2, ncol=3, byrow=FALSE)
summaryMatrix
row.names(summaryMatrix)<- c("Significantly Upregulated", "Significantly Downregulated")
summaryMatrix
colnames(summaryMatrix) <- c("Intertidal", "Subtidal", "Full Model")
summaryMatrix

print(summaryMatrix, type="html")

#PCA for all three models

#Intertidal

vsdInt <- varianceStabilizingTransformation(ddsInt, blind=FALSE)

plotPCA(vsdInt, intgroup=c("health"))

#Subtidal

vsdSub <- varianceStabilizingTransformation(ddsSub, blind=FALSE)

plotPCA(vsdSub, intgroup=c("health"))

#Location

vsdLoc <- varianceStabilizingTransformation(ddsLoc, blind=FALSE)

plotPCA(vsdLoc, intgroup=c("health"))

##Plot counts

#Intertidal

## Check out one of the genes to see interaction between score, health and expression....
dInt <- plotCounts(ddsInt, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
dInt
pInt <- ggplot(dInt, aes(x= score, y=count, color = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
pInt <- pInt + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
pInt

pInt <- ggplot(dInt, aes(x=score, y=count, color=health, group=health)) 
pInt <- pInt +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
pInt

#Subtidal

dSub <- plotCounts(ddsSub, gene="TRINITY_DN42073_c0_g1_TRINITY_DN42073_c0_g1_i1_g.12173_m.12173", intgroup=(c("health","score","location")), returnData=TRUE)
dSub
pSub <- ggplot(dSub, aes(x= score, y=count, color = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
pSub <- pSub + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
pSub

pSub <- ggplot(dSub, aes(x=score, y=count, color=health, group=health)) 
pSub <- pSub +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
pSub

#Full model

dLoc <- plotCounts(ddsLoc, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
dLoc
pLoc <- ggplot(dLoc, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
pLoc <- pLoc + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
pLoc

pLoc <- ggplot(dLoc, aes(x=score, y=count, color=health, group=health)) 
pLoc <- pLoc +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
pLoc

citation("DESeq2")
```



------
<div id='id-section12'/>
### Page 12:

Means of identifying global ancestry: 

- Model based: structure (Bayesian), AdMIXTURE (Maximum Likelihood)
  - admixture is faster, alternates between ancestral and gene frequency matrices
- Nonparametric:
  - multivariate analysis
  - clustering
    - pairwise data matrix (distances or similarities)
    - Program
    - Produce phenograms - type of tree based on overall similarity
  - Ordination methods
    - PCA
    - MDS

Means of identifying local ancestry: each chromosome is a mosaic of segments from different ancestral populations. The goal is to identify the population of origin for each position in the genome.

- HMM
  - Structure (but structure doesn't model LD)
  - LAMP (sliding window, and assigning ancestral populations by PCA)
  - RFMix (discriminate function analysis)

Weaknesses:

- Very reliable for looking at n=2 source populations, beyond that you lose accuracy quickly


- Makes the assumption that we know K (number of ancestral populations) and their allele frequencies
  - we rarely know that
  - Can use simulations to get at this 
  - Greatest impact when ancestor is recent

Applications

- Pharmacogenomics
- Map diseases to genes
  - disease linked genes and variation in disease linked genes

Future research and Challenges

- more SNPs, greater sample sizes
  - improvements on current methods
- personalized medicine
- powerful tools for admixture mapping



Coding Time!

Popgen script three

New vcf file with all 24 individuals

Want to divvy up the vcf by healthy and sick into two different text files so we can compare

```
grep "HH" ssw_healthloc.txt | cut -f1 >~/H_SampleIDs.txt
grep "HS\|SS" ssw_healthloc.txt | cut -f1 >~/S_SampleIDs.txt
```

Now go back to your home directory and check out the text files

```
cd ~/
[lashlock@pbio381 ~]$ cat S_SampleIDs.txt
```



Looks good

Now let's calculate allele frequencies for both files

```
 vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --freq2 --keep H_SampleIDs.txt --out H_AlleleFreqs
 vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --freq2 --keep S_SampleIDs.txt --out S_AlleleFreqs

```

Kept all sites



Now, let's calculate Fst between the healthy and sick individuals

```
vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --weir-fst-pop H_SampleIDs.txt --weir-fst-pop S_SampleIDs.txt --out HvS_Fst
```

Now use winSCP to move the allele frequency files and the Fst file to your computer

Edit the headers on the text files, so R doesn't have trouble reading them into R

R script

```
##Looking at allele frequencies 
##Ecological Genomics
##Lauren Ashlock
##3 20 17

#### read in your text files (remember to edit header first)

H_freq <- read.table("H_AlleleFreqs.frq", header=T)
S_freq <- read.table("S_AlleleFreqs.frq", header=T)

str(H_freq)
str(S_freq)

# Since these files have identical numbers of SNPs in the exact same order, 
# we can concatenate them together into one large dataframe:

All_freq <- merge(H_freq, S_freq, by=c("CHROM", "POS"))

# Looks good, now let's calculate the difference in minor allele frequency
# at each SNP and plot as a histogram

All_freq$diff <- (All_freq$H_ALT - All_freq$S_ALT)

hist(All_freq$diff, breaks=50, col="red", main="Allele frequency difference (H-S)")


# Looks like most loci show little difference (i.e., likely drift), 
# but perhaps a few show very large differences between healthy and sick 
# (drift or selection?)

# How do these highly divergent frequenices compare to Fst at the same SNPs?
fst <- read.table("HvS_Fst.weir.fst", header=T)

All_freq.fst <- merge(All_freq, fst, by=c("CHROM", "POS"))

plot(All_freq.fst$diff, All_freq.fst$WEIR_AND_COCKERHAM_FST, xlab="Allele frequency difference (H-S)", ylab="Fst", main="Healthy vs. Sick SNP divergence")

#As I guess we would expect.... with greater allele frequency difference
#there appears to be a greater Fst
#Each point on this plot is a SNP

# Which are the genes that are showing the highest divergence between Healthy and Sick?
All_freq.fst[which(All_freq.fst$WEIR_AND_COCKERHAM_FST>0.2),]

points(0.2500000,0.2415350,col="red", fill="red", cex=2)

#This modifies a point on the plot you currently have open
```







------
<div id='id-section13'/>
### Page 13:

Start by looking at nucleotide diversity at synonymous sites and the ration of piN/piS

```
cd /data/project_data/snps/reads2snps
head SSW_by24inds.txt.fas
tail SSW_by24inds.txt.fas

```



The next command will takes a long time so set screen and run in background 

```
screen
/data/popgen/dNdSpiNpiS_1.0 -alignment_file=SSW_by24inds.txt.fas -ingroup=sp -out=~/dNdSpiNpiS_output

```

This command calls the program that we will use, what file to use, "ingroup=sp" means we only have 1 species (no outgroup), out is what we want to save the file as

Detach from the screen

```
Ctrl+A+D

```

Let's look at an example while we wait...

```
cat SSW_bamlist.txt.sum

```

Interpreting this file:

- High values mean more homozygosity (vary from 0-1)
- Our population looks like it is randomly mating

```
piS: 0.00585312 [0.005172; 0.006598]
piN: 0.00154546 [0.00133; 0.001782]
ave. piN/piS: 0.264041 [0.223914; 0.310575]

```

Grab the Romiguier file and move to desktop to open in R

```
Rom <- read.csv("Romiguier_nature13685-s3.csv", header=T)

```

Data check

```
str(Rom) 
head(Rom)
```

Now lets make a plot that shows the purifying selection vs effective pop size

```
plot(log(Rom$piS), log(Rom$piNpiS), pch=21, bg="blue", xlab="log Synonymous Nucleotide Diversity (piS)", ylab="log Ratio of Nonysn to Syn Diversity (piN/piS)", main="Purifying Selection vs. Effective Population Size")

```

you can add points to any plot by using this points command; we will add our values to the plot now:

```
points(log(0.00585312), log(0.264041), pch=24, cex=1.5, bg="red") 

```

so we added the PiS value and then the piN/piS then told it what we wanted the point to look at

the sea stars fall in the middle

Now we will add a best fit line; this is just using Rom data

```
reg <- lm(log(Rom$piNpiS) ~ log(Rom$piS)) # Fits a linear regression
abline(reg) # adds the regression line to the plot

```

It would be nice to know which data points are similar species (other echinoderms) to our sea stars

```
echino <- Rom[which(Rom$Phylum=="Echinodermata"),] # subsets the data
points(log(echino$piS), log(echino$piNpiS), pch=21, bg="red") # adds the points

```

They are fairly spread out but there is a cluster around our sea star 

```
legend("bottomleft", cex=1, legend=c("Metazoans", "Echinoderms", "P. ochraceus"), pch=c(21,21,24), col=c("blue", "red", "red"))
```





------
<div id='id-section14'/>
### Page 14:

Info Update: Species divergence with gene flow



- Allopatric speciation - speciation due to physical isolation and subsequent loss of gene flow
- Sympatric speciation - divergence in the presence of gene flow due to diversifying selection
  - selected genes will be divergent
  - neutral alleles appear homogenous

Inferring history of divergence

- Genomic scans
  - Islands of differentiation
    - looking at the distribution of summary statistics that measure differentiation (Fst)
    - genomic region with high Fst is an indicator that it is a region under selection
- Gene vs population trees
  - compare assumed population trees to gene trees
  - compare gene trees
  - D statistic to determine introgression using ABBA-BABA
  - No introgression: D=0 (ABBA=BABA)
  - Introgression: != 0

Limitations

- Throw out data
- Requires many genomes
- Same Fst values can be interpreted in multiple ways

Likelihood/model-based methods

- Allele frequency spectrum
  - uses count data to generate a distribution of allele frequencies
  - neutral/bottleneck negative exponential
  - selective sweep u shape
  - Assumptions
    - SNPs are independent
    - Free recombination among SNPs
    - Mutation rates are equal
  - Limitations
    - computationally challenging
    - lose genomic data by only looking at SNPs
    - Expensive for models with more than three populations
- Genealogy sampling
  - multiple regions, one gene tree
  - From this you can estimate effective population size, mutation rate, and admixture
  - Assumptions
    - Free recombination among genes
    - Complete linkage with loci
    - Mutation rates vary across the genome
    - No recombination since common ancestor
- Likelihood free method
  - Approximate Bayesian Computation
    - simulations under model of interest

Historical Gene Flow and LD Patterns

- Distribution of haplotype lengths
  - If you identify a migrant population, you can measure the time since first migration by looking at the fragments of linked blocks of the genome and their change in length over time
- Approximation of conditional likelihoods
  - Ancestral recombination graphs (ARGs)
    - looks like a gene tree, but instead the tree is based on recombination events
  - Limitations
    - complex
    - difficult to ID correct ARG

NGS advantages and disadvantages

- good estimation of recombination rates
- large area for genome scans
- Ascertainment bias
- Throw out data
- Computationally challenging

Tutorial

- We'll use the **piNpiS** script from Gayral et al. (2013) to run this. We only need a single input file, which is a FASTA formatted sequence file that is output from **reads2snps**, and we'll save the output to our home directories:

```
cd /data/project_data/snps/reads2snps
/data/popgen/dNdSpiNpiS_1.0 -alignment_file=SSW_by24inds.txt.fas -ingroup=sp -out=~/dNdSpiNpiS_output
```



- This command takes a long time to run so before you run , set scree, then run the code, then detach from the screen

```
screen
/data/popgen/dNdSpiNpiS_1.0 -alignment_file=SSW_by24inds.txt.fas -ingroup=sp -out=~/dNdSpiNpiS_output
crtl ad
```

- If you want to reattach...

```
screen -r
crtl ad
```

- While we let this run, let's take a look at an estimate already made from a subsample of our data (only one script per individual)

```
cat SSW_bamlist.txt.sum
```

- From this you get a biological summary for all of your samples
  - take a look at Fit, and piS/piN
  - Fit: Average Fit: -0.0507419 [-0.06817; -0.031933] (positive is high homozygosity)
  - Average piS in focal species: 0.00585312 [0.005172; 0.006598]
    Average piN in focal species: 0.00154546 [0.00133; 0.001782]
    Average piN / average piS: 0.264041 [0.223914; 0.310575]
    - lower diversity in nonsynonymous sites is an indicator of purifying selection removing deleterious alleles caused by new mutations
    - piN/piS is really an estimate of how strongly selection is acting to remove deleterious alleles
    - To compare our values to other metazoans.... download info from Romiguier paper

```
/data/project_data/snps/reads2snps/Romiguier_nature13685-s3.csv
```

- load onto laptop using winscp
- then load it into R
- Ne=365,625 ...

R script

```
# Read in the Romiguier data:
Rom <- read.csv("Romiguier_nature13685-s3.csv", header=T)

# Import OK?
str(Rom) 
head(Rom)

# Looks good!
# Now let's look at how the strength of purifying selection (piN/piS) 
# compares to the size of Ne (piS). We'll plot these on a log scale to 
# linearize the relationship.
plot(log(Rom$piS), log(Rom$piNpiS), pch=21, bg="blue", xlab="log Synonymous Nucleotide Diversity (piS)", ylab="log Ratio of Nonysn to Syn Diversity (piN/piS)", main="Purifying Selection vs. Effective Population Size")

# Now let's add our SSW points to the existing plot and give them a different symbol
points(log(0.00585312), log(0.264041), pch=24, cex=1.5, bg="red") 

# We can also add a regression line to the plot to see how far off the SSW estimates are from expectation
reg <- lm(log(Rom$piNpiS) ~ log(Rom$piS)) # Fits a linear regression
abline(reg) # adds the regression line to the plot

# It would be useful to highlight the other echinoderms in the dataset...do our seastars behave similarly?
echino <- Rom[which(Rom$Phylum=="Echinodermata"),] # subsets the data
points(log(echino$piS), log(echino$piNpiS), pch=21, bg="red") # adds the points

# Lastly, let's add a legend:
legend("bottomleft", cex=1, legend=c("Metazoans", "Echinoderms", "P. ochraceus"), pch=c(21,21,24), col=c("blue", "red", "red"))

# Pisaster seems to be in a group with other echinoderms that have relaxed 
#purifying selection (high piN/piS), given their Ne...Interesting!
#Can we hypothesize why this might be?
```





------
<div id='id-section15'/>
### Page 15:

Pr(G|K,Q,P)

- Pr probability of given genotype (G) such that
- K is the number of populations
- Q is the proportion of the genotype that relates to each populations
- P is the allele frequency

Cross validation

- build a model based on a subset of individuals and cross check your model for the individuals you left out
- Can use these data to perform a cross validation analysis and select the model that has the least amount of error

There is a particular format to run ADMIXTURE

- PGDSpider will convert your files
- Steve already did this...

```
/data/project_data/snps/reads2snps/SSW_tidal.pops

/data/project_data/snps/reads2snps/vcf2admixture_SSW.spid

/data/project_data/snps/reads2snps/vcf2geno.sh
```

Copy all of these files to your home directory using cp 

```
 cp SSW_tidal.pops ~/
[lashlock@pbio381 reads2snps]$ cp vcf2admixture_SSW.spid ~/
[lashlock@pbio381 reads2snps]$ cp vcf2geno.sh ~/
[lashlock@pbio381 reads2snps]$ cd ~/

```

Use vim to open the vcf2geno.sh

```
vim vcf2geno.sh
```



Edit the bash script with your file names

```
java -Xmx512M -jar /data/popgen/PGDSpider_2.0.9.0/PGDSpider2-cli.jar -inputfile ./SSW_all_biallelic.MAF0.02Miss0.8.recode.vcf -inputformat VCF -outputfile ./SSW_all_biallelic.MAF0.02Miss0.8.recode.vcf.geno -outputformat EIGENSOFT -spid ./vcf2admixture_SSW.spid
```

Then run your updated bash script

```
./vcf2geno.sh
```

Copy ADMIX bash script to your home directory and vim to open it and look at the script

```
cp /data/project_data/snps/reads2snps/ADMIX.sh ~/
vim ADMIX.sh

#!/bin/bash

# Run ADMIXTURE to determine the number of genetic clusters in the SNP data,
# and the ancestry proportions of each individual

# Here's the utility of 'for loops'...

for K in {1..10}

do

admixture -C 0.000001 --cv ./SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno $K \
| tee log${K}.out

done

# After the for loop finishes, you can use 'grep' to grab the values of the CV from each separate log file and append them into a new summary text file.

"ADMIX.sh" 22L, 501C        
```





  ```


  

  ```

------
<div id='id-section16'/>
### Page 16:

Concepts:

- Inbreeding produces structured populations
- Selective sweeps change allele frequencies in populations
- Empirical p-values created from distribution of putatively neutral loci are super useful
- Methods - OutFlank

What challenges do outlier detection methods face?

Hoe is LD our friend and foe?



F statistics

- Inbreeding coefficient (probability that any two alleles in a population are identical by descent)

- look at heterozygosity

  - at the level of the individual
  - level of subpopulations
  - level of all populations

- Fst = Htotal population-Hsubpopulation/Htotal population

- Fis = expected heterozygosity - observed heterozygosity/expected heterozygosity

- Fit = Ht-Hi/Ht

  COMPUTING NOTES

  R script 

  ```\
  #OutFLANK
  #ecological genomics
  #4 3 17

  #Install packages
  install.packages("devtools")
  library(devtools)
  source("http://bioconductor.org/biocLite.R")
  biocLite("qvalue")
  install_github("whitlock/OutFLANK")
  install.packages("tcltk")
  install.packages("stringi")
  library("stringi")
  a

  #load these packages
  library(OutFLANK)
  library(vcfR)
  library(adegenet)
  library(tcltk)
  library(tibble)
  library(qvalue)
  ```

# Read in your geno file

  ssw.geno_in <- read.fwf("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno", 

```
                      width=rep(1,24))
```

  ssw.geno <- t(ssw.geno_in)

# read in the metadata

  ssw_meta <- read.table("ssw_healthloc.txt", header=TRUE)
  ssw_meta <- ssw_meta[order(ssw_meta$Individual),]
  ssw_meta$Trajectory[which(ssw_meta$Trajectory=='MM')] = NA

  OF_SNPs <- MakeDiploidFSTMat(ssw.geno,locusNames=seq(1,5317,1) , 

```
                           popNames=ssw_meta$Trajectory)
```

  dim(OF_SNPs)

  head(OF_SNPs)

  OF_out <- OutFLANK(FstDataFrame = OF_SNPs, LeftTrimFraction = 0.05, 

```
                 RightTrimFraction = 0.05, Hmin = 0.1, NumberOfSamples = 3, 
                 qthreshold = 0.1)
```

  OutFLANKResultsPlotter(OF_out, withOutliers = T, NoCorr = T, Hmin = 0.1, 

```
                     binwidth = 0.005, titletext = "Scan for local selection")
```

# find your outliers

  outliers <- which(OF_out$results$OutlierFlag=="TRUE")
  outliers

# we can extract info about the outliers by reading in the vcf file and looking at the annotations

  vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
  vcfann <- as.data.frame(getFIX(vcf1))
  vcfann[outliers,]

```
  Command line 

```

  [lashlock@pbio381 ~]$ cd /data/project_data/snps/reads2snps/
  [lashlock@pbio381 reads2snps]$ ll
  total 436360
  -rwxr-xr-x. 1 srkeller users       501 Mar 26 23:59 ADMIX.sh
  drwxr-xr-x. 3 srkeller users      4096 Mar  8 00:28 bam_file_archive
  drwxr-xr-x. 2 srkeller users      4096 Mar 19 07:45 old_snps_calls
  -rw-r--r--. 1 srkeller users       283 Mar 31 13:32 out.log
  -rw-r--r--. 1 srkeller users     14342 Mar  7 22:16 Romiguier_nature13685-s3.csv
  -rw-r--r--. 1 srkeller users    132925 Mar 22 00:36 SSW_all_biallelic.MAF0.02.Mi                                                                                  ss0.8.recode.vcf.geno
  -rw-r--r--. 1 srkeller users      1636 Mar 19 12:42 SSW_bamlist.txt.sum
  -rw-r--r--. 1 srkeller users      1176 Mar 13 22:59 SSW_by24inds.txt
  -rw-r--r--. 1 srkeller users 414372240 Mar 19 05:21 SSW_by24inds.txt.fas
  -rw-r--r--. 1 srkeller users  32258162 Mar 19 05:21 SSW_by24inds.txt.vcf.gz
  -rwxrwxr-x. 1 srkeller users       324 Mar 19 08:01 ssw_healthloc.txt
  -rw-r--r--. 1 srkeller users       168 Mar 27 00:25 SSW_tidal.pops
  -rwxrwxr-x. 1 srkeller users        87 Mar  6 17:22 unique_inds.txt
  -rwxr-xr-x. 1 srkeller users      1604 Mar 27 00:24 vcf2admixture_SSW.spid
  -rwxr-x--x. 1 srkeller users       223 Mar 27 00:24 vcf2geno.sh
  [lashlock@pbio381 reads2snps]$ cd ~/
  [lashlock@pbio381 ~]$ ll
  total 2928
  -rwxr-xr-x. 1 lashlock users    532 Mar 29 11:12 ADMIX.sh
  -rw-r--r--. 1 lashlock users      0 Mar 29 11:12 chooseK.txt
  drwxr-xr-x. 2 lashlock users   4096 Feb  6 14:16 fastqc_out
  -rw-r--r--. 1 lashlock users 452371 Mar 20 10:40 H_AlleleFreqs.frq
  -rw-r--r--. 1 lashlock users    423 Mar 20 10:40 H_AlleleFreqs.log
  -rw-r--r--. 1 lashlock users     24 Mar  8 11:15 HOneSampPerInd2.txt
  -rw-r--r--. 1 lashlock users     96 Mar  8 11:08 HOneSampPerInd.txt
  -rw-r--r--. 1 lashlock users     24 Mar 20 10:33 H_SampleIDs.txt
  -rw-r--r--. 1 lashlock users    642 Mar 20 10:42 HvS_Fst.log
  -rw-r--r--. 1 lashlock users 114992 Mar 20 10:42 HvS_Fst.weir.fst
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log10.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log1.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log2.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log3.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log4.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log5.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log6.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log7.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log8.out
  -rw-r--r--. 1 lashlock users    674 Mar 29 11:12 log9.out
  drwxr-xr-x. 3 lashlock users   4096 Mar  6 11:07 mydata
  -rw-r--r--. 1 lashlock users 207911 Mar  8 10:37 out.hwe
  -rw-r--r--. 1 lashlock users    300 Mar  8 11:25 out.log
  -rw-r--r--. 1 lashlock users    332 Mar 29 10:58 PGDSpider-cli.log
  -rw-r--r--. 1 lashlock users 483033 Mar 20 10:40 S_AlleleFreqs.frq
  -rw-r--r--. 1 lashlock users    424 Mar 20 10:40 S_AlleleFreqs.log
  drwxr-xr-x. 2 lashlock users   4096 Feb 17 13:36 scripts
  -rw-r--r--. 1 lashlock users     42 Mar  8 11:16 SOneSampPerInd2.txt
  -rw-r--r--. 1 lashlock users    168 Mar  8 11:10 SOneSampPerInd.txt
  -rw-r--r--. 1 lashlock users     42 Mar 20 10:37 S_SampleIDs.txt
  -rw-r--r--. 1 lashlock users    436 Mar  8 10:33 SSW_all_biallelic.MAF0.02.Miss0                                                                                  .8.log
  -rw-r--r--. 1 lashlock users    438 Mar 20 10:27 SSW_all_biallelic.MAF0.02Miss0.                                                                                  8.log
  -rw-r--r--. 1 lashlock users 943565 Mar 20 10:27 SSW_all_biallelic.MAF0.02Miss0.                                                                                  8.recode.vcf
  -rw-r--r--. 1 lashlock users 132925 Mar 29 10:58 SSW_all_biallelic.MAF0.02Miss0.                                                                                  8.recode.vcf.geno
  -rw-r--r--. 1 lashlock users  76352 Mar  8 10:33 SSW_all_biallelic.MAF0.02.Miss0                                                                                  .8.recode.vcf.gz
  -rw-r--r--. 1 lashlock users    312 Mar 29 10:58 SSW.ind
  -rw-r--r--. 1 lashlock users 445407 Mar 29 10:58 SSW.snp
  -rw-r--r--. 1 lashlock users    168 Mar 29 10:38 SSW_tidal.pops
  -rwxr-xr-x. 1 lashlock users   1604 Mar 29 10:39 vcf2admixture_SSW.spid
  -rwxr-x--x. 1 lashlock users    274 Mar 29 10:59 vcf2geno.sh
  [lashlock@pbio381 ~]$ vim SSW_all_biallelic.MAF0.02Miss0.8.recode.vcf
  [lashlock@pbio381 ~]$ cd /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds
  -bash: cd: /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds: Not a directory
  [lashlock@pbio381 ~]$ vim /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds

------
<div id='id-section17'/>
### Page 17:

- Copied original unfiltered vcf file into my home directory.... will try a few different filtering approaches
- filter for biallelic snps

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --min-alleles 2 --max-alleles 2

After filtering, kept 24 out of 24 Individuals
After filtering, kept 57180 out of a possible 7486938 Sites
```

- Filter for a minor allele frequency of 0.02

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --maf 0.02

After filtering, kept 24 out of 24 Individuals
After filtering, kept 3441637 out of a possible 7486938 Sites
```

- Max missing of 20%

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --max-missing 0.8

After filtering, kept 24 out of 24 Individuals
After filtering, kept 353902 out of a possible 7486938 Sites
```

- Max missing of 10%

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --max-missing 0.9

After filtering, kept 24 out of 24 Individuals
After filtering, kept 261702 out of a possible 7486938 Sites
Run Time = 103.00 seconds
```

- Max missing of 15%

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --max-missing 0.85

After filtering, kept 24 out of 24 Individuals
After filtering, kept 307795 out of a possible 7486938 Sites
```

- More stringent minor allele frequency

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --maf 0.04   --maf 0.04

After filtering, kept 24 out of 24 Individuals
After filtering, kept 3434312 out of a possible 7486938 Sites
```

- First data set will filter by
  - maf 0.04
  - biallelic loci
  - max missing of 0.85

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --min-alleles 2 --max-alleles 2 --maf 0.04 --max-missing 0.85

After filtering, kept 24 out of 24 Individuals
After filtering, kept 2005 out of a possible 7486938 Sites
```

- Second set...
  - take out hwe deviants

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --min-alleles 2 --max-alleles 2 --maf 0.04 --max-missing 0.85 --hwe 0.05

After filtering, kept 24 out of 24 Individuals
After filtering, kept 1916 out of a possible 7486938 Sites

```



- Make two separate recode files 

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --remove-indv 37 --remove-indv 38 --min-alleles 2 --max-alleles 2 --maf 0.04 --max-missing 0.85 --recode --out ~/filteredSNPS1.0

After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 1945 out of a possible 7486938 Sites
```

```
[lashlock@pbio381 ~]$ vcftools --vcf SSW_by24inds.txt.vcf --remove-indv 37 --remove-indv 38 --min-alleles 2 --max-alleles 2 --maf 0.04 --max-missing 0.85 --hwe 0.05 --recode --out ~/filteredSNPS2.0

After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 1861 out of a possible 7486938 Sites
```

- Now, time to analyze your vcf files
- Calculate allele frequencies for healthy and sick individuals

```
 vcftools --vcf filteredSNPS1.0.recode.vcf --freq2 --keep H_SampleIDs.txt --out H_alleleFreqs1.0
 vcftools --vcf filteredSNPS1.0.recode.vcf --freq2 --keep S_SampleIDs.txt --out S_alleleFreqs1.0
 vcftools --vcf filteredSNPS2.0.recode.vcf --freq2 --keep H_SampleIDs.txt --out H_alleleFreqs2.0
 vcftools --vcf filteredSNPS2.0.recode.vcf --freq2 --keep S_SampleIDs.txt --out S_alleleFreqs2.0

```

- Calculate Fst for healthy and sick individuals

```
vcftools --vcf filteredSNPS1.0.recode.vcf --weir-fst-pop H_SampleIDs.txt --weir-fst-pop S_SampleIDs.txt --out HvS_Fst_1.0

vcftools --vcf filteredSNPS2.0.recode.vcf --weir-fst-pop H_SampleIDs.txt --weir-fst-pop S_SampleIDs.txt --out HvS_Fst_2.0
```

```
Next step is to recreate your vcf files removing the mm individuals

--remove-indivs ## (corresponds to numbers in text file)
```



------
<div id='id-section18'/>
### Page 18:

R Script for Assignment 3

```
##Ecological Genomics Assignment 3
##Lauren Ashlock

install.packages("vcfR")
install.packages("adegenet")

library(adegenet)
library(vcfR)

#read in vcf files

vcf1.0 <- read.vcfR("filteredSNPS1.0.recode.vcf")

vcf2.0 <- read.vcfR("filteredSNPS2.0.recode.vcf")

gl1.0 <- vcfR2genlight(vcf1.0)

gl2.0 <- vcfR2genlight(vcf2.0)

ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort it by Individual ID

#data check
gl1.0$ind.names
ssw_meta$Individual
ssw_meta$Trajectory
str(gl1.0$other)


gl1.0$pop <- ssw_meta$Location # assign locality info
gl1.0$other <- as.list(ssw_meta$Trajectory) # assign disease status

gl2.0$pop <- ssw_meta$Location # assign locality info
gl2.0$other <- as.list(ssw_meta$Trajectory) # assign disease status

# Now, let's compute the PCA on the SNP genotypes and plot it:
pca1.0 <- glPca(gl1.0, nf=4, parallel=F) # nf = number of PC axes to retain (here, 4)

pca2.0 <- glPca(gl2.0, nf=4, parallel=F) # nf = number of PC axes to retain (here, 4)


# Perhaps we want to show disease status instead of locality:
plot(pca1.0$scores[,1], pca1.0$scores[,2], 
     cex=2, pch=20, col=as.factor(ssw_meta$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on Dataset 1.0 (Freq missing=15%; 1945 SNPs)")
legend("topleft", 
       legend=unique(ssw_meta$Trajectory), 
       pch=20, 
       col=as.factor(unique(ssw_meta$Trajectory)))

plot(pca2.0$scores[,1], pca2.0$scores[,2], 
     cex=2, pch=20, col=as.factor(ssw_meta$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on Dataset 2.0 (Freq missing=15%; 1861 SNPs)")
legend("topleft", 
       legend=unique(ssw_meta$Trajectory), 
       pch=20, 
       col=as.factor(unique(ssw_meta$Trajectory)))



# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca1.0$loadings[,1]),
            threshold=quantile(abs(pca1.0$loadings), 0.999), main="Loading Plot PCA 1.0")

# Get their locus names
gl1.0$loc.names[which(abs(pca1.0$loadings)>quantile(abs(pca1.0$loadings), 0.999))]

# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca2.0$loadings[,1]),
            threshold=quantile(abs(pca2.0$loadings), 0.999), main="Loading Plot PCA 2.0")

# Get their locus names
gl2.0$loc.names[which(abs(pca2.0$loadings)>quantile(abs(pca2.0$loadings), 0.999))]
```



------
<div id='id-section19'/>
### Page 19:

Use 16s Amplicon sequencing to perform metagenomic analyses

- What microbes are present?
- What are they doing?

Microbiome: Assemblage of microbial taxa associated with a host or environment

"Who is there?"

- pathogens
- symbionts (mutualists/commensals)

Holobiont: the host together with its microbiome

How to measure a microbiome?

1. Gene for 16s rRNA (forms small subunit of ribosome)
2. ITS: internal transcribed spacer
   1. spacer between genes that code for subunits of rRNA

Metagenomics: Looks at the genes being expressed by the microbiome. Focuses on function

"What are they doing?"

1. RNAseq or shot-gun DNA-seq

Processing:

1. Sample tissue and extract DNA/RNA
2. PCR target gene (16s - only found in bacteria and archaea) 
3. Barcode amplified fragments and ligate adapters and send for sequencing
4. Cluster sequences based on similarity (97% similarity within OTU)
5. Blast OTUs to database, determine taxonomy

Microbial revolution: most microbes aren't culturable. Only studying microbes that can be cultured... biased. Sequencing has demonstrated the diversity of microbes that have not previously been described. 



Big Questions:

1. Is there a core microbiome shared in common among organisms?
2. How much diversity is unique to individuals, populations, communities?
3. Heritability?
4. Evolution of microbiome? 
5. Adaptive vs maladaptive microbiomes?
6. Can patterns within microbiome predict phenotype?

Coding session

```
login as: lashlock
lashlock@pbio381.uvm.edu's password:
Last login: Wed Apr  5 17:17:57 2017 from ip177051.uvm.edu
[lashlock@pbio381 ~]$ mkdir 16s_analysis
[lashlock@pbio381 ~]$ cd 16s_analysis/
[lashlock@pbio381 16s_analysis]$ cd /data/project_data/16s/
[lashlock@pbio381 16s]$ ll
total 40
drwxr-xr-x. 2 mlloyd users    44 Mar 30 10:45 core_diversity_analysis
drwxr-xr-x. 2 mlloyd users  4096 Mar 29 09:32 data_files
-rw-r--r--. 1 mlloyd users 10516 Mar 27 10:47 map.txt
drwxr-xr-x. 8 mlloyd users  4096 Mar 28 12:46 otu_table
-rw-r--r--. 1 mlloyd users    47 Mar 28 14:42 remove-R1.sh
-rw-r--r--. 1 mlloyd users    40 Mar 28 14:24 remove-underscore.sh
-rw-r--r--. 1 mlloyd users  8733 Mar 29 10:22 R_map.txt
[lashlock@pbio381 16s]$ head map.text
head: cannot open ‘map.text’ for reading: No such file or directory
[lashlock@pbio381 16s]$ head map.txt
#SampleID       BarcodeSequence LinkerPrimerSequence    Day     individual     P                                                                                       henotype        Pheno_num       Final_phenotype Trajectory      Tide    Descript                                                                                       ion
01-5-05                 0       1       Sick    1       Sick    SS      intertid                                                                                       al      1-5/5
01-5-08                 3       1       Sick    3       Sick    SS      intertid                                                                                       al      1-5/8
01-5-11                 6       1       Dead    5       Sick    SS      intertid                                                                                       al      1-5/11
02-5-05                 0       2       Sick    1       Sick    SS      intertid                                                                                       al      2-5/5
02-5-08                 3       2       Sick    2       Sick    SS      intertid                                                                                       al      2-5/8
02-5-11                 6       2       Sick    2       Sick    SS      intertid                                                                                       al      2-5/11
02-5-14                 9       2       Sick    2       Sick    SS      intertid                                                                                       al      2-5/14
02-5-17                 12      2       Sick    2       Sick    SS      intertid                                                                                       al      2-5/17
02-5-20                 15      2       Dead    5       Sick    SS      intertid                                                                                       al      2-5/20
[lashlock@pbio381 16s]$ :set nowrap
-bash: :set: command not found
[lashlock@pbio381 16s]$ vim map.txt
[lashlock@pbio381 16s]$ cd ~/
[lashlock@pbio381 ~]$ cd 16s_analysis
[lashlock@pbio381 16s_analysis]$ ll
total 0
[lashlock@pbio381 16s_analysis]$ validate_mapping_file.py -m /data/project_data/16s/map.txt -o validate_map -p -b
Errors and/or warnings detected in mapping file.  Please check the log and html file for details.
[lashlock@pbio381 16s_analysis]$ ll
total 0
drwxr-xr-x. 2 lashlock users 96 Apr 10 10:23 validate_map
[lashlock@pbio381 16s_analysis]$ cd validate_map/
[lashlock@pbio381 validate_map]$ ll
total 524
-rw-r--r--. 1 lashlock users  10180 Apr 10 10:23 map_corrected.txt
-rw-r--r--. 1 lashlock users 398821 Apr 10 10:23 map.html
-rw-r--r--. 1 lashlock users  68678 Apr 10 10:23 map.log
-rw-r--r--. 1 lashlock users  50732 Apr 10 10:23 overlib.js
[lashlock@pbio381 validate_map]$ multiple_join_paired_ends.py -i /data/project_data/16s/data_files -o ~/16s_analysis/joined --read1_indicator _R1 --read2_indicator _R2
[lashlock@pbio381 validate_map]$ ll
total 524
-rw-r--r--. 1 lashlock users  10180 Apr 10 10:23 map_corrected.txt
-rw-r--r--. 1 lashlock users 398821 Apr 10 10:23 map.html
-rw-r--r--. 1 lashlock users  68678 Apr 10 10:23 map.log
-rw-r--r--. 1 lashlock users  50732 Apr 10 10:23 overlib.js
[lashlock@pbio381 validate_map]$ cd ~/
[lashlock@pbio381 ~]$ cd 16s_analysis/
[lashlock@pbio381 16s_analysis]$ ll
total 4
drwxr-xr-x. 26 lashlock users 4096 Apr 10 10:39 joined
drwxr-xr-x.  2 lashlock users   96 Apr 10 10:23 validate_map
[lashlock@pbio381 16s_analysis]$ cd ~/
[lashlock@pbio381 ~]$ ll
total 1338280
drwxr-xr-x. 4 lashlock users         50 Apr 10 10:33 16s_analysis
-rwxr-xr-x. 1 lashlock users        532 Mar 29 11:12 ADMIX.sh
-rw-r--r--. 1 lashlock users          0 Mar 29 11:12 chooseK.txt
drwxr-xr-x. 2 lashlock users       4096 Feb  6 14:16 fastqc_out
-rw-r--r--. 1 lashlock users        467 Apr  5 17:33 filteredSNPS1.0.log
-rw-r--r--. 1 lashlock users     331705 Apr  5 17:33 filteredSNPS1.0.recode.vcf
-rw-r--r--. 1 lashlock users        479 Apr  5 17:35 filteredSNPS2.0.log
-rw-r--r--. 1 lashlock users     317278 Apr  5 17:35 filteredSNPS2.0.recode.vcf
-rw-r--r--. 1 lashlock users     168185 Apr  4 15:19 H_alleleFreqs1.0.frq
-rw-r--r--. 1 lashlock users        377 Apr  4 15:19 H_alleleFreqs1.0.log
-rw-r--r--. 1 lashlock users     161084 Apr  4 15:21 H_alleleFreqs2.0.frq
-rw-r--r--. 1 lashlock users        377 Apr  4 15:21 H_alleleFreqs2.0.log
-rw-r--r--. 1 lashlock users     452371 Mar 20 10:40 H_AlleleFreqs.frq
-rw-r--r--. 1 lashlock users        423 Mar 20 10:40 H_AlleleFreqs.log
-rw-r--r--. 1 lashlock users         24 Mar  8 11:15 HOneSampPerInd2.txt
-rw-r--r--. 1 lashlock users         96 Mar  8 11:08 HOneSampPerInd.txt
-rw-r--r--. 1 lashlock users         24 Mar 20 10:33 H_SampleIDs.txt
-rw-r--r--. 1 lashlock users        597 Apr  4 15:24 HvS_Fst_1.0.log
-rw-r--r--. 1 lashlock users      50661 Apr  4 15:24 HvS_Fst_1.0.weir.fst
-rw-r--r--. 1 lashlock users        596 Apr  4 15:25 HvS_Fst_2.0.log
-rw-r--r--. 1 lashlock users      46389 Apr  4 15:25 HvS_Fst_2.0.weir.fst
-rw-r--r--. 1 lashlock users        642 Mar 20 10:42 HvS_Fst.log
-rw-r--r--. 1 lashlock users     114992 Mar 20 10:42 HvS_Fst.weir.fst
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log10.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log1.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log2.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log3.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log4.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log5.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log6.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log7.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log8.out
-rw-r--r--. 1 lashlock users        674 Mar 29 11:12 log9.out
drwxr-xr-x. 3 lashlock users       4096 Mar  6 11:07 mydata
-rw-r--r--. 1 lashlock users     207911 Mar  8 10:37 out.hwe
-rw-r--r--. 1 lashlock users        327 Apr  4 12:34 out.log
-rw-r--r--. 1 lashlock users        332 Mar 29 10:58 PGDSpider-cli.log
-rw-r--r--. 1 lashlock users     178095 Apr  4 15:20 S_alleleFreqs1.0.frq
-rw-r--r--. 1 lashlock users        378 Apr  4 15:20 S_alleleFreqs1.0.log
-rw-r--r--. 1 lashlock users     170327 Apr  4 15:22 S_alleleFreqs2.0.frq
-rw-r--r--. 1 lashlock users        378 Apr  4 15:22 S_alleleFreqs2.0.log
-rw-r--r--. 1 lashlock users     483033 Mar 20 10:40 S_AlleleFreqs.frq
-rw-r--r--. 1 lashlock users        424 Mar 20 10:40 S_AlleleFreqs.log
drwxr-xr-x. 2 lashlock users       4096 Feb 17 13:36 scripts
-rw-r--r--. 1 lashlock users         42 Mar  8 11:16 SOneSampPerInd2.txt
-rw-r--r--. 1 lashlock users        168 Mar  8 11:10 SOneSampPerInd.txt
-rw-r--r--. 1 lashlock users         42 Mar 20 10:37 S_SampleIDs.txt
-rw-r--r--. 1 lashlock users        436 Mar  8 10:33 SSW_all_biallelic.MAF0.02.Miss0.8.log
-rw-r--r--. 1 lashlock users        438 Mar 20 10:27 SSW_all_biallelic.MAF0.02Miss0.8.log
-rw-r--r--. 1 lashlock users     943565 Apr  3 10:55 SSW_all_biallelic.MAF0.02Miss0.8.recode.vcf
-rw-r--r--. 1 lashlock users     132925 Mar 29 10:58 SSW_all_biallelic.MAF0.02Miss0.8.recode.vcf.geno
-rw-r--r--. 1 lashlock users      76352 Mar  8 10:33 SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz
-rw-r--r--. 1 lashlock users 1365923459 Apr  4 10:10 SSW_by24inds.txt.vcf
-rw-r--r--. 1 lashlock users        312 Mar 29 10:58 SSW.ind
-rw-r--r--. 1 lashlock users     445407 Mar 29 10:58 SSW.snp
-rw-r--r--. 1 lashlock users        168 Mar 29 10:38 SSW_tidal.pops
-rwxr-xr-x. 1 lashlock users       1604 Mar 29 10:39 vcf2admixture_SSW.spid
-rwxr-x--x. 1 lashlock users        274 Mar 29 10:59 vcf2geno.sh
[lashlock@pbio381 ~]$ cd 16s_analysis/
[lashlock@pbio381 16s_analysis]$ ll
total 4
drwxr-xr-x. 26 lashlock users 4096 Apr 10 10:39 joined
drwxr-xr-x.  2 lashlock users   96 Apr 10 10:23 validate_map
[lashlock@pbio381 16s_analysis]$ cd joined/
[lashlock@pbio381 joined]$ ll
total 8
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:39 02_5-05_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:35 02_5-08_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:34 02_5-11_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:37 02_5-14_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:38 02_5-17_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:35 02_5-20_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:36 04_5-05_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:34 04_5-08_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:37 04_5-11_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:38 04_5-14_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:38 04_5-17_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:38 04_5-20_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:37 10_5-05_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:37 10_5-08_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:35 10_5-11_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:34 10_5-14_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:39 10_5-17_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:39 10_5-20_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:36 27_5-05_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:33 27_5-08_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:34 27_5-11_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:36 27_5-14_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:36 27_5-17_R1
drwxr-xr-x. 2 lashlock users  101 Apr 10 10:35 27_5-20_R1
-rw-r--r--. 1 lashlock users 7728 Apr 10 10:39 log_20170410103335.txt
[lashlock@pbio381 joined]$ bash /data/project_data/16s/remove-underscore.sh
[lashlock@pbio381 joined]$ bash /data/project_data/16s/remove-R1.sh
[lashlock@pbio381 joined]$ cd 16s_analysis/
-bash: cd: 16s_analysis/: No such file or directory
[lashlock@pbio381 joined]$ cd ~/
[lashlock@pbio381 ~]$ cd 16s_analysis/
[lashlock@pbio381 16s_analysis]$ multiple_split_libraries_fastq.py -i ~/16s_analysis/joined -o ~/16s_analysis/filtered -m sampleid_by_file --include_input_dir_path --remove_filepath_in_name  --mapping_indicator ~/16s_analysis/map.txt
^CTraceback (most recent call last):
  File "/usr/bin/multiple_split_libraries_fastq.py", line 219, in <module>
    main()
  File "/usr/bin/multiple_split_libraries_fastq.py", line 216, in main
    close_logger_on_success=True)
  File "/usr/lib/python2.7/site-packages/qiime/workflow/util.py", line 114, in call_commands_serially
    stdout, stderr, return_value = qiime_system_call(e[1])
  File "/usr/lib/python2.7/site-packages/qcli/util.py", line 39, in qcli_system_call
    stdout, stderr = proc.communicate()
  File "/usr/lib64/python2.7/subprocess.py", line 800, in communicate
    return self._communicate(input)
  File "/usr/lib64/python2.7/subprocess.py", line 1401, in _communicate
    stdout, stderr = self._communicate_with_poll(input)
  File "/usr/lib64/python2.7/subprocess.py", line 1455, in _communicate_with_poll
    ready = poller.poll()
KeyboardInterrupt
[lashlock@pbio381 16s_analysis]$ multiple_split_libraries_fastq.py -i ~/16s_analysis/joined -o ~/16s_analysis/filtered -m sampleid_by_file --include_input_dir_path --remove_filepath_in_name  --mapping_indicator /data/project_data/16s/map.txt
[lashlock@pbio381 16s_analysis]$ cd filtered/
[lashlock@pbio381 filtered]$ ll
total 933844
-rw-r--r--. 1 lashlock users      5771 Apr 10 11:01 histograms.txt
-rw-r--r--. 1 lashlock users      6471 Apr 10 10:48 log_20170410104834.txt
-rw-r--r--. 1 lashlock users      6533 Apr 10 11:01 log_20170410105128.txt
-rw-r--r--. 1 lashlock users 956183153 Apr 10 11:01 seqs.fna
-rw-r--r--. 1 lashlock users     36908 Apr 10 11:01 split_library_log.txt
[lashlock@pbio381 filtered]$ head seqs.fna
>10-5-17_0 M02780:138:000000000-ATTB1:1:1101:12235:1197 1:N:0:ACTGAGCG+TATCCTCT orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGGGGCAGCAGCTAAGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCGACACCGCGTGAGTGATGAAGGTCGAGAGGTTGTAAGACTCTTTTATTACCGAAGAATAAGTATTATGCGAAAGTATAATACGATGACATTAAGTGATGAATAAGCCCTGGCTAATTACGTGCCAGCAGCCGCGGTAATACGTAAGGGGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGGGCGCGTAGGCGGATTAGAAGGCATTATGATAAAGGTCATAGCTCAACTATGGTGTATTGTAATGAACCTCTAGTCTTGAGTTTTGGAGAGGGAACTGGAATTCTTGGTGTAAGGGTGAAATCTGTAGATATCAAGAAGAACATCAGTGGCGAAGGCGAGTTCCCAGCCAAAAACTGACGCTGAGGTGCGAAAGCGTGGGTATCGAACAGGATTAGATACCCCAGTAGTC
>10-5-17_1 M02780:138:000000000-ATTB1:1:1101:12229:1216 1:N:0:ACTGAGCG+TATCCTCT orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGGGGCAGCAGCTAAGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCGACACCGCGTGAGTGATGAAGGTCGAGAGGTTGTAAGACTCTTTTATTACCGAGGAATAAGTATTATGCGAAAGTATAATACGATGACATTAAGTGATGAATAAGCCCTGGCTAATTACGTGCCAGCAGCCGCGGTAATACGTAAGGGGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGGGCGCGTAGGCGGATTAGAAGGCATTATGATAAAGGTCATAGCTCAACTATGGTGTATTGTAATGAACCTCTAGTCTTGAGTTTTGGAGAGGGAACTGGAATTCTTGGTGTAAGGGTGAAATCTGTAGATATCAAGAAGAACATCAGTGGCGAAGGCGAGTTCCCAGCCAAAAACTGACGCTGAGGTGCGAAAGCGTGGGTATCGAACAGGATTAGATACCCCAGTAGTC
>10-5-17_2 M02780:138:000000000-ATTB1:1:1101:9797:1235 1:N:0:ACTGAGCG+TATCCTCT orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGCGGCAGCAGCTAAGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCGACACCGCGTGAGTGATGAAGGTCGAGAGGTTATAAAACTCTTTTATTACCGAAGAATAAGTATTATGCGAAAGTATAATACGATGACATTAAGTGATGAATAAGCCCTGGCTAATTACGTGCCAGCAGCCGCGGTAATACGTAAGGGGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGGGCGCGTAGGCGGATTAGAAGGCATTATGATAAAGGTCATAGCTCAACTATGGTGTATTGTAATGAACCTCTAGTCTTGAGTTTTGGAGAGGGAACTGGAATTCTTGGTGTAAGGGTGAAATCTGTAGATATCAAGAAGAACATCGGTGGCGAAGGCGAGTTCCCAGCCAAAAACTGACGCTGAGGTGCGAAAGCGTGGGTATCGAACAGGATTAGATACCCGTGTAGTC
>10-5-17_3 M02780:138:000000000-ATTB1:1:1101:10894:1237 1:N:0:ACTGAGCG+TATCCTCT orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGCGGCAGCAGTGGGGAATCTTGCACAATGGGCGAAAGCCTGATGCAGCCATGCCGCGTGAATGATGAAGGCCTTAGGGTTGTAAAATTCTTTCGCTAGGGATGATAATGACAGTACCTAGTAAAGAAGCACCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGTGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGACTATTAAGTCAGATGTGAAATCCCAAGGCTCAACCTTGGAACTGCATTTGAAACTGGTAGTCTAGAGATCAGGAGAGGTTGGCGGAATACCGAGTGTAGAGGTGAAATTCGTAGATATTCGGTGGAACACCAGTGGCGAAGGCGGCCAACTGGACTGATACTGACGCTGAGGCGCGAAAGTGTGGGGAGCAAACAGGATTAGATACCCGTGTAGTC
>10-5-17_4 M02780:138:000000000-ATTB1:1:1101:14644:1265 1:N:0:ACTGAGCG+TATCCTCT orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGAGGCAGCAGCTAAGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCGACACCGCGTGAGTGATGAAGGTCGAGAGGTTGTAAGACTCTTTTATTACCGAGGAATAAGTATTATGCGAAAGTATAATACGATGACATTAAGTGATGAATAAGCCCTGGCTAATCACGTGCCTGCAGCCGCGGTAATACGTAAGGGGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGGGCGCGTAGGCGGATTAGAAGGCATTATGATAAAGGTCATAGCTCAACTATGGTGTATTGTAATGAACCTCTAGTCTTGAGTTTTGGAGAGGGAACTGGAATTCTTGGTGTAAGGGTGAAATCTGTAGATATCAAGAAGAACATCAGTGGCGAAGGCGAGTTCCCAGCCAAAAACTGACGCTGAGGTGCGAAAGCGTGGGTATCGAACAGGATTAGATACCCGAGTAGTC
[lashlock@pbio381 filtered]$ extract_seqs_by_sample_id.py -i seqs.fna -o test -s 04-5-05
[lashlock@pbio381 filtered]$ head test
>04-5-05_109415 M02780:138:000000000-ATTB1:1:1101:17994:1235 1:N:0:TCCTGAGC+GTAAGGAG orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACAGGCGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGTCAGGAGGAAAGGTTAGTAGTTAATACCTGCTAGCTATGACGTTACTGACAGAAGAAGCACCGGCTAGCTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTACGCGGGCGGTTTGTTAAGCGAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTTCGAACTGGCAAACTAGAGTGTGATAGAGGGTGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCTGAAGGAATTCCGATGGCGAAGGCAGCCACCCGGGTCAACACTGACGCTCATGTACGAAAGCGTGGGGAGCAAACGGGATTAGATACCCTGGTAGTC
>04-5-05_109416 M02780:138:000000000-ATTB1:1:1101:11795:1279 1:N:0:TCCTGAGC+GTAAGGAG orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGCGAGGAGGAAAGGGATGTAGTTAATAACTGCATTCTGTGACGTTACTCGCAGAAGAAGCACCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAGGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGTAGGTGGTTTGTTAAGCAAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTTTGAACTGGCAGGCTAGAGTATTGTAGAGGGTAGTGGAATTTCCAGTGTAGCGGTGAAATGCGTAGAGATTGGAAGGAACATCAGTGGCGAAGGCGGCTACCTGGACAAATACTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCGAGTAGTC
>04-5-05_109417 M02780:138:000000000-ATTB1:1:1101:13929:1361 1:N:0:TCCTGAGC+GTAAGGAG orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGGGGCTGCAGCTAAGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCGACACCGCGTGAGTGATGAAGGTCGAGAGGTTGTAAAACTCTTTTATTACCGAAGAATAAGTATTATGCGAAAGTATAATACGATGACATTAAGTGATGAATAAGCCCTGGCTAATTACGTGCCAGCAGCCGCGGTAATACGTAAGGGGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGGGCGCGTAGGCGGATTAGAAGGCATTATGGTAAAGGTCATAGCTCAACTATGGTGTATTGTAATGAACCTCTAGTCTTGAGTTTTGGAGAGGGAACTGGAATTCTTGGTGTAAGGGTGAAATCTGTAGATATCAAGAAGAACATCAGTGGCGAAGGCGAGTTCCCAGCCAAAAACTGACACTGAGGTGCGAAAGCGTGGGTATCGAACAGGATTAGATACCCCTGTAGTC
>04-5-05_109418 M02780:138:000000000-ATTB1:1:1101:17195:1363 1:N:0:TCCTGAGC+GTAAGGAG orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGTGGCTGCAGTGGGGAATATTGCACAATAGGCGCAAGCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGTCAGGAGGAAAGGTTAGTAGTTAATACATGCTAGCTGTGACGTTACTGACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTACGCAGGCGGTTTGTTAAGCGAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTTCGAACTGGCAAACTAGAGTGTGATAGAGGGTGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCTGAAGGAATACCGATGGCGAAGGCAGCCACCTGGGTCAACACCGACGCTCATGTACGAAAGCGTGGGGAGCAAACGGGATTAGATACCCTTGTAGTC
>04-5-05_109419 M02780:138:000000000-ATTB1:1:1101:8251:1394 1:N:0:TCCTGAGC+GTAAGGAG orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
CCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGTCAGGAGGAAAGGTTAGTAGTTAATACCTGCTAGCTGTGACGTTACTGACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTACGCAGGCGGTTTGTTAAGCGAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTTCGAACTGGCGAACTAGAGTGTGGTAGAGGGTGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCCGAAGGAATACCGATGGCGAAGGCAGCCACCTGGGTCAACACTGACGCTCATGTACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTTGTAGTC
[lashlock@pbio381 filtered]$ rm test
rm: remove regular file ‘test’? y
[lashlock@pbio381 filtered]$
[lashlock@pbio381 filtered]$ pick_open_reference_otus.py -i ~/16s_analysis/filtered/seqs.fna -o ~/16s_analysis/otus  --parallel --jobs_to_start 1
^CTraceback (most recent call last):
  File "/usr/bin/pick_open_reference_otus.py", line 453, in <module>
    main()
  File "/usr/bin/pick_open_reference_otus.py", line 432, in main
    minimum_failure_threshold=minimum_failure_threshold)
  File "/usr/lib/python2.7/site-packages/qiime/workflow/pick_open_reference_otus.py", line 713, in pick_subsampled_open_reference_otus
    close_logger_on_success=False)
  File "/usr/lib/python2.7/site-packages/qiime/workflow/util.py", line 114, in call_commands_serially
    stdout, stderr, return_value = qiime_system_call(e[1])
  File "/usr/lib/python2.7/site-packages/qcli/util.py", line 39, in qcli_system_call
    stdout, stderr = proc.communicate()
  File "/usr/lib64/python2.7/subprocess.py", line 800, in communicate
    return self._communicate(input)
  File "/usr/lib64/python2.7/subprocess.py", line 1401, in _communicate
    stdout, stderr = self._communicate_with_poll(input)
  File "/usr/lib64/python2.7/subprocess.py", line 1455, in _communicate_with_poll
    ready = poller.poll()
KeyboardInterrupt
[lashlock@pbio381 filtered]$ cd ~/
[lashlock@pbio381 ~]$ cd 16s_analysis/
[lashlock@pbio381 16s_analysis]$ ll
total 8
drwxr-xr-x.  2 lashlock users 4096 Apr 10 11:08 filtered
drwxr-xr-x. 26 lashlock users 4096 Apr 10 10:44 joined
drwxr-xr-x.  3 lashlock users   64 Apr 10 11:13 otus
drwxr-xr-x.  2 lashlock users   96 Apr 10 10:23 validate_map
[lashlock@pbio381 16s_analysis]$ rm otus
rm: cannot remove ‘otus’: Is a directory
[lashlock@pbio381 16s_analysis]$ rm -r otus
rm: descend into directory ‘otus’? y
rm: remove regular file ‘otus/log_20170410111307.txt’? y
y
rm: descend into directory ‘otus/step1_otus’? rm: remove directory ‘otus/step1_otus/POTU_5Igq_’? y
rm: remove directory ‘otus/step1_otus’? y
rm: remove directory ‘otus’? y
[lashlock@pbio381 16s_analysis]$
[lashlock@pbio381 16s_analysis]$ biom summarize-table -i /data/project_data/16s/otu_table/otu_table_mc2_w_tax_no_pynast_failures.biom
Num samples: 176
Num observations: 93033
Total count: 8362869
Table density (fraction of non-zero values): 0.028

Counts/sample summary:
 Min: 28412.0
 Max: 77866.0
 Median: 47051.500
 Mean: 47516.301
 Std. dev.: 7637.541
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
24-5-11: 28412.0
07-5-11: 31532.0
04-5-08: 32477.0
07-5-05: 32491.0
38-6-09: 33391.0
17-5-05: 34322.0
27-5-14: 34735.0
14-5-08: 35793.0
26-5-05: 35813.0
18-5-14: 35851.0
06-5-11: 36919.0
31-6-21: 37205.0
38-6-12: 37260.0
17-5-08: 37558.0
03-5-11: 37765.0
11-5-14: 38114.0
13-5-14: 38211.0
04-5-11: 38548.0
37-6-09: 38624.0
27-5-05: 38849.0
17-5-17: 38850.0
11-5-05: 38888.0
05-5-11: 39128.0
08-5-14: 39534.0
12-5-08: 40033.0
05-5-08: 40184.0
21-5-14: 40241.0
10-5-17: 40487.0
07-5-08: 40554.0
25-5-08: 40587.0
16-5-14: 40620.0
04-5-20: 40769.0
06-5-14: 40967.0
37-6-12: 40986.0
18-5-08: 41298.0
21-5-11: 41339.0
13-5-05: 41346.0
04-5-05: 41636.0
02-5-08: 41735.0
10-5-20: 41896.0
13-5-08: 42039.0
19-5-14: 42053.0
11-5-08: 42101.0
02-5-14: 42335.0
23-5-14: 42337.0
27-5-17: 42426.0
22-5-14: 42493.0
28-5-05: 42667.0
21-5-08: 42790.0
15-5-14: 42791.0
03-5-08: 42813.0
15-5-05: 42847.0
29-5-08: 42962.0
02-5-20: 43276.0
17-5-14: 43278.0
18-5-05: 43362.0
33-6-09: 43419.0
06-5-08: 43467.0
26-5-08: 43747.0
01-5-08: 43987.0
09-5-17: 44133.0
16-5-05: 44183.0
29-5-05: 44289.0
31-6-18: 44424.0
15-5-08: 44509.0
09-5-08: 44676.0
12-5-05: 44776.0
06-5-05: 44982.0
04-5-14: 45114.0
31-6-12: 45116.0
28-5-08: 45136.0
18-5-11: 45361.0
32-6-12: 45623.0
35-6-15: 45661.0
19-5-17: 45695.0
20-5-08: 45976.0
11-5-11: 46033.0
10-5-08: 46167.0
10-5-05: 46198.0
08-5-08: 46202.0
37-6-18: 46255.0
36-6-12: 46275.0
20-5-05: 46401.0
14-5-11: 46664.0
35-6-18: 46912.0
09-5-14: 46977.0
19-5-20: 46991.0
16-5-08: 47037.0
02-5-17: 47066.0
04-5-17: 47144.0
34-6-21: 47231.0
09-5-11: 47282.0
17-5-11: 47458.0
31-6-09: 47607.0
31-6-15: 47689.0
14-5-05: 47927.0
13-5-11: 48302.0
24-5-05: 48316.0
32-6-15: 48497.0
24-5-08: 48610.0
24-5-17: 48621.0
02-5-11: 48672.0
29-5-11: 48707.0
28-5-17: 48835.0
23-5-17: 48846.0
29-5-14: 48850.0
01-5-11: 48877.0
10-5-14: 48961.0
15-5-11: 48992.0
22-5-08: 49011.0
21-5-05: 49041.0
16-5-11: 49329.0
23-5-20: 49451.0
19-5-08: 49585.0
35-6-21: 49623.0
08-5-17: 49963.0
32-6-09: 50068.0
28-5-14: 50101.0
09-5-20: 50105.0
05-5-05: 50125.0
21-5-17: 50382.0
33-6-18: 50668.0
08-5-20: 50698.0
15-5-17: 51066.0
20-5-14: 51176.0
33-6-12: 51371.0
32-6-21: 51411.0
36-6-09: 51448.0
23-5-05: 51558.0
10-5-11: 51571.0
34-6-18: 51990.0
19-5-05: 52149.0
33-6-15: 52401.0
08-5-11: 52658.0
35-6-12: 52791.0
35-6-09: 52871.0
25-5-05: 53423.0
20-5-20: 53450.0
23-5-08: 53456.0
28-5-11: 53886.0
26-5-11: 53925.0
38-6-24: 53960.0
20-5-11: 53978.0
38-6-21: 54118.0
32-6-18: 54151.0
34-6-15: 54153.0
37-6-15: 54388.0
34-6-12: 54717.0
02-5-05: 54850.0
24-5-14: 54994.0
36-6-15: 55095.0
01-5-05: 55153.0
38-6-15: 55242.0
27-5-11: 55614.0
08-5-05: 56383.0
15-5-20: 56463.0
20-5-17: 56497.0
09-5-05: 56767.0
35-6-24: 56984.0
03-5-05: 57367.0
34-6-09: 57813.0
22-5-11: 57906.0
27-5-20: 58225.0
21-5-20: 58440.0
24-5-20: 58825.0
22-5-05: 58940.0
27-5-08: 59017.0
33-6-21: 60015.0
34-6-24: 60708.0
36-6-18: 61008.0
31-6-24: 61401.0
23-5-11: 63179.0
38-6-18: 65199.0
33-6-24: 68217.0
32-6-24: 72091.0
37-6-21: 77866.0
[lashlock@pbio381 16s_analysis]$
```

Now we are going to look at the diversity of our samples this takes a long time so we will set a screen



------
<div id='id-section20'/>
### Page 20:

More microbiome stuff

- Filter chimeric sequences
  - first step is to make a fasta file with the chimeric sequences in our data

```
vsearch --uchime_ref /data/project_data/16s/otu_table/rep_set.fna --chimeras ~/16s_analysis/mc2_w_tax_no_pynast_failures_chimeras.fasta --db /usr/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta

```

- Remove chimeras from OTU table
  - we will use the python script below to do this

```
filter_otus_from_otu_table.py -i /data/project_data/16s/otu_table/otu_table_mc2_w_tax_no_pynast_failures.biom -o otu_table_mc2_w_tax_no_pynast_failures_no_chimeras.biom -e ~/16s_analysis/mc2_w_tax_no_pynast_failures_chimeras.fasta

```

- Actually removing files now

```
filter_fasta.py -f /data/project_data/16s/otu_table/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta -o ~/16s_analysis/rep_set_aligned_pfiltered_no_chimeras.fasta -a ~/16s_analysis/mc2_w_tax_no_pynast_failures_chimeras.fasta -n

```



- Set screen and make remake phylogenetic tree without chimeric sequences

```
screen 
make_phylogeny.py -i ~/16s_analysis/rep_set_aligned_pfiltered_no_chimeras.fasta -o ~/16s_analysis/rep_set_no_chimeras.tre
```

Now we will frequency filter our data to get rid of OTUs with a min count of 50 and a min sample of 44

```
filter_otus_from_otu_table.py -i otu_table_mc2_w_tax_no_pynast_failures_no_chimeras.biom -o otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom --min_count 50 --min_samples 44

```

We now have only 1064 OTUs left



Diversity analysis in Qiime (Melanie ran this for us)

```
core_diversity_analyses.py -o core_diversity_filtered -i otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom -m ~/Po/MiSeq/joined/map.txt -t rep_set_no_chimeras.tre -e 20000 -a 8

```

Notes from phyloseq analysis in R script



------
<div id='id-section21'/>
### Page 21:

1. 1. Review of pipeline

      1. Processing data
         1. using fastq files and program Trimmomatic
      2. Transcriptome assembly
         1. fasta file using Trinity
      3. Annotation
         1. fasta file using BLAST
      4. Map reads
         1. sam files using BWA
      5. Differential Gene Expression Analysis
         1. Using counts tables and DESEQ2
      6. Population Genomics
         1. using vcf file analyzing with PCA, DAPC, and Admixture

      We've done the DGE and popgen analyses, and now we want to know biological function...

      - Functional Enrichment Analysis

      Annotation Methods

      - Blast2GO
        - paid
      - Brute force
        - blast and annotate yourself
      - Pipelines
        - Trinotate

      Generic pipeline

      1. .fasta file with genes
      2. use existing databases to find out what these genes are
         1. Using BLAST
            1. Using thresholds (e value cutoff)
         2. Diamond
      3. Blast output
         1. e value (smaller better)
         2. bit score (bigger better)
         3. % identity
         4. length

   ​

Databases

- NCBI NR
  - non redundant protein database
- UniProt
  - proteins and their associated terms


------
<div id='id-section22'/>
### Page 22:

Notes from Phyloseq analysis in R script



------
<div id='id-section23'/>
### Page 23:

Info update: Rarefying microbiome data is not acceptable

Rarefying vs Rarefaction

Past methods

- Proportions of OTUs
  - high rate of false positives
  - does not account for heteroscedasticity
- Rarefying
  - Steps:
    - Find your smallest N reads
    - Discard libraries with fewer than N reads
    - Subsample OTUs with larger than N without replacement
  - What this does:
    - this normalizes the data
    - Removes large discrepancies in the number of reads
  - Problems with rarefying:
    - high rate of false positives
    - requires omission of data
    - reduces statistical power
- Mixture model
  - employing a statistical distribution that includes multiple types of distributions in analysis
  - Not eliminating as much data
  - increase in accuracy
  - distributions: binomial with0-inflated gaussian 
  - accounts for biological variability
  - can use edgeR and DESEQ and phyloseq to emply mixture models 

Working with picrust

1. Make a filtered OTU table with only closed OTUs found in your sample

```
filter_otus_from_otu_table.py -i ~/16s_analysis/otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom -o ~/16s_analysis/closed_otu_table.biom --negate_ids_to_exclude -e /usr/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta 
```

1. How many OTUs do we have?

```
Num samples: 176
Num observations: 259
Total count: 1493357
Table density (fraction of non-zero values): 0.546

Counts/sample summary:
 Min: 152.0
 Max: 32426.0
 Median: 6333.500
 Mean: 8484.983
 Std. dev.: 6937.648
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

```

1. Normalize by copy number of the 16s gene

```
normalize_by_copy_number.py -i ~/16s_analysis/closed_otu_table.biom -o ~/16s_analysis/closed_otu_table_norm.biom
```

1. Predicting metagenomes
   1. output text file

```
predict_metagenomes.py -f -i ~/16s_analysis/closed_otu_table_norm.biom -o ~/16s_analysis/metagenome_predictions.txt -a nsti_per_sample.txt 
```

```
2. output biom file
```

```
predict_metagenomes.py -i ~/16s_analysis/closed_otu_table_norm.biom -o ~/16s_analysis/metagenome_predictions.biom -a nsti_per_sample.txt
```

1. Count the rows in your .txt file to see how many KEGG Orthology terms that were predicted

```
wc -l metagenome_predictions.txt
```

6910 metagenome_predictions.txt

1. Collapse your results to a higher KO heirarchy term

```
categorize_by_function.py -f -i metagenome_predictions.biom -c KEGG_Pathways -l 3 -o metagenome_predictions.L3.txt
```

1. Do this for the biom file as well

```
categorize_by_function.py -i metagenome_predictions.biom -c KEGG_Pathways -l 3 -o metagenome_predictions.L3.biom
```

1. You can run the summarize command again to get an idea of how many terms you have left and then  move the biom file to your computer and run analyses in R

```
library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
library("ggplot2")
theme_set(theme_bw())
library('biom')

x = read_biom("metagenome_predictions.L3.biom")
otumat = as(biom_data(x), "matrix")
OTU = otu_table(otumat, taxa_are_rows=TRUE)


mapping <- import_qiime_sample_data(mapfilename = 'R_map.txt')

phylo <- merge_phyloseq(OTU, mapping)
phylo

###############################################################################
###Test for DE KO terms between individuals that got sick and those that didn't
###############################################################################

final_pheno = phyloseq_to_deseq2(phylo, ~ Final_phenotype)
final_pheno_results = DESeq(final_pheno, test="Wald")
final_pheno_res = results(final_pheno_results)
summary(final_pheno_res)
head(final_pheno_res)

alpha = 0.05
final_pheno_sigtab = final_pheno_res[which(final_pheno_res$padj < alpha), ]
final_pheno_sigtab= cbind(as(final_pheno_sigtab, "data.frame"), as(tax_table(phylo)[rownames(final_pheno_sigtab), ], "matrix"))
head(final_pheno_sigtab)
final_pheno_sigtab
write.table(final_pheno_sigtab, "Final_pheno_L3.txt", sep="\t")
```





------
<div id='id-section24'/>
### Page 24:

Final Project R Script



```
#FinalProj Script

#load in libraries

library("phyloseq")
library("DESeq2")
library(XML)
library("ggplot2")
theme_set(theme_bw())

#Import the OTU table
otutable <- import_biom(BIOMfilename = 'otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom', 
                        treefilename = 'rep_set_no_chimeras.tre', 
                        parseFunction = parse_taxonomy_greengenes)

#The warnings are ok. There is 1 warning for every OTU that doesn't have a taxonomy assignment

#Import the mapping file
mapping <- import_qiime_sample_data(mapfilename = 'R_map.txt')

#Merge the mapping file to the OTU table

phylo <- merge_phyloseq(otutable, mapping)


#Check to make sure the imports worked
phylo

head(otu_table(phylo))

tail(otu_table(phylo))

phylo_subset = subset_samples(phylo, Day == 3 & Phenotype !="Dead")
phylo_int = subset_samples(phylo, Day==3 & Tide=="intertidal")
phylo_int
phylo_sub = subset_samples(phylo, Day==3 & Tide=="subtidal")
phylo_sub
phylo_subset
sample_data(phylo_subset)
sample_data(phylo_int)
sample_data(phylo_sub)

##We want to numbers of each individual to be a factor. Right now it's an integer so we have to change that.
class(sample_data(phylo_subset)$individual)
sample_data(phylo_subset)$individual<-factor(sample_data(phylo_subset)$individual)

##Phyloseq's wrapper to get OTU data into DESeq

loc <- phyloseq_to_deseq2(phylo_subset, ~  Tide)
names(loc)

##Run DESeq. This command takes a bit of time

loc_deseq_test <- DESeq(loc, test="Wald")

str(loc_deseq_test)
##Get results from DESeq

loc_results <- results(loc_deseq_test)

##Have a look at the summary
summary(loc_results)
head(loc_results)


#Subtidal is the reference
# out of 1063 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 73, 6.9% 
# LFC < 0 (down)   : 52, 4.9% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 248, 23% 


##Make table with OTUs padj<0.05
alpha <- 0.05

loc_sigtab <- loc_results[which(loc_results$padj < alpha), ]

loc_sigtab

##Add taxa info to that table

loc_sigtab <- cbind(as(loc_sigtab, "data.frame"), as(tax_table(phylo)[rownames(loc_sigtab), ], "matrix"))
tail(loc_sigtab)

#sort by padj
loc_sigtab <- loc_sigtab[order(resLoc$padj),1:12]
head(loc_sigtab)

loc_sigtab <- loc_sigtab[1:6, 1:12]
loc_sigtab

##Save the table to your desktop

write.table(loc_sigtab, "DE_OTU_intertidal_vs_subtidal.txt", sep="\t")

###Diversity analyses

richness.df <- estimate_richness(phylo_subset, split=TRUE, measures=NULL)

str(richness.df)
richness.df
#make two separate data frames to run richness analyses separately, add on a column denoting tide and then merge back together

int_rich.df <- estimate_richness(phylo_int, split=TRUE, measures=NULL)
int_rich.df$Tide <-as.factor("Intertidal")
str(int_rich.df)
int_rich.df

sub_rich.df <- estimate_richness(phylo_sub, split=TRUE, measures=NULL)
sub_rich.df$Tide <-as.factor("Subtidal")
str(sub_rich.df)
sub_rich.df

richByTide.df <- rbind(int_rich.df, sub_rich.df)

richByTide.df

#removing standard error columns

richByTide.df$se.chao1 <- NULL

richByTide.df$se.ACE <- NULL

richByTide.df

# t1<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$p.value})
# t2<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$estimate})
# t3<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$conf.int[1:2]})
# #t1<-t.test(Fisher~Tide,data=richByTide.df)$p.value
# #t2<-t.test(Fisher~Tide,data=richByTide.df)$estimate
# #t3<-t.test(Fisher~Tide,data=richByTide.df)$conf.int[1:2]
# 
# data.frame(cbind(Pvalue=t1,Estimate=t2,ConfidenceInterval=t3))
# data.frame(cbind(t1,t2,t3))


t1<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$p.value})
t2<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$estimate})
t3<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$conf.int[1:2]})

str(t1)
t2
t3
t1 <- as.data.frame(t1)
t1

t1 <- t(as.data.frame(t1))
t2 <- t(as.data.frame(t2))
t3 <- t(as.data.frame(t3))

t1
t2
t3

SummTable <- rbind(t1,t2,t3)

SummTable

row.names(SummTable) <- c("pvalue", "Mean in Intertidal", "Mean in Subtidal", "LCI", "UCI")

SummTable

write.table(SummTable, "RichnessSummary.txt", sep="\t")

########################################################

##Plot abundances of taxa in different pheno numbers
##First we must rarefy
set.seed(28132)
phyloR = rarefy_even_depth(phylo_subset, sample.size = min(sample_sums(phylo_subset)))

#Then we check to make sure we rarefied corectly.
title = "Sum of reads for each sample, phyloR"
plot(sort(sample_sums(phyloR), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 20000))


##We merge on basis of pheno number so we can make relative abundance graph
#phyloRm = merge_samples(phyloR, "Pheno_num")

##Merge on basis of phenotype
phyloRmD = merge_samples(phyloR, "Tide")

phyloRmD

phyloRt = transform_sample_counts(phyloRmD, function(x) 100 * x/sum(x))


p = plot_bar(phyloRt,x="Tide",facet_grid=NULL, fill="Order")
p + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

###################DGE Analysis###################

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)


conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
colData

colDay3 <- subset(colData, day=="day03")

colDay3
# colInt<-subset(colDay3,location=="int")
# colSub <- subset(colDay3, location=="sub")
# 
# colInt
# colSub
# 
# countDataInt<-countData[, which(colnames(countData) %in% row.names(colInt))]
# countDataSub<-countData[, which(colnames(countData) %in% row.names(colSub))]
# dim(countDataInt)
# dim(countDataSub)

countDataDay3<-countData[, which(colnames(countData) %in% row.names(colDay3))]
head(countDataDay3)


ddsLoc<- DESeqDataSetFromMatrix(countData = countDataDay3, colData = colDay3 ,design = ~ location)

ddsLoc

ddsLoc <- ddsLoc[ rowSums(counts(ddsLoc)) > 100, ]

colData(ddsLoc)$location <- factor(colData(ddsLoc)$location, levels=c("int","sub"))

colData(ddsLoc)

ddsLoc <- DESeq(ddsLoc) 

names(ddsLoc)
str(ddsLoc)

resLoc <- results(ddsLoc)
resLoc



#sort results by padj
resLoc <- resLoc[order(resLoc$padj),]

head(resLoc)

summary(resLoc)

## Merge with normalized count data
resSumm <- merge(as.data.frame(resLoc), as.data.frame(counts(ddsLoc, normalized=TRUE)), by="row.names", sort=FALSE)
names(resSumm)[1] <- "Gene"
resSumm <- head(resSumm)
resSumm

row.names(resSumm) <- c("DN41041_c3_g2", "DN47102_c1_g1", "DN19042_c0_g1", "DN39048_c5_g1", "DN39233_c0_g1", "DN39328_c5_g1")

resSumm <- resSumm[,2:7]

resSumm

## Write results
write.csv(resSumm, file="diffexpr-results.csv")

plotCounts(ddsLoc, gene="TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543", intgroup="location")

vsdata <- vst(ddsLoc, blind=FALSE)

plotPCA(vsdata, intgroup="location")

#plot counts for top 6 genes

normcounts<-data.frame(genes=rownames(counts(ddsLoc, normalized=TRUE)),counts(ddsLoc, normalized=TRUE))
head(normcounts)

#names.to.keep <- as.factor(c("TRINITY_DN41041_c3_g2_TRINITY_DN41041_c3_g2_i1_g.10613_m.10613", "TRINITY_DN47102_c1_g1_TRINITY_DN47102_c1_g1_i3_g.24581_m.24581", "TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710", "TRINITY_DN39048_c5_g1_TRINITY_DN39048_c5_g1_i1_g.8316_m.8316", "TRINITY_DN39233_c0_g1_TRINITY_DN39233_c0_g1_i1_g.8447_m.8447 ", "TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543"))
#length(names.to.keep)
#rows.to.keep<-which(rownames(normcounts) %in% names.to.keep) 
normcounts2 <- normcounts[normcounts$genes %in% as.factor(rownames(head(resLoc))),]

normcounts2$genes <- c("DN41041_c3_g2", "DN47102_c1_g1", "DN19042_c0_g1", "DN39048_c5_g1", "DN39233_c0_g1", "DN39328_c5_g1")
dim(normcounts2)
str(normcounts2)
normcounts2

install.packages("tidyr")

library(tidyr)

countsDF <- gather(normcounts2,individuals,normcounts,I03_5.08_S_2:I38_6.12_H_0)
dim(countsDF)
head(countsDF)

summary(countsDF)

head(conds)
conds$individuals <- as.factor(rownames(conds))
str(conds)
library(dplyr)

topHits <- inner_join(countsDF,conds,by="individuals")
topHits
dim(topHits)
str(topHits)

ggplot(topHits, aes(x=location,y=log(normcounts+1), color=location))+
  geom_boxplot()+facet_grid(.~genes)+theme(legend.position="none")

#####################Annotation of topHits#############################
#Blastx results
#>TRINITY_DN41041_c3_g2_TRINITY_DN41041_c3_g2_i1_g.10613_m.10613 type:5prime_partial len:413 TRINITY_DN41041_c3_g2_i1:2720-1482(-)
#CCTAATGAGGTGCTTGTGCGCATGGACCCAAGTGTTAAGCAGCCCAAGTGTGTCAAGAAAGAGTTTCCTCCATGCGTTTACAGAGATCGGGCCTGTCCCAATTATCGAGTTGTGGACACATTGACGGATAACATTGAGCATCGTCAGTTAGCTGCTAGCCATTACGTCATGAAGAGGACTCGAACCTGTAACATGAGCCAGGCATACGATGATGCATATATGCCACTTTACAATTACCTCAAAGGAGCAAACAGCGAGTCTGTGACTATGAGCCCTACAGCTCCGGACATGATCCTGGTCCACAAATCAACAGACAACCCACCCGAGGGGTGTGACTTTGCTTACAGTCTGTTGTTTTACCTTCCCGTTGACAGACAGCTCAGTGCCCCAGCACCGACTGAGGATGGTGTTGAACTCATCAAGCAAGAAGCCAAGGAAGCGTTTGTCATTACCTTCACAGATTCAACAACGGACCGAGTACTTGCAGCCAAGGAGACCGAGCTACACAATGAACTGGACAAACGTGAACTCTGTTATAGCACCGATGAATTCTTCTTTGCAAGCTATGATGCATCGTGGAGACCAGAGCCCCATCGCAAAGAGATTTGGATCCCAGAAGCAATGTGCACTACACCAAGTGTCCAAGGTCCTGGTTACACCATCCTGGACAATGGATGTGGGGAAGAGGTCAACTGTCCATCTTATCAGAAGGGTCAATCCCACGGGACATTCGAGGAGAGGAACTACGAGGAGAGCATGTGGATTTGTAAGACGACCGTGAGCTGTAATGTTGAGCAGGCCTTCAGTCAGTCCATCCTGCCGCTGTACCAGTATTTCACCAGCAGCTATGACATCCAACCTGTAGCCAGACCTATCATTTCCTATATGAGAATTGCCGACTTACTGAGCACGGAGTGCAACAAGGAGATTAAGACCTGTGCCTACCTACCAGAAACTGCTTTCAGACTGCCCCTGTTGACCCCACCTAGTCAATTGACATTGCTGAGCGTTCCCAGTGATGACCATCCATGGCAGGCAGCTTATGTAACCGCCCTCCAAGGTCCACCGACAGTCACCAACATCCAGTCAAGTATTACAACACTGCTGCAACAGATCATGGATGCAGACTTGTTCCACAGTGGCCGTGTGTTTGTGGCTCGTTACACCCTTCCGATGGAGGATGGACAGCAATACATTGAAGTCGGAGCTCTGACTGACACTATATCTTATGACGTCTAA
#Can promote mitochondrial permeability transition and facilitate necrotic cell death under different types of stress conditions
#Heme binding protein 2
#human protein

#"Induction of necrotic cell death and mitochondrial permeabilization by heme binding protein 2/SOUL."
#Szigeti A., Bellyei S., Gasz B., Boronkai A., Hocsak E., Minik O., Bognar Z., Varbiro G., Sumegi B., Gallyas F. Jr.
#FEBS Lett. 580:6447-6454(2006) [PubMed] [Europe PMC] [Abstract]

#>TRINITY_DN47102_c1_g1_TRINITY_DN47102_c1_g1_i3_g.24581_m.24581 type:5prime_partial len:127 TRINITY_DN47102_c1_g1_i3:420-40(-)
#GATCGATGTGCACCATCTTCGTACCGGAAATGGCCTGCAGGCCGGGACACTCAGAATGCCACTCAAAAGCATCTTCAGCGGGATTTTTGTACACAGAGCTGCAAAAAAGAAAAAAGAAGAAGAAACGAAATAAAATACAAGGAGGGTCAAATCACACCCATAACAACATCGGATTTGCACAAGTTTTCCACCACTGGAAAGCCACGATTCAGGGCTCTCACCTCGTTGTGGATCCGTGTCTCTGACCATTTTCTAGGCGGGGCGAGCGGCCCCGGAATAGGGCTACAACCACCGGCAAATTTGCCCACGCGCCGGGGCCGAGCAACAACAAAGGCAAGACATGTCCATCCGAAGAATGATTTGAAACAACACCCGGGCTAA
#No significant similarity found

#>TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710 type:internal len:168 TRINITY_DN19042_c0_g1_i1:2-502(+)
#TACACTATGGTTGATGCGGTAGCTAAACACTCCTTCAAAGCAACTAGCGAAGATGAGCTGTCCTTCGAGAAAGGAGATACAATCGTTGTTAACGTTTTAGAGAACTCGGAGCAGCACTGGTACAAAGCCTCTTTGAACGGCAAAAAAGGATTAGTCCCAGCAAACTACATTGAAATAAAACCTTACCCTTTCTTACATGGAAACATCTCAAGAGACGGCGCTGTTTGGAAGCTTAACGACCAGCCAGACGGAGCTTTCCTTATCCGTAGAAGTGAGAGTGACATCAGCGAAATCAGTTACTCACTGTCTGTCAAGTACGGTGATGGAGTACAGCACTTCAAGATCTTAACAGATGGCAGCTACAAATACTTCCTGTGGGTAGTGAAGTTTAGATCGCTAAATGAACTGGTGGAATACTACAAGACGTCCTCTGTTAGCCGGACACAAACAATTTTATTGAAGGATATGGAAGTGGAGCACGAAGTAGCGGTGGCAGAATTT
#Adapter protein that provides a critical link between cell surface growth factor receptors and the Ras signaling pathway.
#Isoform 2 does not bind to phosphorylated epidermal growth factor receptor (EGFR) but inhibits EGF-induced transactivation of a RAS-responsive element. Isoform 2 acts 
#as a dominant negative protein over GRB2 and by suppressing proliferative signals, may trigger active programmed cell death.
#Associate with host-virus interaction
#human protein


#>TRINITY_DN39048_c5_g1_TRINITY_DN39048_c5_g1_i1_g.8316_m.8316 type:5prime_partial len:107 TRINITY_DN39048_c5_g1_i1:361-41(-)
#TTTCGAGGGCCGGCGGGAGCTCGCCGGACGCCGCCGGAGACGCGGCGCTCTACGGGGCCGCAGCCCCTAGCTCCGGACAAGCCGATTCCAGGGACCCCGCCCCTTACCAAGAAAAGAAAACTCTTCCCGGGGCCCCGGCCGGCGTCTGCGGGCTCGCTCGCTTTACAGCTCGGACCGCCACGAAGACGGAGGGCTCCGCCTCCGGGTTCGGGAATATTGACCCGATTCCCTTTCGCCCGGGGAGCGAGGCCCCCGCTCACATGGAGAGGGGCGCCCGCCGCCGCTCGGAGCGGAACTTCCCTGGGGCTTAGGATCTACTGA
#No Blastx results
#Blast Asterias forbesii sea star 5.8S rRNA gene

#>TRINITY_DN39233_c0_g1_TRINITY_DN39233_c0_g1_i1_g.8447_m.8447 type:internal len:150 TRINITY_DN39233_c0_g1_i1:2-448(+)
#GGACGCGACGGATCAACTGGTTCTCCTGGTAATAAGGGATTCTCTGGTGATGTTGGACGAATTGGAATCCCAGGACTTCCTGGAGCCCCAGGTGACCGTGGTGGTGTGGGACAATCAGGTATTTCTGGAGCACCTGGTGCCCAAGGTTCATCAGGCAAGCGAGGAATTTCTGGAAAAGCCGGATCTAGAGGATCTCCAGGTCCACAAGGACCAGCTGGGTCTACAGGTCTGAGCGGTGAGCCGGGTGATATTGGAAGATCAGGAAGGACTGGTCCCCCTGGACCTAAAGGAAACCCGGGAGAGAGAGGAGGCTCTGGTCAACGAGGAGCTCAGGGTATTCCAGGAGGACAGGGACTAGCAGGAGGTGCCGGAGCCCGTGGGGAGAAGGGACCAAAAGGCTCCATCGGTGAGGCTGGAGTCCAGGGTGGCTCTGGCCCAACTGGACGA
#Human
#Collagen alpha-1(XI) chain
#May play an important role in fibrillogenesis by controlling lateral growth of collagen II fibrils.



#>TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543 type:internal len:103 TRINITY_DN39328_c5_g1_i1:3-308(+)
#TCCTCGAGGACCGAGTGGCGGAGAAGGGTTCCATGTGAACAGCAGTTGTACATGGGTCAGTCGATCCTAAGCCCCAGGGAAGTTCCGCTCCGAGCGGAGGCGGGCGCCCCTCTCCATGTGAGCGGGGGCCTCGCTCCCCGGGCGAAAGGGAATCGGGTCAATATTCCCGAACCCGGAGGCGGAGCCCTCCGTCTTCGTGGCGGTCCGAGCTGTAAAGCGAGCGAGCCCGCAGACGCCGGCCGGGGCCCCGGGAAGAGTTTTCTTCTCTTGGTAAGGGGCGGGATCCCTGGAATCGGCGTGTCCGGA
#BlastX results: uncharacterized protein
#blast results 5.8s rRNA


##################Plotcounts for DGE####################


par(mfrow=c(2,3))

plotCounts(ddsLoc, gene="TRINITY_DN41041_c3_g2_TRINITY_DN41041_c3_g2_i1_g.10613_m.10613", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN47102_c1_g1_TRINITY_DN47102_c1_g1_i3_g.24581_m.24581", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN39048_c5_g1_TRINITY_DN39048_c5_g1_i1_g.8316_m.8316", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN39233_c0_g1_TRINITY_DN39233_c0_g1_i1_g.8447_m.8447", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543", intgroup="location")
```



------
<div id='id-section25'/>
### Page 25:
------
<div id='id-section26'/>
### Page 26:
------
<div id='id-section27'/>
### Page 27:
------
<div id='id-section28'/>
### Page 28:
------
<div id='id-section29'/>
### Page 29:
------
<div id='id-section30'/>
### Page 30:
------
<div id='id-section31'/>
### Page 31:
------
<div id='id-section32'/>
### Page 32:
------
<div id='id-section33'/>
### Page 33:
------
<div id='id-section34'/>
### Page 34:
------
<div id='id-section35'/>
### Page 35:
------
<div id='id-section36'/>
### Page 36:
------
<div id='id-section37'/>
### Page 37:
------
<div id='id-section38'/>
### Page 38:
------
<div id='id-section39'/>
### Page 39:
------
<div id='id-section40'/>
### Page 40:
------
<div id='id-section41'/>
### Page 41:
------
<div id='id-section42'/>
### Page 42:
------
<div id='id-section43'/>
### Page 43:
------
<div id='id-section44'/>
### Page 44:
------
<div id='id-section45'/>
### Page 45:
------
<div id='id-section46'/>
### Page 46:
------
<div id='id-section47'/>
### Page 47:
------
<div id='id-section48'/>
### Page 48:
------
<div id='id-section49'/>
### Page 49:
------
<div id='id-section50'/>
### Page 50:
------
<div id='id-section51'/>
### Page 51:
------
<div id='id-section52'/>
### Page 52:
------
<div id='id-section53'/>
### Page 53:
------
<div id='id-section54'/>
### Page 54:
------
<div id='id-section55'/>
### Page 55:
------
<div id='id-section56'/>
### Page 56:
------
<div id='id-section57'/>
### Page 57:
------
<div id='id-section58'/>
### Page 58:
------
<div id='id-section59'/>
### Page 59:
------
<div id='id-section60'/>
### Page 60:

------




