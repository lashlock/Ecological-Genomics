# 2017 Ecological Genomics Course

### Author: Lauren Ashlock     


## Overall Description of notebook      
Notes from class material,and class project will populate this notebook. 

## Date started: (2017-02-01)   
## Date end:   ongoing    


### Table of contents for 60 entries (Format is *Page: Date(with year-month-day). Title*)        
* [Page 1: 2017-02-01](#id-section1). Sequencing strategies applied to biological questions
* [Page 2: 2017-02-06](#id-section2).  RNAseq   
* [Page 3: 2017-02-08](#id-section3). 
* [Page 4: 2017-02-13](#id-section4) . 
* [Page 5: 2017-02-15](#id-section5). 
* [Page 6: 2017-02-22](#id-section6). DESEQ 2 Tutorial
* [Page 7: 2017-02-27](#id-section7). Scott Edwards and Differential Expression Analysis
* [Page 8: 2017-03-01](#id-section8). Differential expression - Catch up Day
* [Page 9:2017-03-06](#id-section9). Population Genomics
* [Page 10:2017-03-08](#id-section10).Effective population size
* [Page 11:](#id-section11).
* [Page 12:](#id-section12).
* [Page 13:](#id-section13).
* [Page 14:](#id-section14).
* [Page 15:](#id-section15).
* [Page 16:](#id-section16).
* [Page 17:](#id-section17).
* [Page 18:](#id-section18).
* [Page 19:](#id-section19).
* [Page 20:](#id-section20).
* [Page 21:](#id-section21).
* [Page 22:](#id-section22).
* [Page 23:](#id-section23).
* [Page 24:](#id-section24).
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



------
<div id='id-section5'/>
### Page 5:   

**Info Update: SNPs and Population Genomics**

- SNP data - expressed sequences

1.    Tissue
      1. breadth of tissue from different developmental stages
         1. This controls for exon skipping
2.    Pool samples and create your sequence libraries (30-100Million p.e. long reads)
3.    Process raw sequence data
      1. Important for SNP detection
4.    Digital normalization
      1. Remove high coverage reads and associated errors
         1. Loss of quantitative info
5.    Assemble clean p.e. reads
6.    Prune assembled transcripts
      1. Reduce DNA contaminations, noncoding DNA, and gene fragments
7.    Assembly evaluation using either a reference genome or conserved genes in other eukaryotic organisms
8.    SNP detection
      1. Software: constant patterns of sequence variation
         1. sequence errors - hopefully your software will eliminate reads of low frequency
         2. Errors can also be found in homogeneous regions... amplified by PCR... can filter for these
         3. Artifacts caused by InDels - filter SNP clusters near indels, quality scores
9.    SNP validation - primers
      1. Use Sanger sequencing or mass spec to quality control a portion of your sequence data
10.    Applications
11.    Differences in population structure
12.    How natural selection is acting on particular loci
13.    Methods for Applications
14.    Outlier - for a given locus, whats the level of differentiation compared to differences across the genome? Using Fst
15.    Non outlier - Tests high Fst loci for other features associated with selection
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
------
<div id='id-section12'/>
### Page 12:
------
<div id='id-section13'/>
### Page 13:
------
<div id='id-section14'/>
### Page 14:
------
<div id='id-section15'/>
### Page 15:
------
<div id='id-section16'/>
### Page 16:
------
<div id='id-section17'/>
### Page 17:
------
<div id='id-section18'/>
### Page 18:
------
<div id='id-section19'/>
### Page 19:
------
<div id='id-section20'/>
### Page 20:
------
<div id='id-section21'/>
### Page 21:
------
<div id='id-section22'/>
### Page 22:
------
<div id='id-section23'/>
### Page 23:
------
<div id='id-section24'/>
### Page 24:
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


