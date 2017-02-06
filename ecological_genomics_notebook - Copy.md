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
* [Page 6:](#id-section6).
* [Page 7:](#id-section7).
* [Page 8:](#id-section8).
* [Page 9:](#id-section9).
* [Page 10:](#id-section10).
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
   - ​

   ​

   ​

------
<div id='id-section3'/>
### Page 3:   
------
<div id='id-section4'/>
### Page 4:   
------
<div id='id-section5'/>
### Page 5:   
------
<div id='id-section6'/>
### Page 6:
------
<div id='id-section7'/>
### Page 7:
------
<div id='id-section8'/>
### Page 8:
------
<div id='id-section9'/>
### Page 9:
------
<div id='id-section10'/>
### Page 10:
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


