# ONT_Illumina_transcriptome_assembly-annotation

## **ASSEMBLY PROCEDURE**

Assembling RNA-seq data from Illumina reads and Direct RNA Nanopore sequencing reads with a reference genome.

Outline of the process:

**1) Quality Control:** Trim adapters, filter low-quality reads

- For Illumina reads use ``FastQC`` before and after Trimming adapters (``trimmgalore``)
- For Direct RNA Nanopore reads (``Porechop or Guppy``) 
  
**2) Alignment to Reference Genome:**

First, we must create separate alignments with Illumina and direct ONT reads.

  **Short-Reads Alignment:**
  
  **A)** Index the reference genome using ``HISAT2`` or another appropriate indexing tool for short reads
  
    hisat2-build reference_genome.fasta reference_index

  **B)** Align Illumina reads to the reference genome
    
    hisat2 -x reference_index -1 illumina_read1.fastq -2 illumina_read2.fastq -S illumina_alignment.sam 
  
  **Long-Reads Alignment:**
  
  **C)** Align Direct RNA Nanopore reads to the reference genome using the long-read aligner ``minimap2``.
    
    minimap2 -ax splice -uf reference_genome.fasta nanopore_reads.fastq | samtools sort -o nanopore_alignment.bam   

  **D)** Convert SAM to BAM and Sort:
    
    samtools view -bS illumina_alignment.sam > illumina_alignment.bam samtools sort -o sorted_illumina_alignment.bam illumina_alignment.bam 

# Index the reference genome or transcriptome flair index -g reference.fasta # Align the long reads to the reference flair align -r reads.fastq -i index_dir -v # Quantify isoform abundance flair quantify -r reads.bam -i index_dir -q 
Another tool you can consider is Sqanti (from the IsoSeq toolkit), which is specifically designed for analyzing isoform diversity and quantification in long-read sequencing data

**3) Transcriptome Assembly:**

I assembled the transcriptome using StringTie because it has a "-mix" option, which allows the assembly of both Illumina and Direct RNA Nanopore sequencing reads simultaneously. This step generates a GTF file containing predicted transcripts. You can include a reference GTF file to 
    
    stringtie sorted_illumina_alignment.bam flair_aligned_ONT.bam --mix -p 24 -G reference.gtf -o transcriptome_output.gtf

**3) Transcriptome Evaluation:**

I used ``BUSCO`` to asses transcriptome completeness.

-	Visualize the assembled transcripts along with the reference genome using genome browsers like IGV or UCSC Genome Browser.
  
## **ANNOTATION**

-	Annotate the assembled transcripts using tools like gffcompare or PASA.

# IN PROCESS

Quantification and Differential Expression Analysis:
-	Quantify gene expression using tools like featureCounts or Salmon.
-	Perform differential expression analysis using tools like DESeq2 or edgeR
-	Isoform Identification:
Use tools like FLAIR, SQANTI, or TAMA to classify and identify isoforms. These tools leverage long-read information to identify full-length isoforms.
- Compare with Reference Annotation:
Compare the assembled transcriptome with a reference annotation (if available) using tools like gffcompare.
- Isoform Detection and Quantification:
After aligning the reads, you can use tools designed for isoform detection and quantification. One popular tool for this purpose is Flair. Flair combines information from splice junctions and poly(A) sites to identify and quantify full-length isoforms. You can find it here: Flair GitHub repository.
#flair align -r <genome_annotation.gtf> -b <input_bam_file1> <input_bam_file2> ... -o <output_directory>

