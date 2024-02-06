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

Index the reference genome or transcriptome flair index -g reference.fasta # Align the long reads to the reference flair align -r reads.fastq -i index_dir -v # Quantify isoform abundance flair quantify -r reads.bam -i index_dir -q 
Another tool you can consider is Sqanti (from the IsoSeq toolkit), which is specifically designed for analyzing isoform diversity and quantification in long-read sequencing data

**3) Transcriptome Assembly:**

I assembled the transcriptome using StringTie because it has a "-mix" option, which allows the assembly of both Illumina and Direct RNA Nanopore sequencing reads simultaneously. This step generates a GTF file containing predicted transcripts. You can include a reference GTF file to keep a relationship with the genome and make the annotation: 
    
    stringtie sorted_illumina_alignment.bam flair_aligned_ONT.bam --mix -p 24 -G reference.gtf -o transcriptome_output.gtf

**4) Creating transcriptome:** 

I used ``gffread`` read to create a FASTA file containing the transcripts base in the genome Fasta file and the GTF generated from the previous step.

		> ./gffread/gffread -w transcriptome_output.gtf -g <genome.fasta> <assemblied_transcriptome.fasta>

**5) Transcriptome Evaluation:**

I used ``BUSCO`` to asses transcriptome completeness.

		> busco -i Ahem_transcripts.fa -l metazoa_odb10 -o transcriptome_stringtie -m tran -f




sed '/>/!s/U/T/g' /ibex/scratch/projects/c2078/Villegmb/20230913_P4U2_2/Ahem_transcript_rename.fasta > /ibex/scratch/pro
jects/c2078/Villegmb/20230913_P4U2_2/Ahem_transcript_DNA.fasta

export BUSCO_CONFIG_FILE="/path/to/myconfig.ini"

busco -i /ibex/scratch/projects/c2078/Villegmb/20230913_P4U2_2/Ahem_transcript_DNA.fasta -l metazoa_odb10 -o transcripto
me_AHEM -m tran -f -c 24
/ibex/tmp/c2078/Heat_stress_analysis/scripts/transcripts_without_predicted_proteins.tsv


/ibex/tmp/c2078/Heat_stress_analysis/scripts/Ahem_transcripts_with_reference.fa.transdecoder_dir/longest_orfs.pep



-	Visualize the assembled transcripts and the reference genome using genome browsers like IGV or UCSC Genome Browser.
  
## **ANNOTATION**

**A)** In this part, TransDecoder identifies possible Open Reading Frames (ORF) in our transcripts and produces a multiFasta file with the longest possible proteins (.pep extension). Then, this file is used to find homologies with protein domains in the PFAM database and annotated proteins in the Swiss-Prot, TrEMBL, and NR databases.

For more information about TransDecoder https://github.com/TransDecoder/TransDecoder/wiki

By default, TransDecoder.LongOrfs will identify ORFs that are at least 100 amino acids long. 

**Step 1:** Extract the long open reading frames
	
 	> TransDecoder.LongOrfs -t <assemblied_transcriptome.fasta>

TransDecoder produces a file in a new directory "./<assemblied_transcriptome.fasta>.transdecoder_dir/", where you can find the file "longest_orfs.pep" that we are going to use in the next step. 

**Step 2:** Identify ORFs with homology 

To search for homologies, we need to download and decompress the files to create the databases. Let's start with the databases for Blastp searches: 

_________________________________________________________________
**Databases**

**2.1)** Use these links to download the files, and decompress them:

**SwissProt:** ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

**TrEMBL:** ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

**2.2)** Rename the files and use ``ncbi-blast+`` (I used version 2.13.0) to create the blast databases:

		> mv uniprot_sprot.fasta sport
		> mv uniprot_trembl.fasta trembl
		> makeblastdb -in sprot -dbtype prot -title "UniProt/Swiss-Prot (Jan 2024)" -parse_seqids
		> makeblastdb -in trembl -dbtype prot -title "UniProt/TrEMBL (Jan 2024)" -parse_seqids

_________________________________________________________________

**NR Database from NCBI**

This database is composed of 83 databases. You can use wget "link" and tar loops to decompress the files from:

ftp://ftp.ncbi.nlm.nih.gov/blast/db/ and look for files in the pattern of nr.*.tar.gz.


or you can use the next command from ``ncbi-blast+``:

		> update_blastdb.pl --decompress nr

It checks all the files and updates them. However, while I was running this code, it stopped several times. So, I recommend running it repeatedly until you can see all the decompressed files in the folder. To date (February 2024), you should see up to 83 sets of files. These files are already formatted as a database, so you can skip the "makeblatdb" step you did for the other two databases. Just ensure you run the blast in the same folder where you stored the databases.

Now we can run Blastp:

		> blast -query /.../<assemblied_transcriptome.fasta>.transdecoder_dir/longest_orfs.pep \
		> -db sprot -max_target_seqs 20 \
		> -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6_sprot
_________________________________________________________________

		> blastp -query /.../<assemblied_transcriptome.fasta>.transdecoder_dir/longest_orfs.pep \
		> -db trembl -max_target_seqs 20 \
		> -outfmt 6 -evalue 1e-5 -num_threads 50 > blastp.outfmt6_trembl

_________________________________________________________________

		> blastp -query /.../<assemblied_transcriptome.fasta>.transdecoder_dir/longest_orfs.pep \
		> -db nr -max_target_seqs 20 \
		> -outfmt 6 -evalue 1e-5 -num_threads 50 > blastp.outfmt6_nr

_________________________________________________________________

**PFAM** 

For downloading the database: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 

Pfam checks for homologies with protein domain families; we use HMMER (http://hmmer.org/) to run the search.

		> export HMMERDB=/ibex/tmp/c2078/Heat_stress_analysis/scripts/Pfam-A.hmm
		> hmmsearch --cpu 24 -E 1e-10 --domtblout pfam.domtblout $HMMERDB /ibex/tmp/c2078/Heat_stress_analysis/scripts/Ahem_transcripts_with_reference.fa.transdecoder_dir/longest_orfs.pep

**Step 3:** Predict the likely coding regions

		> TransDecoder.Predict -t Ahem_transcripts_with_reference.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

From this step, I obtained a new Fasta containing the proteins that showed homology to Pfam and Swiss-prot entries. This Fasta file is going to be used then to find GO terms.  
Usando ese Blastp y el PFAM 

A ver: 
PFAM me da una lista de proteinas con homologia a algo
Blasp Swiss prot otra homologia a otra cosa

Pueden haber proteinas que coinciden
Pueden haber algunas que hacen match con PFAM y no SP
y vice-versa


Con el final PEP luego de predicted es que debería repetir busquedas para encontrar GO 

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



I checked using R if there where predicted_proteins from transcoder (including PFAM and Blastp against Swiss-Prot) for all the transcripts. Then I filtered the "longest PEP proteins fasta", to run blast ensembl and NR against this.  

awk -F'\t' '{print $2}' transcripts_without_predicted_proteins.tsv > output_column_2.txt
https://github.com/meiyang12/Genome-annotation-pipeline?tab=readme-ov-file


gffread -E input.gff3 -T -o output.gtf

#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J gffread_gff_to_gtf
#SBATCH -o gffread_gff_to_gtf.%J.out
#SBATCH -e gffread_gff_to_gtf.%J.err
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=24

./gffread/gffread -E /ibex/tmp/c2078/Heat_stress_analysis/scripts/Ahem_transcripts_with_reference.fa.transdecoder.gff3 -T -o Ahem_transcripts_with_reference.fa.transdecoder_output.gtf

gffcompare/0.12.7



cat blastp.outfmt6_trembl_* > merged_output.outfmt6





