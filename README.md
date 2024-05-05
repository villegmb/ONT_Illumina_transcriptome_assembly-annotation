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
  
    > hisat2-build reference_genome.fasta reference_index

  **B)** Align Illumina reads to the reference genome
    
    > hisat2 -x reference_index -1 illumina_read1.fastq -2 illumina_read2.fastq -S illumina_alignment.sam

**C)** Convert SAM to BAM, Sort and Index:
 
    > samtools view -bS illumina_alignment.sam > illumina_alignment.bam 
    > samtools sort -o sorted_illumina_alignment.bam illumina_alignment.bam
		> samtools index sorted_illumina_alignment.bam
	
  **Long-Reads Alignment:**
  
  **A)** Align Direct RNA Nanopore reads to the reference genome using the long-read aligner ``flair``, sort and index:

		flair uses ``minimap2 -ax splice`` which is recommended in Stringtie's manual	
	
   	> flair align -g genome.fa -r <reads.fq>|<reads.fa> --threads 24
		> samtools sort -o sorted_flair.aligned.bam flair.aligned.bam
		> samtools index sorted_flair.aligned.bam

**3) Transcriptome Assembly:**

I assembled the transcriptome using StringTie because it has a "-mix" option, which allows the assembly of both Illumina and Direct RNA Nanopore sequencing reads simultaneously. This step generates a GTF file containing predicted transcripts. You can include a reference GTF file to keep a relationship with the genome and make the annotation: 
    
    > stringtie sorted_illumina_alignment.bam sorted_flair.aligned.bam --mix -p 24 -G reference.gtf -o transcriptome_annotation_output.gtf

**4) Creating transcriptome:** 

NOTE: If you are using TransDecoder to predict protein regions skip this step and use the TransDecoder utility to produce the fasta with the transcripts (ANNOTATION -> TRANSDECODER. A - Step 0)

``gffread`` read to create a FASTA file containing the transcripts base in the genome Fasta file and the GTF generated from the previous step.

	> ./gffread/gffread -w transcriptome_annotation_output.gtf -g <genome.fasta> <assemblied_transcriptome.fasta> 

**5) Transcriptome Evaluation:**

I used ``BUSCO`` to asses transcriptome completeness.

	> busco -i Ahem_transcripts.fa -l metazoa_odb10 -o transcriptome_stringtie -m tran -f
		
	> sed '/>/!s/U/T/g' /ibex/scratch/projects/c2078/Villegmb/20230913_P4U2_2/Ahem_transcript_rename.fasta > Ahem_transcript_DNA.fasta

	> busco -i /ibex/scratch/projects/c2078/Villegmb/20230913_P4U2_2/Ahem_transcript_DNA.fasta -l metazoa_odb10 -o transcriptome_AHEM -m tran -f -c 24
 
 **6) Transcriptome visualization:**

I visualized the assembled transcripts, the GTF and the reference genome using IGV genome browser
 
## **ANNOTATION** 

### TRANSDECODER

**A)** In this part, TransDecoder identifies possible Open Reading Frames (ORF) in our transcripts and produces a multiFasta file with the longest possible proteins (.pep extension). Then, this file is used to find homologies with protein domains in the PFAM database and annotated proteins in the Swiss-Prot, TrEMBL, and NR databases.

For more information about TransDecoder https://github.com/TransDecoder/TransDecoder/wiki

**Step 0:**

	> ./TransDecoder/util/gtf_genome_to_cdna_fasta.pl transcriptome_annotation_output.gtf reference_genome.fasta > transcripts.fasta

 You must also convert the GTF to GFF3, as TransDecoder uses this format.

 	> ./TransDecoder/util/gtf_to_alignment_gff3.pl transcriptome_annotation_output.gtf > transcripts.gff3

**Step 1:** Extract the long open reading frames

	By default, TransDecoder.LongOrfs will identify ORFs that are at least 100 amino acids long. 

 	> TransDecoder.LongOrfs -t transcripts.fasta

TransDecoder produces a file in a new directory "./<transcripts.fasta>.transdecoder_dir/", where you can find the file "longest_orfs.pep" that we are going to use in the next step. 

**Step 2:** Identify ORFs with homology 

To search for homologies, we need to download and decompress the files to create the databases. Let's start with the databases for Blastp searches: 

_________________________________________________________________
**Databases**

**2.1)** Use these links to download the files and decompress them:

**SwissProt:** ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

**TrEMBL:** ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

**2.2)** Rename the files and use ``ncbi-blast+`` (I used version 2.13.0) to create the blast databases:

		> mv uniprot_sprot.fasta sprot
		> mv uniprot_trembl.fasta trembl
		> makeblastdb -in sprot -dbtype prot -title "UniProt/Swiss-Prot (Jan 2024)" -parse_seqids
		> makeblastdb -in trembl -dbtype prot -title "UniProt/TrEMBL (Jan 2024)" -parse_seqids

_________________________________________________________________

**NR Database from NCBI**

This database is composed of 83 databases. You can use wget "link" and tar loops to decompress the files from:

ftp://ftp.ncbi.nlm.nih.gov/blast/db/ and look for files in the pattern of nr.*.tar.gz.


or you can use the next command from ``ncbi-blast+``:

		> update_blastdb.pl --decompress nr

It checks all the files and updates them. However, while I was running this code, it stopped several times. So, I recommend running it repeatedly until you can see all the decompressed files in the folder. To date (February 2024), you should see up to 83 sets of files. These files are already formatted as a database, so you can skip the "makeblatdb" step you did for the other two databases. Please ensure you run the blast in the same folder where you stored the databases.

Now we can run Blastp:

		for transdecoder:
		> blast -query /.../<transcripts.fasta>.transdecoder_dir/longest_orfs.pep \
		> -db sprot -max_target_seqs 1 \
		> -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6_1_sprot
	
		> blast -query /.../<transcripts.fasta>.transdecoder_dir/longest_orfs.pep \
		> -db sprot -max_target_seqs 20 \
		> -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6_20_sprot
_________________________________________________________________

		> blastp -query /.../<transcripts.fasta>.transdecoder_dir/longest_orfs.pep \
		> -db trembl -max_target_seqs 20 \
		> -outfmt 6 -evalue 1e-5 -num_threads 50 > blastp.outfmt6_trembl

_________________________________________________________________

		> blastp -query /.../<transcripts.fasta>.transdecoder_dir/longest_orfs.pep \
		> -db nr -max_target_seqs 20 \
		> -outfmt 6 -evalue 1e-5 -num_threads 50 > blastp.outfmt6_nr

_________________________________________________________________

**PFAM** 

For downloading the database: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 

Pfam checks for homologies with protein domain families; we use HMMER (http://hmmer.org/) to run the search.

		> export HMMERDB=/.../Pfam-A.hmm
		> hmmsearch --cpu 24 -E 1e-10 --domtblout pfam.domtblout $HMMERDB /.../<transcripts.fasta>.transdecoder_dir/longest_orfs.pep

**Step 3:** Predict the likely coding regions

		> TransDecoder.Predict -t <transcripts.fasta> --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6_1_sprot

From this step, I obtained a new Fasta containing the proteins that showed homology to Pfam and Swiss-prot entries. I also get a new GFF3 with the annotation information for the predicted proteins, then I am going to produce the final GFF3 with genomic coordinates:

		> ./TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl \
		> ./<transcripts.fasta>.transdecoder.gff3 \
		> transcripts.gff3 \
		> <transcripts.fasta> > <transcripts.fasta>.transdecoder.genome_coordinates.gff3

 Then, you can use ´´gffread´´ (0.12.7) to get a new fasta with the predicted transcripts using the final GFF3:

	> ./gffread/gffread <transcripts.fasta>.transdecoder.genome_coordinates.gff3 -g reference_genome.fasta -w <new_transcripts_using_transdecoder_final_GFF3.fasta>

or you can convert GFF3 to GTF, and use the GTF to produce the final fasta file using TransDecoder utilities.


### **GO TERMS ANNOTATION** 

Adapted from https://github.com/lyijin/annotating_proteomes/

**Step 4:** Prepare GO Terms files

Use wget to download the information needed to match proteins with GO terms:
- GO annotation file: http://www.geneontology.org/gene-associations/goa_uniprot_all.gaf.gz
- GO term hierarchy: http://purl.obolibrary.org/obo/go/go-basic.obo

Run the shell script parse_gp_assoc.sh to produce goa_uniprot_all.parsed.gaf and goa_uniprot_all.unique_ids.txt. Take note of the directories where you store those files and modify the paths in the next codes: 

- Modify line 64 of create_go_annots_sprot_trembl.py to where you kept goa_uniprot_all.parsed.gaf.
- Modify line 53 of get_top_hit_with_amigo_annot.py to where you kept goa_uniprot_all.unique_ids.txt.
- Modify line 35 of parse_go_obo.py to where you kept go-basic.obo.

NOTE: In Lyijin's pipeline (https://github.com/lyijin/annotating_proteomes/tree/master) at this point he does the step of "Parsing the XML outputs", but as I filtered my proteins directly when I ran blastp using "-evalue 1e-5" and setting the output to be tabular (Output 6), then I skipped this step.

Now, to guarantee I keep just the proteins with a GO term, I run the following codes:

	> get_top_hit_with_amigo_annot.py blastp.outfmt6_sprot > spis_vs_sprot.tGO.tsv

	> get_top_hit_with_amigo_annot.py blastp.outfmt6_trembl > spis_vs_trembl.tGO.tsv

I must clarify that I divided my files into many small pieces to save some time. Otherwise, this annotation process would have taken more than a month. For example, I divided the file "longest_orfs.pep" into 100 different files to run the search against trEMBL and NR databases because it takes more than a week (even when I was using 50 threads). I can access a supercomputer with different nodes to run many processes in parallel. Then, I merged the output into one file called "blastp.outfmt6_*_complete. For running the "get_top_hit_with_amigo_annot.py" I divided the "blastp.outfmt6_*_complete," taking care that all the proteins from the same transcript were in the same file, and then I assigned the go terms for these files.    

Follow his pipeline to finish the annotation process.


Now, the transcriptome is annotated and ready to use for your analysis. 

# ANALYSIS

Some examples of analysis I did:

**1)** Differential Gene Expression Analysis


**2)** Differential Isoforms Usage and Alternative Splicing

https://github.com/villegmb/ONT_Transcriptome_analysis

**3)** Identification of mRNA modifications

https://github.com/villegmb/Methylation_ONT_transcriptome
