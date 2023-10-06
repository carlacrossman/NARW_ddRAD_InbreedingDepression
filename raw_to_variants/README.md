## From raw reads to genotype calls

- *01_prelim_multiqc.sh*	Run multiqc on the fastq from sequencing facility
- *02_ddRAD_trim_array.sh*	Trim raw sequence reads
- *03_rw_align_genome.sh*	Prepare reference genome
- *04_rw_map_bwa_mem.sh*	Map to reference assembly and sort .bam files
- *05_qualimap.sh*	Run qualimap to assess mapping statistics
- *06a_supp_addrg.sh*	Add readgroup header to bam files
- *06b_index_bam.sh*	Index bam files
- *07a_freebayes_region.sh*	Variant calling with freebayes
- *07b_freebayes_regions_complex.sh*	Variant calling with freebayes with additonal flag to limit haplotypes called in complex regions 
- *08a_stacks_ref_gt0.01.sh*	Running the stacks_ref pipeline with alphas of 0.01
- *08b_stacks_ref_gt0.001.sh*	Running the stacks_ref pipeline with alphas of 0.001


![Fig1_workflow](https://github.com/carlacrossman/NARW_ddRAD_InbreedingDepression/blob/main/raw_to_variants/ddRAD_preprocessing_streamlined.png)