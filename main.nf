#!/usr/bin/env nextflow

// Enable DSL2 syntax (Required for modern Nextflow scripts and Project 2)
nextflow.enable.dsl=2

/*
 * ====================================================================================
 *  Define Process Modules
 * ====================================================================================
 */

// 1. STAR_INDEX: Generate a genome index for STAR aligner
// This process runs only once. It takes the FASTA and GTF files as input.
process STAR_INDEX {
    label 'process_high'
    container 'ghcr.io/bf528/star:latest'
    
    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    """
    # Create the directory required by STAR to store index files
    mkdir star_index
    
    # Run STAR genomeGenerate mode with all available cores
    STAR --runThreadN ${task.cpus} \\
         --runMode genomeGenerate \\
         --genomeDir star_index \\
         --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile ${gtf}
    """
}

// 2. FASTQC: Perform quality control on raw sequencing reads
// Runs in parallel for each paired-end sample.
process FASTQC {
    label 'process_low'
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    // Emit both zip and html files for MultiQC and user inspection
    tuple val(sample_id), path("*.zip"), emit: zips
    path "*.html", emit: htmls

    script:
    """
    # Run fastqc on both R1 and R2
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}

// 3. GTF_SYMBOL: Parse GTF to create a Gene ID to Symbol mapping file
process GTF_SYMBOL {
    label 'process_low'
    container 'ghcr.io/bf528/pandas:latest' // Using pandas container for python env
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path parse_script
    path gtf

    output:
    path "gene_id_to_symbol.csv", emit: mapping

    script:
    """
    # Execute the python script you wrote in bin/ directory
    python ${parse_script} --gtf ${gtf} --out gene_id_to_symbol.csv
    """
}

// 4. STAR_ALIGN: Align reads to the reference genome
// Must wait for STAR_INDEX to complete.
process STAR_ALIGN {
    label 'process_high'
    container 'ghcr.io/bf528/star:latest'
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path gtf

    output:
    // Emit BAM file for VERSE quantification
    tuple val(sample_id), path("${sample_id}Aligned.out.bam"), emit: bam
    // Emit the final log file for MultiQC
    tuple val(sample_id), path("${sample_id}.Log.final.out"), emit: log

    script:
    """
    # Run STAR alignment. Redirect standard error (2>) to create the log file required for QC
    STAR --runThreadN ${task.cpus} \\
         --genomeDir ${index} \\
         --readFilesIn ${reads[0]} ${reads[1]} \\
         --outFileNamePrefix ${sample_id} \\
         --outSAMtype BAM Unsorted \\
         2> ${sample_id}.Log.final.out
    """
}

// 5. VERSE: Quantify alignments to gene level features
process VERSE_COUNT {
    label 'process_medium'
    container 'ghcr.io/bf528/verse:latest'
    publishDir "${params.outdir}/verse", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    // Output format is required to be .exon.txt
    path "${sample_id}.exon.txt", emit: counts

    script:
    """
    # Run verse with strandedness (-S) or other default options
    verse -S -a ${gtf} -o ${sample_id} ${bam}
    """
}

// 6. MERGE_COUNTS: Aggregate individual VERSE counts into a single matrix
// Executes only after ALL VERSE processes have successfully completed.
process MERGE_COUNTS {
    label 'process_low'
    container 'ghcr.io/bf528/pandas:latest'
    publishDir "${params.outdir}/matrix", mode: 'copy'

    input:
    path merge_script
    path count_files // A list containing all .exon.txt files

    output:
    path "raw_counts_matrix.csv", emit: matrix

    script:
    """
    # Execute python script to merge all provided count files
    python ${merge_script} --out raw_counts_matrix.csv ${count_files}
    """
}

// 7. MULTIQC: Aggregate QC metrics from FastQC and STAR into an HTML report
process MULTIQC {
    label 'process_low'
    container 'ghcr.io/bf528/multiqc:latest'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files // A collection of FastQC zips and STAR logs

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    # Run multiqc on current directory (.) and force overwrite (-f)
    multiqc -f .
    """
}


/*
 * ====================================================================================
 *  Main Workflow Execution Logic
 * ====================================================================================
 */

workflow {
    // ---------------------------------------------------------
    // 1. Create Input Channels
    // ---------------------------------------------------------
    // Pair R1 and R2 fastq files based on common prefix
    align_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    
    // Create value channels for reference files (they are reused multiple times)
    genome_ch = file(params.genome)
    gtf_ch    = file(params.gtf)

    // ---------------------------------------------------------
    // 2. Execute Independent/Upstream Processes
    // ---------------------------------------------------------
    // Run FastQC on raw reads
    FASTQC(align_ch)

    // Run GTF parsing script
    script_gtf_ch = file('bin/gtf_to_symbol.py')
    GTF_SYMBOL(script_gtf_ch, gtf_ch)

    // Generate STAR index
    STAR_INDEX(genome_ch, gtf_ch)

    // ---------------------------------------------------------
    // 3. Execute Dependent/Downstream Processes
    // ---------------------------------------------------------
    // Align reads. Wait for STAR_INDEX to emit its output.
    // The inputs match the definition: tuple(reads), path(index), path(gtf)
    ALIGN_HARD(align_ch, STAR_INDEX.out.index, gtf_ch)

    // Run VERSE quantification on the BAM files produced by ALIGN_HARD
    VERSE_COUNT(ALIGN_HARD.out.bam, gtf_ch)

    // Merge counts. Use .collect() to gather all scatter outputs into a single list
    script_merge_ch = file('bin/merge_counts.py')
    MERGE_COUNTS(script_merge_ch, VERSE_COUNT.out.counts.collect())

    // ---------------------------------------------------------
    // 4. Aggregate QC Reports
    // ---------------------------------------------------------
    // Strip sample_id from the tuples and extract only the actual files
    fastqc_files = FASTQC.out.zips.map { sample_id, files -> files }
    star_logs    = ALIGN_HARD.out.log.map { sample_id, log -> log }

    // Mix the two streams of files and collect them into a single list for MultiQC
    multiqc_input = fastqc_files.mix(star_logs).collect()
    MULTIQC(multiqc_input)
}
