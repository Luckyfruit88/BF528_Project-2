nextflow.enable.dsl=2

params.reads = '/projectnb/bf528/materials/project-2-rnaseq/subsampled_files/*_R{1,2}.subset.fastq.gz'
params.outdir = 'results'
params.gtf = "/projectnb/bf528/materials/project-2-rnaseq/refs/gencode.v45.primary_assembly.annotation.gtf"
params.genome = "/projectnb/bf528/materials/project-2-rnaseq/refs/GRCh38.primary_assembly.genome.fa"

include { FASTQC }       from './modules/fastqc'
include { GTF_SYMBOL }   from './modules/gtf2symbol'
include { ALIGN_HARD }   from './modules/align_hard'
include { STAR_INDEX }   from './modules/star_index'
include { MULTIQC }      from './modules/multiqc'
include { VERSE_COUNT }  from './modules/verse'
include { MERGE_COUNTS } from './modules/merge_counts'

workflow {
    // 1. Create paired reads channel
    align_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    fastqc_channel = align_ch.transpose()

    // 2. FastQC quality control
    FASTQC(fastqc_channel)

    // 3. GTF annotation conversion
    script_ch = Channel.value(file('bin/gtf_to_symbol.py'))
    gtf_ch    = Channel.value(params.gtf)
    GTF_SYMBOL(script_ch, gtf_ch)

    // 4. STAR index generation
    genome_ch = Channel.value(params.genome)
    STAR_INDEX(genome_ch, gtf_ch)

    // 5. Sequence alignment
    align_input = align_ch
        .combine(STAR_INDEX.out.index)
        .combine(gtf_ch)
        .map { sample, reads, star_index, gtf -> 
            tuple(sample, reads, star_index, gtf)
        }
    ALIGN_HARD(align_input)

    // 6. VERSE gene counting (simplified because bam is now a tuple)
    verse_input = ALIGN_HARD.out.bam.combine(gtf_ch)
    VERSE_COUNT(verse_input)

    // 7. Merge all VERSE counts into single matrix
    merge_script_ch = Channel.value(file('bin/merge_counts.py'))
    all_counts = VERSE_COUNT.out.counts
        .map { sample, file -> file }
        .collect()
    MERGE_COUNTS(merge_script_ch, all_counts)

    // 8. MultiQC summary report
    multiqc_input = FASTQC.out.zip
        .map { sample, file -> file }
        .mix(ALIGN_HARD.out.log)
        .collect()
    MULTIQC(multiqc_input)
}
