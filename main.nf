#!/usr/bin/env nextflow

log.info """\
         POLISH MY READS! - N F   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/*
 * Specify two channels with read pairs.
 * read_pairs_ch --> use for fastqc
 * read_pairs2_ch --> use as input for superdeduper or directly into quality trim
 */


Channel
    .fromFilePairs( params.reads, flat:true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch } 


process fastqc {
    
/*
 * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
 * Note that the fastqc results are stored in a channel but that this channel will be unused (for now)
 */

    tag "FASTQC on $sample_id"
    publishDir "${params.outdir}/01.fastqc", mode:'copy'

    input:
    set val(sample_id), file(r1), file(r2) from read_pairs_ch

    output:
    file("${sample_id}")
    set val(sample_id), file("$sample_id/${sample_id}.stats") into fastqc_stats_ch

    script:
    """
    mkdir $sample_id

    fastqc --threads ${task.cpus} --outdir $sample_id $r1 $r2

    seqkit stats $r1 $r2 > $sample_id/${sample_id}.stats
    """
}

if (!params.skip_dedup) {
    
    process dedup {
    
    /*
     * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
     */

        label 'RAM_high'

        tag "Superdeduper on $sample_id"
        publishDir "${params.outdir}/02.dedup/$sample_id", mode:'copy'

        input:
        set val(sample_id), file(r1), file(r2) from read_pairs2_ch

        output:
        set val(sample_id), \
        file("PE_R1.fastq.gz"), \
        file("PE_R2.fastq.gz") into dedup_ch
        set val(sample_id), file("${sample_id}_dedup.stats") into dedup_stats_ch

        script:
        """
        hts_SuperDeduper -1 $r1 -2 $r2 -f PE

        seqkit stats PE_R1.fastq.gz PE_R2.fastq.gz > ${sample_id}_dedup.stats
        """
    }

    process trim_adapt_dedup {
        
    /*
     * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
     */

        tag "Trimmomatic on $sample_id"
        publishDir "${params.outdir}/03.trim_adapt/$sample_id", pattern: '*_{R1,R2,U}.fastq.gz', mode:'copy'
        publishDir "${params.outdir}/03.trim_adapt/$sample_id", pattern: '*.stats', mode:'copy', overwrite: false

        input:
        set val(sample_id), file(r1), file(r2) from dedup_ch 

        output:
        set val(sample_id), \
        file("${sample_id}_adapt_R1.fastq.gz"), \
        file("${sample_id}_adapt_R2.fastq.gz"), \
        file("${sample_id}_adapt_U.fastq.gz") into trim_adapt_ch
        set val(sample_id), file("${sample_id}_adapt.stats") into trim_adapt_stats_ch

        script:
        """

        trimmomatic PE \
        -threads ${task.cpus} \
        $r1 $r2 \
        ${sample_id}_adapt_R1.fastq.gz \
        ${sample_id}_1_u.fastq.gz \
        ${sample_id}_adapt_R2.fastq.gz \
        ${sample_id}_2_u.fastq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10:8:TRUE

        cat ${sample_id}_1_u.fastq.gz ${sample_id}_2_u.fastq.gz > ${sample_id}_adapt_U.fastq.gz

        seqkit stats ${sample_id}_adapt_R1.fastq.gz \
        ${sample_id}_adapt_R2.fastq.gz \
        ${sample_id}_adapt_U.fastq.gz > ${sample_id}_adapt.stats

        """
    }

} else {

    process trim_adapt {
        
    /*
     * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
     */

        tag "Trimmomatic on $sample_id"
        publishDir "${params.outdir}/02.trim_adapt/$sample_id", pattern: '*_{R1,R2,U}.fastq.gz', mode:'copy'
        publishDir "${params.outdir}/02.trim_adapt/$sample_id", pattern: '*.stats', mode:'copy', overwrite: false
   
        input:
        set val(sample_id), file(r1), file(r2) from read_pairs2_ch

        output:
        set val(sample_id), \
        file("${sample_id}_adapt_R1.fastq.gz"), \
        file("${sample_id}_adapt_R2.fastq.gz"), \
        file("${sample_id}_adapt_U.fastq.gz") into trim_adapt_ch
        set val(sample_id), file("${sample_id}_adapt.stats") into trim_adapt_stats_ch

        script:
        """
        trimmomatic PE \
        -threads ${task.cpus} \
        $r1 $r2 \
        ${sample_id}_R1_adapt.fastq.gz \
        ${sample_id}_1_u.fastq.gz \
        ${sample_id}_R2_adapt.fastq.gz \
        ${sample_id}_2_u.fastq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10:8:TRUE

        cat ${sample_id}_1_u.fastq.gz ${sample_id}_2_u.fastq.gz > ${sample_id}_adapt_U.fastq.gz

        seqkit stats ${sample_id}_adapt_R1.fastq.gz \
        ${sample_id}_adapt_R2.fastq.gz \
        ${sample_id}_adapt_U.fastq.gz > ${sample_id}_adapt.stats

        """
    }
}

if (!params.skip_dedup){
    merge_path = "${params.outdir}/04.merged_reads/"
} else {
    merge_path = "${params.outdir}/03.merged_reads/"
}


process merge {
    
/*
 * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
 */

    tag "Pear on $sample_id"

    publishDir "$merge_path/$sample_id", pattern: '*.unassembled.forward.fastq.gz', saveAs: {"${sample_id}_merge_R1.fastq.gz"}, mode:'copy'
    publishDir "$merge_path/$sample_id", pattern: '*.unassembled.reverse.fastq.gz', saveAs: {"${sample_id}_merge_R2.fastq.gz"}, mode:'copy', overwrite: false
    publishDir "$merge_path/$sample_id", pattern: '*_U.fastq.gz', mode:'copy', overwrite: false
    publishDir "$merge_path/$sample_id", pattern: '*.stats', mode:'copy', overwrite: false

    input:
    set val(sample_id), \
    file(r1_trimAdapt), \
    file(r2_trimAdapt), \
    file(u_trimAdapt) from trim_adapt_ch

    output:
    set val(sample_id), \
    file("${sample_id}.unassembled.forward.fastq.gz"), \
    file("${sample_id}.unassembled.reverse.fastq.gz"), \
    file("${sample_id}_merge_U.fastq.gz") into merge_ch
    set val(sample_id), file("${sample_id}_merge.stats") into merge_stats_ch

    script:
    """
    pear \
    -f $r1_trimAdapt \
    -r $r2_trimAdapt \
    -o ${sample_id} \
    -p ${params.pvalue} \
    -v ${params.min_overlap} \
    -n ${params.min_read_len} \
    -k \
    -j ${task.cpus}

    pigz -p ${task.cpus} -v -r ${sample_id}.unassembled.forward.fastq
    pigz -p ${task.cpus} -v -r ${sample_id}.unassembled.reverse.fastq
    pigz -p ${task.cpus} -v -r ${sample_id}.assembled.fastq

    cat $u_trimAdapt ${sample_id}.assembled.fastq.gz > ${sample_id}_merge_U.fastq.gz

    seqkit stats ${sample_id}.unassembled.forward.fastq.gz \
    ${sample_id}.unassembled.reverse.fastq.gz \
    ${sample_id}_merge_U.fastq.gz > ${sample_id}_merge.stats

    """
}

if (!params.skip_dedup){
    qual_path = "${params.outdir}/05.trim_qual/"
} else {
    qual_path = "${params.outdir}/04.trim_qual/"
}


process trim_qual {
    
/*
 * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
 */

    tag "Trimmomatic on $sample_id"

    publishDir "$qual_path/$sample_id", pattern: '*_{R1,R2,U}.fastq.gz', mode:'copy'
    publishDir "$qual_path/$sample_id", pattern: '*.stats', mode:'copy', overwrite: false

    input:
    set val(sample_id), file(r1_merge), file(r2_merge), file(u_merge) from merge_ch

    output:
    set val(sample_id), \
    file("${sample_id}_qual_R1.fastq.gz"), \
    file("${sample_id}_qual_R2.fastq.gz"), \
    file("${sample_id}_qual_U.fastq.gz") into qual_ch
    set val(sample_id), file("${sample_id}_qual.stats") into qual_stats_ch

    script:
    """
    trimmomatic PE \
    -threads ${task.cpus} \
    $r1_merge $r2_merge \
    ${sample_id}_qual_R1.fastq.gz \
    ${sample_id}_1_u.fastq.gz \
    ${sample_id}_qual_R2.fastq.gz \
    ${sample_id}_2_u.fastq.gz \
    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:${params.min_read_len}

    trimmomatic SE \
    -threads ${task.cpus} \
    $u_merge \
    ${sample_id}_unpaired.fastq.gz \
    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:${params.min_read_len}

    cat ${sample_id}_1_u.fastq.gz \
    ${sample_id}_2_u.fastq.gz \
    ${sample_id}_unpaired.fastq.gz > \
    ${sample_id}_qual_U.fastq.gz

    seqkit stats ${sample_id}_qual_R1.fastq.gz \
    ${sample_id}_qual_R2.fastq.gz \
    ${sample_id}_qual_U.fastq.gz > ${sample_id}_qual.stats
    """
}

if (!params.skip_dedup){
    complex_path = "${params.outdir}/06.low_complex/"
} else {
    complex_path = "${params.outdir}/05.low_complex/"
}

process low_complex {
    
/*
 * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
 */

    tag "Remove low-complexity reads from $sample_id"

    publishDir "$complex_path/$sample_id", pattern: '*_{R1,R2,U}.fastq.gz', mode:'copy'
    publishDir "$complex_path/$sample_id", pattern: '*.stats', mode:'copy', overwrite: false

    input:
    set val(sample_id), file(r1_qual), file(r2_qual), file(u_qual) from qual_ch

    output:
    set val(sample_id), \
    file("${sample_id}_R1.fastq.gz"), \
    file("${sample_id}_R2.fastq.gz"), \
    file("${sample_id}_U.fastq.gz") into low_complex_ch
    set val(sample_id), file("${sample_id}_complex.stats") into low_complex_stats_ch

    when: 
    !params.skip_low_complex

    script:
    """
    remove_low_complex.py \
    -1 $r1_qual \
    -2 $r2_qual\
    -u $u_qual\
    -c ${params.lc_cut_off_val}\
    -p ${sample_id}

    pigz -p ${task.cpus} -v -r ${sample_id}_R1.fastq
    pigz -p ${task.cpus} -v -r ${sample_id}_R2.fastq
    pigz -p ${task.cpus} -v -r ${sample_id}_U.fastq

    seqkit stats ${sample_id}_R1.fastq.gz \
    ${sample_id}_R2.fastq.gz \
    ${sample_id}_U.fastq.gz > ${sample_id}_complex.stats
    """
}

if (!params.skip_dedup){

    process summarize_stats_withDedup {
    
    /*
     * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
     */

        tag "Gather all read stats for $sample_id"

        publishDir "${params.outdir}/07.polishing_stats/", mode:'copy'

        input:
        set val(sample_id), \
        file(fqc_stats_fn), \
        file(dedup_stats_fn), \
        file(adapt_stats_fn), \
        file(merge_stats_fn), \
        file(qual_stats_fn), \
        file(complex_stats_fn) from fastqc_stats_ch.join(dedup_stats_ch).join(trim_adapt_stats_ch).join(merge_stats_ch).join(qual_stats_ch).join(low_complex_stats_ch)

        output:
        file("${sample_id}_polishing_table.tsv")
        file("${sample_id}_polishing_overview.png")

        script:
        """
        echo fqc_stats_fn is $fqc_stats_fn
        echo dedup_stats_fn is $dedup_stats_fn
        echo adapt_stats_fn is $adapt_stats_fn
        echo merge_stats_fn is $merge_stats_fn
        echo qual_stats_fn is $qual_stats_fn
        echo complex_stats_fn is $complex_stats_fn

        parse_visualize_stats.py \
        -f $fqc_stats_fn \
        -d $dedup_stats_fn \
        -a $adapt_stats_fn \
        -m $merge_stats_fn \
        -q $qual_stats_fn \
        -c $complex_stats_fn \
        -p ${sample_id}

        """
    }

} else {

    process summarize_stats {
    
    /*
     * First specify the SLURM reqs for MfN (this should later be replaced with a config file and label)
     */

        tag "Gather all read stats for $sample_id"

        publishDir "${params.outdir}/06.polishing_stats/", mode:'copy'

        input:
        set val(sample_id), \
        file(fqc_stats_fn), \
        file(adapt_stats_fn), \
        file(merge_stats_fn), \
        file(qual_stats_fn), \
        file(complex_stats_fn) from fastqc_stats_ch.join(trim_adapt_stats_ch).join(merge_stats_ch).join(qual_stats_ch).join(low_complex_stats_ch)

        output:
        file("${sample_id}_polishing_table.tsv")
        file("${sample_id}_polishing_overview.png")

        script:
        """
        parse_visualize_stats.py \
        -f $fqc_stats_fn \
        -a $adapt_stats_fn \
        -m $merge_stats_fn \
        -q $qual_stats_fn \
        -c $complex_stats_fn \
        -p ${sample_id}

        """
    }
}

