process GATHER {
    tag { mode }
    publishDir "${params.outdir}/logs/gather", mode: 'copy', pattern: '*.log'
    container params.container

    input:
    val trigger
    val mode

    output:
    path "gather_${mode}.done", emit: done

    script:
    """
    log_file=gather_${mode}.log
    : > \$log_file
    Rscript ${projectDir}/bin/07_gatherResults.R ${mode} &>> \$log_file
    touch gather_${mode}.done
    """
}