process HARMONIZE {
    tag { study }
    publishDir "${params.outdir}/logs/harmonize", mode: 'copy', pattern: '*.log'
    container params.container

    input:
    val study

    output:
    path "harmonize_${study}.done", emit: done

    script:
    """
    log_file=harmonize_${study}.log
    : > \$log_file
    Rscript ${projectDir}/bin/04_harmonization.R ${params.mode} ${study} &>> \$log_file
    touch harmonize_${study}.done
    """
}