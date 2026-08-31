process ESTIMATE {
    tag { "${study}_chunk${chunk}_IS${is_idx}" }
    publishDir "${params.outdir}/logs/estimate", mode: 'copy', pattern: '*.log'
    container params.container

    input:
    tuple val(chunk), val(study), val(is_idx)

    output:
    path "estimate_${study}_chunk${chunk}_IS${is_idx}.done", emit: done

    script:
    def rcmd =
        'source("~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R"); ' +
        'cfg <- load_mr_config("' + params.mode + '"); ' +
        'cat(match("' + study + '", names(cfg$gwas)))'

    def study_index_cmd = 'Rscript -e ' + "'" + rcmd.replace("'", "'\"'\"'") + "'"

    """
    log_file=estimate_${study}_chunk${chunk}_IS${is_idx}.log
    : > \$log_file
    study_index=\$(${study_index_cmd})
    Rscript ${projectDir}/bin/05_estimation.R ${chunk} \$study_index ${is_idx} ${params.mode} ${params.n_chunks} &>> \$log_file
    touch estimate_${study}_chunk${chunk}_IS${is_idx}.done
    """
}