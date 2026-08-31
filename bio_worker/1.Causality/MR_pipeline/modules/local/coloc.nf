process COLOC {
    tag { "${study}_chunk${chunk}" }
    publishDir "${params.outdir}/logs/coloc", mode: 'copy', pattern: '*.log'
    container params.container

    input:
    tuple val(chunk), val(study)

    output:
    path "coloc_${study}_chunk${chunk}.done", emit: done

    script:
    def rcmd =
        'source("~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R"); ' +
        'cfg <- load_mr_config("' + params.mode + '"); ' +
        'cat(match("' + study + '", names(cfg$gwas)))'

    def study_index_cmd = 'Rscript -e ' + "'" + rcmd.replace("'", "'\"'\"'") + "'"

    """
    log_file=coloc_${study}_chunk${chunk}.log
    : > \$log_file
    study_index=\$(${study_index_cmd})
    Rscript ${projectDir}/bin/06_colocalization.R ${chunk} \$study_index ${params.mode} ${params.n_chunks} &>> \$log_file
    touch coloc_${study}_chunk${chunk}.done
    """
}