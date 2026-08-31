nextflow.enable.dsl=2

include { HARMONIZE } from './modules/local/harmonize'
include { ESTIMATE  } from './modules/local/estimate'
include { COLOC     } from './modules/local/coloc'

if( !params.containsKey('mode') )          params.mode = 'pQTL'
if( !params.containsKey('run_harmonize') ) params.run_harmonize = true
if( !params.containsKey('run_estimate') )  params.run_estimate = true
if( !params.containsKey('run_coloc') )     params.run_coloc = true
if( !params.containsKey('n_is') )          params.n_is = 3
if( !params.containsKey('n_chunks') )      params.n_chunks = 10

workflow {

    def studies = getStudies(params.mode)
    def nchunks = params.n_chunks as Integer

    println "Studies: ${studies}"
    println "nchunks: ${nchunks}"
    println "run_harmonize=${params.run_harmonize} run_estimate=${params.run_estimate} run_coloc=${params.run_coloc}"

    studies_ch = Channel.fromList(studies)

    estimate_jobs_ch = Channel
        .fromList(studies)
        .combine(Channel.from(1..nchunks))
        .combine(Channel.from(1..params.n_is))
        .map { study, chunk, is_idx -> tuple(chunk, study, is_idx) }

    coloc_jobs_ch = Channel
        .fromList(studies)
        .combine(Channel.from(1..nchunks))
        .map { study, chunk -> tuple(chunk, study) }

    if (params.run_harmonize) {
        HARMONIZE(studies_ch)
    }

    if (params.run_estimate) {
        ESTIMATE(estimate_jobs_ch)
    }

    if (params.run_coloc) {
        COLOC(coloc_jobs_ch)
    }
}

def getStudies(mode) {
    def rcmd =
        'source("~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R"); ' +
        'cfg <- load_mr_config("' + mode + '"); ' +
        'cat(paste(names(cfg$gwas), collapse="\\n"))'

    def cmd = ['bash', '-lc', 'Rscript -e ' + shellQuote(rcmd)]
    def proc = cmd.execute()
    proc.waitFor()

    if (proc.exitValue() != 0) {
        throw new RuntimeException("Failed to read studies for mode=${mode}: ${proc.err.text}")
    }

    return proc.in.text.readLines().findAll { it?.trim() }
}

def shellQuote(s) {
    return "'" + s.replace("'", "'\"'\"'") + "'"
}