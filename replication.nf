#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.security.MessageDigest

def quoted(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
def sha256(path) {
    def digest=MessageDigest.getInstance('SHA-256')
    path.toFile().withInputStream { stream ->
        byte[] buffer=new byte[1048576]
        int size
        while((size=stream.read(buffer))!=-1) digest.update(buffer,0,size)
    }
    digest.digest().encodeHex().toString()
}
def localPath(value) {
    def path=java.nio.file.Paths.get(value.toString())
    if(path.isAbsolute() || path.normalize().toString()!=value.toString() || value.toString().contains('..'))
        error "Expected normalized project-relative path: ${value}"
    def full=projectDir.resolve(path).toRealPath()
    if(!full.startsWith(projectDir.toRealPath())) error "Path escapes project: ${value}"
    full
}
def verify(meta) {
    def manifestPath=localPath(meta.evaluation_freeze)
    def frozen=new JsonSlurper().parse(manifestPath.toFile())
    def phase=meta.synthetic ? 'SYNTHETIC_REPLICATION_TEST_ONLY' : 'REPLICATION_IMPLEMENTATION_FROZEN_PRE_EXTERNAL'
    if(frozen.phase!=phase || !(frozen.artifacts instanceof Map) || frozen.artifacts.isEmpty())
        error 'Replication freeze phase/schema mismatch'
    def required=[meta.counts,meta.design,meta.clinical,meta.stat_config,
        meta.pre_effect_freeze,meta.discovery_freeze,
        'replication.nf','replication.config','scripts/fit_donor_de.R',
        'scripts/donor_de_lib.R','scripts/evaluate_replication.py',
        'scripts/statistical_runtime_env.sh','requirements/statistical-runtime.explicit.txt']
    required.each { if(!frozen.artifacts.containsKey(it)) error "Freeze omitted ${it}" }
    def paths=frozen.artifacts.keySet().sort().collect { key ->
        def path=localPath(key)
        if(sha256(path)!=frozen.artifacts[key]) error "Frozen artifact changed: ${key}"
        path
    }
    if(meta.synthetic) {
        [meta.counts,meta.design,meta.clinical,meta.pre_effect_freeze,meta.discovery_freeze].each {
            if(!it.startsWith('work/synthetic/')) error 'Synthetic inputs must be isolated under work/synthetic/'
        }
    } else {
        if(meta.counts!='data/processed/heca/counts_stable.tsv.gz' ||
           meta.design!='reports/final_clinical_design_heca_primary.tsv' ||
           meta.clinical!='manifests/final_heca_clinical_selected.tsv' ||
           meta.pre_effect_freeze!='receipts/statistical_pre_effect_freeze.json' ||
           meta.discovery_freeze!='receipts/discovery_rank_freeze.json')
            error 'Production replication must use the locked HECA and discovery inputs'
    }
    [manifestPath,paths,sha256(manifestPath)]
}
def guard(meta, expected) {
    """
    python - <<'PY'
import hashlib,json,pathlib
def digest(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
manifest=${JsonOutput.toJson(meta.evaluation_freeze)}
if digest(manifest)!=${JsonOutput.toJson(expected)}:raise SystemExit('Replication manifest changed')
for path,checksum in json.loads(pathlib.Path(manifest).read_text())['artifacts'].items():
    if digest(path)!=checksum:raise SystemExit('Frozen artifact changed: '+path)
PY
    """
}
process HECA_DONOR_DE {
    tag 'heca'
    publishDir "${params.outdir}",mode:'copy',overwrite:false,failOnError:true
    input:
    val meta
    path frozen_manifest
    path frozen_artifacts,stageAs:'frozen_inputs/artifact????'
    val frozen_sha
    output:
    path 'de',emit:de
    script:
    def runtime=params.runtime_mode=='native' ? "source ${quoted(projectDir.resolve('scripts/statistical_runtime_env.sh'))}" : ''
    def synthetic=meta.synthetic ? '--synthetic-test true' : ''
    """
    ${runtime}
    output="\$PWD/de"
    cd ${quoted(projectDir)}
    ${guard(meta,frozen_sha)}
    Rscript scripts/fit_donor_de.R --counts ${quoted(meta.counts)} --design ${quoted(meta.design)} \
      --clinical ${quoted(meta.clinical)} --cohort heca --config ${quoted(meta.stat_config)} \
      --pre-effect-freeze ${quoted(meta.pre_effect_freeze)} --discovery-freeze ${quoted(meta.discovery_freeze)} \
      --output-dir "\$output" ${synthetic}
    """
}
process EVALUATE_REPLICATION {
    tag 'fixed_lists_and_nulls'
    publishDir "${params.outdir}",mode:'copy',overwrite:false,failOnError:true
    input:
    val meta
    path frozen_manifest
    path frozen_artifacts,stageAs:'frozen_inputs/artifact????'
    val frozen_sha
    path de
    output:
    path 'evaluation',emit:evaluation
    script:
    def runtime=params.runtime_mode=='native' ? "source ${quoted(projectDir.resolve('scripts/statistical_runtime_env.sh'))}" : ''
    def synthetic=meta.synthetic ? '--synthetic-test' : ''
    """
    ${runtime}
    de_input="\$PWD/${de}"
    output="\$PWD/evaluation"
    cd ${quoted(projectDir)}
    ${guard(meta,frozen_sha)}
    python scripts/evaluate_replication.py --discovery-freeze ${quoted(meta.discovery_freeze)} \
      --external-dir "\$de_input" --output-dir "\$output" ${synthetic}
    """
}
workflow {
    if(!params.evaluation_freeze || !params.outdir) error '--evaluation_freeze and --outdir required'
    if(!(params.synthetic_test instanceof Boolean)) error '--synthetic_test must be boolean'
    def out=java.nio.file.Paths.get(params.outdir.toString()).toAbsolutePath().normalize()
    if(!out.startsWith(projectDir.toAbsolutePath().normalize()) || out==projectDir.toAbsolutePath().normalize())
        error 'Output must stay within project'
    if(java.nio.file.Files.exists(out)) error 'Output must be new; preserve previous runs'
    def work=java.nio.file.Paths.get(workflow.workDir.toString()).toAbsolutePath().normalize()
    if(!work.startsWith(projectDir.toAbsolutePath().normalize())) error 'Work directory must stay within project'
    if(params.runtime_mode=='docker' && !(params.container_image==~ /.+@sha256:[0-9a-f]{64}/))
        error 'Docker execution requires an immutable image digest'
    def meta=[counts:params.counts,design:params.design,clinical:params.clinical,
        stat_config:params.stat_config,pre_effect_freeze:params.pre_effect_freeze,
        discovery_freeze:params.discovery_freeze,evaluation_freeze:params.evaluation_freeze,
        synthetic:params.synthetic_test]
    def verified=verify(meta)
    def gate=Channel.value(meta)
    def manifest=Channel.value(verified[0])
    def artifacts=Channel.value(verified[1])
    def digest=Channel.value(verified[2])
    HECA_DONOR_DE(gate,manifest,artifacts,digest)
    EVALUATE_REPLICATION(gate,manifest,artifacts,digest,HECA_DONOR_DE.out.de)
}
