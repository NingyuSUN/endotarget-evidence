#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.security.MessageDigest

def quoted(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
def sha256(path) {
    def digest = MessageDigest.getInstance('SHA-256')
    path.toFile().withInputStream { stream ->
        byte[] buffer = new byte[1048576]
        int size
        while ((size = stream.read(buffer)) != -1) digest.update(buffer, 0, size)
    }
    digest.digest().encodeHex().toString()
}
def projectPath(value) {
    def path = java.nio.file.Paths.get(value.toString())
    if (path.isAbsolute() || path.normalize().toString() != value.toString() || value.toString().contains('..'))
        error "Inputs must use normalized project-relative paths: ${value}"
    def full = projectDir.resolve(path).toRealPath()
    if (!full.startsWith(projectDir.toRealPath())) error "Input escapes project root: ${value}"
    full
}
def verifyManifest(meta) {
    def path = projectPath(meta.freeze)
    def frozen = new JsonSlurper().parse(path.toFile())
    def phase = meta.synthetic ? 'SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE' : 'STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT'
    if (frozen.phase != phase || !(frozen.artifacts instanceof Map) || frozen.artifacts.isEmpty())
        error "Freeze phase/schema does not authorize this run"
    def required = [
        meta.counts, meta.design, meta.clinical, meta.edges, meta.stat_config, meta.rank_config,
        'main.nf', 'nextflow.config', 'scripts/fit_donor_de.R', 'scripts/donor_de_lib.R',
        'scripts/rank_targets.py', 'scripts/sample_matched_nulls.py', 'scripts/statistical_runtime_env.sh',
        'requirements/statistical-runtime.explicit.txt', 'requirements/statistical-runtime.yml',
        'receipts/runtime_statistics_install.json', 'receipts/runtime_statistics_R_smoke.json'
    ]
    required.each { item ->
        if (!frozen.artifacts.containsKey(item)) error "Required input/code/runtime omitted from freeze: ${item}"
    }
    def paths = frozen.artifacts.keySet().sort().collect { key ->
        if (!(frozen.artifacts[key] ==~ /[0-9a-f]{64}/)) error "Invalid SHA256 in freeze: ${key}"
        def artifact = projectPath(key)
        if (sha256(artifact) != frozen.artifacts[key]) error "Frozen artifact changed: ${key}"
        artifact
    }
    if (meta.synthetic) {
        [meta.counts, meta.design, meta.clinical, meta.edges].each {
            if (!it.startsWith('work/synthetic/')) error "Synthetic inputs must be under work/synthetic/"
        }
    } else {
        if (meta.counts != 'data/processed/discovery/counts_stable.tsv.gz' ||
            meta.design != 'reports/final_clinical_design_discovery_primary.tsv' ||
            meta.clinical != 'manifests/final_discovery_clinical_selected.tsv' ||
            meta.edges != 'data/reference/string/direct_evidence_gene_edges.tsv.gz')
            error "Production discovery requires the adopted primary inputs"
    }
    [path, paths, sha256(path)]
}

// Recheck original project paths inside each task, after switching to project cwd.
// Deep path inputs + the manifest digest prevent stale -resume cache reuse.
def taskFreezeGuard(meta, expected_sha) {
    def literal = JsonOutput.toJson(meta.freeze)
    def expected = JsonOutput.toJson(expected_sha)
    """
    python - <<'PY'
import hashlib,json,pathlib
def digest(path):
    value=hashlib.sha256()
    with pathlib.Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1048576),b''):
            value.update(block)
    return value.hexdigest()
freeze_path=${literal}
if digest(freeze_path)!=${expected}:
    raise SystemExit('Freeze manifest changed after workflow preflight')
frozen=json.loads(pathlib.Path(freeze_path).read_text())
for name,expected in frozen['artifacts'].items():
    if digest(name)!=expected:
        raise SystemExit('Frozen artifact changed inside task: '+name)
PY
    """
}

process DONOR_DE {
    tag 'discovery'
    publishDir "${params.outdir}", mode: 'copy', overwrite: false, failOnError: true
    input:
    val meta
    path freeze_manifest
    path frozen_inputs, stageAs: 'frozen_inputs/artifact????'
    val freeze_sha
    output:
    path 'de', emit: de
    script:
    def runtime = params.runtime_mode == 'native' ? "source ${quoted(projectDir.resolve('scripts/statistical_runtime_env.sh'))}" : ''
    def testflag = meta.synthetic ? '--synthetic-test true' : ''
    """
    ${runtime}
    output="\$PWD/de"
    cd ${quoted(projectDir)}
    ${taskFreezeGuard(meta, freeze_sha)}
    Rscript scripts/fit_donor_de.R --counts ${quoted(meta.counts)} --design ${quoted(meta.design)} \
      --clinical ${quoted(meta.clinical)} --cohort discovery --config ${quoted(meta.stat_config)} \
      --pre-effect-freeze ${quoted(meta.freeze)} --output-dir "\$output" ${testflag}
    """
}

process RANK_TARGETS {
    tag 'B0_B1_G1'
    publishDir "${params.outdir}", mode: 'copy', overwrite: false, failOnError: true
    input:
    val meta
    path freeze_manifest
    path frozen_inputs, stageAs: 'frozen_inputs/artifact????'
    val freeze_sha
    path de
    output:
    path 'ranking', emit: ranking
    script:
    def runtime = params.runtime_mode == 'native' ? "source ${quoted(projectDir.resolve('scripts/statistical_runtime_env.sh'))}" : ''
    def testflag = meta.synthetic ? '--synthetic-test' : ''
    """
    ${runtime}
    de_input="\$PWD/${de}"
    output="\$PWD/ranking"
    cd ${quoted(projectDir)}
    ${taskFreezeGuard(meta, freeze_sha)}
    python scripts/rank_targets.py --de "\$de_input/de_fitted.tsv" --de-all "\$de_input/de_all.tsv" \
      --edges ${quoted(meta.edges)} --config ${quoted(meta.rank_config)} --output-dir "\$output" ${testflag}
    """
}

process MATCHED_NULLS {
    tag 'matched_nulls'
    publishDir "${params.outdir}", mode: 'copy', overwrite: false, failOnError: true
    input:
    val meta
    path freeze_manifest
    path frozen_inputs, stageAs: 'frozen_inputs/artifact????'
    val freeze_sha
    path de
    path ranking
    output:
    path 'nulls', emit: nulls
    script:
    def runtime = params.runtime_mode == 'native' ? "source ${quoted(projectDir.resolve('scripts/statistical_runtime_env.sh'))}" : ''
    def testflag = meta.synthetic ? '--synthetic-test' : ''
    """
    ${runtime}
    de_input="\$PWD/${de}"
    rank_input="\$PWD/${ranking}"
    output="\$PWD/nulls"
    cd ${quoted(projectDir)}
    ${taskFreezeGuard(meta, freeze_sha)}
    python scripts/sample_matched_nulls.py --de "\$de_input/de_fitted.tsv" --ranking-dir "\$rank_input" \
      --config ${quoted(meta.stat_config)} --output-dir "\$output" ${testflag}
    """
}

workflow {
    if (!params.pre_effect_freeze) error '--pre_effect_freeze is required; no disease effects are authorized by defaults'
    if (!params.outdir) error '--outdir is required; choose a new empty project directory'
    if (!(params.synthetic_test instanceof Boolean)) error '--synthetic_test must be a boolean'
    if (params.runtime_mode == 'docker' && !(params.container_image ==~ /.+@sha256:[0-9a-f]{64}/))
        error 'Docker requires an immutable --container_image IMAGE@sha256:DIGEST'
    def out = java.nio.file.Paths.get(params.outdir.toString()).toAbsolutePath().normalize()
    if (!out.startsWith(projectDir.toAbsolutePath().normalize()) || out == projectDir.toAbsolutePath().normalize())
        error '--outdir must be inside this project and cannot be the project root'
    if (java.nio.file.Files.exists(out)) error '--outdir must not exist; preserve prior outputs, including on -resume'
    def work = java.nio.file.Paths.get(workflow.workDir.toString()).toAbsolutePath().normalize()
    if (!work.startsWith(projectDir.toAbsolutePath().normalize())) error 'Nextflow workDir must stay inside this project'
    def meta = [counts:params.counts, design:params.design, clinical:params.clinical, edges:params.edges,
                stat_config:params.stat_config, rank_config:params.rank_config,
                freeze:params.pre_effect_freeze, synthetic:params.synthetic_test]
    def verified = verifyManifest(meta)
    def gate = Channel.value(meta)
    def manifest = Channel.value(verified[0])
    def bundle = Channel.value(verified[1])
    def digest = Channel.value(verified[2])
    DONOR_DE(gate, manifest, bundle, digest)
    RANK_TARGETS(gate, manifest, bundle, digest, DONOR_DE.out.de)
    MATCHED_NULLS(gate, manifest, bundle, digest, DONOR_DE.out.de, RANK_TARGETS.out.ranking)
}
