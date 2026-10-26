#!/usr/bin/env bash
# Phase positive-control feasibility check (protocol step 1; metadata only).
#
# Reads the clinical and cell-label tables already on disk, applies the frozen donor
# rules in docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md, writes candidate ledgers and
# numeric designs in the v1 design-file format, checks rank and residual degrees of
# freedom, and reads archive sizes from the public servers with HEAD requests.
# It reads no expression values, downloads no data and fits no model.
# Requires bash, gawk, sha256sum, curl, git, gzip. Optional operator-reported
# environment: HOST_RAM_GIB, WSL_MEMORY_GB, WSL_SWAP_GB (recorded, not measured).
set -euo pipefail
export LC_ALL=C
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CFG=config/phase_positive_control.json
D_V1=manifests/final_discovery_clinical_selected.tsv
D_ALL=manifests/discovery_clinical_candidates.tsv
H_V1=manifests/final_heca_clinical_selected.tsv
H_OBS=data/metadata/heca_h5/obs_selected.tsv
H_DON=data/metadata/heca_author/metadata-donors.csv
H_LIB=data/metadata/heca_author/metadata-libraries.csv
D_COUNTS=data/processed/discovery/counts_stable.tsv.gz
H_PROV=data/processed/heca/provenance.json

D_LEDGER=manifests/phase_positive_control_discovery_candidates.tsv
H_LEDGER=manifests/phase_positive_control_heca_candidates.tsv
D_DESIGN=reports/phase_positive_control_design_discovery_primary.tsv
D_DESIGN_FB=reports/phase_positive_control_design_discovery_fallback.tsv
H_DESIGN=reports/phase_positive_control_design_heca_primary.tsv
H_DESIGN_FB=reports/phase_positive_control_design_heca_fallback.tsv
DL=manifests/phase_positive_control_heca_download_manifest.tsv
RECEIPT=receipts/phase_positive_control_feasibility.json

GEO_RDS_URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266265/suppl/GSE266265_EndoAtlas.rds
GEO_RDS_RECORDED_BYTES=41729517596
BIOSTUDIES=https://www.ebi.ac.uk/biostudies/files/E-MTAB-14039
HEAD_TIMEOUT=${HEAD_TIMEOUT:-60}
MIN_CELLS=50
MIN_DF=5
MIN_D=8
MIN_H=5

for f in "$CFG" "$D_V1" "$D_ALL" "$H_V1" "$H_OBS" "$H_DON" "$H_LIB" "$D_COUNTS" "$H_PROV"; do
  [ -s "$f" ] || { echo "missing input: $f" >&2; exit 2; }
done
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

sha() { sha256sum "$1" | cut -d' ' -f1; }
ids_in_config() { grep -o -E "\"$1\"([^:]|\$)" "$CFG" | grep -o -E "\"$1\"" | tr -d '"' | sort -u; }

# ---- 1. Discovery candidate ledger ---------------------------------------------
{
  printf 'donor_id\tphase\tsub_phase\tdisease\tage\tcapture_batch\tsequencing_batch\tadenomyosis\thormone_original\tminor_exclusion\tdisease_inconsistency\tprimary_lineage_cells\tcounts_available\tsource\n'
  gawk -F'\t' 'NR>1 {printf "%s\tproliferative\tproliferative\t%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\t%s\tyes\tv1_primary_ledger\n", $2,$4,$5,$6,$7,$8,$14,$16,$19}' "$D_V1"
  gawk -F'\t' 'NR>1 && $17 ~ /^secretory_/ && $13=="FALSE" && $20=="0" && ($25=="" || $25=="False") && $24 ~ /^[0-9.]+$/ {printf "%s\tsecretory\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tno\tmetadata_candidate_rule\n", $1,$17,$6,$24,$18,$19,$11,$13,$20,$25}' "$D_ALL"
} > "$D_LEDGER"
D_EXCLUDED=$(gawk -F'\t' 'NR>1 && $17 ~ /^secretory_/ && !($13=="FALSE" && $20=="0" && ($25=="" || $25=="False") && $24 ~ /^[0-9.]+$/) {printf "%s%s[%s;minor=%s;detail=%s]", (n++?";":""), $1, $17, $20, $21}' "$D_ALL")
D_PERIOVULATORY=$(gawk -F'\t' 'NR>1 && $17=="periovulatory"' "$D_ALL" | wc -l | tr -d ' ')
D_N=$(gawk -F'\t' 'NR>1' "$D_LEDGER" | wc -l | tr -d ' ')
D_SEC=$(gawk -F'\t' 'NR>1 && $2=="secretory"' "$D_LEDGER" | wc -l | tr -d ' ')
D_PRO=$(gawk -F'\t' 'NR>1 && $2=="proliferative"' "$D_LEDGER" | wc -l | tr -d ' ')
D_SEC_ENDO=$(gawk -F'\t' 'NR>1 && $2=="secretory" && $4=="ENDO"' "$D_LEDGER" | wc -l | tr -d ' ')
D_CFG_MATCH=false
if [ "$(ids_in_config 'sample[0-9]+')" = "$(gawk -F'\t' 'NR>1{print $1}' "$D_LEDGER" | sort -u)" ]; then D_CFG_MATCH=true; fi

# ---- 2. HECA candidate ledger ------------------------------------------------------
gawk -F'\t' -v OFS='\t' 'BEGIN{split("eStromal eStromal_MMPs eStromal_cycling dStromal_early dStromal_mid dStromal_late Fibroblast_basalis",a," "); for(i in a)P[a[i]]=1}
NR>1 { s=$1; if(!(s in g)){g[s]=$3; grp[s]=$5; age[s]=$6; horm[s]=$7; bs[s]=$8; st[s]=$9; bio[s]=$10; tis[s]=$11} tot[s]++; if($19=="no" && ($14 in P)) p[s]++ }
END{ for(s in g) print s, g[s], grp[s], age[s], horm[s], bs[s], st[s], bio[s], tis[s], tot[s], p[s]+0 }' "$H_OBS" | sort > "$TMP/heca_samples.tsv"
H_OBS_ROWS=$(gawk 'END{print NR-1}' "$H_OBS")

gawk -F'\t' -v OFS='\t' -v LIB="$H_LIB" -v DON="$H_DON" -v MIN="$MIN_CELLS" '
BEGIN{
  CSV="([^,]*)|(\"[^\"]+\")"
  while((getline line < LIB)>0){ n=patsplit(line,f,CSV); if(f[1]=="FolderName"){for(i=1;i<=n;i++)lh[f[i]]=i; continue} lane[f[1]]=f[lh["Lane"]]; qc[f[1]]=f[lh["QC_single sample analysis"]]; cr[f[1]]=f[lh["QC_CellRanger"]] }
  close(LIB)
  while((getline line < DON)>0){ n=patsplit(line,f,CSV); if(f[1]=="Donor"){for(i=1;i<=n;i++)dh[f[i]]=i; continue} d=f[1]; if(d=="") continue; cnt[d]++; path[d]=f[dh["Endometrial pathology_A-F-H-P"]]; dage[d]=f[dh["Age"]] }
  close(DON)
}
{ s=$1; d=$2; grp=$3; age=$4; horm=$5; bs=$6; st=$7; bio=$8; tis=$9; prim=$11
  if(!(bs=="Secretory" && st!="Hormones" && st!="Menstrual" && st!="Unknown" && horm=="nan" && bio=="Superficial_Eutopic" && tis=="eutopic_endometrium")) next
  lib=s; sub(/_[^_]+$/,"",lib)
  row[s]=d OFS lib OFS grp OFS age OFS st OFS prim OFS horm
  if(prim+0<MIN){ reason[s]="lt" MIN "_primary_nuclei"; next }
  if(cnt[d]>1){ reason[s]="duplicate_clinical_rows"; next }
  if(!(d in best) || prim+0>bp[d] || (prim+0==bp[d] && s<best[d])){ if(d in best) reason[best[d]]="repeat_smaller_or_later"; best[d]=s; bp[d]=prim+0 } else reason[s]="repeat_smaller_or_later"
}
END{ for(s in row){ split(row[s],r,OFS); d=r[1]; lib=r[2]; grp=r[3]; age=r[4]; st=r[5]; prim=r[6]; horm=r[7]
  dis=(grp=="Endo_Superficial")?"ENDO":"non-ENDO"
  sel=(best[d]==s)?"True":"False"; rs=(sel=="True")?"":reason[s]
  tok=path[d]; gsub(/"/,"",tok); ad=(tok ~ /(^|[^A-Za-z])A([^A-Za-z]|$)/)?1:0
  agec=(dage[d]+0==age+0)?"":"age_mismatch_obs_vs_clinical"
  if(agec!="") rs=(rs==""?agec:rs ";" agec)
  q=qc[lib]; gsub(/"/,"",q); c=cr[lib]; gsub(/"/,"",c)
  print d, s, lib, "secretory", st, dis, age+0, lane[lib], ad, tok, prim, horm, q, c, "no", sel, rs, "metadata_candidate_rule"
}}' "$TMP/heca_samples.tsv" | sort -t$'\t' -k16,16r -k6,6 -k1,1 > "$TMP/heca_secretory.tsv"

{
  printf 'donor_id\tsample_id\tlibrary_id\tphase\tstage_fine\tdisease\tage\tlibrary_group_Lane\tadenomyosis\tpathology_original\tprimary_nuclei\thormone_h5_original\tlibrary_QC_original\tlibrary_CellRanger_QC_original\tcounts_available\tselected_primary\texclusion_reason\tsource\n'
  gawk -F'\t' -v OFS='\t' 'NR>1 {print $1,$2,$3,"proliferative",$8,$4,$6+0,$28,$20,$19,$37,$12,$32,$33,"yes","True","","v1_primary_ledger"}' "$H_V1"
  cat "$TMP/heca_secretory.tsv"
} > "$H_LEDGER"
H_N=$(gawk -F'\t' 'NR>1 && $16=="True"' "$H_LEDGER" | wc -l | tr -d ' ')
H_SEC=$(gawk -F'\t' 'NR>1 && $16=="True" && $4=="secretory"' "$H_LEDGER" | wc -l | tr -d ' ')
H_PRO=$(gawk -F'\t' 'NR>1 && $16=="True" && $4=="proliferative"' "$H_LEDGER" | wc -l | tr -d ' ')
H_SEC_ENDO=$(gawk -F'\t' 'NR>1 && $16=="True" && $4=="secretory" && $6=="ENDO"' "$H_LEDGER" | wc -l | tr -d ' ')
H_SEC_SAMPLES=$(gawk -F'\t' 'NR>1 && $4=="secretory"' "$H_LEDGER" | wc -l | tr -d ' ')
H_NOT_SELECTED=$(gawk -F'\t' 'NR>1 && $4=="secretory" && $16=="False" {printf "%s%s:%s", (n++?";":""), $2, $17}' "$H_LEDGER")
H_AGE_MISMATCH=$(gawk -F'\t' 'NR>1 && $17 ~ /age_mismatch/' "$H_LEDGER" | wc -l | tr -d ' ')
H_CFG_MATCH=false
if [ "$(ids_in_config '(EX|FX)[0-9]{4}')" = "$(gawk -F'\t' 'NR>1 && $16=="True"{print $1}' "$H_LEDGER" | sort -u)" ]; then H_CFG_MATCH=true; fi

# ---- 3. Numeric designs and rank checks ------------------------------------------------
# build_design LEDGER OUT id_col phase_col disease_col age_col batch_col batch_prefix adeno_col sel_col
build_design() {
  gawk -F'\t' -v OFS='\t' -v idc="$3" -v pc="$4" -v dc="$5" -v ac="$6" -v bc="$7" -v bp="$8" -v adc="$9" -v sc="${10}" '
  NR>1 && (sc==0 || $sc=="True") { n++; id[n]=$idc; ph[n]=($pc=="secretory")?1:0; ds[n]=($dc=="ENDO")?1:0; ag[n]=$ac+0; sum+=$ac; ad[n]=$adc+0; if(bc>0){ b[n]=$bc; lv[$bc]=1 } }
  END{
    mean=sum/n; m=0
    for(l in lv){ m++; L[m]=l; k=l; gsub(/[^0-9]/,"",k); K[m]=(k==""?-1:k+0) }
    for(i=1;i<=m;i++) for(j=i+1;j<=m;j++) if(K[j]<K[i] || (K[j]==K[i] && L[j]<L[i])){t=L[i];L[i]=L[j];L[j]=t;t=K[i];K[i]=K[j];K[j]=t}
    ref=L[1]; adconst=1; for(i=2;i<=n;i++) if(ad[i]!=ad[1]) adconst=0
    hdr="donor_id" OFS "intercept" OFS "SEC_minus_PRO" OFS "ENDO_minus_non_ENDO" OFS "age_centered_years"
    for(k=2;k<=m;k++) hdr=hdr OFS bp ":" L[k] "_vs_" ref
    if(!adconst) hdr=hdr OFS "adenomyosis"
    print hdr
    for(i=1;i<=n;i++){ line=id[i] OFS "1.0" OFS sprintf("%.1f",ph[i]) OFS sprintf("%.1f",ds[i]) OFS sprintf("%.16g",ag[i]-mean)
      for(k=2;k<=m;k++) line=line OFS ((b[i]==L[k])?"1.0":"0.0")
      if(!adconst) line=line OFS sprintf("%.1f",ad[i]); print line }
    printf "%.16g\t%d\t%s\t%d\n", mean, m, (m>0?ref:"none"), adconst > "/dev/stderr"
  }' "$1" > "$2"
}
design_rank() {  # prints n, p, rank, residual_df; Gaussian elimination with partial pivoting
  gawk -F'\t' 'NR==1{p=NF-1; next} {n++; for(j=1;j<=p;j++) M[n,j]=$(j+1)+0}
  END{ mx=0; for(i=1;i<=n;i++) for(j=1;j<=p;j++){v=M[i,j]; if(v<0)v=-v; if(v>mx)mx=v} tol=1e-7*mx; r=0
    for(j=1;j<=p && r<n;j++){ piv=0; maxv=tol; for(i=r+1;i<=n;i++){v=M[i,j]; if(v<0)v=-v; if(v>maxv){maxv=v; piv=i}} if(!piv) continue; r++
      if(piv!=r) for(k=1;k<=p;k++){t=M[r,k]; M[r,k]=M[piv,k]; M[piv,k]=t}
      for(i=r+1;i<=n;i++){f=M[i,j]/M[r,j]; if(f!=0) for(k=j;k<=p;k++) M[i,k]-=f*M[r,k]} }
    printf "%d\t%d\t%d\t%d\n", n, p, r, n-r }' "$1"
}
crosstab() { gawk -F'\t' -v pc="$2" -v bc="$3" -v sc="$4" 'NR>1 && (sc==0 || $sc=="True") {c[$bc" "$pc]++} END{for(k in c) print k, c[k]}' "$1" | sort | tr '\n' ';'; }

build_design "$D_LEDGER" "$D_DESIGN"    1 2 4 5 6 capture_batch    8 0 2> "$TMP/d1.meta"
build_design "$D_LEDGER" "$D_DESIGN_FB" 1 2 4 5 7 sequencing_batch 8 0 2> "$TMP/d2.meta"
build_design "$H_LEDGER" "$H_DESIGN"    1 4 6 7 8 library_group_Lane 9 16 2> "$TMP/h1.meta"
build_design "$H_LEDGER" "$H_DESIGN_FB" 1 4 6 7 0 none 9 16 2> "$TMP/h2.meta"
read -r D1_N D1_P D1_R D1_DF < <(design_rank "$D_DESIGN")
read -r D2_N D2_P D2_R D2_DF < <(design_rank "$D_DESIGN_FB")
read -r H1_N H1_P H1_R H1_DF < <(design_rank "$H_DESIGN")
read -r H2_N H2_P H2_R H2_DF < <(design_rank "$H_DESIGN_FB")
D_XTAB=$(crosstab "$D_LEDGER" 2 6 0)
H_XTAB=$(crosstab "$H_LEDGER" 4 8 16)
decide() { # n p rank df min_group_pass -> level or NOT_ESTIMABLE
  local n=$1 p=$2 r=$3 df=$4 fbn=$5 fbp=$6 fbr=$7 fbdf=$8 gp=$9
  if [ "$gp" != "true" ]; then echo NOT_ESTIMABLE; return; fi
  if [ "$r" -eq "$p" ] && [ "$df" -ge "$MIN_DF" ]; then echo primary; return; fi
  if [ "$fbr" -eq "$fbp" ] && [ "$fbdf" -ge "$MIN_DF" ]; then echo fallback; return; fi
  echo NOT_ESTIMABLE
}
D_GROUPS_PASS=false; [ "$D_SEC" -ge "$MIN_D" ] && [ "$D_PRO" -ge "$MIN_D" ] && D_GROUPS_PASS=true
H_GROUPS_PASS=false; [ "$H_SEC" -ge "$MIN_H" ] && [ "$H_PRO" -ge "$MIN_H" ] && H_GROUPS_PASS=true
D_DECISION=$(decide "$D1_N" "$D1_P" "$D1_R" "$D1_DF" "$D2_N" "$D2_P" "$D2_R" "$D2_DF" "$D_GROUPS_PASS")
H_DECISION=$(decide "$H1_N" "$H1_P" "$H1_R" "$H1_DF" "$H2_N" "$H2_P" "$H2_R" "$H2_DF" "$H_GROUPS_PASS")

# ---- 4. Download manifests (HEAD requests only) ------------------------------------------
head_size() { curl -sIL --max-time "$HEAD_TIMEOUT" "$1" | tr -d '\r' | gawk 'tolower($1) ~ /^http\// {code=$2} tolower($1)=="content-length:" {len=$2} END{printf "%s\t%s\n", (code==""?"NA":code), (len==""?"NA":len)}'; }
{
  printf 'donor_id\tsample_id\turl\thttp_status\tbytes\n'
  gawk -F'\t' 'NR>1 && $4=="secretory" && $16=="True" {print $1"\t"$2}' "$H_LEDGER" | sort | while IFS=$'\t' read -r d s; do
    u="$BIOSTUDIES/$s.tar.gz"; r=$(head_size "$u"); printf '%s\t%s\t%s\t%s\n' "$d" "$s" "$u" "$r"
  done
} > "$DL"
DL_N=$(gawk -F'\t' 'NR>1' "$DL" | wc -l | tr -d ' ')
DL_OK=$(gawk -F'\t' 'NR>1 && $4=="200" && $5 ~ /^[0-9]+$/' "$DL" | wc -l | tr -d ' ')
DL_BYTES=$(gawk -F'\t' 'NR>1 && $5 ~ /^[0-9]+$/ {s+=$5} END{printf "%d", s}' "$DL")
read -r GEO_STATUS GEO_BYTES < <(head_size "$GEO_RDS_URL")
GEO_MATCH=false; [ "$GEO_BYTES" = "$GEO_RDS_RECORDED_BYTES" ] && GEO_MATCH=true
DISK_AVAIL_KB=$(df -kP . | gawk 'NR==2{print $4}')
# (zcat || true): head closes the pipe early; without this, pipefail turns SIGPIPE into a fatal status.
D_COUNT_COLS=$( (zcat "$D_COUNTS" || true) | head -1 | gawk -F'\t' '{print NF-1}')
D_SEC_IN_COUNTS=$( (zcat "$D_COUNTS" || true) | head -1 | tr '\t' '\n' | { grep -c -x -F -f <(gawk -F'\t' 'NR>1 && $2=="secretory"{print $1}' "$D_LEDGER") || true; })
H_SEC_IN_PROV=$(grep -c -E "\"($(gawk -F'\t' 'NR>1 && $4=="secretory" && $16=="True"{printf "%s%s", (n++?"|":""), $1}' "$H_LEDGER"))\"" "$H_PROV" || true)

# ---- 5. Receipt ----------------------------------------------------------------------------
json_kv_sha() { local n=$# i=0; for f in "$@"; do i=$((i+1)); printf '    "%s": "%s"%s\n' "$f" "$(sha "$f")" "$([ "$i" -lt "$n" ] && echo ',')"; done; }
TS=$(date -u +%Y-%m-%dT%H:%M:%S+00:00)
{
cat <<EOF
{
  "schema_version": 1,
  "phase": "PHASE_POSITIVE_CONTROL_FEASIBILITY_METADATA_ONLY",
  "created_utc": "$TS",
  "git_head": "$(git rev-parse HEAD)",
  "protocol": "docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md",
  "script": "scripts/phase_positive_control_feasibility.sh",
  "inputs_sha256": {
$(json_kv_sha "$CFG" "$D_V1" "$D_ALL" "$H_V1" "$H_OBS" "$H_DON" "$H_LIB" "$D_COUNTS" "$H_PROV")
  },
  "outputs_sha256": {
$(json_kv_sha "$D_LEDGER" "$H_LEDGER" "$D_DESIGN" "$D_DESIGN_FB" "$H_DESIGN" "$H_DESIGN_FB" "$DL")
  },
  "discovery": {
    "donors_total": $D_N, "proliferative_v1": $D_PRO, "secretory_candidates": $D_SEC, "secretory_ENDO": $D_SEC_ENDO, "secretory_non_ENDO": $((D_SEC-D_SEC_ENDO)),
    "periovulatory_excluded_as_ambiguous": $D_PERIOVULATORY,
    "secretory_excluded_by_rule": "$D_EXCLUDED",
    "candidate_ids_match_frozen_config": $D_CFG_MATCH,
    "group_minimum": $MIN_D, "group_minimum_pass": $D_GROUPS_PASS,
    "primary_design": {"file": "$D_DESIGN", "n": $D1_N, "p": $D1_P, "rank": $D1_R, "residual_df": $D1_DF, "age_mean_levels_ref_adenoconst": "$(cat "$TMP/d1.meta")"},
    "fallback_design": {"file": "$D_DESIGN_FB", "n": $D2_N, "p": $D2_P, "rank": $D2_R, "residual_df": $D2_DF, "age_mean_levels_ref_adenoconst": "$(cat "$TMP/d2.meta")"},
    "phase_by_capture_batch": "$D_XTAB",
    "design_decision": "$D_DECISION",
    "count_matrix_donor_columns": $D_COUNT_COLS,
    "secretory_donors_present_in_count_matrix": ${D_SEC_IN_COUNTS:-0},
    "count_source_required": {"url": "$GEO_RDS_URL", "recorded_bytes": $GEO_RDS_RECORDED_BYTES, "head_status": "$GEO_STATUS", "head_bytes": "$GEO_BYTES", "head_matches_recorded": $GEO_MATCH}
  },
  "heca": {
    "obs_rows_scanned": $H_OBS_ROWS,
    "donors_total_selected": $H_N, "proliferative_v1": $H_PRO, "secretory_selected": $H_SEC, "secretory_ENDO": $H_SEC_ENDO, "secretory_non_ENDO": $((H_SEC-H_SEC_ENDO)),
    "secretory_samples_considered": $H_SEC_SAMPLES,
    "secretory_samples_not_selected": "$H_NOT_SELECTED",
    "age_mismatch_rows": $H_AGE_MISMATCH,
    "candidate_ids_match_frozen_config": $H_CFG_MATCH,
    "group_minimum": $MIN_H, "group_minimum_pass": $H_GROUPS_PASS,
    "primary_design": {"file": "$H_DESIGN", "n": $H1_N, "p": $H1_P, "rank": $H1_R, "residual_df": $H1_DF, "age_mean_levels_ref_adenoconst": "$(cat "$TMP/h1.meta")"},
    "fallback_design": {"file": "$H_DESIGN_FB", "n": $H2_N, "p": $H2_P, "rank": $H2_R, "residual_df": $H2_DF, "age_mean_levels_ref_adenoconst": "$(cat "$TMP/h2.meta")"},
    "phase_by_library_group_Lane": "$H_XTAB",
    "design_decision": "$H_DECISION",
    "secretory_donors_present_in_count_provenance": ${H_SEC_IN_PROV:-0},
    "download_manifest": {"file": "$DL", "archives": $DL_N, "head_200_with_size": $DL_OK, "total_bytes": $DL_BYTES}
  },
  "resources": {
    "repository_volume_available_kib": $DISK_AVAIL_KB,
    "operator_reported_host_ram_gib": "${HOST_RAM_GIB:-NA}",
    "operator_reported_wsl_memory_gb": "${WSL_MEMORY_GB:-NA}",
    "operator_reported_wsl_swap_gb": "${WSL_SWAP_GB:-NA}"
  },
  "expression_values_read": false,
  "downloads_started": false,
  "models_fitted": 0,
  "v1_artifacts_modified": false
}
EOF
} > "$RECEIPT"

echo "discovery: n=$D_N (pro $D_PRO / sec $D_SEC; ENDO among sec $D_SEC_ENDO); config match $D_CFG_MATCH; primary rank $D1_R/$D1_P df $D1_DF; fallback rank $D2_R/$D2_P df $D2_DF; decision $D_DECISION"
echo "discovery excluded secretory by rule: $D_EXCLUDED"
echo "discovery phase x batch: $D_XTAB"
echo "discovery secretory donors already in count matrix: ${D_SEC_IN_COUNTS:-0} of $D_SEC; GEO RDS HEAD $GEO_STATUS $GEO_BYTES (recorded $GEO_RDS_RECORDED_BYTES, match $GEO_MATCH)"
echo "heca: selected n=$H_N (pro $H_PRO / sec $H_SEC; ENDO among sec $H_SEC_ENDO); samples considered $H_SEC_SAMPLES; config match $H_CFG_MATCH; age mismatches $H_AGE_MISMATCH"
echo "heca not selected: $H_NOT_SELECTED"
echo "heca primary rank $H1_R/$H1_P df $H1_DF; fallback rank $H2_R/$H2_P df $H2_DF; decision $H_DECISION"
echo "heca phase x Lane: $H_XTAB"
echo "heca download: $DL_OK/$DL_N archives sized; total $DL_BYTES bytes; secretory donors in count provenance: ${H_SEC_IN_PROV:-0}"
echo "receipt: $RECEIPT"
