#!/bin/bash
# Ranged, resumable download of the GEO object into DIR with PARTS concurrent streams,
# then concatenation and SHA-256 verification against the author-route record.
# Each part resumes from its current size; stalls below 20 KB/s for 120 s are retried.
set -u
URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266265/suppl/GSE266265_EndoAtlas.rds
SIZE=41729517596
EXPECT=44a832cdbb7b43ebd3b4a2b6e375c52b08a8059229cc60228ab20a29249cf94f
DIR=${1:?target directory}; PARTS=${2:-4}
mkdir -p "$DIR" && cd "$DIR" || exit 2
avail=$(df -kP . | awk 'NR==2{printf "%.0f", $4*1024}'); need=$((2*SIZE + 20*1024*1024*1024))
if [ "$avail" -lt "$need" ]; then echo "insufficient space: avail=$avail need=$need"; exit 2; fi
chunk=$(( (SIZE + PARTS - 1) / PARTS )); chunk=$(( (chunk + 4194303) / 4194304 * 4194304 ))
start_utc=$(date -u +%FT%TZ); t0=$(date +%s)
echo "start $start_utc parts=$PARTS chunk=$chunk avail=$avail"
fetch_part() {
  local k=$1 s=$((k*chunk)) e=$(( (k+1)*chunk - 1 )); [ "$e" -ge "$SIZE" ] && e=$((SIZE-1))
  local want=$((e-s+1)) f=part.$k tries=0 have
  : >> "$f"
  while :; do
    have=$(stat -c %s "$f")
    [ "$have" -eq "$want" ] && break
    if [ "$have" -gt "$want" ]; then echo "part $k overrun $have>$want"; return 3; fi
    tries=$((tries+1)); if [ "$tries" -gt 60 ]; then echo "part $k: too many retries"; return 4; fi
    curl -sS --fail --connect-timeout 30 --speed-limit 20000 --speed-time 120 -r $((s+have))-$e "$URL" >> "$f" || { echo "part $k retry $tries at $(stat -c %s "$f") $(date -u +%T)"; sleep 20; }
  done
  echo "part $k complete: $(stat -c %s "$f") bytes ($s-$e) $(date -u +%T)"
}
pids=(); for k in $(seq 0 $((PARTS-1))); do fetch_part "$k" & pids+=($!); done
rc=0; for pid in "${pids[@]}"; do wait "$pid" || rc=1; done
[ $rc -eq 0 ] || { echo "a part failed"; exit 4; }
total=0; for k in $(seq 0 $((PARTS-1))); do total=$((total + $(stat -c %s part.$k))); done
[ "$total" -eq "$SIZE" ] || { echo "size mismatch $total != $SIZE"; exit 5; }
t1=$(date +%s); echo "download_seconds=$((t1-t0)) rate_MBps=$(awk -v b=$SIZE -v s=$((t1-t0)) 'BEGIN{printf "%.2f", b/1e6/(s>0?s:1)}')"
cat $(for k in $(seq 0 $((PARTS-1))); do printf 'part.%d ' "$k"; done) > GSE266265_EndoAtlas.rds
t2=$(date +%s); echo "concat_seconds=$((t2-t1))"
got=$(sha256sum GSE266265_EndoAtlas.rds | cut -d' ' -f1); t3=$(date +%s)
echo "sha256=$got expected=$EXPECT sha_seconds=$((t3-t2))"
match=false; [ "$got" = "$EXPECT" ] && match=true
printf '{\n  "url": "%s",\n  "bytes": %s,\n  "expected_bytes": %s,\n  "sha256": "%s",\n  "expected_sha256": "%s",\n  "sha256_match": %s,\n  "parts": %s,\n  "start_utc": "%s",\n  "end_utc": "%s",\n  "download_seconds": %s,\n  "concat_seconds": %s,\n  "sha_seconds": %s\n}\n' "$URL" "$(stat -c %s GSE266265_EndoAtlas.rds)" "$SIZE" "$got" "$EXPECT" "$match" "$PARTS" "$start_utc" "$(date -u +%FT%TZ)" "$((t1-t0))" "$((t2-t1))" "$((t3-t2))" > download_receipt.json
cat download_receipt.json
if $match; then rm -f part.*; echo "SHA256_MATCH"; else echo "SHA256_MISMATCH"; exit 6; fi
