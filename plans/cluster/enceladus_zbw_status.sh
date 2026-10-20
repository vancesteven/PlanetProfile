#!/usr/bin/env bash
# Progress / health check for a running Enceladus (zb, w) shard build
# (enceladus_zbw_local.sh or the SLURM array). Read-only.
#
#   plans/cluster/enceladus_zbw_status.sh [SHARD_DIR]
#
# SHARD_DIR defaults to $SHARD_DIR, else ~/scratch/enc_zbw_shards.

SHARD_DIR="${1:-${SHARD_DIR:-$HOME/scratch/enc_zbw_shards}}"
N_SHARDS="${N_SHARDS:-11}"
cd "$SHARD_DIR" 2>/dev/null || { echo "no shard dir: $SHARD_DIR" >&2; exit 1; }

echo "shard dir : $SHARD_DIR   ($(date))"
echo "commit    : $(cat BUILD_COMMIT 2>/dev/null || echo 'BUILD_COMMIT missing')"
n_meta=$(ls shard_*_of_*.meta.json 2>/dev/null | wc -l | tr -d ' ')
echo "finished  : $n_meta/$N_SHARDS shards with .meta.json"
for m in shard_*_of_*.meta.json; do [[ -f "$m" ]] && echo "            $m"; done
echo "running   : $(pgrep -f enceladus_zbw_shard_build | wc -l | tr -d ' ') builder processes"
[[ -f DONE ]] && echo "DONE      : $(cat DONE)" || echo "DONE      : not yet"

# Health: any of these means stop and look before merging.
tb=$(grep -l Traceback shard_*.log 2>/dev/null)
echo "traceback : ${tb:-none}"
q=$(ls ./*QUARANTINED* 2>/dev/null)
echo "C2 quarantine dumps: ${q:-none}"

echo "--- log tails ---"
for f in shard_*.log; do
  [[ -f "$f" ]] || continue
  echo "$f: $(tail -c 200 "$f" | tr '\r\n' '  ' | tail -c 110)"
done
