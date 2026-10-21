#!/usr/bin/env bash
# Tally WHY each ocean node of an Enceladus (zb, w) shard build is None
# (r7 rebuild adjudication condition C1). The cache records only the reason
# category ('build_exception' lumps every exception); the shard logs hold
# the actual exception per node. C1 requires every None to be the
# bracket / same-sign "target outside achievable zb" failure, with no I/O,
# cache-path or EOS errors. Read-only on the shard dir.
#
#   plans/cluster/enceladus_zbw_none_reasons.sh [SHARD_DIR]
#
# Writes the full per-node list and the tally to
# validation_reports/enceladus_isostasy/r7_rebuild_none_reasons.txt so the
# result is in the repo for review.

SHARD_DIR="${1:-${SHARD_DIR:-$HOME/scratch/enc_zbw_shards}}"
EXPECT="${EXPECT:-121}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO/validation_reports/enceladus_isostasy/r7_rebuild_none_reasons.txt"

# Builder lines look like:
#   × (zb=12.0000, w=4.9239) → None — ValueError: <message>
lines=$(cat "$SHARD_DIR"/shard_*.log 2>/dev/null | grep -a ') → None — ')
n=$(printf '%s\n' "$lines" | grep -c . )

{
  echo "# r7 rebuild None reasons — $(date)"
  echo "# shard dir: $SHARD_DIR   build commit: $(cat "$SHARD_DIR/BUILD_COMMIT" 2>/dev/null)"
  echo "# None lines in logs: $n (cache says $EXPECT)"
  echo
  echo "## Tally by exception type + message (numbers replaced by #)"
  printf '%s\n' "$lines" | sed -E 's/.*→ None — //; s/[-+]?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?/#/g' \
    | sort | uniq -c | sort -rn
  echo
  echo "## Per node"
  printf '%s\n' "$lines" | sed -E 's/.*(\(zb=[^)]*\)) → None — /\1 /'
} | tee "$OUT"

echo
if [[ "$n" -ne "$EXPECT" ]]; then
  echo "MISMATCH: $n None lines in logs vs $EXPECT None nodes in the cache" >&2
  exit 1
fi
echo "wrote $OUT"
