#!/bin/sh
# Deep run: per-row stats plus row periods sharded across cores.
#   ./deep.sh [rows] [out_dir] [shards]
set -e
cd "$(dirname "$0")"
ROWS=${1:-1000000}
OUT=${2:-out_1m}
SHARDS=${3:-12}
mkdir -p build "$OUT"
cc -O3 -march=native -Wall -o build/rule30 rule30.c

LAST=$((ROWS - 1))
build/rule30 stats "$ROWS" "$OUT" 1000 10000 100000 250000 500000 750000 "$LAST" 2> "$OUT/stats.log" &
i=0
while [ $i -lt "$SHARDS" ]; do
  build/rule30 period "$ROWS" "$OUT" $i "$SHARDS" &
  i=$((i + 1))
done
wait
echo "done: $ROWS rows in $OUT"
