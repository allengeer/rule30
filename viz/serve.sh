#!/bin/sh
# Starts the Rule 30 viewer on http://127.0.0.1:${1:-8730}/ (builds the C library first if needed).
cd "$(dirname "$0")" || exit 1
[ -f build/libr30.dylib ] || make >/dev/null || exit 1
mkdir -p /Volumes/BLIPS/rule30/viz/tmp
TMPDIR=/Volumes/BLIPS/rule30/viz/tmp exec ../.venv/bin/python server.py "${1:-8730}"
