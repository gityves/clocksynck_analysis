#!/usr/bin/env bash
# usage: make_provenance.sh <dir> "<title>" "<source url>" "<note>"
set -euo pipefail
d=$1; { echo "# PROVENANCE: $2"; echo; echo "- Source: $3"; echo "- Downloaded: 2026-11-12 (public, open access; processed files only, nothing executed)"; echo "- Note: $4"; echo; echo "| file | bytes | sha256 |"; echo "|---|---|---|";
cd "$d"; for f in $(ls -A | grep -v -e PROVENANCE.md -e '^\.dl_done$'); do echo "| $f | $(stat -c %s "$f") | $(sha256sum "$f" | cut -d' ' -f1) |"; done; } > "$d/PROVENANCE.md"
