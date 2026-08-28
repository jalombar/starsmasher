#!/bin/sh
# One-off rebuild of the documentation.  For editing, prefer ./livedocs.sh,
# which rebuilds by itself whenever a source file is saved.
#
# Unlike livedocs.sh this builds with -W, so a bad cross-reference or a
# malformed directive fails here rather than becoming a dead link.  Run it once
# before committing.
set -e
cd "$(dirname "$0")"
VENV=/home/jalombar/research/docs-venv/bin
# Regenerate the sph.input reference from init.f first -- it is not hand-written.
python3 generate_reference.py
"$VENV/sphinx-build" -b html -W source build
echo
echo "Rebuilt, with warnings treated as errors."
