#!/bin/sh
# Watch docs/source and rebuild whenever anything is saved, then reload any
#
# Anything the pre-build step writes into source/ must be ignored here, or the
# write is seen as an edit and the build loops.  generate_reference.py also
# leaves unchanged files alone, so the two guards are independent.
# browser tab showing the page.  Leave this running while editing .rst files.
#
# Warnings are NOT errors here: a half-finished directive should show you the
# rest of the page rather than blanking it.  Run ./rebuild.sh before committing
# to get the strict check.
cd "$(dirname "$0")"
VENV=/home/jalombar/research/docs-venv/bin
exec "$VENV/sphinx-autobuild" source build \
    --host 127.0.0.1 --port 45645 \
    --pre-build "python3 generate_reference.py" \
    --watch ../parallel_bleeding_edge/src/init.f \
    --ignore "*/reference/sph_input.rst" \
    --ignore "*/reference/settings_by_script.rst" \
    --open-browser
