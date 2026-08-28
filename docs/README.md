# StarSmasher documentation

Sphinx sources for the StarSmasher documentation.

## Building locally

```
python3 -m venv env
env/bin/pip install -r docs/requirements.txt
env/bin/sphinx-build -b html docs/source docs/build
```

Then open `docs/build/index.html`.

## The sph.input reference is generated

`docs/source/reference/sph_input.rst` is produced from
`parallel_bleeding_edge/src/init.f` by `docs/generate_reference.py`, so the
listed settings, defaults and descriptions cannot drift from the code. After
changing the namelist or its defaults in `init.f`, run:

```
python3 docs/generate_reference.py
```

and commit the regenerated page. The CI workflow fails if you forget.
