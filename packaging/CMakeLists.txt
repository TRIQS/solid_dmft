# Configure the version of packaging files
configure_file(conda/meta.yaml.in conda/meta.yaml)
configure_file(pypi/MANIFEST.in pypi/MANIFEST.in COPYONLY)
configure_file(pypi/pyproject.toml.in pypi/pyproject.toml)
configure_file(solid_dmft-foss-2021b.eb.in solid_dmft-foss-2021b.eb)


