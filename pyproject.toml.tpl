[build-system]
requires = [ "setuptools" ]
build-backend = "setuptools.build_meta"

[project]
name = "kimvv"
version = "0.1.6"
dynamic = ["dependencies", "readme"]
authors = [
    { name = "ilia Nikiforov <nikif002@umn.edu>" },
    { name = "Eric Fuemmeler <efuemmel@umn.edu>" },
]
maintainers = [
    { name = "ilia Nikiforov", email = "nikif002@umn.edu" },
]
description = "OpenKIM material property computations for arbitrary crystal structures as Python classes"
requires-python = ">=3.8"
classifiers=[
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: Common Development and Distribution License 1.0 (CDDL-1.0)",
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Physics",
]
keywords = [ "interatomic potential", "openkim", "crystal genome",  ]

[project.urls]
Homepage = "https://github.com/openkim/kimvv"
Issues = "https://github.com/openkim/kimvv/issues"

[tool.setuptools.dynamic]
dependencies = {file = ["requirements.txt"]}
readme = {file = ["README.rst"]}

[tool.setuptools.packages.find]
include = [ "kimvv*" ]
