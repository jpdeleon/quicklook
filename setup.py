#!/usr/bin/env python
from setuptools import setup


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="tql",
    packages=["tql"],
    version="1.0",
    author="Jerome de Leon",
    author_email="jpdeleon@astron.s.u-tokyo.ac.jp",
    url="https://github.com/jpdeleon/tql",
    license=["GNU GPLv3"],
    description="TESS QuickLook plot generator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    scripts=[
        "scripts/tql",
        "scripts/rank_tls",
    ],
    keywords=["TESS", "exoplanets", "stars"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python",
    ],
    install_requires=[
        "pre-commit",
        "black",
        "flake8",
        "toml",
        "jupytext",
        "nbsphinx",
    ],
    extras_require={
        "tests": [
            "nose",
            "pytest",
            "pytest-dependency",
            "pytest-env",
            "pytest-cov",
        ],
    },
)
