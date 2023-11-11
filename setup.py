from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="quicklook",
    version="1.0",
    description="Quicklook lightcurve plot generator",
    long_description=long_description,
    author="Jerome de Leon",
    author_email="jpdeleon@g.ecc.u-tokyo.ac.jp",
    url="https://github.com/jpdeleon/quicklook",
    packages=find_packages(),
    include_package_data=True,
    scripts=["scripts/ql", "scripts/read_tls"],
)
