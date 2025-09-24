from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

# Read the dependencies from requirements.txt
with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

# see pyproject.toml for details such as version and requirements
setup(
    name="quicklook-package",
    version="1.3.2",
    description="Quicklook lightcurve plot generator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jpdeleon/quicklook",
    packages=find_packages(),
    install_requires=install_requires,
    include_package_data=True,
    scripts=["scripts/ql", "scripts/read_tls"],
)
