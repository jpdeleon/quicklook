from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

# Read the dependencies from requirements.txt
with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="quicklook-package",
    version="1.4",
    description="Quicklook lightcurve plot generator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jpdeleon/quicklook",
    packages=find_packages(),
    include_package_data=True,
    install_requires=install_requires,
    extras_require={
        "gui": [
            "flask>=2.3",
            "matplotlib",  # needed for plotting in GUI
        ],
    },
    scripts=["scripts/ql", "scripts/read_tls"],
    entry_points={
        "console_scripts": [
            "ql-gui=app.app:main",  # app/app.py must have a main() function
        ],
    },
    python_requires=">=3.8",
)
