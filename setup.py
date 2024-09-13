from setuptools import setup, find_packages

# Read the dependencies from requirements.txt
with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

# see pyproject.toml for details such as version
setup(
    name="quicklook",  
    description="Quicklook lightcurve plot generator",
    url="https://github.com/jpdeleon/quicklook", 
    packages=find_packages(), 
    include_package_data=True,  
    scripts=["scripts/ql", "scripts/read_tls"], 
    install_requires=install_requires,
    project_urls={  # Additional URLs related to your project
        "Bug Tracker": "https://github.com/jpdeleon/quicklook/issues",
        "Documentation": "https://github.com/jpdeleon/quicklook#readme",
        "Source Code": "https://github.com/jpdeleon/quicklook",
    },
)