"""
Proteotyping setup.

Online documentation:
    https://proteotyping.readthedocs.org

Online repository:
    https://bitbucket.org/chalmersmathbioinformatics/proteotyping
"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Read long description from README file
with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="proteotyping",

    # This project aims to follow semantic versioning (www.semver.org)
    version="1.0.0",

    description="Proteotyping",
    long_description=long_description,

    url="https://bitbucket.org/chalmersmathbioinformatics/proteotyping",

    author="Fredrik Boulund",
    author_email="fredrik.boulund@chalmers.se",

    license="BSD",

    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],

    keywords="proteomics proteotyping antibiotic resistance",

    packages=find_packages(exclude=["docs", "tests"]),

    install_requires=[
        "ete3",
        "xlsxwriter",
        "requests",
        "six",
        "scipy",
        "numpy",
    ],

    extras_require={
        "dev": ["twine", "wheel"],
        "docs": ["sphinx", "sphinx-autobuild", "sphinx-rtd-theme"]
    },

    entry_points={
        "console_scripts": [
            "annotation_db = proteotyping.annotation_db:main",
            "antibiotic_resistance = proteotyping.antibiotic_resistance:main",
            "construct_resfinder_db = proteotyping.construct_resfinder_db:main",
            "taxonomic_composition = proteotyping.taxonomic_composition:main",
            "taxref_db = proteotyping.taxref_db:main"
        ]
    }
)

