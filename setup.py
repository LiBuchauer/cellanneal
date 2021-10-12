from setuptools import find_packages
from setuptools import setup


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="cellanneal",
    version="0.1.0",
    url="https://github.com/libuchauer/cellanneal",
    license='MIT',

    author="Lisa Buchauer",
    author_email="lisa.buchauer@posteo.de",

    description="A user-friendly application for deconvolving bulk RNAseq samples.",
    long_description=long_description,
    packages=find_packages(),

    # package dependencies
    install_requires=[
            'numpy',
            'scipy',
            'matplotlib',
            'pandas',
            'seaborn',
            'xlrd',
            'openpyxl'],

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],

    entry_points = {
        'console_scripts': [
            'cellanneal = cellanneal.__main__:main'
            ]
        }
)
