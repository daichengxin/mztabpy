# mztabpy

[![Python application](https://github.com/bigbio/mztabpy/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/mztabpy/actions/workflows/python-app.yml)

Python library to handle mztab files. mzTab is a tab-delimited file format created by HUPO-PSI containing protein/peptide quantification and identification data. 


## Input
#### mztab_convert
mzTab
#### hdf5_search
HDF5

## Output
#### mztab_convert
Metadata, protein, peptide, and PSM subtables(`.csv`) or HDF5 file(`.hdf5`) with the four parts of information.
#### hdf5_search
A dataframe that can be filtered or already filtered based on the condition.

## Usage
#### mztab_convert
`python mztabpy_click.py mztab_convert --mztab_path {mztab_path}
 --directory {result folder}
 --type {result type}
 --section {section of the mzTab}
 --removemeta {True or False}`
#### hdf5_search
`python mztabpy_click.py hdf5_search --hdf5_path {hdf5_path}
 --subtable {section of the mzTab}
 --where {filtering condition}`


## Parameters
#### mztab_convert
-   --mztab_path: The path to mzTab
-   --directory: Folder to result files. Default "./"
-   --type: Result type(`"tsv"`, `"hdf5"` or `"all"`). Default "all"
-   --section: Indicates the data section of the mzTab that is required. `"all"`, `"protein"`, `"peptide"` or `"psm"`.Default "all"
-   --removemeta: Whether to remove `metadata`. Default False
#### hdf5_search
-   --hdf5_path: Path to HDF5
-   --section: Indicates the data section of the mzTab that is required. `"protein"`, `"peptide"` or `"psm"`.
-   --where: The filtering condition of the corresponding chunk is expressed as the key-value pair in one string, e.g. `"accession:P45464,sequence:TIQQGFEAAK"`, default None
