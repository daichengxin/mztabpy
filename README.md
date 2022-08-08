# mztabpy
Python library to handle mzTab files
## Input
mzTab
## Output
Protein, peptide, and PSM subtables(.csv) or HDF5 file(.hdf5) with three parts of information.
## Usage
`python mztabpy_click.py mztab_parse --mztab_path {mztab_path}
 --subtable {True or False}
 --hdf5 {True or False}`
 ## Parameters
 -   --mztab_path: The path to mzTab
 -   --subtable: Whether to convert mzTab into subtables(.csv)
 -  --hdf5: Whether to convert mzTab into HDF5(.hdf5)
