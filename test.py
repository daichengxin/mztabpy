from mztabpy import MzTabPy

# mztab storage
mztab = MzTabPy("test/diatest.sdrf_openms_design_out.mztab", single_cache_size=50*1024*1024)
mztab.storage(to_tsv=True, to_hdf5=True)

# Read HDF5 and filter
result = MzTabPy.loadHDF5("test/diatest.sdrf_openms_design_out.hdf5", subtable='peptide', where={'accession': 'P68066', 'sequence': 'AANDDLLNSFWLLDSEKGEAR'})
print(result)
