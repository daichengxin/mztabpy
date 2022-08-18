from mzTabpy import mzTabpy

# mztab storage
mztab = mzTabpy("test/PXD026600.sdrf_openms_design_out.mztab", single_cache_size=50*1024*1024)
mztab.storage(to_csv=True, to_hdf5=True)

# Read HDF5 and filter
result = mzTabpy.loadHDF5("test/PXD026600.sdrf_openms_design_out.hdf5", subtable='peptide', where={'accession': 'P60240', 'sequence': 'NGLDLILSGDTGSSTISLLK'})
print(result)