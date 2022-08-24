from mztabpy import MzTabPy

import unittest

class TestCases(unittest.TestCase):

    def test_something(self):
        self.assertEqual(True, True)

    def test_mztabpy_to_hdf5(self):
        # mztab storage
        mztab = MzTabPy("../testdata/diatest.sdrf_openms_design_out.mztab", single_cache_size=50 * 1024 * 1024)
        mztab.storage(to_tsv=True, to_hdf5=True)

        # Read HDF5 and filter
        result = MzTabPy.loadHDF5("tests/diatest.sdrf_openms_design_out.hdf5", subtable='peptide',
                                  where={'accession': 'P68066', 'sequence': 'AANDDLLNSFWLLDSEKGEAR'})
        print(result)

if __name__ == '__main__':
    unittest.main()


