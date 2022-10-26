from mztabpy import MzTabPy
from diannconvert import DiannConvert
from mztabmerge import MzTabMerge

import unittest


class TestCases(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)

    def test_mztabpy_to_hdf5(self):
        # mztab storage
        mztab = MzTabPy(
            "testdata/split&convert/diatest.sdrf_openms_design_out.mztab",
            single_cache_size=50 * 1024 * 1024,
            result_folder="./",
        )
        mztab.storage(type="all", section="all", removemeta=True)

        # Read HDF5 and filter
        result = MzTabPy.loadHDF5(
            "testdata/diatest.sdrf_openms_design_out.hdf5",
            subtable="peptide",
            where={"accession": "P68066", "sequence": "AANDDLLNSFWLLDSEKGEAR"},
        )
        print(result)

    def test_diannconvert(self):
        DiannConvert(
            "testdata/diannconvert/",
            dia_params="20;ppm;10;ppm;Trypsin;Carbamidomethyl (C);Oxidation (M)",
            charge=3,
            missed_cleavages=1,
            qvalue_threshold=0.05,
            out="./",
        )

    def test_mztabmerge(self):
        MzTabMerge(
            mztab1="testdata/merge/out1.mztab",
            mztab2="testdata/merge/out2.mztab",
            single_cache_size=50 * 1024 * 1024,
            result_folder="./",
        )


if __name__ == "__main__":
    unittest.main()
