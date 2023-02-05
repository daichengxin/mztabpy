import os
import pandas as pd
import dask.dataframe as dd
import numpy as np
from collections import OrderedDict
import re
import logging as log
from pyopenms import AASequence, FASTAFile, ModificationsDB
from copy import copy

pd.set_option("mode.chained_assignment", None)
log.basicConfig(format="%(levelname)s:%(message)s", level=log.DEBUG)
# os.environ["NUMEXPR_MAX_THREADS"] = "16"


class MzTabPy:
    """
    The MzTabPy class is used to separate mzTab into metadata, protein, peptide, and PSM. It can convert mzTab to
    4 split tables(.tsv) or HDF5(.hdf5) depending on the options. The four subtables will be cached as dataframes
    if the given cache size allows.

    :param mztab_path: The path to mzTab
    :type mztab_path: str
    :param single_cache_size: Single cache size, default 500*1024*1024
    :type single_cache_size: int
    :param result_folder: Folder to result files, default './'
    :type result_folder: str
    :param section_tag: Whether to keep section tags, default False
    :type section_tag: bool
    """

    def __init__(
        self,
        mztab_path,
        single_cache_size=500 * 1024 * 1024,
        result_folder="./",
        section_tag=False,
        df_index=True,
    ):
        self.mztab_path = mztab_path
        self.chunk_size = single_cache_size
        self.folder_path = result_folder
        self.sectag = section_tag
        self.df_index = df_index
        file_name = os.path.split(mztab_path)[1]
        self.basename = os.path.splitext(file_name)[0]
        self.file_bytes = os.path.getsize(mztab_path)
        self.datatype = dict()
        self.chunks_num = int(self.file_bytes / self.chunk_size) + 1
        self.psm_compression = None if self.file_bytes < 1 * pow(1024, 3) else "gzip"

    def loadmzTab(self):
        """If the computer has enough cache, read mzTab as four global dataframes('metadata', 'protein', 'peptide' and 'psm')."""
        self.meta_dict = OrderedDict()
        with open(self.mztab_path, "r") as f:
            chunk = f.read(self.file_bytes)
            (self.protein, self.peptide, self.psm) = self.loadChunk(chunk)
            self.metadata = pd.DataFrame(self.meta_dict, index=[0])
        log.info("Four subtables have been cached!")

        return self.meta_dict, self.protein, self.peptide, self.psm

    def loadChunk(self, chunk):
        """This function processes data blocks from streams, classifies them into metadata, protein, peptide and
            PSM, and caches them as dataframes.

        :param chunk: A block of data
        :type chunk: str
        """
        data = chunk.split("\n")
        if not chunk.endswith("\n"):
            self.row_remains = data[-1]
            data = data[0:-1]
        else:
            self.row_remains = ""
        meta_dict = dict()
        pro_data, pep_data, psm_data = [], [], []
        for i in data:
            i = i.strip("\n")
            row_list = i.split("\t")
            row_list = [cast_value(i) for i in row_list]
            if row_list[0] == "MTD":
                meta_dict[row_list[1]] = row_list[2]
            elif row_list[0] == "PRH":
                log.info(
                    "Metadata processing complete. Start processing the protein subtable..."
                )
                self.pro_cols = row_list[1:]
                self.classify_datatype(self.pro_cols)
            elif row_list[0] == "PRT":
                pro_data.append(row_list[1:])
            elif row_list[0] == "PEH":
                log.info(
                    "Protein subtable processing complete. Start processing the peptide subtable..."
                )
                self.pep_cols = row_list[1:]
                self.classify_datatype(self.pep_cols)
            elif row_list[0] == "PEP":
                pep_data.append(row_list[1:])
            elif row_list[0] == "PSH":
                log.info(
                    "Peptide subtable processing complete. Start processing the psm subtable..."
                )
                self.psm_cols = row_list[1:]
                self.classify_datatype(self.psm_cols)
            elif row_list[0] == "PSM":
                psm_data.append(row_list[1:])
            else:
                continue

        for i in meta_dict.keys():
            self.meta_dict.update(
                {
                    i: ";".join(list(meta_dict[i]))
                    if type(meta_dict[i]) == tuple
                    else meta_dict[i]
                }
            )

        if len(pro_data) > 0:
            protein_table = pd.DataFrame(pro_data, columns=self.pro_cols)
            del pro_data
            if any(protein_table.columns.duplicated()):
                protein_table = protein_table.loc[
                    :, ~protein_table.columns.duplicated()
                ]
        else:
            protein_table = pd.DataFrame()

        if len(pep_data) > 0:
            peptide_table = pd.DataFrame(pep_data, columns=self.pep_cols)
            del pep_data
            if any(peptide_table.columns.duplicated()):
                peptide_table = peptide_table.loc[
                    :, ~peptide_table.columns.duplicated()
                ]
        else:
            peptide_table = pd.DataFrame()

        if len(psm_data) > 0:
            psm_table = pd.DataFrame(psm_data, columns=self.psm_cols)
            del psm_data
            if any(psm_table.columns.duplicated()):
                psm_table = psm_table.loc[:, ~psm_table.columns.duplicated()]
            if "opt_global_cv_MS:1002217_decoy_peptide" not in psm_table.columns.values:
                psm_table["opt_global_cv_MS:1002217_decoy_peptide"] = psm_table.apply(
                    lambda x: 1
                    if any(i in x["accession"] for i in ["DECOY_", "decoy"])
                    else 0,
                    axis=1,
                )
            psm_table[["start", "end"]] = psm_table[["start", "end"]].astype(str)
        else:
            psm_table = pd.DataFrame()

        return protein_table, peptide_table, psm_table

    def storage(self, type="all", section="all", removemeta=False):
        """Store mzTab as tsv subtables or HDF5.

        :param type: Result type("tsv", "hdf5" or "all"). Default "all"
        :type type: str
        :param section: Indicates the data section of the mzTab that is required. "all", "protein", "peptide" or "psm".
        Default "all"
        :type section: str
        :param removemeta: Whether to remove metadata. Default False
        :type removemeta: bool
        """
        if type == "tsv" or type == "all":
            self.to_tsv = True
            self.meta_path = self.folder_path + self.basename + "_metadata_openms.tsv"
            self.pro_path = self.folder_path + self.basename + "_protein_openms.tsv"
            self.pep_path = self.folder_path + self.basename + "_peptide_openms.tsv"
            self.psm_path = (
                self.folder_path
                + self.basename
                + (
                    "_psm_openms.tsv"
                    if self.file_bytes < 1 * pow(1024, 3)
                    else "_openms.psms.mzTab.gz"
                )
            )
        else:
            self.to_tsv = False

        if type == "hdf5" or type == "all":
            self.to_hdf5 = True
            self.h5_path = self.folder_path + self.basename + ".hdf5"
        else:
            self.to_hdf5 = False

        log.info("Start converting...")
        self.stream_storage(section, removemeta)

    def stream_storage(self, section, removemeta):
        """This function streams the mzTab and stores it in four tsv subtables (Metadata, protein, peptide, PSM)
            or one HDF5 containing these four parts.

        :param section: Indicates the data section of the mzTab that is required. "all", "protein", "peptide" or "psm"
        :type section: str
        :param removemeta: Whether to remove metadata
        :type removemeta: bool
        """
        self.meta_dict = OrderedDict()
        self.chunk_index = 0
        pro_chunk, pep_chunk, psm_chunk = [], [], []
        self.chunk_dict = dict.fromkeys(["protein", "peptide", "psm"], "")
        self.chunk_dict.update(
            {"mzTab_file_path": self.mztab_path, "mztab_file_bytes": self.file_bytes}
        )
        self.row_remains = ""
        log.info("Chunk Size: {} MB".format(self.chunk_size / pow(1024, 2)))
        with open(self.mztab_path, "r") as f:
            while True:
                chunk = self.row_remains + f.read(self.chunk_size)
                if not chunk:
                    break
                start = time()
                (pro_df, pep_df, psm_df) = self.loadChunk(chunk)
                print(f"Process chunk {self.chunk_index} cost {time() - start} secs")
                if len(pro_df) > 0 and (section == "all" or section == "protein"):
                    if self.to_tsv:
                        if os.path.exists(self.pro_path):
                            pro_df.to_csv(
                                self.pro_path,
                                mode="a",
                                sep="\t",
                                index=False,
                                header=False,
                            )
                        else:
                            pro_df.to_csv(
                                self.pro_path,
                                mode="w",
                                sep="\t",
                                index=False,
                                header=True,
                            )
                    if self.to_hdf5:
                        pro_chunk.append(f"chunk_{self.chunk_index}_protein")
                        if "protein_cols" not in self.chunk_dict:
                            self.chunk_dict.update(
                                {"protein_cols": ";".join(pro_df.columns.values)}
                            )
                        pro_df.fillna(np.nan, inplace=True)
                        pro_df.to_hdf(
                            self.h5_path,
                            mode="a",
                            key=f"chunk_{self.chunk_index}_protein",
                            format="f",
                            index=False,
                            complevel=9,
                            complib="zlib",
                        )
                    del pro_df

                if len(pep_df) > 0 and (section == "all" or section == "peptide"):
                    if self.to_tsv:
                        if os.path.exists(self.pep_path):
                            pep_df.to_csv(
                                self.pep_path,
                                mode="a",
                                sep="\t",
                                index=False,
                                header=False,
                            )
                        else:
                            pep_df.to_csv(
                                self.pep_path,
                                mode="w",
                                sep="\t",
                                index=False,
                                header=True,
                            )
                    if self.to_hdf5:
                        pep_chunk.append(f"chunk_{self.chunk_index}_peptide")
                        if "peptide_cols" not in self.chunk_dict:
                            self.chunk_dict.update(
                                {"peptide_cols": ";".join(pep_df.columns.values)}
                            )
                        pep_df.fillna(np.nan, inplace=True)
                        pep_df.to_hdf(
                            self.h5_path,
                            mode="a",
                            key=f"chunk_{self.chunk_index}_peptide",
                            format="f",
                            index=False,
                            complevel=9,
                            complib="zlib",
                        )
                    del pep_df

                if len(psm_df) > 0 and (section == "all" or section == "psm"):
                    if self.to_tsv:
                        if os.path.exists(self.psm_path):
                            psm_df.to_csv(
                                self.psm_path,
                                mode="a",
                                sep="\t",
                                index=False,
                                header=False,
                                compression=self.psm_compression,
                            )
                        else:
                            psm_df.to_csv(
                                self.psm_path,
                                mode="w",
                                sep="\t",
                                index=False,
                                header=True,
                                compression=self.psm_compression,
                            )
                    if self.to_hdf5:
                        psm_chunk.append(f"chunk_{self.chunk_index}_psm")
                        if "psm_cols" not in self.chunk_dict:
                            self.chunk_dict.update(
                                {"psm_cols": ";".join(psm_df.columns.values)}
                            )
                        psm_df.fillna(np.nan, inplace=True)
                        psm_df.to_hdf(
                            self.h5_path,
                            mode="a",
                            key=f"chunk_{self.chunk_index}_psm",
                            format="f",
                            complevel=9,
                            complib="zlib",
                        )
                    del psm_df

                self.chunk_index += 1

        meta_df = pd.DataFrame(self.meta_dict, index=[0])
        self.chunk_dict.update(
            {
                "protein": ";".join(pro_chunk),
                "peptide": ";".join(pep_chunk),
                "psm": ";".join(psm_chunk),
            }
        )

        chunks_info = pd.DataFrame(self.chunk_dict, index=[0])

        if self.to_tsv:
            if not removemeta:
                meta_df.to_csv(
                    self.meta_path, mode="a", sep="\t", index=False, header=True
                )
            log.info("mzTab separation and storage complete!")
        if self.to_hdf5:
            chunks_info.to_hdf(
                self.h5_path,
                mode="a",
                key="CHUNKS_INFO",
                format="f",
                complevel=9,
                complib="zlib",
            )
            if not removemeta:
                meta_df.to_hdf(
                    self.h5_path,
                    mode="a",
                    key="meta",
                    format="f",
                    complevel=9,
                    complib="zlib",
                )
            log.info("mzTab binary storage complete!")

    def loadHDF5(h5, subtable, where=None, tsv=False):
        """Load HDF5 into dataframe

        :param h5: Path to HDF5
        :type h5: str
        :param subtable: Name of the subtable, which contains 'metadata', 'protein', 'peptide' and 'psm'
        :type subtable: str
        :param where: The filtering condition of the corresponding chunk is expressed as the key-value pair in the dictionary,
            e.g. {'accession': 'P45464', 'sequence': 'TIQQGFEAAK'}, default None
        :type where: dict
        :param tsv: Whether to save the result as tsv, default False
        :type tsv: bool
        :return result: Results filtered by criteria in the target subtable
        :rtype result: dataframe
        """
        chunks_info = pd.read_hdf(h5, key="CHUNKS_INFO")
        chunks = (
            chunks_info[subtable].values[0].split(";")
            if subtable in chunks_info
            else []
        )
        result = pd.DataFrame(
            columns=chunks_info[subtable + "_cols"].values[0].split(";")
        )

        for i in chunks:
            df = pd.read_hdf(h5, key=i)
            if where:
                for j in where.keys():
                    df = df[df[j] == where[j]]
            result = pd.concat([result, df])

        return result

    def classify_datatype(self, columns):
        """Classify data types according to column names

        :param columns: Columns of target dataframe
        :type columns: list
        """
        for col in columns:
            if any(
                [
                    i in col
                    for i in [
                        "score",
                        "abundance",
                        "mass_to_charge",
                        "study_variable",
                        "q-value",
                    ]
                ]
            ) or col in ["protein_coverage", "retention_time"]:
                self.datatype.update({col: float})
            elif any([i in col for i in ["ID", "_decoy_", "index"]]) or col in [
                "opt_global_nr_found_peptides",
                "unique",
                "charge",
                "start",
                "end",
            ]:
                self.datatype.update({col: int})
            else:
                self.datatype.update({col: str})

    def resetdatatype(self, df):
        """Reset data type in target dataframe

        :param df: Target dataframe
        :type df: dataframe
        """
        for col in df.columns:
            if self.datatype[col] != "str":
                df.loc[:, col] = df.loc[:, col].astype(
                    self.datatype[col], errors="ignore"
                )


class DiannConvert:
    """DiannConvert is a tool for DiaNN's mass spectrometry-based protein identification and quantitative result
    conversion. It can output three proteomics standard format files: Triqler, MSstats and mzTab

    :param FolderPath: DiannConvert specifies the folder where the required file resides. The folder contains the
     DiaNN main report, experimental design file, protein sequence FASTA file, version file of DiaNN and mzml_info
     TSVs
    :type FolderPath: str
    :param dia_params: A string contains DIA parameters, e.g. "20;ppm;10;ppm;Trypsin;Carbamidomethyl (C);Oxidation (M)"
    :type dia_params: str
    :param charge: The charge assigned by DiaNN(max_precursor_charge)
    :type charge: int
    :param missed_cleavages: Allowed missed cleavages assigned by DiaNN
    :type missed_cleavages: int
    :param qvalue_threshold: Threshold for filtering q value
    :type qvalue_threshold: float
    :param processors: Number of used processors, defaults to 20
    :type processors: int
    :param threads_per_processor: Number of threads used per processor, defaults to 8
    :type threads_per_processor: int
    :param out: Path to out directory, defaults to "./"
    :type out: str
    :param block_size: Chunk size, defaults to 500e6
    :type block_size: int
    """

    def __init__(
        self,
        FolderPath,
        dia_params,
        charge,
        missed_cleavages,
        qvalue_threshold,
        processors=20,
        threads_per_processor=8,
        out="./",
        block_size=500e6,
    ):
        """Init"""
        pathdict = {
            key: []
            for key in ["report", "exp_design", "fasta", "diann_version", "mzml_info"]
        }
        self.folder = FolderPath
        fileslist = os.listdir(FolderPath)
        if not FolderPath.endswith("/"):
            FolderPath = FolderPath + "/"
        for i in fileslist:
            if i.endswith("report.tsv"):
                pathdict["report"].append(i)
            elif i.endswith("design.tsv"):
                pathdict["exp_design"].append(i)
            elif i.endswith(".fasta") or i.endswith(".fa"):
                pathdict["fasta"].append(i)
            elif i.endswith("versions.yml"):
                pathdict["diann_version"].append(i)
            elif i.endswith("mzml_info.tsv"):
                pathdict["mzml_info"].append(i)
            else:
                pass

        for item in pathdict.items():
            if len(item[1]) == 0:
                log.error(f"can't find {item[0]}, please check your file folder!")

            if item[0] != "mzml_info" and len(item[1]) > 1:
                log.error(
                    f"{item[0]} is duplicate, check whether the file is redundant or change the file name!"
                )

        if len(pathdict["diann_version"]) == 1:
            with open(FolderPath + pathdict["diann_version"][0]) as f:
                for line in f:
                    if "DIA-NN" in line:
                        diann_version_id = line.rstrip("\n").split(": ")[1]
                        break
            f.close()
        else:
            log.error("Can't find 'version.yml', diann_version_id defaults to '1.8.1'")
            diann_version_id = "1.8.1"

        if diann_version_id != "1.8.1":
            raise ValueError(
                "DiannConvert is currently available only for DiaNN v1.8.1"
            )

        self.diann_report = FolderPath + pathdict["report"][0]
        self.exp_design = FolderPath + pathdict["exp_design"][0]
        self.fasta = FolderPath + pathdict["fasta"][0]
        self.mzml_info = pathdict["mzml_info"]
        self.msstats_path = (
            out
            + os.path.splitext(os.path.basename(self.exp_design))[0]
            + "_msstats_in.csv"
        )
        self.triqler_path = (
            out
            + os.path.splitext(os.path.basename(self.exp_design))[0]
            + "_triqler_in.tsv"
        )
        self.mztab_path = (
            out + os.path.splitext(os.path.basename(self.exp_design))[0] + "_out.mzTab"
        )

        self.dia_params = dia_params
        self.charge = charge
        self.missed_cleavages = missed_cleavages
        self.qvalue_threshold = qvalue_threshold
        self.diann_version = diann_version_id
        self.block_size = block_size
        self.processors = processors
        self.threads_per_processor = threads_per_processor

    def convert(self):
        """This function is designed to convert the DiaNN output into three standard formats: MSstats, Triqler and mzTab. These documents are
        used for quality control and downstream analysis.
        """
        log.info("Start loading...")
        remain_cols = [
            "File.Name",
            "Run",
            "Protein.Group",
            "Protein.Names",
            "Protein.Ids",
            "First.Protein.Description",
            "PG.MaxLFQ",
            "RT.Start",
            "Global.Q.Value",
            "Lib.Q.Value",
            "PEP",
            "Precursor.Normalised",
            "Precursor.Id",
            "Q.Value",
            "Modified.Sequence",
            "Stripped.Sequence",
            "Precursor.Charge",
            "Precursor.Quantity",
            "Global.PG.Q.Value",
        ]
        self.report = dd.read_table(
            self.diann_report,
            blocksize=self.block_size,
            assume_missing=True,
            usecols=remain_cols,
        )
        # filter based on qvalue parameter for downstream analysiss
        self.report = self.report[self.report["Q.Value"] < self.qvalue_threshold]

        with open(self.exp_design, "r") as f:
            data = f.readlines()
            empty_row = data.index("\n")
            f_table = [i.replace("\n", "").split("\t") for i in data[1:empty_row]]
            f_header = data[0].replace("\n", "").split("\t")
            f_table = pd.DataFrame(f_table, columns=f_header)
            f_table.loc[:, "run"] = f_table.apply(
                lambda x: os.path.splitext(os.path.basename(x["Spectra_Filepath"]))[0],
                axis=1,
            )

            s_table = [i.replace("\n", "").split("\t") for i in data[empty_row + 1 :]][
                1:
            ]
            s_header = data[empty_row + 1].replace("\n", "").split("\t")
            s_DataFrame = pd.DataFrame(s_table, columns=s_header)

        # Convert to MSstats
        out_msstats = self.report[
            [
                "Protein.Names",
                "Modified.Sequence",
                "Precursor.Charge",
                "Precursor.Quantity",
                "File.Name",
                "Run",
            ]
        ]
        out_msstats.columns = [
            "ProteinName",
            "PeptideSequence",
            "PrecursorCharge",
            "Intensity",
            "Reference",
            "Run",
        ]
        out_msstats = out_msstats[out_msstats["Intensity"] != 0]
        out_msstats["PeptideSequence"] = out_msstats.apply(
            lambda x: AASequence.fromString(x["PeptideSequence"]).toString(),
            axis=1,
            meta=(None, "str"),
        )
        out_msstats["FragmentIon"] = "NA"
        out_msstats["ProductCharge"] = "0"
        out_msstats["IsotopeLabelType"] = "L"
        out_msstats["Reference"] = out_msstats.apply(
            lambda x: os.path.basename(x["Reference"]), axis=1, meta=(None, "str")
        )

        out_msstats[["Fraction", "BioReplicate", "Condition"]] = out_msstats.apply(
            lambda x: self.query_expdesign_value(x["Run"], f_table, s_DataFrame),
            axis=1,
            result_type="expand",
            meta={0: "int64", 1: "int64", 2: "str"},
        )

        out_msstats.to_csv(self.msstats_path, single_file=True, index=False)

        # Convert to Triqler
        out_triqler = pd.DataFrame()
        out_triqler = out_msstats[
            [
                "ProteinName",
                "PeptideSequence",
                "PrecursorCharge",
                "Intensity",
                "Run",
                "Condition",
            ]
        ]
        del out_msstats
        out_triqler.columns = [
            "proteins",
            "peptide",
            "charge",
            "intensity",
            "run",
            "condition",
        ]
        out_triqler = out_triqler[out_triqler["intensity"] != 0]

        out_triqler["searchScore"] = self.report["Q.Value"]
        out_triqler["searchScore"] = 1 - out_triqler["searchScore"]
        out_triqler.to_csv(self.triqler_path, sep="\t", single_file=True, index=False)
        del out_triqler

        # Convert to mzTab
        if self.diann_version == "1.8.1":
            self.fasta_df = pd.DataFrame()
            entries = []
            f = FASTAFile()
            f.load(self.fasta, entries)
            line = 0
            for e in entries:
                self.fasta_df.loc[line, "id"] = e.identifier
                self.fasta_df.loc[line, "seq"] = e.sequence
                self.fasta_df.loc[line, "len"] = len(e.sequence)
                line += 1

            self.index_ref = f_table
            self.index_ref.loc[:, "ms_run"] = self.index_ref.apply(
                lambda x: x["Fraction_Group"], axis=1
            )
            self.index_ref.loc[:, "study_variable"] = self.index_ref.apply(
                lambda x: x["Sample"], axis=1
            )
            self.index_ref.loc[:, "ms_run"] = self.index_ref.loc[:, "ms_run"].astype(
                "int"
            )
            self.index_ref.loc[:, "study_variable"] = self.index_ref.loc[
                :, "study_variable"
            ].astype("int")

            self.max_study_variable = self.index_ref["study_variable"].max()
            self.max_ms_run = self.index_ref["ms_run"].max()
            self.report[["ms_run", "study_variable"]] = self.report.apply(
                lambda x: self.add_info(x["Run"]),
                axis=1,
                result_type="expand",
                meta={0: "int64", 1: "int64"},
            )
            self.report["Calculate.Precursor.Mz"] = self.report.apply(
                lambda x: calculate_mz(x["Stripped.Sequence"], x["Precursor.Charge"]),
                axis=1,
                meta=(None, "float64"),
            )

            PRH_cols = (
                [
                    "PRH",
                    "accession",
                    "description",
                    "database",
                    "taxid",
                    "species",
                    "database_version",
                    "search_engine",
                ]
                + [
                    f"protein_abundance_assay[{ms_run}]"
                    for ms_run in range(1, self.max_ms_run + 1)
                ]
                + np.array(
                    [
                        [
                            f"protein_abundance_study_variable[{study_variable}]",
                            f"protein_abundance_stdev_study_variable[{study_variable}]",
                            f"protein_abundance_std_error_study_variable[{study_variable}]",
                        ]
                        for study_variable in range(1, self.max_study_variable + 1)
                    ]
                )
                .flatten()
                .tolist()
                + [
                    "best_search_engine_score[1]",
                    "opt_global_result_type",
                    "ambiguity_members",
                    "modifications",
                    "protein_converage",
                    "opt_global_Posterior_Probability_score",
                    "opt_global_nr_found_peptides",
                    "opt_global_cv_PRIDE:0000303_decoy_hit",
                ]
            )

            PEH_cols = (
                [
                    "PEH",
                    "accession",
                    "database",
                    "sequence",
                    "charge",
                    "retention_time",
                    "mass_to_charge",
                    "database_version",
                    "search_engine",
                    "retention_time_window",
                    "modifications",
                    "unique",
                ]
                + [
                    f"search_engine_score[1]_ms_run[{ms_run}]"
                    for ms_run in range(1, self.max_ms_run + 1)
                ]
                + ["best_search_engine_score[1]"]
                + np.array(
                    [
                        [
                            f"peptide_abundance_study_variable[{study_variable}]",
                            f"peptide_abundance_stdev_study_variable[{study_variable}]",
                            f"peptide_abundance_std_error_study_variable[{study_variable}]",
                            f"opt_global_mass_to_charge_study_variable[{study_variable}]",
                            f"opt_global_retention_time_study_variable[{study_variable}]",
                        ]
                        for study_variable in range(1, self.max_study_variable + 1)
                    ]
                )
                .flatten()
                .tolist()
                + [
                    "opt_global_cv_MS:1000889_peptidoform_sequence",
                    "opt_global_cv_MS:1002217_decoy_peptide",
                    "opt_global_feature_id",
                    "opt_global_q-value",
                    "opt_global_SpecEValue_score",
                ]
            )

            PSH_cols = [
                "PSH",
                "sequence",
                "accession",
                "search_engine_score[1]",
                "retention_time",
                "charge",
                "calc_mass_to_charge",
                "exp_mass_to_charge",
                "PSM_ID",
                "unique",
                "database",
                "database_version",
                "search_engine",
                "pre",
                "post",
                "start",
                "end",
                "modifications",
                "spectra_ref",
                "opt_global_cv_MS:1000889_peptidoform_sequence",
                "opt_global_SpecEValue_score",
                "opt_global_q-value",
                "opt_global_q-value_score",
                "opt_global_spectrum_reference",
                "opt_global_cv_MS:1002217_decoy_peptide",
                "opt_global_feature_id",
                "opt_global_map_index",
            ]

            (MTD, self.database) = self.mztab_MTD()
            with open(
                self.mztab_path,
                "a",
                newline="",
            ) as f:
                MTD.to_csv(f, mode="a", sep="\t", index=0, header=0)
                f.write("\n")
                f.write("\t".join(PRH_cols) + "\n")

            self.report.groupby("Protein.Ids").apply(self.mztab_PRH_dask).compute(
                n_workers=self.processors, threads_per_worker=self.threads_per_processor
            )
            with open(
                self.mztab_path,
                "a",
                newline="",
            ) as f:
                f.write("\n")
                f.write("\t".join(PEH_cols) + "\n")

            self.report.groupby(["Modified.Sequence", "Precursor.Charge"]).apply(
                self.mztab_PEH_dask
            ).compute(
                n_workers=self.processors, threads_per_worker=self.threads_per_processor
            )

            with open(
                self.mztab_path,
                "a",
                newline="",
            ) as f:
                f.write("\n")
                f.write("\t".join(PSH_cols) + "\n")

            self.report.groupby("File.Name").apply(self.mztab_PSH_dask).compute(
                n_workers=self.processors, threads_per_worker=self.threads_per_processor
            )

    def mztab_MTD(self):
        """Construct MTD sub-table.

        :return: MTD sub-table
        :rtype: pandas.core.frame.DataFrame
        """
        dia_params_list = self.dia_params.split(";")
        dia_params_list = ["null" if i == "" else i for i in dia_params_list]
        FragmentMassTolerance = dia_params_list[0]
        FragmentMassToleranceUnit = dia_params_list[1]
        PrecursorMassTolerance = dia_params_list[2]
        PrecursorMassToleranceUnit = dia_params_list[3]
        Enzyme = dia_params_list[4]
        FixedModifications = dia_params_list[5]
        VariableModifications = dia_params_list[6]
        out_mztab_MTD = pd.DataFrame()
        out_mztab_MTD.loc[1, "mzTab-version"] = "1.0.0"
        out_mztab_MTD.loc[1, "mzTab-mode"] = "Summary"
        out_mztab_MTD.loc[1, "mzTab-type"] = "Quantification"
        out_mztab_MTD.loc[1, "title"] = "ConsensusMap export from OpenMS"
        out_mztab_MTD.loc[1, "description"] = "OpenMS export from consensusXML"
        out_mztab_MTD.loc[
            1, "protein_search_engine_score[1]"
        ] = "[, , DiaNN Global.PG.Q.Value, ]"
        out_mztab_MTD.loc[
            1, "peptide_search_engine_score[1]"
        ] = "[, , DiaNN Q.Value (minimum of the respective precursor q-values), ]"
        out_mztab_MTD.loc[
            1, "psm_search_engine_score[1]"
        ] = "[MS, MS:MS:1001869, protein-level q-value, ]"
        out_mztab_MTD.loc[
            1, "software[1]"
        ] = "[MS, MS:1003253, DiaNN, Release (v1.8.1)]"
        out_mztab_MTD.loc[1, "software[1]-setting[1]"] = self.fasta
        out_mztab_MTD.loc[1, "software[1]-setting[2]"] = "db_version:null"
        out_mztab_MTD.loc[1, "software[1]-setting[3]"] = (
            "fragment_mass_tolerance:" + FragmentMassTolerance
        )
        out_mztab_MTD.loc[1, "software[1]-setting[4]"] = (
            "fragment_mass_tolerance_unit:" + FragmentMassToleranceUnit
        )
        out_mztab_MTD.loc[1, "software[1]-setting[5]"] = (
            "precursor_mass_tolerance:" + PrecursorMassTolerance
        )
        out_mztab_MTD.loc[1, "software[1]-setting[6]"] = (
            "precursor_mass_tolerance_unit:" + PrecursorMassToleranceUnit
        )
        out_mztab_MTD.loc[1, "software[1]-setting[7]"] = "enzyme:" + Enzyme
        out_mztab_MTD.loc[1, "software[1]-setting[8]"] = "enzyme_term_specificity:full"
        out_mztab_MTD.loc[1, "software[1]-setting[9]"] = "charges:" + str(self.charge)
        out_mztab_MTD.loc[1, "software[1]-setting[10]"] = "missed_cleavages:" + str(
            self.missed_cleavages
        )
        out_mztab_MTD.loc[1, "software[1]-setting[11]"] = (
            "fixed_modifications:" + FixedModifications
        )
        out_mztab_MTD.loc[1, "software[1]-setting[12]"] = (
            "variable_modifications:" + VariableModifications
        )

        (fixed_mods, variable_mods, fix_flag, var_flag) = self.MTD_mod_info(
            FixedModifications, VariableModifications
        )
        if fix_flag == 1:
            for i in range(1, len(fixed_mods) + 1):
                out_mztab_MTD.loc[1, "fixed_mod[" + str(i) + "]"] = fixed_mods[i - 1][0]
                out_mztab_MTD.loc[1, "fixed_mod[" + str(i) + "]-site"] = fixed_mods[
                    i - 1
                ][1]
                out_mztab_MTD.loc[1, "fixed_mod[" + str(i) + "]-position"] = "Anywhere"
        else:
            out_mztab_MTD.loc[1, "fixed_mod[1]"] = fixed_mods[0]

        if var_flag == 1:
            for i in range(1, len(variable_mods) + 1):
                out_mztab_MTD.loc[1, "variable_mod[" + str(i) + "]"] = variable_mods[
                    i - 1
                ][0]
                out_mztab_MTD.loc[
                    1, "variable_mod[" + str(i) + "]-site"
                ] = variable_mods[i - 1][1]
                out_mztab_MTD.loc[
                    1, "variable_mod[" + str(i) + "]-position"
                ] = "Anywhere"
        else:
            out_mztab_MTD.loc[1, "variable_mod[1]"] = variable_mods[0]

        out_mztab_MTD.loc[
            1, "quantification_method"
        ] = "[MS, MS:1001834, LC-MS label-free quantitation analysis, ]"
        out_mztab_MTD.loc[1, "protein-quantification_unit"] = "[, , Abundance, ]"
        out_mztab_MTD.loc[1, "peptide-quantification_unit"] = "[, , Abundance, ]"

        for i in range(1, max(self.index_ref["ms_run"]) + 1):
            out_mztab_MTD.loc[
                1, "ms_run[" + str(i) + "]-format"
            ] = "[MS, MS:1000584, mzML file, ]"
            out_mztab_MTD.loc[1, "ms_run[" + str(i) + "]-location"] = (
                "file://"
                + self.index_ref[self.index_ref["ms_run"] == i][
                    "Spectra_Filepath"
                ].values[0]
            )
            out_mztab_MTD.loc[
                1, "ms_run[" + str(i) + "]-id_format"
            ] = "[MS, MS:1000777, spectrum identifier nativeID format, ]"

        for i in range(1, max(self.index_ref["ms_run"]) + 1):
            out_mztab_MTD.loc[
                1, "assay[" + str(i) + "]-quantification_reagent"
            ] = "[MS, MS:1002038, unlabeled sample, ]"
            out_mztab_MTD.loc[1, "assay[" + str(i) + "]-ms_run_ref"] = (
                "ms_run[" + str(i) + "]"
            )

        for i in range(1, max(self.index_ref["study_variable"]) + 1):
            study_variable = []
            for j in list(
                self.index_ref[self.index_ref["study_variable"] == i]["ms_run"].values
            ):
                study_variable.append("assay[" + str(j) + "]")
            out_mztab_MTD.loc[
                1, "study_variable[" + str(i) + "]-assay_refs"
            ] = ",".join(study_variable)
            out_mztab_MTD.loc[
                1, "study_variable[" + str(i) + "]-description"
            ] = "no description given"

        # Transpose out_mztab_MTD
        out_mztab_MTD_T = pd.DataFrame(
            out_mztab_MTD.values.T,
            index=list(out_mztab_MTD.columns),
            columns=out_mztab_MTD.index,
        )
        out_mztab_MTD_T.insert(0, "key", out_mztab_MTD_T.index)
        out_mztab_MTD_T.insert(0, "tag", "MTD")
        database = os.path.basename(self.fasta.split(".")[-2])

        return out_mztab_MTD_T, database

    def mztab_PRH_dask(self, pro_group):
        """Group the main report into a single task by protein and write the mztab line by line

        :param pro_group: Data frame of a single protein
        :type pro_group: dataframe
        """
        null_col = [
            "taxid",
            "species",
            "database_version",
            "search_engine",
            "opt_global_Posterior_Probability_score",
            "opt_global_nr_found_peptides",
            "opt_global_cv_PRIDE:0000303_decoy_hit",
        ]
        protein_df = pd.DataFrame()
        protein_df.loc[1, "PRH"] = "PRT"
        protein_df.loc[1, "accession"] = pro_group["Protein.Ids"].values[0]
        if protein_df.loc[1, "accession"] == "foo":
            return
        protein_df.loc[1, "description"] = pro_group[
            "First.Protein.Description"
        ].values[0]
        protein_df.loc[1, "database"] = self.database
        for col in null_col[0:4]:
            protein_df.loc[1, col] = "null"

        for ms_run in range(1, self.max_ms_run + 1):
            target = pro_group[pro_group["ms_run"] == ms_run]
            protein_df.loc[1, f"protein_abundance_assay[{ms_run}]"] = (
                target["PG.MaxLFQ"].values[0] if len(target) > 0 else "null"
            )

        for study_variable in range(1, self.max_study_variable + 1):
            target = pro_group[pro_group["study_variable"] == study_variable]
            protein_df.loc[1, f"protein_abundance_study_variable[{study_variable}]"] = (
                target["PG.MaxLFQ"].values[0] if len(target) > 0 else "null"
            )
            protein_df.loc[
                1, f"protein_abundance_stdev_study_variable[{study_variable}]"
            ] = "null"
            protein_df.loc[
                1, f"protein_abundance_std_error_study_variable[{study_variable}]"
            ] = "null"

        protein_df.loc[1, "best_search_engine_score[1]"] = pro_group[
            "Global.PG.Q.Value"
        ].min()

        protein_df.loc[1, "opt_global_result_type"] = (
            "indistinguishable_protein_group"
            if ";" in pro_group["Protein.Ids"].values[0]
            else "single_protein"
        )
        protein_df.loc[1, "ambiguity_members"] = (
            protein_df.loc[1, "accession"]
            if protein_df.loc[1, "opt_global_result_type"]
            == "indistinguishable_protein_group"
            else "null"
        )
        protein_details_df = protein_df[
            protein_df["opt_global_result_type"] == "indistinguishable_protein_group"
        ]
        if len(protein_details_df) > 0:
            prh_series = (
                protein_details_df["accession"]
                .str.split(";", expand=True)
                .stack()
                .reset_index(level=1, drop=True)
            )
            prh_series.name = "accession"
            protein_details_df = (
                protein_details_df.drop("accession", axis=1)
                .join(prh_series)
                .reset_index()
                .drop(columns="index")
            )
            protein_details_df["opt_global_result_type"] = "protein_details"
            protein_df = pd.concat([protein_df, protein_details_df]).reset_index(
                drop=True
            )

        protein_df["accession"] = protein_df.apply(
            lambda x: x["accession"].split(";")[0],
            axis=1,
        )
        protein_df[["modifications", "protein_converage"]] = protein_df.apply(
            lambda x: self.protein_calculate(
                list(set(pro_group["Stripped.Sequence"])),
                list(set(pro_group["Modified.Sequence"])),
                x["accession"],
            ),
            axis=1,
            result_type="expand",
        )
        for col in null_col[4:]:
            protein_df[col] = "null"
        protein_df.to_csv(self.mztab_path, mode="a", sep="\t", header=0, index=0)

    def mztab_PEH_dask(self, pep_group):
        """The main report is grouped into single precursors by sequence, modification, and charge and written to mztab line by line

        :param pep_group: Data frame of a single precursor
        :type pep_group: dataframe
        """
        if "foo" in pep_group.values:
            return
        null_col = [
            "database_version",
            "search_engine",
            "retention_time_window",
            "opt_global_feature_id",
        ]
        peptide_df = pd.DataFrame()
        peptide_df.loc[1, "PEH"] = "PEP"
        peptide_df.loc[1, "accession"] = pep_group["Protein.Ids"].values[0]
        peptide_df.loc[1, "database"] = self.database
        peptide_df.loc[1, "sequence"] = pep_group["Stripped.Sequence"].values[0]
        peptide_df.loc[1, "charge"] = pep_group["Precursor.Charge"].values[0]
        peptide_df.loc[1, "retention_time"] = pep_group["RT.Start"].mean()
        peptide_df.loc[1, "mass_to_charge"] = pep_group["Calculate.Precursor.Mz"].mean()
        for col in null_col[:-1]:
            peptide_df.loc[1, col] = "null"
        peptide_df.loc[1, "modifications"] = self.find_modification(
            pep_group["Modified.Sequence"].values[0]
        )
        peptide_df.loc[1, "unique"] = 0 if ";" in peptide_df.loc[1, "accession"] else 1

        for ms_run in range(1, self.max_ms_run + 1):
            target = pep_group[pep_group["ms_run"] == ms_run]
            peptide_df.loc[1, f"search_engine_score[1]_ms_run[{ms_run}]"] = (
                target["Q.Value"].min() if len(target) > 0 else "null"
            )
        peptide_df.loc[1, "best_search_engine_score[1]"] = pep_group["Q.Value"].min()
        for study_variable in range(1, self.max_study_variable + 1):
            target = pep_group[pep_group["study_variable"] == study_variable]
            peptide_df.loc[1, f"peptide_abundance_study_variable[{study_variable}]"] = (
                target["Precursor.Normalised"].mean() if len(target) > 0 else "null"
            )
            null_study = [
                f"peptide_abundance_stdev_study_variable[{study_variable}]",
                f"peptide_abundance_std_error_study_variable[{study_variable}]",
                f"opt_global_mass_to_charge_study_variable[{study_variable}]",
                f"opt_global_retention_time_study_variable[{study_variable}]",
            ]
            peptide_df.loc[1, null_study] = "null"
        peptide_df.loc[
            1, "opt_global_cv_MS:1000889_peptidoform_sequence"
        ] = AASequence.fromString(pep_group["Modified.Sequence"].values[0]).toString()
        peptide_df.loc[1, "opt_global_cv_MS:1002217_decoy_peptide"] = "0"
        peptide_df.loc[1, null_col[-1]] = "null"
        peptide_df.loc[1, "opt_global_q-value"] = pep_group["Global.Q.Value"].values[0]
        peptide_df.loc[1, "opt_global_SpecEValue_score"] = pep_group[
            "Lib.Q.Value"
        ].values[0]
        peptide_df.to_csv(self.mztab_path, sep="\t", mode="a", header=0, index=0)

    def mztab_PSH_dask(self, psm_group):
        """Group the main report into different mass spectrum files and write it to mztab

        :param psm_group: psm from the same ms_run
        :type psm_group: datafrme
        """
        if "foo" in psm_group.values:
            return
        psm_df = pd.DataFrame()
        file = self.folder + psm_group["Run"].values[0] + "_mzml_info.tsv"
        target = pd.read_csv(file, sep="\t")
        psm_group.sort_values(by="RT.Start", inplace=True)
        target = target[["Retention_Time", "SpectrumID", "Exp_Mass_To_Charge"]]
        target.columns = [
            "RT.Start",
            "opt_global_spectrum_reference",
            "exp_mass_to_charge",
        ]
        # TODO seconds returned from precursor.getRT()
        target.loc[:, "RT.Start"] = target.apply(lambda x: x["RT.Start"] / 60, axis=1)
        psm_df = pd.concat(
            [
                psm_df,
                pd.merge_asof(psm_group, target, on="RT.Start", direction="nearest"),
            ]
        )
        ## Score at PSM level: Q.Value
        psm_df = psm_df[
            [
                "Stripped.Sequence",
                "Protein.Ids",
                "Q.Value",
                "RT.Start",
                "Precursor.Charge",
                "Calculate.Precursor.Mz",
                "exp_mass_to_charge",
                "Modified.Sequence",
                "PEP",
                "Global.Q.Value",
                "Global.Q.Value",
                "opt_global_spectrum_reference",
                "ms_run",
            ]
        ]
        psm_df.columns = [
            "sequence",
            "accession",
            "search_engine_score[1]",
            "retention_time",
            "charge",
            "calc_mass_to_charge",
            "exp_mass_to_charge",
            "opt_global_cv_MS:1000889_peptidoform_sequence",
            "opt_global_SpecEValue_score",
            "opt_global_q-value",
            "opt_global_q-value_score",
            "opt_global_spectrum_reference",
            "ms_run",
        ]

        psm_df.loc[:, "opt_global_cv_MS:1002217_decoy_peptide"] = "0"
        psm_df.loc[:, "PSM_ID"] = "null"
        psm_df.loc[:, "unique"] = psm_df.apply(
            lambda x: "0" if ";" in str(x["accession"]) else "1",
            axis=1,
        )
        psm_df.loc[:, "database"] = self.database

        null_col = [
            "database_version",
            "search_engine",
            "pre",
            "post",
            "start",
            "end",
            "opt_global_feature_id",
            "opt_global_map_index",
        ]
        for i in null_col:
            psm_df.loc[:, i] = "null"

        psm_df.loc[:, "modifications"] = psm_df.apply(
            lambda x: self.find_modification(
                x["opt_global_cv_MS:1000889_peptidoform_sequence"]
            ),
            axis=1,
        )

        psm_df.loc[:, "spectra_ref"] = psm_df.apply(
            lambda x: "ms_run[{}]:".format(x["ms_run"])
            + x["opt_global_spectrum_reference"],
            axis=1,
        )

        psm_df.loc[:, "opt_global_cv_MS:1000889_peptidoform_sequence"] = psm_df.apply(
            lambda x: AASequence.fromString(
                x["opt_global_cv_MS:1000889_peptidoform_sequence"]
            ).toString(),
            axis=1,
        )

        psm_df["PSH"] = "PSM"
        psm_df.drop(["ms_run"], axis=1, inplace=True)
        psm_df.fillna("null")
        new_cols = (
            ["PSH"]
            + [
                col
                for col in psm_df.columns
                if not col.startswith("opt_") and col != "PSH"
            ]
            + [col for col in psm_df.columns if col.startswith("opt_")]
        )
        psm_df = psm_df[new_cols]
        psm_df.to_csv(self.mztab_path, mode="a", sep="\t", header=0, index=0)

    def find_modification(self, peptide):
        """Identify the modification site based on the peptide containing modifications.

        :param peptide: Sequences of peptides
        :type peptide: str
        :return: Modification sites
        :rtype: str
        """
        peptide = str(peptide)
        pattern = re.compile(r"\((.*?)\)")
        original_mods = re.findall(pattern, peptide)
        peptide = re.sub(r"\(.*?\)", ".", peptide)
        position = [i.start() for i in re.finditer(r"\.", peptide)]
        for j in range(1, len(position)):
            position[j] -= j

        for k in range(0, len(original_mods)):
            original_mods[k] = str(position[k]) + "-" + original_mods[k].upper()

        original_mods = (
            ",".join(str(i) for i in original_mods)
            if len(original_mods) > 0
            else "null"
        )

        return original_mods

    def protein_calculate(self, peptide_list, modified_seq_list, target):
        """Calculate protein coverage and find locations of modifications in protein sequence

        :param peptide_list: List contains peptides from one protein
        :type peptide_list: list
        :param modified_seq_list: List contains modified peptides from one protein
        :type modified_seq_list: list
        :param target: Protein accession
        :type target: str
        :return: protein coverage and modifications
        :rtype: tuple
        """

        resultlist = []
        resultdict = dict()
        ref = self.fasta_df[self.fasta_df["id"].str.contains(target)]["seq"].values[0]

        for i in peptide_list:
            resultdict[i] = list()
            result = re.finditer(i, ref)
            if result:
                for j in result:
                    pos = [j.span()[0], j.span()[1] - 1]
                    resultlist.append(pos)
                resultdict[i].append(pos[0])

        # Sort and merge the interval list
        resultlist.sort()
        l, r = 0, 1
        while r < len(resultlist):
            x1, y1 = resultlist[l][0], resultlist[l][1]
            x2, y2 = resultlist[r][0], resultlist[r][1]
            if x2 > y1:
                l += 1
                r += 1
            else:
                resultlist[l] = [x1, max(y1, y2)]
                resultlist.pop(r)

        coverage_length = np.array([i[1] - i[0] + 1 for i in resultlist]).sum()
        protein_coverage = format(coverage_length / len(ref), ".3f")

        mod_list = list()
        for i in modified_seq_list:
            if "(" in i:
                pattern = re.compile(r"\((.*?)\)")
                original_mods = pattern.findall(i)
                original_seq = re.sub("\\(.*?\\)", "", i)
                peptide = re.sub(r"\(.*?\)", ".", i)
                position = [i.start() for i in re.finditer(r"\.", peptide)]
                for j in range(1, len(position)):
                    position[j] -= j

                for k in range(0, len(original_mods)):
                    for l in resultdict[original_seq]:
                        mod_list.append(
                            str(position[k] + l - 1) + "-" + original_mods[k].upper()
                        )
            else:
                continue
        modifications = ";".join(list(set(mod_list))) if len(mod_list) > 0 else "null"

        return modifications, protein_coverage

    def query_expdesign_value(self, reference, f_table, s_table):
        """By matching the "Run" column in f_table or the "Sample" column in s_table, this function returns a tuple containing Fraction,
        BioReplicate and Condition.

        :param reference: The value of "Run" column in out_msstats
        :type reference: str
        :param f_table: A table contains experiment settings(search engine settings etc.)
        :type f_table: pandas.core.frame.DataFrame
        :param s_table: A table contains experimental design
        :type s_table: pandas.core.frame.DataFrame
        :return: A tuple contains Fraction, BioReplicate and Condition
        :rtype: tuple
        """
        query_reference = f_table[f_table["run"] == reference]
        Fraction = query_reference["Fraction"].values[0]
        row = s_table[s_table["Sample"] == query_reference["Sample"].values[0]]
        BioReplicate = row["MSstats_BioReplicate"].values[0]
        Condition = row["MSstats_Condition"].values[0]

        return Fraction, BioReplicate, Condition

    def MTD_mod_info(self, fix_mod, var_mod):
        """Convert fixed and variable modifications to the format required by the MTD sub-table.

        :param fix_mod: Fixed modifications from DIA parameter list
        :type fix_mod: str
        :param var_mod: Variable modifications from DIA parameter list
        :type var_mod: str
        :return: A tuple contains fixed and variable modifications, and flags indicating whether they are null
        :rtype: tuple
        """
        var_ptm = []
        fix_ptm = []
        mods_db = ModificationsDB()

        if fix_mod != "null":
            fix_flag = 1
            for mod in fix_mod.split(","):
                mod_obj = mods_db.getModification(mod)
                mod_name = mod_obj.getId()
                mod_accession = mod_obj.getUniModAccession()
                site = mod_obj.getOrigin()
                fix_ptm.append(
                    (
                        "[UNIMOD, " + mod_accession.upper() + ", " + mod_name + ", ]",
                        site,
                    )
                )
        else:
            fix_flag = 0
            fix_ptm.append("[MS, MS:1002453, No fixed modifications searched, ]")

        if var_mod != "null":
            var_flag = 1
            for mod in var_mod.split(","):
                mod_obj = mods_db.getModification(mod)
                mod_name = mod_obj.getId()
                mod_accession = mod_obj.getUniModAccession()
                site = mod_obj.getOrigin()
                var_ptm.append(
                    (
                        "[UNIMOD, " + mod_accession.upper() + ", " + mod_name + ", ]",
                        site,
                    )
                )
        else:
            var_flag = 0
            var_ptm.append("[MS, MS:1002454, No variable modifications searched, ]")

        return fix_ptm, var_ptm, fix_flag, var_flag

    def add_info(self, target):
        """On the basis of f_table, two columns "ms_run" and "study_variable" are added for matching.

        :param target: The value of "Run" column in f_table
        :type target: str
        :return: A tuple contains ms_run and study_variable
        :rtype: tuple
        """
        match = self.index_ref[self.index_ref["run"] == target]
        ms_run = match["ms_run"].values[0]
        study_variable = match["study_variable"].values[0]

        return ms_run, study_variable


class MzTabMerge:
    """
    The MzTabMerge class is used to merge two mzTabs into one, it streams and merges the original mztab, but
    requires that the mztab to be merged must be cacheable.

    :param mztab1: Path to the original mztab
    :type mztab1: str
    :param mztab2: Path to the mztab to be merged
    :type mztab2: str
    :param single_cache_size: Single cache size, default 500*1024*1024
    :type single_cache_size: int
    :param result_folder: Folder to result files, default './'
    :type result_folder: str
    """

    def __init__(
        self, mztab1, mztab2, single_cache_size=800 * 1024 * 1024, result_folder="./"
    ):
        if os.path.getsize(mztab2) > single_cache_size:
            log.error(
                "The size of mztab2 must be able to be cached, please increase single cache size!"
            )
            return
        self.mztab_path1, self.mztab_path2 = mztab1, mztab2
        log.info(f"Start loading {self.mztab_path2}")
        (meta_dict2, pro2, pep2, psm2) = MzTabPy(
            mztab_path=self.mztab_path2, section_tag=True, df_index=False
        ).loadmzTab()
        self.merge_dict, self.meta_dict = {
            "meta": meta_dict2,
            "pro": pro2,
            "pep": pep2,
            "psm": psm2,
        }, dict()
        self.psmid = 0
        self.chunk_size = single_cache_size
        self.file_bytes = os.path.getsize(mztab1)
        self.file_path = result_folder + "merge.mztab"
        if self.file_bytes < self.chunk_size:
            log.info(f"Start loading {self.mztab_path1}")
            (meta_dict1, pro1, pep1, psm1) = MzTabPy(
                mztab_path=self.mztab_path1, section_tag=True, df_index=False
            ).loadmzTab()
            merge_metadict = self.meta_merge(meta_dict1, meta_dict2)
            merge_meta = pd.DataFrame(merge_metadict, index=[0])
            meta = pd.DataFrame(
                merge_meta.values.T,
                index=list(merge_meta.columns),
                columns=merge_meta.index,
            )
            meta.insert(0, "key", meta.index)
            meta.insert(0, "tag", "MTD")
            merge_pro = self.protein_merge(pro1, pro2)
            merge_pep = self.peptide_merge(pep1, pep2)
            map_list = psm1.loc[:, "opt_global_map_index"].to_list()
            max_map_index = max([int(i) for i in map_list if not "null" in map_list])

            def reindex(x):
                return int(x) + max_map_index

            psm2.loc[:, "opt_global_map_index"] = psm2.loc[
                :, "opt_global_map_index"
            ].apply(reindex)
            merge_psm = self.psm_merge(psm1, psm2)

            meta.loc["", :] = ""
            merge_pro.loc[len(merge_pro) + 1, :] = ""
            merge_pep.loc[len(merge_pep) + 1, :] = ""
            with open(self.file_path, "a", newline="") as f:
                meta.to_csv(f, mode="a", sep="\t", index=False, header=False)
                merge_pro.to_csv(f, mode="a", sep="\t", index=False, header=True)
                merge_pep.to_csv(f, mode="a", sep="\t", index=False, header=True)
                merge_psm.to_csv(f, mode="a", sep="\t", index=False, header=True)

        else:
            log.warning(
                "No subtable is cached because the single cache is smaller than file bytes!"
            )
            (
                self.pro_cache,
                self.pro_merge_status,
                self.pep_cache,
                self.pep_merge_status,
            ) = (False, False, False, False)
            self.chunkDict = {"meta": [], "pro": [], "pep": [], "psm": []}
            self.chunkID = 0
            self.row_remains = ""
            log.info(f"Start loading {self.mztab_path1}")
            with open(self.mztab_path1, "r") as f:
                while True:
                    chunk = self.row_remains + f.read(self.chunk_size)
                    if not chunk:
                        break
                    self.loadMergeChunk(chunk)
            f.close()

    def select_non_missing(self, values):
        if values[0] == np.nan or values[0] == "null":
            values[0] = values[1]
        else:
            pass

        return values

    def max_or_null(self, x, y):
        data_list = [i for i in [x, y] if i not in [np.nan, None, "", "null"]]
        if len(data_list) == 0:
            return None
        elif len((data_list)) == 1:
            return data_list[0]
        else:
            return max(float(x), float(y))

    def min_or_null(self, x, y):
        data_list = [i for i in [x, y] if i not in [np.nan, None, "", "null"]]
        if len(data_list) == 0:
            return None
        elif len((data_list)) == 1:
            return data_list[0]
        else:
            return min(float(x), float(y))

    def paste_or_null(self, x, y):
        paste_list = [i for i in [x, y] if i not in [np.nan, None, "", "null"]]
        if len(paste_list) == 0:
            return None
        elif len((paste_list)) == 1:
            return paste_list[0]
        else:
            return "|".join([str(x), str(y)])

    def parse_extral_peps(self, x, y, field):
        if x == "PEP":
            if y != "PEP":
                cols = np.array(field)[1, :]
            else:
                pass
        else:
            if y == "PEP":
                cols = np.array(field)[:, 1]
            else:
                log.warning("Some errors occur in peptide merging!")
        return cols

    def pro_recalculate(self, merged, fields):
        # if not set(sum(fields, [])) <= set(merged.columns):
        #     return(merged)
        # recalculate
        x_cols = [col for col in merged.columns if col.endswith("_x")]
        y_cols = [col for col in merged.columns if col.endswith("_y")]
        for i in fields:
            merged[i] = merged.apply(lambda x: self.select_non_missing(x[i]), axis=1)
        if "best_search_engine_score[1]_x" in merged.columns:
            merged.loc[:, "best_search_engine_score[1]_x"] = merged.apply(
                lambda x: self.min_or_null(
                    x["best_search_engine_score[1]_x"],
                    x["best_search_engine_score[1]_y"],
                ),
                axis=1,
            )

        if "protein_coverage_x" in merged.columns:
            merged.loc[:, "protein_coverage_x"] = merged.apply(
                lambda x: self.max_or_null(
                    x["protein_coverage_x"], x["protein_coverage_y"]
                ),
                axis=1,
            )

        if "opt_global_nr_found_peptides_x" in merged.columns:
            merged.loc[:, "opt_global_nr_found_peptides_x"] = merged.apply(
                lambda x: self.max_or_null(
                    x["opt_global_nr_found_peptides_x"],
                    x["opt_global_nr_found_peptides_y"],
                ),
                axis=1,
            )

        if "opt_global_Posterior_Probability_score_x" in merged.columns:
            merged.loc[:, "opt_global_Posterior_Probability_score_x"] = merged.apply(
                lambda x: self.min_or_null(
                    x["opt_global_Posterior_Probability_score_x"],
                    x["opt_global_Posterior_Probability_score_y"],
                ),
                axis=1,
            )

        # strip columns from pro2
        pro2_variable_cols = np.array(fields)[:, 1]
        combine_columns = [
            item for item in merged.columns if item not in pro2_variable_cols
        ]
        combine_columns = [x.strip() for x in combine_columns if x.strip() != ""]
        merged = merged[combine_columns]

        merged_cols = dict()
        for i in merged.columns:
            if i.endswith("_x"):
                merged_cols.update({i: i.rstrip("_x")})

        merged = merged.rename(columns=merged_cols)
        return merged

    def pep_recalculate(self, merged, fields, mz_study_variable):
        # recalculate
        for i in fields:
            merged[i] = merged.apply(lambda x: self.select_non_missing(x[i]), axis=1)

        if "best_search_engine_score[1]_x" in merged.columns:
            merged.loc[:, "best_search_engine_score[1]_x"] = merged.apply(
                lambda x: self.min_or_null(
                    x["best_search_engine_score[1]_x"],
                    x["best_search_engine_score[1]_y"],
                ),
                axis=1,
            )

        if "opt_global_Posterior_Error_Probability_score_x" in merged.columns:
            merged.loc[
                :, "opt_global_Posterior_Error_Probability_score_x"
            ] = merged.apply(
                lambda x: self.min_or_null(
                    x["opt_global_Posterior_Error_Probability_score_x"],
                    x["opt_global_Posterior_Error_Probability_score_y"],
                ),
                axis=1,
            )

        if "opt_global_q-value_x" in merged.columns:
            merged.loc[:, "opt_global_q-value_x"] = merged.apply(
                lambda x: self.min_or_null(
                    x["opt_global_q-value_x"], x["opt_global_q-value_y"]
                ),
                axis=1,
            )

        if "spectra_ref_x" in merged.columns:
            merged["spectra_ref_x"] = merged.apply(
                lambda x: self.paste_or_null(x["spectra_ref_x"], x["spectra_ref_y"]),
                axis=1,
            )

        if "mass_to_charge_x" in merged.columns:
            if len(mz_study_variable) > 0:
                for i in mz_study_variable:
                    merged.loc[:, i] = merged.apply(
                        lambda x: float(x[i]) if not x[i] in ["null", None] else None,
                        axis=1,
                    )
                merged.loc[:, "mass_to_charge_x"] = merged[mz_study_variable].mean()

        ## TODO apply the extral peptides (less or more)
        # Ignore the missing '_y' section because the '_x' section is always truncated as the final result
        unconsensus_peps = merged[merged["PEH_x"] != "PEP"]
        merged = merged[merged["PEH_x"] == "PEP"]
        unconsensus_cols = [
            col for col in unconsensus_peps.columns if not col.endswith("_x")
        ]
        unconsensus = unconsensus_peps[unconsensus_cols]
        unconsensus.columns = [i.rstrip("_y") for i in unconsensus.columns]

        pep2_variable_cols = np.array(fields)[:, 1]
        combine_columns = [
            item for item in merged.columns if item not in pep2_variable_cols
        ]
        merged = merged[combine_columns]
        merged.columns = [i.rstrip("_x") for i in merged.columns]
        merged = pd.concat([merged, unconsensus[merged.columns]])

        return merged

    def fields_count(self, df, variable_fields):
        count = dict.fromkeys(variable_fields, 0)
        for i in df.columns:
            for j in variable_fields:
                if i.startswith(j):
                    count[j] += 1
        for i in list(count.keys()):
            if count[i] == 0:
                del count[i]
                variable_fields.remove(i)
        return count, variable_fields

    def adjust_spectra_ref(self, value):
        if str(value).startswith("ms_run"):
            ms_run = re.findall(r"[[](.*?)[]]", value)[0]
            if ms_run in self.filed_change_dict["ms_run"].keys():
                new_spectra_ref = value.replace(
                    ms_run.join(["[", "]"]),
                    self.filed_change_dict["ms_run"][ms_run].join(["[", "]"]),
                )
                return new_spectra_ref
        else:
            pass

    def loadMergeChunk(self, chunk):
        """This function processes data blocks from streams, classifies them into metadata, protein, peptide and
            PSM, and caches them as dataframes.

        :param chunk: A block of data
        :type chunk: str
        """
        data = chunk.split("\n")
        if not chunk.endswith("\n"):
            self.row_remains = data[-1]
            data = data[0:-1]
        else:
            self.row_remains = ""
        meta_dict = dict()
        pro_data, pep_data, psm_data = [], [], []
        for i in data:
            i = i.strip("\n")
            row_list = i.split("\t")
            if row_list[0] == "MTD":
                meta_dict[row_list[1]] = row_list[2]
            elif row_list[0] == "PRH":
                log.info(
                    "Metadata processing complete. Start processing the protein subtable..."
                )
                self.pro_cols = row_list
                self.meta_cache = True
            elif row_list[0] == "PRT":
                pro_data.append(row_list)
            elif row_list[0] == "PEH":
                log.info(
                    "Protein subtable processing complete. Start processing the peptide subtable..."
                )
                self.pep_cols = row_list
                self.pro_cache = True
            elif row_list[0] == "PEP":
                pep_data.append(row_list)
            elif row_list[0] == "PSH":
                log.info(
                    "Peptide subtable processing complete. Start processing the psm subtable..."
                )
                self.psm_cols = row_list
                self.pep_cache = True
            elif row_list[0] == "PSM":
                psm_data.append(row_list)
            else:
                continue

        for i in meta_dict.keys():
            meta_dict.update(
                {
                    i: ";".join(list(meta_dict[i]))
                    if type(meta_dict[i]) == tuple
                    else meta_dict[i]
                }
            )
        self.meta_dict.update(self.meta_merge(meta_dict, self.merge_dict["meta"]))

        if self.meta_cache and self.chunkID == 0:
            merge_meta = pd.DataFrame(self.meta_dict, index=[0])
            meta = pd.DataFrame(
                merge_meta.values.T,
                index=list(merge_meta.columns),
                columns=merge_meta.index,
            )
            meta.insert(0, "key", meta.index)
            meta.insert(0, "tag", "MTD")
            meta.loc["", :] = ""
            with open(self.file_path, "a", newline="") as f:
                meta.to_csv(f, mode="a", sep="\t", index=False, header=False)

        if len(pro_data) > 0:
            self.chunkDict["pro"].append(self.chunkID)
            protein_table = pd.DataFrame(pro_data, columns=self.pro_cols, dtype="str")
            if not self.pro_cache:
                self.pro_differ_chunks = True
            sub_pro = self.merge_dict["pro"][
                self.merge_dict["pro"]["accession"].isin(protein_table.columns)
            ]
            if len(sub_pro) == 0:
                sub_pro = self.merge_dict["pro"]
                self.pro_merge_status = True
            else:
                self.merge_dict["pro"] = self.merge_dict["pro"][
                    ~self.merge_dict["pro"]["accession"].isin(protein_table.columns)
                ]
            merge_pro = self.protein_merge(protein_table, sub_pro)

            if self.pro_merge_status:
                merge_pro.loc[len(merge_pro) + 1, :] = ""

            merge_pro.fillna("null", inplace=True)
            if len(self.chunkDict["pro"]) == 1:
                pro_header = True
            else:
                pro_header = False
            with open(self.file_path, "a", newline="") as f:
                merge_pro.to_csv(f, mode="a", sep="\t", index=False, header=pro_header)

            del protein_table, sub_pro, merge_pro

        if len(pep_data) > 0:
            self.chunkDict["pep"].append(self.chunkID)
            peptide_table = pd.DataFrame(pep_data, columns=self.pep_cols, dtype="str")
            if not self.pep_cache:
                self.pep_differ_chunks = True
            sub_pep = self.merge_dict["pep"][
                self.merge_dict["pep"]["accession"].isin(peptide_table.columns)
            ]
            if len(sub_pep) == 0:
                sub_pep = self.merge_dict["pep"]
                self.pep_merge_status = True
            else:
                self.merge_dict["pep"] = self.merge_dict["pep"][
                    ~self.merge_dict["pep"]["accession"].isin(peptide_table.columns)
                ]
            merge_pep = self.peptide_merge(peptide_table, sub_pep)

            if self.pep_merge_status:
                merge_pep.loc[len(merge_pep) + 1, :] = ""

            merge_pep.fillna("null", inplace=True)
            if len(self.chunkDict["pep"]) == 1:
                pep_header = True
            else:
                pep_header = False
            with open(self.file_path, "a", newline="") as f:
                merge_pep.to_csv(f, mode="a", sep="\t", index=False, header=pep_header)

            del peptide_table, sub_pep, merge_pep

        if len(psm_data) > 0:
            self.chunkDict["psm"].append(self.chunkID)
            psm_table = pd.DataFrame(psm_data, columns=self.psm_cols, dtype="str")
            if "opt_global_cv_MS:1002217_decoy_peptide" not in psm_table.columns.values:
                psm_table["opt_global_cv_MS:1002217_decoy_peptide"] = psm_table.apply(
                    lambda x: 1
                    if any(i in x["accession"] for i in ["DECOY_", "decoy"])
                    else 0,
                    axis=1,
                )

            sub_psm = self.merge_dict["psm"][
                self.merge_dict["psm"]["accession"].isin(psm_table.columns)
            ]
            if len(sub_psm) == 0:
                sub_psm = self.merge_dict["psm"]
            else:
                self.merge_dict["psm"] = self.merge_dict["psm"][
                    ~self.merge_dict["psm"]["accession"].isin(psm_table.columns)
                ]

            merge_psm = self.psm_merge(psm_table, sub_psm)
            merge_psm.reset_index(drop=True, inplace=True)
            self.psmid += merge_psm["PSM_ID"].max() + 1
            if len(self.chunkDict["psm"]) == 1:
                psm_header = True
            else:
                merge_psm.loc[:, "PSM_ID"] = merge_psm.loc[:, "PSM_ID"] + self.psmid
                psm_header = False
            ## TODO abandon "opt_global_map_index" if mztab can't be cached
            merge_psm.loc[:, "opt_global_map_index"] = "null"
            merge_psm.fillna("null", inplace=True)
            with open(self.file_path, "a", newline="") as f:
                merge_psm.to_csv(f, mode="a", sep="\t", index=False, header=psm_header)

            del psm_table, sub_psm, merge_psm

        self.chunkID += 1

    def meta_merge(self, meta1, meta2):
        # Warning values clash and missing keys
        prefix_list = ["ms_run", "assay", "study_variable"]
        meta_fix, meta_ms_run, meta_assay, meta_study_variable = (
            OrderedDict(),
            OrderedDict(),
            OrderedDict(),
            OrderedDict(),
        )
        self.filed_change_dict = dict.fromkeys(prefix_list, dict())
        miss = [item for item in meta1.keys() if item not in meta2.keys()]
        for i in miss:
            log.warning("Metadata entry for field '" + i + "' is missing")

        fix_columns = list()
        for i in meta1.keys():
            if i.startswith("ms_run"):
                meta_ms_run.update({i: meta1[i]})
            elif i.startswith("assay"):
                meta_assay.update({i: meta1[i]})
            elif i.startswith("study_variable"):
                meta_study_variable.update({i: meta1[i]})
            else:
                fix_columns.append(i)
        diff_columns = [item for item in fix_columns if meta1[item] != meta2[item]]
        for i in diff_columns:
            log.warning("Metadata entries for field '" + i + "' differ")

        for i in fix_columns:
            meta_fix.update({i: meta1[i]})

        variable_fields = {
            "ms_run": "location",
            "assay": "ms_run_ref",
            "study_variable": "assay_refs",
        }

        def fields_count(meta, variable_fields):
            count = dict.fromkeys(list(variable_fields.keys()), 0)
            for i in list(meta.keys()):
                for j in variable_fields:
                    if i.startswith(j) and i.endswith(variable_fields[j]):
                        count[j] += 1
            return count

        fields_count1 = fields_count(meta1, variable_fields)

        for i in list(meta2.keys()):
            if i.startswith("ms_run"):
                index = re.findall(r"[[](.*?)[]]", i)[0]
                new_run = str(int(index) + fields_count1["ms_run"])
                meta_ms_run.update({i.replace(index, new_run): meta2[i]})
                self.filed_change_dict["ms_run"].update({index: new_run})

            elif i.startswith("assay"):
                index = re.findall(r"[[](.*?)[]]", i)[0]
                assay_value_index = re.findall(r"[[](.*?)[]]", meta2[i])
                if i.endswith(variable_fields["assay"]):
                    for j in assay_value_index:
                        meta2[i] = meta2[i].replace(
                            j, str(int(j) + fields_count1["ms_run"])
                        )
                new_assay = str(int(index) + fields_count1["assay"])
                meta_assay.update({i.replace(index, new_assay): meta2[i]})
                self.filed_change_dict["assay"].update({index: new_assay})

            elif i.startswith("study_variable"):
                index = re.findall(r"[[](.*?)[]]", i)[0]
                study_variable_index = re.findall(r"[[](.*?)[]]", meta2[i])
                if i.endswith(variable_fields["study_variable"]):
                    for j in study_variable_index:
                        meta2[i] = meta2[i].replace(
                            j, str(int(j) + fields_count1["assay"])
                        )
                new_study_variable = str(int(index) + fields_count1["study_variable"])
                meta_study_variable.update(
                    {i.replace(index, new_study_variable): meta2[i]}
                )
                self.filed_change_dict["study_variable"].update(
                    {index: new_study_variable}
                )

        meta = meta_fix.copy()
        meta.update(meta_ms_run)
        meta.update(meta_assay)
        meta.update(meta_study_variable)

        return meta

    def protein_merge(self, pro1, pro2):
        variable_fields = [
            "protein_abundance_assay",
            "protein_abundance_study_variable",
            "protein_abundance_stdev_study_variable",
            "protein_abundance_std_error_study_variable",
            "search_engine_score[1]_ms_run",
        ]

        (fields_count1, variable_fields) = self.fields_count(pro1, variable_fields)
        update_cols = dict()

        for i in pro2.columns:
            for j in variable_fields:
                if i.startswith(j + "[") and i.endswith("]"):
                    index = re.findall(r"[[](.*?)[]]", i)[0]
                    update_cols.update(
                        {i: i.replace(index, str(int(index) + fields_count1[j]))}
                    )

        pro2 = pro2.rename(columns=update_cols)
        merge = pd.merge(pro1, pro2, how="outer", left_index=True, right_index=True)
        duplicate_cols = list()
        for i in merge.columns:
            if i.endswith("_x"):
                duplicate_cols.append([i, i.rstrip("_x") + "_y"])
        # recalculate
        merge = self.pro_recalculate(merge, fields=duplicate_cols)
        merge.sort_values(by="opt_global_result_type", inplace=True, ascending=False)
        merge.reset_index(drop=True, inplace=True)
        merge.fillna("null", inplace=True)

        return merge

    def peptide_merge(self, pep1, pep2):
        # reindex
        variable_fields = [
            "search_engine_score[1]_ms_run",
            "peptide_abundance_study_variable",
            "peptide_abundance_stdev_study_variable",
            "peptide_abundance_std_error_study_variable",
            "opt_global_mass_to_charge_study_variable",
            "opt_global_retention_time_study_variable",
        ]

        (fields_count1, variable_fields) = self.fields_count(pep1, variable_fields)
        pep2["spectra_ref"] = pep2["spectra_ref"].apply(
            lambda x: self.adjust_spectra_ref(x)
        )
        update_cols = dict()

        for i in pep2.columns:
            for j in variable_fields:
                if i.startswith(j + "[") and i.endswith("]"):
                    index = re.findall(r"[[](.*?)[]]", i)[-1]
                    new_col = str(int(index) + fields_count1[j]).join([j + "[", "]"])
                    update_cols.update({i: new_col})

        pep2 = pep2.rename(columns=update_cols)
        merge = pd.merge(
            pep1, pep2, on=["sequence", "modifications", "charge"], how="outer"
        )
        duplicate_cols = list()
        for i in merge.columns:
            if i.endswith("_x"):
                duplicate_cols.append([i, i.rstrip("_x") + "_y"])
        ## TODO: Need to recalculate mass-to-charge
        mz_cols = [
            col
            for col in merge.columns
            if col.startswith("opt_global_mass_to_charge_study_variable")
        ]
        merge = self.pep_recalculate(
            merge, fields=duplicate_cols, mz_study_variable=mz_cols
        )
        merge.sort_values(by="accession", inplace=True, ascending=False)
        new_cols = [col for col in merge.columns if not col.startswith("opt_")] + [
            col for col in merge.columns if col.startswith("opt_")
        ]
        merge = merge[new_cols]
        merge.reset_index(drop=True, inplace=True)
        merge.fillna("null", inplace=True)

        return merge

    def psm_merge(self, psm1, psm2):
        psm2["spectra_ref"] = psm2["spectra_ref"].apply(
            lambda x: self.adjust_spectra_ref(x)
        )
        psm1, psm2 = (
            psm1[[item for item in psm1.columns if item != ""]],
            psm2[[item for item in psm2.columns if item != ""]],
        )
        merge = pd.concat([psm1, psm2])
        merge.sort_values(by="accession", inplace=True, ascending=False)
        merge.reset_index(drop=True, inplace=True)

        merge.loc[:, "PSM_ID"] = merge.index
        merge.fillna("null", inplace=True)
        return merge


def cast_value(value):
    if value in ["", "null"]:
        value = None
    else:
        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass
    return value


def calculate_mz(seq, charge):
    """Remove unknown aminoacids and calculate mz
    :param seq: Sequences of peptides
    :type seq: str
    :param charge: charge of peptides
    :type seq: str
    :return: mz
    :rtype: float or NoneType
    """
    ref = "ARNDBCEQZGHILKMFPSTWYV"
    seq = "".join([i for i in seq if i in ref])
    if charge == "":
        return None
    else:
        return AASequence.fromString(seq).getMZ(int(charge))
