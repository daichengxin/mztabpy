import os
import pandas as pd
from collections import OrderedDict
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

class MzTabPy:
    '''
    The MzTabPy class is used to separate mzTab into metadata, protein, peptide, and PSM. It can convert mzTab to 4 split
    tables(.tsv) or HDF5(.hdf5) depending on the options. The four subtables will be cached as dataframes if the given
    cache size allows.
            
    :param mztab_path: The path to mzTab
    :type mztab_path: str
    :param single_cache_size: Single cache size, default 500*1024*1024
    :type single_cache_size: int
    :param result_folder: Folder to result files, default './' 
    :type result_folder: bool
    '''
    def __init__(self, mztab_path, single_cache_size = 500*1024*1024, result_folder = "./"):
        self.mztab_path = mztab_path
        self.chunk_size = single_cache_size
        self.folder_path = result_folder
        file_name = os.path.split(mztab_path)[1]
        self.basename = os.path.splitext(file_name)[0]
        self.file_bytes = os.path.getsize(mztab_path)
        self.chunks_num = int(self.file_bytes / self.chunk_size) + 1
        self.psm_compression = None if self.file_bytes < 1 * pow(1024, 3) else "gzip"
        if self.file_bytes < self.chunk_size:
            logging.info("Start loading...")
            self.loadmzTab()
            logging.info("Four subtables have been cached!")
        else:
            logging.warning("No subtable is cached because the single cache is smaller than file bytes!")


    def loadmzTab(self):
        '''If the computer has enough cache, read mzTab as four global dataframes('metadata', 'protein', 'peptide' and 'psm').
        '''
        self.meta_dict = OrderedDict()
        with open(self.mztab_path, 'r') as f:
            chunk = f.read(self.file_bytes)
            (self.protein, self.peptide, self.psm) = self.loadChunk(chunk)
            self.metadata = pd.DataFrame(self.meta_dict, index = [0])


    def loadChunk(self, chunk):
        '''This function processes data blocks from streams, classifies them into metadata, protein, peptide and 
            PSM, and caches them as dataframes.
        
        :param chunk: A block of data
        :type chunk: str
        '''
        data = chunk.split('\n')
        meta_dict = dict()
        pro_data, pep_data, psm_data = [], [], []
        for i in data:
            i = i.strip("\n")
            row_list = i.split('\t')
            if row_list[0] == "MTD":
                meta_dict[row_list[1]] = row_list[2]
            elif row_list[0] == "PRH":
                logging.info("Metadata processing complete. Start processing the protein subtable...")
                self.pro_cols = row_list
            elif row_list[0] == "PRT":
                pro_data.append(row_list)
            elif row_list[0] == "PEH":
                logging.info("Protein subtable processing complete. Start processing the peptide subtable...")
                self.pep_cols = row_list
            elif row_list[0] == "PEP":
                pep_data.append(row_list)
            elif row_list[0] == "PSH":
                logging.info("Peptide subtable processing complete. Start processing the psm subtable...")
                self.psm_cols = row_list
            elif row_list[0] == "PSM":
                psm_data.append(row_list)
            else:
                continue

        for i in meta_dict.keys():
            self.meta_dict.update({i: ';'.join(list(meta_dict[i])) if type(meta_dict[i]) == tuple else meta_dict[i]})

        if len(pro_data) > 0:
            protein_table = pd.DataFrame(pro_data, columns=self.pro_cols, dtype='str')
            protein_table.drop(columns='PRH', inplace=True)
            protein_table.set_index("accession", inplace=True, drop=False)
        else: 
            protein_table = pd.DataFrame()

        if len(pep_data) > 0:
            peptide_table = pd.DataFrame(pep_data, columns=self.pep_cols, dtype='str')
            peptide_table.drop(columns='PEH', inplace=True)
            peptide_table.set_index("accession", inplace=True, drop=False)
        else: 
            peptide_table = pd.DataFrame()

        if len(psm_data) > 0:
            psm_table = pd.DataFrame(psm_data, columns=self.psm_cols, dtype='str')
            psm_table.drop(columns='PSH', inplace=True)
            psm_table.set_index("accession", inplace=True, drop=False)
            if "opt_global_cv_MS:1002217_decoy_peptide" not in psm_table.columns.values:
                psm_table['opt_global_cv_MS:1002217_decoy_peptide'] = psm_table.apply(
                    lambda x: 1 if any(i in x['accession'] for i in ["DECOY_", "decoy"]) else 0, axis=1) 
        else: 
            psm_table = pd.DataFrame()
        
        protein_table.fillna("null", inplace = True)
        peptide_table.fillna("null", inplace = True)
        psm_table.fillna("null", inplace = True)

        return protein_table, peptide_table, psm_table
    

    def storage(self, type = "all", section = "all", removemeta = False):
        '''Store mzTab as tsv subtables or HDF5.

        :param type: Result type("tsv", "hdf5" or "all"). Default "all"
        :type type: str
        :param section: Indicates the data section of the mzTab that is required. "all", "protein", "peptide" or "psm".
        Default "all"
        :type section: str
        :param removemeta: Whether to remove metadata. Default False
        :type removemeta: bool
        '''
        if type == "tsv" or type == "all":
            self.to_tsv = True 
            self.meta_path = self.folder_path + self.basename + "_metadata_openms.tsv"
            self.pro_path = self.folder_path + self.basename + "_protein_openms.tsv"
            self.pep_path = self.folder_path + self.basename + "_peptide_openms.tsv"
            self.psm_path = self.folder_path + self.basename + \
                ("_psm_openms.tsv" if self.file_bytes < 1 * pow(1024, 3) else "_openms.psms.mzTab.gz")
        else:
            self.to_tsv = False

        if type == "hdf5" or type == "all":
            self.to_hdf5 = True
            self.h5_path = self.folder_path + self.basename + ".hdf5"
        else: 
            self.to_hdf5 = False

        logging.info("Start converting...")
        self.stream_storage(section, removemeta)


    def stream_storage(self, section, removemeta):
        '''This function streams the mzTab and stores it in four tsv subtables (Metadata, protein, peptide, PSM) 
            or one HDF5 containing these four parts.
        
        :param section: Indicates the data section of the mzTab that is required. "all", "protein", "peptide" or "psm"
        :type section: str
        :param removemeta: Whether to remove metadata
        :type removemeta: bool
        '''
        self.meta_dict = OrderedDict()
        chunk_index = 0
        pro_chunk, pep_chunk, psm_chunk = [], [], []
        chunk_dict = dict.fromkeys(['protein', 'peptide', 'psm'], '')
        chunk_dict.update({'mzTab_file_path': self.mztab_path, 'mztab_file_bytes': self.file_bytes})

        with open(self.mztab_path, 'r') as f:
            while True:
                chunk = f.read(self.chunk_size)
                if not chunk:
                    break

                (pro_df, pep_df, psm_df) = self.loadChunk(chunk)
                if self.to_tsv:
                    if section == "all" or section == "protein":
                        if os.path.exists(self.pro_path):
                            pro_df.to_csv(self.pro_path, mode="a", sep = '\t', index = False, header = False)
                        else:
                            pro_df.to_csv(self.pro_path, mode="a", sep = '\t', index = False, header = True)
                    
                    if section == "all" or section == "peptide":
                        if os.path.exists(self.pep_path):
                            pep_df.to_csv(self.pep_path, mode="a", sep = '\t', index = False, header = False)
                        else:
                            pep_df.to_csv(self.pep_path, mode="a", sep = '\t', index = False, header = True)

                    if section == "all" or section == "psm":
                        if os.path.exists(self.psm_path):
                            psm_df.to_csv(self.psm_path, mode="a", sep = '\t', index = False, header = False, compression=self.psm_compression)
                        else:
                            psm_df.to_csv(self.psm_path, mode="a", sep = '\t', index = False, header = True, compression=self.psm_compression)
                            
                if self.to_hdf5:
                    if len(pro_df) != 0 and (section == "all" or section == "protein"):
                        pro_df.to_hdf(self.h5_path, mode='a', key=f'Chunk{chunk_index}_protein', format='t', complevel=9, complib='zlib')
                        pro_chunk.append(f'Chunk{chunk_index}_protein')
                        if 'protein_cols' not in chunk_dict:
                            chunk_dict.update({'protein_cols': ';'.join(pro_df.columns.values)})

                    if len(pep_df) != 0 and (section == "all" or section == "peptide"):
                        pep_df.to_hdf(self.h5_path, mode='a', key=f'Chunk{chunk_index}_peptide', format='t', complevel=9, complib='zlib')
                        pep_chunk.append(f'Chunk{chunk_index}_peptide')
                        if 'peptide_cols' not in chunk_dict:
                            chunk_dict.update({'peptide_cols': ';'.join(pep_df.columns.values)})

                    if len(psm_df) != 0 and (section == "all" or section == "psm"): 
                        psm_df.to_hdf(self.h5_path, mode='a', key=f'Chunk{chunk_index}_psm', format='t', complevel=9, complib='zlib')
                        psm_chunk.append(f'Chunk{chunk_index}_psm')
                        if 'psm_cols' not in chunk_dict:
                            chunk_dict.update({'psm_cols': ';'.join(psm_df.columns.values)})
                chunk_index += 1

        meta_df = pd.DataFrame(self.meta_dict, index = [0])
        chunk_dict.update({'protein': ';'.join(pro_chunk),
                           'peptide': ';'.join(pep_chunk),
                           'psm': ';'.join(psm_chunk)})

        chunks_info = pd.DataFrame(chunk_dict, index = [0])

        if self.to_tsv:
            if not removemeta:
                meta_df.to_csv(self.meta_path, mode="a", sep = '\t', index = False, header = True)
            logging.info("mzTab separation and storage complete!")
        if self.to_hdf5:
            chunks_info.to_hdf(self.h5_path, mode='a', key='CHUNKS_INFO', format='t', complevel=9, complib='zlib')
            if not removemeta:
                meta_df.to_hdf(self.h5_path, mode='a', key='meta', format='t', complevel=9, complib='zlib')
            logging.info("mzTab binary storage complete!")


    def loadHDF5(h5, subtable, where=None):
        ''' Load HDF5 into dataframe

        :param h5: Path to HDF5
        :type h5: str
        :param subtable: Name of the subtable, which contains 'metadata', 'protein', 'peptide' and 'psm'
        :type subtable: str
        :param where: The filtering condition of the corresponding chunk is expressed as the key-value pair in the dictionary,
            e.g. {'accession': 'P45464', 'sequence': 'TIQQGFEAAK'}, default None
        :type where: dict
        :return result: Results filtered by criteria in the target subtable
        :rtype result: dataframe
        '''  
        chunks_info = pd.read_hdf(h5, key='CHUNKS_INFO')
        chunks = chunks_info[subtable].values[0].split(';') if subtable in chunks_info else []

        result = pd.DataFrame(columns=chunks_info[subtable + '_cols'].values[0].split(';'))
        if where:
            for i in chunks:
                df = pd.read_hdf(h5, key=i)
                if 'accession' in where.keys() and where['accession'] not in df.index:
                    continue
                else:
                    for j in where.keys():
                        df = df[df[j] == where[j]]
                    target = df
                result = pd.concat([result, target])
            return result

        else:
            if os.path.getsize(h5) < 100*pow(1024, 2):
                for i in chunks:
                    df = pd.read_hdf(h5, key=i)
                    result = pd.concat([result, df])
                return result
            else:
                logging.warning("No result has been cached because the single cache is too small!")

        


