from collections import OrderedDict
from copy import copy
from mztabpy import MzTabPy
import pandas as pd
import numpy as np
import re
import os
import logging

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
    def __init__(self, mztab1, mztab2, single_cache_size = 800*1024*1024, result_folder = "./"):
        if os.path.getsize(mztab2) > single_cache_size:
            logging.error("The size of mztab2 must be able to be cached, please increase single cache size!")
            return
        self.mztab_path1, self.mztab_path2 = mztab1, mztab2
        logging.info(f"Start loading {self.mztab_path2}")
        (meta_dict2, pro2, pep2, psm2) = MzTabPy(mztab_path = self.mztab_path2, section_tag = True, df_index = False).loadmzTab()
        self.merge_dict, self.meta_dict = {"meta": meta_dict2, "pro": pro2, "pep": pep2, "psm": psm2}, dict()
        self.psmid = 0
        self.chunk_size = single_cache_size
        self.file_bytes = os.path.getsize(mztab1)
        self.file_path = result_folder + "merge.mztab"
        if self.file_bytes < self.chunk_size:
            logging.info(f"Start loading {self.mztab_path1}")
            (meta_dict1, pro1, pep1, psm1) = MzTabPy(mztab_path = self.mztab_path1, section_tag = True, df_index = False).loadmzTab()
            merge_metadict = self.meta_merge(meta_dict1, meta_dict2)
            merge_meta = pd.DataFrame(merge_metadict, index = [0])
            meta = pd.DataFrame(merge_meta.values.T, index = list(merge_meta.columns), columns = merge_meta.index)
            meta.insert(0, "key", meta.index)
            meta.insert(0, "tag", "MTD")
            merge_pro = self.protein_merge(pro1, pro2)
            merge_pep = self.peptide_merge(pep1, pep2)
            map_list = psm1.loc[:, "opt_global_map_index"].to_list()
            max_map_index = max([int(i) for i in map_list if not "null" in map_list])

            def reindex(x):
                return int(x) + max_map_index

            psm2.loc[:, "opt_global_map_index"] = psm2.loc[:, "opt_global_map_index"].apply(reindex)
            merge_psm = self.psm_merge(psm1, psm2)

            meta.loc["", :] = ""
            merge_pro.loc[len(merge_pro) + 1, :] = ""
            merge_pep.loc[len(merge_pep) + 1, :] = ""
            with open(self.file_path, "a", newline = "") as f:
                meta.to_csv(f, mode="a", sep = '\t', index = False, header = False)
                merge_pro.to_csv(f, mode="a", sep = '\t', index = False, header = True)
                merge_pep.to_csv(f, mode="a", sep = '\t', index = False, header = True)
                merge_psm.to_csv(f, mode="a", sep = '\t', index = False, header = True)

        else:
            logging.warning("No subtable is cached because the single cache is smaller than file bytes!")
            self.pro_cache, self.pro_merge_status, self.pep_cache, self.pep_merge_status = False, False, False, False
            self.chunkDict = {"meta": [], "pro": [], "pep": [], "psm": []}
            self.chunkID = 0
            self.row_remains = ''
            logging.info(f"Start loading {self.mztab_path1}")
            with open(self.mztab_path1, 'r') as f:
                while True:
                    chunk = self.row_remains + f.read(self.chunk_size)
                    if not chunk:
                        break
                    self.loadMergeChunk(chunk)
            f.close()


    def select_non_missing(self, values):
        if values[0] == np.nan or values[0] == "null":
            values[0] = values[1]
        else:pass

        return values


    def max_or_null(self, x, y):
        data_list = [i for i in [x, y] if i not in [np.nan, None, '', 'null']]
        if len(data_list) == 0: return None
        elif len((data_list)) == 1: return data_list[0]
        else:
            return max(float(x), float(y))
    

    def min_or_null(self, x, y):
        data_list = [i for i in [x, y] if i not in [np.nan, None, '', 'null']]
        if len(data_list) == 0: return None
        elif len((data_list)) == 1: return data_list[0]
        else:
            return min(float(x), float(y))


    def paste_or_null(self, x, y):
        paste_list = [i for i in [x, y] if i not in [np.nan, None, '', 'null']]
        if len(paste_list) == 0: return None
        elif len((paste_list)) == 1: return paste_list[0]
        else:
            return "|".join([str(x), str(y)])
            

    def parse_extral_peps(self, x, y, field):
        if x == "PEP":
            if y != "PEP":
                cols = np.array(field)[1, :]
            else:pass
        else:
            if y == "PEP":
                cols = np.array(field)[:, 1]
            else:
                logging.warning("Some errors occur in peptide merging!")
        return cols


    def pro_recalculate(self, merged, fields):
        # if not set(sum(fields, [])) <= set(merged.columns): 
        #     return(merged)
        # recalculate
        x_cols = [col for col in merged.columns if col.endswith("_x")]
        y_cols = [col for col in merged.columns if col.endswith("_y")]
        for i in fields:
            merged[i] = merged.apply(lambda x: self.select_non_missing(x[i]), axis = 1)
        if "best_search_engine_score[1]_x" in merged.columns:
            merged.loc[:, "best_search_engine_score[1]_x"] = merged.apply(lambda x: 
                self.min_or_null(x["best_search_engine_score[1]_x"], x["best_search_engine_score[1]_y"]), axis=1)

        if "protein_coverage_x" in merged.columns:
            merged.loc[:, "protein_coverage_x"] = merged.apply(lambda x: 
                self.max_or_null(x["protein_coverage_x"], x["protein_coverage_y"]), axis=1)

        if "opt_global_nr_found_peptides_x" in merged.columns:
            merged.loc[:, "opt_global_nr_found_peptides_x"] = merged.apply(lambda x: 
                self.max_or_null(x["opt_global_nr_found_peptides_x"], x["opt_global_nr_found_peptides_y"]), axis=1)

        if "opt_global_Posterior_Probability_score_x" in merged.columns:
            merged.loc[:, "opt_global_Posterior_Probability_score_x"] = merged.apply(lambda x: 
                self.min_or_null(x["opt_global_Posterior_Probability_score_x"], x["opt_global_Posterior_Probability_score_y"]), axis=1)

        # strip columns from pro2
        pro2_variable_cols = np.array(fields)[:, 1]
        combine_columns = [item for item in merged.columns if item not in pro2_variable_cols]
        combine_columns = [x.strip() for x in combine_columns if x.strip()!='']
        merged = merged[combine_columns]

        merged_cols = dict()
        for i in merged.columns:
            if i.endswith("_x"):
                merged_cols.update({i: i.rstrip("_x")})
        
        merged = merged.rename(columns = merged_cols)
        return merged


    def pep_recalculate(self, merged, fields, mz_study_variable):
        # recalculate
        for i in fields:
            merged[i] = merged.apply(lambda x: self.select_non_missing(x[i]), axis = 1)

        if "best_search_engine_score[1]_x" in merged.columns:
            merged.loc[:, "best_search_engine_score[1]_x"] = merged.apply(lambda x: 
                self.min_or_null(x["best_search_engine_score[1]_x"], x["best_search_engine_score[1]_y"]), axis=1)

        if "opt_global_Posterior_Error_Probability_score_x" in merged.columns:
            merged.loc[:, "opt_global_Posterior_Error_Probability_score_x"] = merged.apply(lambda x: 
                self.min_or_null(x["opt_global_Posterior_Error_Probability_score_x"], x["opt_global_Posterior_Error_Probability_score_y"]), axis=1)

        if "opt_global_q-value_x" in merged.columns:
            merged.loc[:, "opt_global_q-value_x"] = merged.apply(lambda x: 
                self.min_or_null(x["opt_global_q-value_x"], x["opt_global_q-value_y"]), axis=1)

        if "spectra_ref_x" in merged.columns:
            merged["spectra_ref_x"] = merged.apply(lambda x: 
                self.paste_or_null(x["spectra_ref_x"], x["spectra_ref_y"]), axis=1)
        
        if "mass_to_charge_x" in merged.columns:
            if len(mz_study_variable) > 0:
                for i in mz_study_variable:
                    merged.loc[:, i] = merged.apply(lambda x: float(x[i]) if not x[i] in ["null", None] else None, axis = 1)
                merged.loc[:, "mass_to_charge_x"] = merged[mz_study_variable].mean()

        ## TODO apply the extral peptides (less or more)
        # Ignore the missing '_y' section because the '_x' section is always truncated as the final result
        unconsensus_peps = merged[merged["PEH_x"] != "PEP"]
        merged = merged[merged["PEH_x"] == "PEP"]
        unconsensus_cols = [col for col in unconsensus_peps.columns if not col.endswith("_x")]
        unconsensus = unconsensus_peps[unconsensus_cols]
        unconsensus.columns = [i.rstrip("_y") for i in unconsensus.columns]

        pep2_variable_cols = np.array(fields)[:, 1]
        combine_columns = [item for item in merged.columns if item not in pep2_variable_cols]
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
                new_spectra_ref = value.replace(ms_run.join(["[", "]"]), self.filed_change_dict["ms_run"][ms_run].join(["[", "]"]))
                return new_spectra_ref
        else:pass


    def loadMergeChunk(self, chunk):
        '''This function processes data blocks from streams, classifies them into metadata, protein, peptide and 
            PSM, and caches them as dataframes.
        
        :param chunk: A block of data
        :type chunk: str
        '''
        data = chunk.split('\n')
        if not chunk.endswith('\n'):
            self.row_remains = data[-1]
            data = data[0:-1]
        else:
            self.row_remains = ""
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
                self.meta_cache = True
            elif row_list[0] == "PRT":
                pro_data.append(row_list)
            elif row_list[0] == "PEH":
                logging.info("Protein subtable processing complete. Start processing the peptide subtable...")
                self.pep_cols = row_list
                self.pro_cache = True
            elif row_list[0] == "PEP":
                pep_data.append(row_list)
            elif row_list[0] == "PSH":
                logging.info("Peptide subtable processing complete. Start processing the psm subtable...")
                self.psm_cols = row_list
                self.pep_cache = True
            elif row_list[0] == "PSM":
                psm_data.append(row_list)
            else:
                continue
        
        for i in meta_dict.keys():
            meta_dict.update({i: ';'.join(list(meta_dict[i])) if type(meta_dict[i]) == tuple else meta_dict[i]})
        self.meta_dict.update(self.meta_merge(meta_dict, self.merge_dict["meta"]))

        if self.meta_cache and self.chunkID == 0:
            merge_meta = pd.DataFrame(self.meta_dict, index = [0])
            meta = pd.DataFrame(merge_meta.values.T, index = list(merge_meta.columns), columns = merge_meta.index)
            meta.insert(0, "key", meta.index)
            meta.insert(0, "tag", "MTD")
            meta.loc["", :] = ""
            with open(self.file_path, "a", newline = "") as f:
                meta.to_csv(f, mode="a", sep = '\t', index = False, header = False)

        if len(pro_data) > 0:
            self.chunkDict["pro"].append(self.chunkID)
            protein_table = pd.DataFrame(pro_data, columns=self.pro_cols, dtype='str')
            if not self.pro_cache:
                self.pro_differ_chunks = True
            sub_pro = self.merge_dict["pro"][self.merge_dict["pro"]["accession"].isin(protein_table.columns)]
            if len(sub_pro) == 0:
                sub_pro = self.merge_dict["pro"]
                self.pro_merge_status = True
            else:
                self.merge_dict["pro"] = self.merge_dict["pro"][~self.merge_dict["pro"]["accession"].isin(protein_table.columns)]
            merge_pro = self.protein_merge(protein_table, sub_pro)

            if self.pro_merge_status:
                merge_pro.loc[len(merge_pro) + 1, :] = ""
                
            merge_pro.fillna("null", inplace = True)
            if len(self.chunkDict["pro"]) == 1: 
                pro_header = True
            else:
                pro_header = False
            with open(self.file_path, "a", newline = "") as f:
                merge_pro.to_csv(f, mode="a", sep = '\t', index = False, header = pro_header)

            del protein_table, sub_pro, merge_pro

        if len(pep_data) > 0:
            self.chunkDict["pep"].append(self.chunkID)
            peptide_table = pd.DataFrame(pep_data, columns=self.pep_cols, dtype='str')
            if not self.pep_cache:
                self.pep_differ_chunks = True
            sub_pep = self.merge_dict["pep"][self.merge_dict["pep"]["accession"].isin(peptide_table.columns)]
            if len(sub_pep) == 0:
                sub_pep = self.merge_dict["pep"]
                self.pep_merge_status = True
            else:
                self.merge_dict["pep"] = self.merge_dict["pep"][~self.merge_dict["pep"]["accession"].isin(peptide_table.columns)]
            merge_pep = self.peptide_merge(peptide_table, sub_pep)

            if self.pep_merge_status:
                merge_pep.loc[len(merge_pep) + 1, :] = ""

            merge_pep.fillna("null", inplace = True)
            if len(self.chunkDict["pep"]) == 1: 
                pep_header = True
            else:
                pep_header = False
            with open(self.file_path, "a", newline = "") as f:
                merge_pep.to_csv(f, mode="a", sep = '\t', index = False, header = pep_header)

            del peptide_table, sub_pep, merge_pep

        if len(psm_data) > 0:
            self.chunkDict["psm"].append(self.chunkID)
            psm_table = pd.DataFrame(psm_data, columns=self.psm_cols, dtype='str')
            if "opt_global_cv_MS:1002217_decoy_peptide" not in psm_table.columns.values:
                psm_table['opt_global_cv_MS:1002217_decoy_peptide'] = psm_table.apply(
                    lambda x: 1 if any(i in x['accession'] for i in ["DECOY_", "decoy"]) else 0, axis=1) 

            sub_psm = self.merge_dict["psm"][self.merge_dict["psm"]["accession"].isin(psm_table.columns)]
            if len(sub_psm) == 0:
                sub_psm = self.merge_dict["psm"]
            else:
                self.merge_dict["psm"] = self.merge_dict["psm"][~self.merge_dict["psm"]["accession"].isin(psm_table.columns)]

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
            merge_psm.fillna("null", inplace = True)
            with open(self.file_path, "a", newline = "") as f:
                merge_psm.to_csv(f, mode="a", sep = '\t', index = False, header = psm_header)

            del psm_table, sub_psm, merge_psm
     
        self.chunkID += 1


    def meta_merge(self, meta1, meta2):
        # Warning values clash and missing keys
        prefix_list = ["ms_run", "assay", "study_variable"]
        meta_fix, meta_ms_run, meta_assay, meta_study_variable = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()
        self.filed_change_dict = dict.fromkeys(prefix_list, dict())
        miss = [item for item in meta1.keys() if item not in meta2.keys()]
        for i in miss:
            logging.warning("Metadata entry for field '" + i + "' is missing")

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
            logging.warning("Metadata entries for field '" + i + "' differ")

        for i in fix_columns:
            meta_fix.update({i: meta1[i]})
        
        variable_fields={"ms_run": "location", "assay": "ms_run_ref", "study_variable": "assay_refs"}

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
                        meta2[i] = meta2[i].replace(j, str(int(j) + fields_count1["ms_run"]))
                new_assay = str(int(index) + fields_count1["assay"])
                meta_assay.update({i.replace(index, new_assay): meta2[i]})
                self.filed_change_dict["assay"].update({index: new_assay})

            elif i.startswith("study_variable"):
                index = re.findall(r"[[](.*?)[]]", i)[0]
                study_variable_index = re.findall(r"[[](.*?)[]]", meta2[i])
                if i.endswith(variable_fields["study_variable"]):
                    for j in study_variable_index:
                        meta2[i] = meta2[i].replace(j, str(int(j) + fields_count1["assay"]))
                new_study_variable = str(int(index) + fields_count1["study_variable"])
                meta_study_variable.update({i.replace(index, new_study_variable): meta2[i]})
                self.filed_change_dict["study_variable"].update({index: new_study_variable})

        meta = meta_fix.copy()
        meta.update(meta_ms_run)
        meta.update(meta_assay)
        meta.update(meta_study_variable)
        
        return meta


    def protein_merge(self, pro1, pro2):
        variable_fields = ["protein_abundance_assay", "protein_abundance_study_variable", "protein_abundance_stdev_study_variable",
                        "protein_abundance_std_error_study_variable", "search_engine_score[1]_ms_run"]

        (fields_count1, variable_fields) = self.fields_count(pro1, variable_fields)
        update_cols = dict()

        for i in pro2.columns:
            for j in variable_fields:
                if i.startswith(j+"[") and i.endswith("]"):
                    index = re.findall(r"[[](.*?)[]]", i)[0]
                    update_cols.update({i: i.replace(index, str(int(index) + fields_count1[j]))})

        pro2 = pro2.rename(columns = update_cols)
        merge = pd.merge(pro1, pro2, how="outer", left_index=True, right_index=True)
        duplicate_cols = list()
        for i in merge.columns:
            if i.endswith("_x"):
                duplicate_cols.append([i, i.rstrip("_x") + "_y"])    
        # recalculate
        merge = self.pro_recalculate(merge, fields=duplicate_cols)
        merge.sort_values(by="opt_global_result_type" , inplace=True, ascending=False) 
        merge.reset_index(drop=True, inplace=True)
        merge.fillna("null", inplace = True)

        return merge
        

    def peptide_merge(self, pep1, pep2):
        # reindex
        variable_fields = ["search_engine_score[1]_ms_run", "peptide_abundance_study_variable", "peptide_abundance_stdev_study_variable",
                        "peptide_abundance_std_error_study_variable", "opt_global_mass_to_charge_study_variable", "opt_global_retention_time_study_variable"]

        (fields_count1, variable_fields) = self.fields_count(pep1, variable_fields)
        pep2["spectra_ref"] = pep2["spectra_ref"].apply(lambda x: self.adjust_spectra_ref(x))
        update_cols = dict()

        for i in pep2.columns:
            for j in variable_fields:
                if i.startswith(j+"[") and i.endswith("]"):
                    index = re.findall(r"[[](.*?)[]]", i)[-1]
                    new_col = str(int(index) + fields_count1[j]).join([j+"[", "]"])
                    update_cols.update({i: new_col})

        pep2 = pep2.rename(columns = update_cols)
        merge = pd.merge(pep1, pep2, on=["sequence", "modifications", "charge"], how="outer")
        duplicate_cols = list()
        for i in merge.columns:
            if i.endswith("_x"):
                duplicate_cols.append([i, i.rstrip("_x") + "_y"]) 
        ## TODO: Need to recalculate mass-to-charge
        mz_cols = [col for col in merge.columns if col.startswith("opt_global_mass_to_charge_study_variable")]
        merge = self.pep_recalculate(merge, fields=duplicate_cols, mz_study_variable=mz_cols)
        merge.sort_values(by="accession" , inplace=True, ascending=False) 
        new_cols = [col for col in merge.columns if not col.startswith("opt_")] + [col for col in merge.columns if col.startswith("opt_")]
        merge = merge[new_cols]
        merge.reset_index(drop=True, inplace=True)
        merge.fillna("null", inplace = True)
        
        return merge
        

    def psm_merge(self, psm1, psm2):
        psm2["spectra_ref"] = psm2["spectra_ref"].apply(lambda x: self.adjust_spectra_ref(x))
        psm1, psm2 = psm1[[item for item in psm1.columns if item != ""]], psm2[[item for item in psm2.columns if item != ""]]
        merge = pd.concat([psm1, psm2])
        merge.sort_values(by="accession" , inplace=True, ascending=False) 
        merge.reset_index(drop=True, inplace=True)

        merge.loc[:, "PSM_ID"] = merge.index
        merge.fillna("null", inplace = True)
        return merge


