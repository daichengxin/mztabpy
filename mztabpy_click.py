import os
import pandas as pd
from pyteomics import mztab
import click

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@click.command("mztab_parse")
@click.option("--mztab_path", "-p")
@click.option("--subtable", "-s")
@click.option("--hdf5", "-h5")
@click.option("--fillna", "-fa")
@click.option("--filter", "-fr")
@click.pass_context
def mztab_parse(ctx, mztab_path, subtable = False, hdf5 = False, fillna = False, filter = False):
    '''This script is used to separate mzTab into metadata, protein, peptide, and PSM. Metadata is orderedDict, and the rest is dataframe. 
        It converts mzTab to a split table(.csv) or HDF5(.hdf5) depending on the options.
    
    :param mztab_path: The path to mzTab
    :type mztab_path: str
    :param subtable: Whether to convert mzTab into subtables(.csv)
    :type subtable: bool
    :param hdf5: Whether to convert mzTab into HDF5
    :type hdf5: bool
    :param fillna: Whether to replace null values with "null"
    :type fillna: bool
    :param filter: Whether to remove empty columns
    :type filter: bool
    '''
    # Print file info
    folder_path = os.path.split(mztab_path)[0]
    file_name = os.path.split(mztab_path)[1]
    basename = os.path.splitext(file_name)[0]
    file_bytes = os.path.getsize(mztab_path)
    print("FILE INFO:"+"\n\t"+f"File Name: {file_name}"+"\n\t"+f"File Bytes: {file_bytes} B")

    mztab_data = mztab.MzTab(mztab_path)
    meta_dict = mztab_data.metadata
    protein_table =  mztab_data.protein_table
    peptide_table =  mztab_data.peptide_table
    psm_table = mztab_data.spectrum_match_table
    for i in meta_dict.keys():
        meta_dict.update({i: ';'.join(list(meta_dict[i])) if type(meta_dict[i]) == tuple else meta_dict[i]})
    meta_table = pd.DataFrame(meta_dict, index = [0])

    if fillna and not hdf5:
        protein_table.fillna("null", inplace = True)
        peptide_table.fillna("null", inplace = True)
        psm_table.fillna("null", inplace = True)
    
    if filter:
        protein_table = protein_table.dropna(axis=1, how="all")
        peptide_table = peptide_table.dropna(axis=1, how="all")
        psm_table = psm_table.dropna(axis=1, how="all")

    if subtable:
        meta_table.to_csv(f"./{basename}_meta.csv", index = False)
        protein_table.to_csv(f"./{basename}_protein.csv", index = False)
        peptide_table.to_csv(f"./{basename}_peptide.csv", index = False)
        if file_bytes/pow(1024, 3) >= 2:
            psm_table.to_csv(f"./{basename}_psm.csv.gz", index = False, compression = 'gzip')
        else:
            psm_table.to_csv(f"./{basename}_psm.csv", index = False)
        print("INFO: Subtables completed!")
    
    if hdf5:
        h5_name = f"./{basename}_mztab.hdf5"
        meta_table.to_hdf(h5_name, mode='a', key='meta', format = 't', complevel=9, complib='zlib')
        protein_table.to_hdf(h5_name, mode='a', key='protein', format = 't', complevel=9, complib='zlib')
        peptide_table.to_hdf(h5_name, mode='a', key='peptide', format = 't', complevel=9, complib='zlib')
        psm_table.to_hdf(h5_name, mode='a', key='psm', format = 't', complevel=9, complib='zlib')
        print("INFO: HDF5 completed!")

cli.add_command(mztab_parse)

if __name__ == "__main__":
    cli()
