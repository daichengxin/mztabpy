import pandas as pd
from mztabpy import MzTabPy
import click

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@click.command("mztab_convert")
@click.option("--mztab_path", "-p")
@click.option("--directory", "-d", default = "./")
@click.option("--type", "-t", default = "all")
@click.option("--section", "-s", default = "all")
@click.option("--removemeta", "-r", default = False)
@click.pass_context
def mztab_convert(ctx, mztab_path, directory, type, section, removemeta):
    '''This script is used to separate mzTab into metadata, protein, peptide, and PSM. Metadata is orderedDict,
    and the rest is dataframe. It converts mzTab to a split table(.csv) or HDF5(.hdf5) depending on the options.
    
    :param mztab_path: The path to mzTab
    :type mztab_path: str
    :param directory: Folder to result files. Default "./"
    :type directory: bool
    :param type: Result type("tsv", "hdf5" or "all"). Default "all"
    :type type: str
    :param section: Indicates the data section of the mzTab that is required. "all", "protein", "peptide" or "psm".
    Default "all"
    :type section: str
    :param removemeta: Whether to remove metadata. Default False
    :type removemeta: bool
    '''
    result_type, data_section, meta_flag = type, section, removemeta
    mztab = MzTabPy(mztab_path, single_cache_size=50 * 1024 * 1024, result_folder=directory)
    mztab.storage(type=result_type, section=data_section, removemeta=meta_flag)

cli.add_command(mztab_convert)


@click.command("hdf5_search")
@click.option("--hdf5_path", "-p")
@click.option("--section", "-s")
@click.option("--where", "-w", default = None)
@click.pass_context
def hdf5_search(ctx, hdf5_path, section, where):
    '''This script is used to Load HDF5 into dataframe
    
    :param hdf5_path: Path to HDF5
    :type hdf5_path: str
    :param section: Indicates the data section of the mzTab that is required. "protein", "peptide" or "psm".
    :type section: str
    :param where: The filtering condition of the corresponding chunk is expressed as the key-value pair in one string,
        e.g. 'accession:P45464,sequence:TIQQGFEAAK', default None
    :type where: str
    '''
    # pd.set_option("expand_frame_repr", False)
    condition = dict()
    if where:
        for i in where.replace(" ", "").split(","):
            condition.update({i.split(":")[0] : i.split(":")[1]})

    result = MzTabPy.loadHDF5(hdf5_path, section, condition)
    print(result)

cli.add_command(hdf5_search)

if __name__ == "__main__":
    cli()
