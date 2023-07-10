import pandas as pd
from mztabpy import MzTabPy, DiannConvert, MzTabMerge
import click

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@click.command("mztab_convert")
@click.option("--mztab_path", "-p")
@click.option("--directory", "-d", default="./")
@click.option("--type", "-t", default="all")
@click.option("--section", "-s", default="all")
@click.option("--removemeta", "-r", default=False)
@click.pass_context
def mztab_convert(ctx, mztab_path, directory, type, section, removemeta):
    """This script is used to separate mzTab into metadata, protein, peptide, and PSM. Metadata is orderedDict,
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
    """
    result_type, data_section, meta_flag = type, section, removemeta
    mztab = MzTabPy(
        mztab_path, single_cache_size=50 * 1024 * 1024, result_folder=directory
    )
    mztab.storage(type=result_type, section=data_section, removemeta=meta_flag)


cli.add_command(mztab_convert)


@click.command("hdf5_search")
@click.option("--hdf5_path", "-p")
@click.option("--section", "-s")
@click.option("--where", "-w", default=None)
@click.pass_context
def hdf5_search(ctx, hdf5_path, section, where):
    """This script is used to Load HDF5 into dataframe

    :param hdf5_path: Path to HDF5
    :type hdf5_path: str
    :param section: Indicates the data section of the mzTab that is required. "protein", "peptide" or "psm".
    :type section: str
    :param where: The filtering condition of the corresponding chunk is expressed as the key-value pair in one string,
        e.g. 'accession:P45464,sequence:TIQQGFEAAK', default None
    :type where: str
    """
    # pd.set_option("expand_frame_repr", False)
    condition = dict()
    if where:
        for i in where.replace(" ", "").split(","):
            condition.update({i.split(":")[0]: i.split(":")[1]})

    result = MzTabPy.loadHDF5(hdf5_path, section, condition)
    print(result)


cli.add_command(hdf5_search)

# TODO: diannconvert
@click.command("diannconvert")
@click.option("--directory", "-d")
@click.option("--diannparams", "-p")
@click.option("--charge", "-c", type=int)
@click.option("--missed_cleavages", "-m", type=int)
@click.option("--qvalue_threshold", "-q", type=float)
@click.option("--processors", "-pr", default=20)
@click.option("--threads_per_processor", "-t", default=8)
@click.option("--out", "-o", default="./")
@click.option("--block_size", "-b", default=500e6)
@click.pass_context
def diannconvert(
    ctx,
    directory,
    diannparams,
    charge,
    missed_cleavages,
    qvalue_threshold,
    processors,
    threads_per_processor,
    out,
    block_size,
):
    """This function is designed to convert the DIA-NN output into three standard formats: MSstats, Triqler and mzTab. These documents are
    used for quality control and downstream analysis.

    :param directory: DiannConvert specifies the folder where the required file resides. The folder contains
        the Dia-NN main report, experimental design file, mzml_info TSVs and protein sequence FASTA file
    :type directory: str
    :param diannparams: A list contains DIA parameters
    :type diannparams: list
    :param charge: The charge assigned by DIA-NN(max_precursor_charge)
    :type charge: int
    :param missed_cleavages: Allowed missed cleavages assigned by DIA-NN
    :type missed_cleavages: int
    :param qvalue_threshold: Threshold for filtering q value
    :type qvalue_threshold: float
    :param processors: Number of used processors, defaults to 20
    :type processors: int
    :param threads_per_processor: Number of threads used per processor, defaults to 8
    :type threads_per_processor: int
    :param out: Path to out directory, defaults to "./"
    :type out: str
    :param block_size: _description_, defaults to 500e6
    :type block_size: int
    """

    DiannConvert(
        directory,
        diannparams,
        charge,
        missed_cleavages,
        qvalue_threshold,
        processors,
        threads_per_processor,
        out,
        block_size,
    ).convert()


cli.add_command(diannconvert)


# TODO: mztabmerge
@click.command("mztabmerge")
@click.option("--mztab1", "-m1")
@click.option("--mztab2", "-m2")
@click.option("--single_cache_size", "-s", default=800 * 1024 * 1024)
@click.option("--out", "-o", default="./")
@click.pass_context
def mztabmerge(ctx, mztab1, mztab2, single_cache_size, out):
    """This script is used to merge two mzTabs into one

    :param mztab1: Path to the original mztab
    :type mztab1: str
    :param mztab2: Path to the mztab to be merged
    :type mztab2: str
    :param single_cache_size: Single cache size, default 500*1024*1024
    :type single_cache_size: int
    :param out: Folder to result files, default './'
    :type out: str
    """
    MzTabMerge(mztab1, mztab2, single_cache_size, out)


cli.add_command(mztabmerge)

if __name__ == "__main__":
    cli()
