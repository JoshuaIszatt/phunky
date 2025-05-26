import os
import json
import subprocess
import gzip
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import logging
import logging.config
import importlib.resources


# _____________________________________________________BASE


def configure_log(location=None, configuration=None):
    # Default logging settings if needed
    if configuration is None:
        configuration = importlib.resources.files("phunky") / "logging.json"
    # Read logging configuration
    with open(str(configuration), "r") as f:
        config = json.load(f)
    # Set the log file location
    logfile = 'phunky.log'
    if location is None:
        location = str(importlib.resources.files('phunky'))
    # Create logfile location and file if it does not exist
    os.makedirs(location, exist_ok=True)
    logfile = os.path.join(location, logfile)
    # Update the log file path in the logging configuration
    if 'handlers' in config and 'file' in config['handlers']:
        config['handlers']['file']['filename'] = logfile
    # Configure and set first message
    logging.config.dictConfig(config)
    logger = logging.getLogger(__name__)
    logger.info(f"Logging to {logfile}")
    logger.info(f"Log configuration: {str(config)}")
    return logger


# _____________________________________________________BIO


def gzip_file(file_in):
    """Compresses a file using gzip and returns the path of the compressed file.

    :param file_in: The path to the file to be compressed.
    :return: The path to the compressed file.
    :raises Exception: Raises if the file to be compressed is not found.
    """
    file_out = os.path.abspath(f'{file_in}.gz')
    try:
        with open(file_in, 'rb') as f_in:
            with gzip.open(file_out, 'wb') as f_out:
                content = f_in.read()
                f_out.write(content)
        return file_out
    except FileNotFoundError as e:
        raise Exception(e)


def convert_bam_to_fastq(bam_file, output_file):
    """Converts a BAM file into a FASTQ file using reformat.sh.

    :param bam_file: The path to the BAM file to be converted.
    :type bam_file: str
    :param output_file: The path to the output FASTQ file.
    :type output_file: str
    :raises Exception: Raises if the conversion from BAM to FASTQ process fails.
    """
    command = [
        'reformat.sh',
        f'in={bam_file}',
        f'out={output_file}'
    ]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        raise Exception(f"BAM to FQ conversion failed: {e}")


def porechop_abi(input_fq, output_fq):
    """Trims adaptors from the FASTQ file and saves trimmed reads.

    :param input_fq: The path to the FASTQ reads file to be trimmed.
    :type input_fq: str
    :param output_fq: The path to the output FASTQ reads file. 
    :type output_fq: str
    :raises Exception: Raises if the adaptor trimming process fails.
    """
    command = [
        'porechop_abi',
        '-abi',
        '-i', input_fq,
        '-o', output_fq
    ]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        raise Exception(f"porechop_abi failed: {e}")


def filtlong(reads_fastq_gz, output_fq, minlen=1000,
             target_bases=20000000, keep_percent=90):
    """Filters a FASTQ using Filtlong based on a minimum length, percentage and target bases and outputs filtered FASTQ.

    Parameters can be set using JSON config file, otherwise default parameters are used.

    :param reads_fastq_gz: The path to the FASTQ reads (.gz) file to be filtered.
    :type reads_fastq_gz: str
    :param output_fq: The path to the filtered FASTQ reads (.gz) file.
    :type output_fq: str
    :param minlen: Minimum read length, defaults to 1000.
    :type minlen: int, optional
    :param target_bases: Approximate total number of bases to retain, defaults to 20000000.
    :type target_bases: int, optional
    :param keep_percent: Percent of highest quality reads to retain, defaults to 90.
    :type keep_percent: int, optional
    :raises Exception: Raises if the Filtlong process fails.
    """
    command = [
        'filtlong', reads_fastq_gz,
        '--min_length', str(minlen),
        '--keep_percent', str(keep_percent),
        '--target_bases', str(target_bases),
        '--mean_q_weight', str(10)
    ]
    try:
        process = subprocess.run(command, check=True,
                                 capture_output=True)
        out_str = process.stdout.decode('utf-8')
        with open(output_fq, 'w') as f_out:
            f_out.write(out_str)
    except Exception as e:
        raise Exception(f"FiltLong failed: {e}")


def nanoplot(reads_fastq_gz, output_directory, barcode=None):
    """Generates read quality and length statistics using NanoPlot and returns the number of bases.

    :param reads_fastq_gz: Path to the input FASTQ (.gz) file.
    :type reads_fastq_gz: str
    :param output_directory: Directory where NanoPlot output files will be saved.
    :type output_directory: str
    :param barcode: Optional barcode prefix for plot filenames, defaults to None.
    :type barcode: str, optional
    :raises Exception: Raises if NanoPlot failes to execute or output files cannot be read.
    :return: Total number of bases reported in NanoStats.txt.
    :rtype: int
    """
    command = [
        'NanoPlot',
        '--fastq', reads_fastq_gz,
        '-o', output_directory,
        '--tsv_stats',
        '--format', 'png'
    ]
    if barcode is not None:
        command.extend(['-p', barcode])
    try:
        subprocess.run(command, check=True)
        stats = os.path.join(output_directory, 'NanoStats.txt')
        df = pd.read_csv(stats, sep='\t')
        number_of_bases = df.loc[df['Metrics'] == 'number_of_bases', 'dataset'].values[0]
        return number_of_bases
    except Exception as e:
        raise Exception(f"NanoPlot failed: {e}")


def flye_assembly(reads_fastq, output_directory, threads=8,
                  raise_on_fail=True):
    """Runs Flye assembly on the FASTQ reads and returns the path to the assembled contigs.
 
    :param reads_fastq: Path to the input FASTQ reads file.
    :type reads_fastq: str
    :param output_directory: Directory where Flye output files will be saved.
    :type output_directory: str
    :param threads: Number of CPU threads to be used for assembly, defaults to 8.
    :type threads: int, optional
    :param raise_on_fail: Whether to raise an execption if the assembler fails, defaults to True.
    :type raise_on_fail: bool, optional
    :raises Exception: Raises if Flye assembly fails and raise_on_fail is True.
    :return: Path to the assembled contigs FASTA file.
    :rtype: str
    """
    command = [
        'flye',
        '--nano-hq', reads_fastq,
        '-o', output_directory,
        '--threads', str(threads),
        '--iterations', str(5)
    ]
    try:
        subprocess.run(command, check=True)
        contigs = os.path.join(output_directory, 'assembly.fasta')
        if os.path.exists(contigs):
            print('Flye successful')
        return contigs
    except Exception as e:
        if raise_on_fail:
            raise Exception(f"Flye assembly failed: {e}")
        else:
            return None


def checkv(contigs, output_directory):
    """Runs CheckV on assembled contigs and outputs the results to the output directory.

    :param contigs: Path to the contigs FASTA file.
    :type contigs: str
    :param output_directory: Path to the output directory where CheckV results will be saved.
    :type output_directory: str
    :raises Exception: Raises if CheckV process fails.
    """
    command = [
        "checkv", "end_to_end",
        f"{contigs}",
        f"{output_directory}"
    ]
    try:
        subprocess.run(command, check=True)
        print("CheckV successful")
    except subprocess.CalledProcessError:
        raise Exception("CheckV failed")


def read_mapping(contigs_fasta, reads, output_directory, ram_mb=20000, mapped_sam=False):
    """Maps reads to contigs using BBMap and generates coverage statistics.

    :param contigs_fasta: Path to the contigs FASTA file.
    :type contigs_fasta: str
    :param reads: Path to the FASTQ reads file to be mapped.
    :type reads: str
    :param output_directory: Path to the output directory where results will be saved.
    :type output_directory: str
    :param ram_mb: Amount of RAM to allocate for BBMap (mb), defaults to 20000.
    :type ram_mb: int, optional
    :param mapped_sam: Option to output mapped reads in SAM format, defaults to False.
    :type mapped_sam: bool, optional
    :raises Exception: Raises if the read mapping process fails.
    :return: Paths to the base coverage, coverage statistics, scaffold statistics, and optionally the mapped SAM file.
    :rtype: str, str, str, bool
    """
    covstats = os.path.join(output_directory, "covstats.tsv")
    basecov = os.path.join(output_directory, "basecov.tsv")
    scafstats = os.path.join(output_directory, "scafstats.tsv")
    qhist = os.path.join(output_directory, "qhist.tsv")
    aqhist = os.path.join(output_directory, "aqhist.tsv")
    lhist = os.path.join(output_directory, "lhist.tsv")
    gchist = os.path.join(output_directory, "gchist.tsv")
    command = [
        "bbmap.sh",
        f"-Xmx{ram_mb}m",
        f"ref={contigs_fasta}",
        f"in={reads}",
        f"covstats={covstats}",
        f"basecov={basecov}",
        f"scafstats={scafstats}",
        f"qhist={qhist}",
        f"aqhist={aqhist}",
        f"lhist={lhist}",
        f"gchist={gchist}",
        "nodisk",
        'fastareadlen=600'
    ]
    if mapped_sam:
        mapped = os.path.join(output_directory, "mapped.sam")
        command.append(f"out={mapped}")
    else:
        mapped = False
    try:
        subprocess.run(command, check=True)
        return basecov, covstats, scafstats, mapped
    except Exception as e:
        raise Exception(f"Read mapping failed {e}")


def extract_contig(contigs_fasta, header, output_file, rename=None):
    """Extracts a specific contig from a multi-FASTA file based on the header.

    :param contigs_fasta: Path to the multi-FASTA file containing contigs.
    :type contigs_fasta: str
    :param header: Header (name) of the contig to be extracted.
    :type header: str
    :param output_file: Path to the output file where the extracted contig will be saved.
    :type output_file: str
    :param rename: Optional new name for the contig in the output file, defaults to None.
    :type rename: str, optional
    :raises Exception: Raises if the contig cannot be extracted or written to the output file.
    """
    with open(contigs_fasta, 'r') as handle:
        entries = SeqIO.parse(handle, 'fasta')
        with open(output_file, 'w') as textfile:
            for entry in entries:
                if header in entry.id:
                    if rename:
                        entry.id = rename
                    try:
                        SeqIO.write(entry, textfile, 'fasta')
                        break
                    except Exception as e:
                        raise Exception(f"Could not extract {output_file}: {e}")


def extract_contig_header(fasta_file):
    """This function extracts the header (name) and size of the largest contig from the multi-FASTA file.

    :param fasta_file: Path to the multi-FASTA file.
    :type fasta_file: str
    :return: A tuple with the header (name) and the size (length) of the largest contig.
    :rtype: str
    """
    seq_id = None
    size = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_length = len(record.seq)
        if seq_length > size:
            seq_id = record.id
            size = seq_length
    return seq_id, size


def generate_coverage_graph(header, basecov, output_directory):
    """Generates a per-base coverage graph for a specific contig based on the base coverage file.


    :param header: The header (name) of the contig for which the coverage graph is generated.
    :type header: str
    :param basecov: Path to the base coverage file (basecov.tsv).
    :type basecov: str
    :param output_directory: Directory where the coverage graph image will be saved.
    :type output_directory: str
    """
    headers = ["ID", "Pos", "Coverage"]
    df = pd.read_csv(basecov, sep='\t', comment='#', names=headers)
    coverage = df[df['ID'].str.contains(header)]
    mean_cov = df['Coverage'].mean()
    # Plot
    x_values = coverage['Pos']
    y_values = coverage['Coverage']
    plt.figure(figsize=(15, 8))
    plt.plot(x_values,
             y_values,
             marker=',',
             markersize=0.1,
             linestyle='-',
             color='b')
    plt.title(f"Per base coverage for {header} (Mean coverage: {mean_cov})")
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.grid(True)
    outfile = os.path.join(output_directory, f"{header}.png")
    plt.savefig(outfile, dpi=300)
