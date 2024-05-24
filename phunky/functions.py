import os
import subprocess
import gzip
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd


###_____________________________________________________FUNCTIONS


def gzip_file(file_in):
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
    except Exception as e:
        raise Exception(f"NanoPlot failed: {e}")


def flye_assembly(reads_fastq, output_directory, threads=8):
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
        raise Exception(f"Flye assembly failed: {e}")


def checkv(contigs, output_directory):
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


def extract_contig_header(fasta_file):
    """
    This function takes a path to a multi-FASTA file and returns
    the header (name) and size of the largest contig in the file.

    :param fasta_file: Path to the multi-FASTA file.
    :return: A tuple with the header (name) and the size (length) of the largest contig.
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
    headers = ["ID", "Pos", "Coverage"]
    df = pd.read_csv(basecov, sep='\t', comment='#', names=headers)
    coverage = df[df['ID'].str.contains(header)]
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
    plt.title(f"Per base coverage for {header}")
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.grid(True)
    outfile = os.path.join(output_directory, f"{header}.png")
    plt.savefig(outfile, dpi=300)
