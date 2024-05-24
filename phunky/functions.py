import os
import subprocess
import gzip


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
        raise e

def convert_bam_to_fastq(bam_file, output_file):
    command = [
        'reformat.sh',
        f'in={bam_file}', 
        f'out={output_file}'
    ]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"BAM to FQ conversion failed: {e}")
        raise


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
        print(f"porechop_abi failed: {e}")
        raise


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
        print(f"FiltLong failed: {e}")
        raise


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
        print(f"NanoPlot failed: {e}")
        raise


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
        print(f"Flye assembly failed: {e}")
        raise


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
        print("CheckV failed")
