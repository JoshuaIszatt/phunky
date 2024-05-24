import os
from phunky.functions import (
    convert_bam_to_fastq,
    porechop_abi,
    gzip_file,
    filtlong,
    nanoplot,
    flye_assembly,
    checkv
)


###_____________________________________________________PIPELINES


def phage_assembly_pipeline(input_file, output_dir):
    # Create output location
    name = os.path.basename(str(input_file)[:-4])
    out = os.path.join(output_dir, name)
    os.makedirs(out, exist_ok=False)

    # Convert
    fq_raw = os.path.join(output_dir, f'{name}_raw.fastq')
    convert_bam_to_fastq(input_file, fq_raw)

    # Remove adapters
    fq_trim = os.path.join(output_dir, f'{name}_trimmed.fastq')
    porechop_abi(fq_raw, fq_trim)

    # gzip file
    fq_trim_gz = gzip_file(fq_trim)

    # Filter
    fq_filt = os.path.join(output_dir, f'{name}_filtered.fastq')
    filtlong(fq_trim_gz, fq_filt)

    # Reads QC
    outdir = os.path.join(output_dir, 'nanoplot_raw')
    nanoplot(fq_raw, outdir)

    outdir = os.path.join(output_dir, 'nanoplot_filtered')
    nanoplot(fq_filt, outdir)

    # Genome assembly
    outdir = os.path.join(output_dir, 'Flye_assembly')
    contigs = flye_assembly(fq_filt, outdir)

    # CheckV
    outdir = os.path.join(output_dir, 'CheckV')
    checkv(contigs, outdir)


def bacterial_assembly_pipeline(input_file, output_dir):
    # Create output location
    name = os.path.basename(str(input_file)[:-4])
    out = os.path.join(output_dir, name)
    os.makedirs(out, exist_ok=False)

    # Convert
    fq_raw = os.path.join(output_dir, f'{name}_raw.fastq')
    convert_bam_to_fastq(input_file, fq_raw)

    # Remove adapters
    fq_trim = os.path.join(output_dir, f'{name}_trimmed.fastq')
    porechop_abi(fq_raw, fq_trim)

    # gzip file
    fq_trim_gz = gzip_file(fq_trim)

    # Filter
    fq_filt = os.path.join(output_dir, f'{name}_filtered.fastq')
    filtlong(fq_trim_gz, fq_filt)

    # Reads QC
    outdir = os.path.join(output_dir, 'nanoplot_raw')
    nanoplot(fq_raw, outdir)

    outdir = os.path.join(output_dir, 'nanoplot_filtered')
    nanoplot(fq_filt, outdir)

    # Genome assembly
    outdir = os.path.join(output_dir, 'Flye_assembly')
    contigs = flye_assembly(fq_filt, outdir)

    # CheckM
    print(contigs)


###_____________________________________________________BATCHES


def batch_phage_assembly_pipeline(input_dir, output_dir):
    for file in os.listdir(input_dir):
        try:
            phage_assembly_pipeline(file, output_dir)
        except Exception as e:
            print(f"ERROR {file}")
            continue


def batch_bacterial_assembly_pipeline(input_dir, output_dir):
    for file in os.listdir(input_dir):
        try:
            bacterial_assembly_pipeline(file, output_dir)
        except Exception as e:
            print(f"ERROR {file}")
            continue
