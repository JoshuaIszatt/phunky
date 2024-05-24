import os
from .functions import (
    convert_bam_to_fastq,
    porechop_abi,
    gzip_file,
    filtlong,
    nanoplot,
    flye_assembly,
    checkv,
    read_mapping,
    extract_contig_header,
    generate_coverage_graph
)


###_____________________________________________________PIPELINES


def phage_assembly_pipeline(input_file, output_dir):
    # Create output location
    name = os.path.basename(str(input_file)[:-4])
    out = os.path.join(output_dir, name)
    os.makedirs(out, exist_ok=False)

    # Convert
    fq_raw = os.path.join(out, f'{name}_raw.fastq')
    convert_bam_to_fastq(input_file, fq_raw)

    # Remove adapters
    fq_trim = os.path.join(out, f'{name}_trimmed.fastq')
    porechop_abi(fq_raw, fq_trim)

    # gzip file
    fq_trim_gz = gzip_file(fq_trim)

    # Filter
    fq_filt = os.path.join(out, f'{name}_filtered.fastq')
    filtlong(fq_trim_gz, fq_filt)

    # Reads QC
    outdir = os.path.join(out, 'nanoplot_raw')
    nanoplot(fq_raw, outdir)

    outdir = os.path.join(out, 'nanoplot_filtered')
    nanoplot(fq_filt, outdir)

    # Genome assembly
    outdir = os.path.join(out, 'Flye_assembly')
    contigs = flye_assembly(fq_filt, outdir)

    # CheckV
    if os.getenv('CHECKVDB'):
        outdir = os.path.join(out, 'CheckV')
        checkv(contigs, outdir)

    # Read mapping
    fa_filt = os.path.join(out, f'{name}_filtered.fasta')
    convert_bam_to_fastq(fq_filt, fa_filt)

    outdir = os.path.join(out, 'Read_mapping')
    basecov = read_mapping(
        contigs_fasta=contigs,
        reads=fa_filt,
        output_directory=outdir
    )[0]

    # Using basecov.tsv and header to generate coverage graph
    header = extract_contig_header(contigs)[0]
    generate_coverage_graph(
        header=header,
        basecov=basecov,
        output_directory=out)


def bacterial_assembly_pipeline(input_file, output_dir):
    # Create output location
    name = os.path.basename(str(input_file)[:-4])
    out = os.path.join(output_dir, name)
    os.makedirs(out, exist_ok=False)

    # Convert
    fq_raw = os.path.join(out, f'{name}_raw.fastq')
    convert_bam_to_fastq(input_file, fq_raw)

    # Remove adapters
    fq_trim = os.path.join(out, f'{name}_trimmed.fastq')
    porechop_abi(fq_raw, fq_trim)

    # gzip file
    fq_trim_gz = gzip_file(fq_trim)

    # Filter
    fq_filt = os.path.join(out, f'{name}_filtered.fastq')
    filtlong(fq_trim_gz, fq_filt, target_bases=500000000)

    # Reads QC
    outdir = os.path.join(out, 'nanoplot_raw')
    nanoplot(fq_raw, outdir)

    outdir = os.path.join(out, 'nanoplot_filtered')
    nanoplot(fq_filt, outdir)

    # Genome assembly
    outdir = os.path.join(out, 'Flye_assembly')
    contigs = flye_assembly(fq_filt, outdir)

    # Read mapping
    outdir = os.path.join(out, 'Read_mapping')
    basecov = read_mapping(
        contigs_fasta=contigs,
        reads=fq_filt,
        output_directory=out
    )[0]

    # Using basecov.tsv and header to generate coverage graph
    header = extract_contig_header(contigs)[0]
    generate_coverage_graph(
        header=header,
        basecov=basecov,
        output_directory=out)


###_____________________________________________________BATCHES


def batch_phage_assembly_pipeline(input_dir, output_dir):
    for file in os.listdir(input_dir):
        if not file.endswith('.bam'):
            print(f"Skipping {file}")
            continue
        path = os.path.join(input_dir, file)
        try:
            phage_assembly_pipeline(path, output_dir)
        except Exception as e:
            print(f"ERROR {e}")
            continue


def batch_bacterial_assembly_pipeline(input_dir, output_dir):
    for file in os.listdir(input_dir):
        if not file.endswith('.bam'):
            print(f"Skipping {file}")
        path = os.path.join(input_dir, file)
        try:
            bacterial_assembly_pipeline(path, output_dir)
        except Exception as e:
            print(f"ERROR {e}")
            continue
