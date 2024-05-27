import os
import sys
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


# _____________________________________________________PIPELINES


def assembly_pipeline(input_file, output_dir, isolate='phage'):
    # Check if isolate value is allowed
    if isolate not in ['phage', 'bacterial', 'fungal']:
        raise ValueError("Isolate must be one of 'phage', 'bacterial', or 'fungal' ")

    # Create output location
    extensions = ['.bam', '.fastq', '.fastq.gz']
    basename = os.path.basename(str(input_file))
    name = out = None
    for extension in extensions:
        if basename.endswith(extension):
            name = basename[:-len(extension)]
            out = os.path.join(output_dir, name)
            os.makedirs(out, exist_ok=False)
            break
        else:
            raise Exception("File type not accepted")

    # Ensure variables are set
    if out is None or name is None:
        raise Exception("Could not determine output location or filename")

    # Convert if required
    if input_file.endswith('.bam'):
        fq_raw = os.path.join(out, f'{name}_raw.fastq')
        convert_bam_to_fastq(input_file, fq_raw)
    elif input_file.endswith('.fastq.gz'):
        fq_raw = os.path.join(out, f'{name}_raw.fastq')
        convert_bam_to_fastq(input_file, fq_raw)
    else:
        fq_raw = input_file

    # Remove adapters
    fq_trim = os.path.join(out, f'{name}_trimmed.fastq')
    porechop_abi(fq_raw, fq_trim)

    # gzip file
    fq_trim_gz = gzip_file(fq_trim)

    # Filter
    target = 20000000
    if isolate == 'phage':
        target = 20000000
    elif isolate == 'bacterial':
        target = 500000000
    elif isolate == 'fungal':
        target = 10000000000  # Unsure of this one, check!
    fq_filt = os.path.join(out, f'{name}_filtered.fastq')
    filtlong(fq_trim_gz, fq_filt, target)

    # Reads QC
    outdir = os.path.join(out, 'nanoplot_raw')
    nanoplot(fq_raw, outdir)

    outdir = os.path.join(out, 'nanoplot_trimmed')
    nanoplot(fq_trim, outdir)

    outdir = os.path.join(out, 'nanoplot_filtered')
    filt_bases = nanoplot(fq_filt, outdir)

    # Genome assembly
    outdir = os.path.join(out, 'Flye_assembly')
    if int(filt_bases) > 1000000:
        print("Using filtered reads for assembly")
        read_type = 'filtered'
        contigs = flye_assembly(fq_filt, outdir)
    else:
        print("Using trimmed reads for assembly")
        read_type = 'trimmed'
        contigs = flye_assembly(fq_trim, outdir)

    # Read mapping
    fa_filt = os.path.join(out, f'{name}_{read_type}.fasta')
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

    # CheckV
    if os.getenv('CHECKVDB'):
        outdir = os.path.join(out, 'CheckV')
        checkv(contigs, outdir)

# todo: Create a summary file

# _____________________________________________________BATCHES


def batch_phage_assembly_pipeline(input_dir, output_dir):
    for file in os.listdir(input_dir):
        path = os.path.join(input_dir, file)
        try:
            assembly_pipeline(path, output_dir, 'phage')
        except Exception as e:
            print(f"ERROR {e}")
            continue


def batch_bacterial_assembly_pipeline(input_dir, output_dir):
    for file in os.listdir(input_dir):
        path = os.path.join(input_dir, file)
        try:
            assembly_pipeline(path, output_dir, 'bacterial')
        except Exception as e:
            print(f"ERROR {e}")
            continue
