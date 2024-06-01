import os
from .functions import (
    configure_log,
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


def assembly_pipeline(input_file, output_dir, isolate='unknown',
                      logfile_location=None, logfile_configuration=None):
    # Setting logfile
    logger = configure_log(
        location=logfile_location,
        configuration=logfile_configuration
    )
    logger.info(f'Beginning assembly pipeline: {os.path.basename(input_file)}')

    # Check if isolate value is allowed
    if isolate not in ['phage', 'bacterial', 'fungal', 'unknown']:
        raise ValueError("Isolate must be: 'phage', 'bacterial', 'fungal' or 'unknown' (Skip filtering)")

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

    # Raw QC
    outdir = os.path.join(out, 'nanoplot_raw')
    nanoplot(fq_raw, outdir)

    # Trimmed QC
    outdir = os.path.join(out, 'nanoplot_trimmed')
    nanoplot(fq_trim, outdir)

    # gzip file
    fq_trim_gz = gzip_file(fq_trim)

    # Filter
    target = True
    if isolate == 'phage':
        target = 30000000
    elif isolate == 'bacterial':
        target = 500000000
    elif isolate == 'fungal':
        target = 10000000000  # Unsure of this one, check!
    elif isolate == 'unknown':
        target = False

    fq_filt = False
    if target:
        fq_filt = os.path.join(out, f'{name}_filtered.fastq')
        filtlong(fq_trim_gz, fq_filt, target)

        # Filtered QC
        outdir = os.path.join(out, 'nanoplot_filtered')
        nanoplot(fq_filt, outdir)

    # Genome assembly
    outdir = os.path.join(out, 'Flye_assembly')
    read_type = None
    contigs = False

    if fq_filt:
        print("Using filtered reads for assembly")
        read_type = 'filtered'
        contigs = flye_assembly(fq_filt, outdir, raise_on_fail=False)

    if not contigs:
        print("Filtered assembly failed. Using trimmed reads for assembly")
        read_type = 'trimmed'
        contigs = flye_assembly(fq_trim, outdir, raise_on_fail=False)

    if not contigs:
        print("Trimmed reads assembly failed. Using raw reads for assembly")
        read_type = 'raw'
        contigs = flye_assembly(fq_raw, outdir)

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


# _____________________________________________________BATCHES


def batch_assembly_pipeline(input_dir, output_dir, isolate=None,
                            logfile_location=None, logfile_configuration=None):
    # Setting logfile
    logger = configure_log(
        location=logfile_location,
        configuration=logfile_configuration
    )

    # Check inputs
    if not os.path.isdir(input_dir):
        e = f'Input directory {input_dir} is not a directory'
        logger.error(e)
        raise ValueError(e)
    else:
        e = f'Beginning assembly pipeline: input_dir: {input_dir} output_dir: {output_dir}'
        logger.info(e)

    # Phage default isolate
    if isolate is None:
        logger.warning('Isolate type not specified, defaulting to phage')
        isolate = 'phage'

    # Batch pipeline
    input_dir = os.path.expanduser(input_dir)
    output_dir = os.path.expanduser(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        path = os.path.join(input_dir, file)
        try:
            assembly_pipeline(
                input_file=path,
                output_dir=output_dir,
                logfile_location=logfile_location,
                logfile_configuration=logfile_configuration,
                isolate=isolate
            )
        except Exception as e:
            print(f"ERROR {e}")
            continue
