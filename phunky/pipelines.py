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
    extract_contig,
    generate_coverage_graph
)


# _____________________________________________________PIPELINES


def assembly_pipeline(input_file, output_dir, isolate='phage',
                      logger=None, logfile_location=None, logfile_configuration=None):
    # Setting logger if None has been passed
    if logger is None:
        if logfile_location is None:
            logfile_location = output_dir
        logger = configure_log(
            location=logfile_location,
            configuration=logfile_configuration
        )

    # Attempt to use logger
    try:
        logger.info(f'Beginning assembly pipeline: {os.path.basename(input_file)}')
    except Exception as e:
        logger.error(e)
        raise Exception(f"Logging error: {e}")

    # Check if isolate value is allowed
    if isolate == 'phage':
        target = 30000000  # Approx 100 X coverage for a 300kb genome
    elif isolate == 'bacterial':
        target = 500000000  # Approx 100 X coverage for a 5mb genome
    elif isolate == 'fungal':
        target = 5000000000  # Approx 100 X coverage for a 50mb genome
    elif isolate == 'unknown':
        target = 10000000000
    elif isinstance(isolate, (int, float)):
        target = int(isolate)
    else:
        raise ValueError("Isolate can be: 'phage', 'bacterial', 'fungal', 'unknown', or a numeric value")

    # Log
    logger.info(f"Target bases set to {target}")

    # Create output location
    extensions = ['.bam', '.fastq', '.fastq.gz']
    basename = os.path.basename(str(input_file))
    logger.debug(f'Basename: {basename}')
    name = out = None
    for extension in extensions:
        if basename.endswith(extension):
            name = basename[:-len(extension)]
            logger.debug(f'Name: {name}')
            out = os.path.join(output_dir, name)
            logger.info(f'File type ({extension}) accepted: output_directory: {out}')
            os.makedirs(out, exist_ok=False)
            break
    else:
        logger.error(f"File type rejected: {basename}")
        raise Exception("File type rejected")

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
    logger.info(f"Path to raw data: {fq_raw}")

    # Remove adapters
    logger.info("Trimming adaptors from raw reads")
    fq_trim = os.path.join(out, f'{name}_trimmed.fastq')
    porechop_abi(fq_raw, fq_trim)
    logger.info(f"Path to trimmed data: {fq_trim}")

    # Raw QC
    outdir = os.path.join(out, 'nanoplot_raw')
    nanoplot(fq_raw, outdir)
    logger.info("Raw reads QC complete")

    # Trimmed QC
    outdir = os.path.join(out, 'nanoplot_trimmed')
    nanoplot(fq_trim, outdir)
    logger.info("Trimmed reads QC complete")

    # gzip file
    fq_trim_gz = gzip_file(fq_trim)

    # Filter
    logger.info(f"Filtering data to {target}bp")
    fq_filt = os.path.join(out, f'{name}_filtered.fastq')
    filtlong(fq_trim_gz, fq_filt,
             target_bases=target)
    logger.info(f"Path to filtered data: {fq_filt}")

    # Filtered QC
    outdir = os.path.join(out, 'nanoplot_filtered')
    nanoplot(fq_filt, outdir)
    logger.info("Filtered reads QC complete")

    # Genome assembly
    logger.info("Beginning genome assembly...")
    outdir = os.path.join(out, 'Flye_assembly')
    read_type = None
    contigs = False

    if fq_filt:
        logger.info("Using filtered reads for assembly")
        read_type = 'filtered'
        contigs = flye_assembly(fq_filt, outdir, raise_on_fail=False)

    if not contigs:
        logger.info("Filtered assembly failed. Using trimmed reads for assembly")
        read_type = 'trimmed'
        contigs = flye_assembly(fq_trim, outdir, raise_on_fail=False)

    if not contigs:
        logger.info("Trimmed reads assembly failed. Using raw reads for assembly")
        read_type = 'raw'
        contigs = flye_assembly(fq_raw, outdir)

    # Read mapping
    logger.info("Mapping reads to the assembly")
    fa_filt = os.path.join(out, f'{name}_{read_type}.fasta')
    convert_bam_to_fastq(fq_filt, fa_filt)
    outdir = os.path.join(out, 'Read_mapping')
    basecov = read_mapping(
        contigs_fasta=contigs,
        reads=fa_filt,
        output_directory=outdir
    )[0]

    # Using basecov.tsv and header to generate coverage graph
    logger.info("Generating coverage graph")
    header, seq_len = extract_contig_header(contigs)
    generate_coverage_graph(
        header=header,
        basecov=basecov,
        output_directory=out)

    # Extracting contig
    name = f"{name}_{seq_len}bp"
    file = os.path.join(out, name)
    extract_contig(
        contigs_fasta=contigs,
        header=header,
        output_file=file,
        rename=name
    )

    # CheckV
    if os.getenv('CHECKVDB'):
        logger.info(f"CHECKVDB variable detected, running analysis: {os.getenv('CHECKVDB')}")
        if not os.path.isdir(os.getenv('CHECKVDB')):
            logger.error(f"CheckV database variable does not exist: {os.getenv('CHECKVDB')}")
        else:
            try:
                outdir = os.path.join(out, 'CheckV')
                checkv(contigs, outdir)
            except Exception as e:
                logger.warning(e)
    else:
        logger.info(f"No Checkv database detected, skipping...")

    # Finish pipeline
    logger.info(f"Putative genome: {name}")
    logger.info(f"Pipeline complete:")

# _____________________________________________________BATCHES


def batch_assembly_pipeline(input_dir, output_dir, isolate=None,
                            logger=None, logfile_location=None, logfile_configuration=None):
    # Batch pipeline
    input_dir = os.path.expanduser(input_dir)
    output_dir = os.path.expanduser(output_dir)

    # Setting logger if None has been passed
    if logger is None:
        if logfile_location is None:
            logfile_location = output_dir
        logger = configure_log(
            location=logfile_location,
            configuration=logfile_configuration
        )

    # Setting in/out to log
    logger.info(f"input_dir: {input_dir}")
    logger.info(f"output_dir: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    count = 0

    # Check inputs
    if not os.path.isdir(input_dir):
        e = f'Input directory {input_dir} is not a directory'
        logger.error(e)
        raise ValueError(e)
    else:
        e = f'Batch assembly pipeline:'
        logger.info(e)

    # Phage default isolate
    if isolate is None:
        logger.warning('Isolate type not specified, defaulting to phage')
        isolate = 'phage'
    else:
        logger.info(f'Isolate type set to: {isolate}')

    # Processing
    for file in os.listdir(input_dir):
        count = count+1
        logger.info(f"Running input file #{count}")
        path = os.path.join(input_dir, file)
        try:
            assembly_pipeline(
                input_file=path,
                output_dir=output_dir,
                logger=logger,
                isolate=isolate
            )
        except Exception as e:
            logger.error(f"Pipeline failure ({file}): {e}")
            continue

    # Finish batch pipeline
    logger.info("Batch pipeline complete")
