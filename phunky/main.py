import os
import argparse
from phunky.pipelines import (
    phage_assembly_pipeline,
    bacterial_assembly_pipeline,
    batch_phage_assembly_pipeline,
    batch_bacterial_assembly_pipeline
)


###_____________________________________________________ARGS


def parse_args():
    parser = argparse.ArgumentParser(description="Phage and Bacterial Assembly Pipeline")
    parser.add_argument("input_file", help="Path to input BAM file or directory")
    parser.add_argument("output_dir", help="Path to output directory")
    parser.add_argument(
        "--pipeline",
        choices=["phage", "bacterial"],
        required=True,
        help="Choose 'phage' or 'bacterial' pipeline",
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Process all BAM files in the input directory",
    )
    return parser.parse_args()


###_____________________________________________________RUN


if __name__ == "__main__":
    args = parse_args()

    if args.pipeline == "phage":
        phage_assembly_pipeline(args.input_file, args.output_dir)
    elif args.pipeline == "bacterial":
        bacterial_assembly_pipeline(args.input_file, args.output_dir)
    else:
        print("Invalid pipeline choice. Use 'phage' or 'bacterial'.")
