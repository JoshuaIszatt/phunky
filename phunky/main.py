import os
import sys
import argparse
from phunky.pipelines import (
    assembly_pipeline,
    batch_assembly_pipeline
)


# _____________________________________________________ARGS


def parse_args():
    parser = argparse.ArgumentParser(description="Phage and Bacterial Assembly Pipeline")
    parser.add_argument("-i", "--input_file", required=True, 
                        help="Path to input BAM file or directory")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Path to output directory")
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


# _____________________________________________________RUN


def main():
    args = parse_args()
    if args.batch:
        if args.pipeline == "phage":
            batch_assembly_pipeline(args.input_file, args.output_dir)
        elif args.pipeline == "bacterial":
            batch_assembly_pipeline(args.input_file, args.output_dir)
        else:
            print("Invalid pipeline choice. Use 'phage' or 'bacterial'.")
            sys.exit(0)
        print("Running a Phunky pipeline")
    else:
        if args.pipeline == "phage":
            print("Running Phunky pipelines")
            assembly_pipeline(args.input_file, args.output_dir)
        elif args.pipeline == "bacterial":
            print("Running Phunky pipelines")
            assembly_pipeline(args.input_file, args.output_dir)
        else:
            print("Invalid pipeline choice. Use 'phage' or 'bacterial'.")


if __name__ == "__main__":
    main()
