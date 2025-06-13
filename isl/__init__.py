import argparse
from lib.ComputeDescriptors import ComputeDescriptors


def main():
    parser = argparse.ArgumentParser(description="ISL command-line tool")
    parser.add_argument("input_file", type=str, help="Path to input *csv file.")
    parser.add_argument("output_directory", type=str, help="Path to output directory")
    parser.add_argument(
        "--output_file_name",
        type=str,
        default="output_descr.csv",
        help="Optional output file name",
    )
    args = parser.parse_args()

    ComputeDescriptors(args.input_file, args.output_directory, args.output_file_name)()
