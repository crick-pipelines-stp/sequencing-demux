#!/usr/bin/env python3

"""
Checks input samplesheet for common errors
"""

import sys
import argparse
import platform
import shutil

HEADER = ["sample_id", "sample_name", "group", "user", "project_id", "barcode"]


def check_samplesheet(samplesheet_path, output_path):
    """
    Checks input samplesheet for common errors
    """

    with open(samplesheet_path, "r", encoding="UTF-8") as fin:
        # Read header
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        # Check header length
        if len(HEADER) != len(header):
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        # Check header content
        if HEADER != header:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        # Cycle through samplesheet line by line
        line_no = 2
        for line in fin:
            data_line = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check if its just a blank line so we dont error
            if line.strip() == "":
                continue

            # Unpack values
            # sample_id, group, user, project_id, barcode = data_line

            # Check valid number of columns per row
            if len(data_line) != len(HEADER):
                print(f"Invalid number of columns (found {len(data_line)} should be {len(HEADER)})! - line no. {line_no}")
                sys.exit(1)

            # Check number of filled columns
            num_cols = len([x for x in data_line if x])
            if num_cols != len(HEADER):
                print(f"Invalid number of filled columns (found {num_cols} should be {len(HEADER)})! - line no. {line_no}")
                sys.exit(1)

            # Loop admin
            line_no = line_no + 1

    # Write sample sheet as a valid file
    shutil.copy(samplesheet_path, output_path)


def dump_versions(process_name):
    """
    Output versions file
    """
    with open("versions.yml", "w", encoding="UTF-8") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")


def main(args):
    """
    Main input function
    """

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--sample", default="!{sample}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args(args)

    print(args)

    # Run main function
    check_samplesheet(args.sample, args.output)

    # Dump versions
    dump_versions(args.process_name)


if __name__ == "__main__":
    main(sys.argv[1:])
