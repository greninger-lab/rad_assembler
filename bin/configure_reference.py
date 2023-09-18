#!/usr/bin/env python

import argparse
import logging
import sys
import os
import json
from pathlib import Path
from Bio import SeqIO


logger = logging.getLogger()

def generate_fasta(genbank_ref):
    ref_name = os.path.basename(os.path.realpath(genbank_ref))
    genbank_ref_fasta_name = str(ref_name).replace('.gb', '.fasta')
    SeqIO.convert(genbank_ref, 'genbank', genbank_ref_fasta_name, 'fasta')
    return genbank_ref_fasta_name

def extract_regions(region_map, genbank_ref):
    ref_name = generate_fasta(genbank_ref)
    seq_records = SeqIO.parse(ref_name, "fasta")
    ref_base_name = ref_name.replace(".fasta", "")
    with open(region_map, encoding='utf-8') as region_file:
        region_info = json.load(region_file)
        for record in seq_records:
            for region in region_info["region_order"]:
                start = region_info["regions"][region["region"]][0]
                end = region_info["regions"][region["region"]][1]
                region_fasta_name = "region_" + region["region"] + "_" + ref_base_name + ".fasta"
                with open(region_fasta_name, 'w') as out:
                    out.write(">region_" + region["region"] + "_" + ref_base_name + "\n")
                    out.write(str(record.seq[start:end]))

def build_reference(region_map, genbank_ref, sample_name):
    ref_name = os.path.basename(os.path.realpath(genbank_ref))
    genbank_ref_name = str(ref_name).replace('.gb', '')
    new_ref_sequence = ">" + sample_name + "_consensus\n"
    with open(region_map, encoding='utf-8') as region_map_file:
        region_info = json.load(region_map_file)
        for reg in region_info["region_order"]:
            path_to_region_file = sample_name + "_" + str(reg["region"]) + "_" + genbank_ref_name + "_consensus.fasta"
            if os.path.isfile(path_to_region_file):
                with open(path_to_region_file, "r") as region_fasta:
                    seq_records = SeqIO.parse(region_fasta, "fasta")
                    for record in seq_records:
                        if reg["reverse"]:
                            sequence = str(record.seq.reverse_complement())
                        else:
                            sequence = str(record.seq)               
                        new_ref_sequence += sequence
        new_ref_file_name = sample_name + "_new_ref_consensus.fasta"
        with open(new_ref_file_name, "w") as new_ref_file:
            new_ref_file.write(new_ref_sequence + "\n")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Set genbank reference and extract regions for region based reference building.",
        epilog="Example: python configure_reference.py -r HSV1_NC_001806.json NC_001806.2.gb",
    )
    parser.add_argument(
        "genbank_ref",
        metavar="GENBANK_REF",
        type=Path,
        help="Genbank eferene (.gb) file.",
    )
    parser.add_argument(
        "-r",
        "--region_map",
        nargs='?',
        const='arg_was_not_given',
        help="A json file defining regions to extract for individual scaffold mapping.",
        default=None,
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        nargs='?',
        help="Name of sample (required if using --build_reference)",
        default=None,
    )    
    parser.add_argument(
        "--build_reference", 
        action="store_true",
        help="Concatenate regions into new reference"
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.genbank_ref.is_file():
        logger.error(f"The given input file {args.genbank_ref} was not found!")
        sys.exit(2)
                                                                                                                           
    if args.region_map is None and args.genbank_ref:
        generate_fasta(args.genbank_ref)
    elif args.region_map == 'arg_was_not_given':
        print('Option given, but no command-line argument: "-r, --region_map"')
    elif args.region_map and args.genbank_ref and not args.build_reference:
        extract_regions(args.region_map, args.genbank_ref)
    elif args.region_map and args.genbank_ref and args.build_reference and not args.sample_name:
        print("--sample_name SAMPLE_NAME must be provided when using --build_reference")
    elif args.region_map and args.genbank_ref and args.build_reference and args.sample_name:        
        build_reference(args.region_map, args.genbank_ref, args.sample_name)


if __name__ == "__main__":
    sys.exit(main())
