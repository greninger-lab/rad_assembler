import pandas as pd
import numpy as np
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, os
import subprocess
import logging
from enum import Enum


# 4 command line args are required - scaffolds_file, sample_id, fastq_1, fastq_2 
scaffolds_file = sys.argv[1]
meta_id = sys.argv[2]
fastq_1 = sys.argv[3]
fastq_2 = sys.argv[4]

class ExtendDirection(Enum):
    LEFT = 1
    RIGHT = 2

def extend_contig_map(contig_record, extend_direction):
    base_name = meta_id + "_" + contig_record.name
    fasta_name = base_name + ".fasta"
    base_name_extended = base_name + "_extended"
    contig_length = len(contig_record.seq)
    with open(fasta_name, "w") as contig_out:
        SeqIO.write(contig_record, "fasta")
        sam_out = open(base_name_extended + ".sam", "w")
        subprocess.run(["mkdir", "bowtie2"])        
        subprocess.run(["bowtie2-build", fasta_name, "bowtie2/" + base_name])
        subprocess.run(["bowtie2", "-x", "bowtie2/" + base_name, "-1 " + fastq_1 + " -2 " + fastq_2], stdout=sam_out)
        subprocess.run(["samtools", "view", "-bS", "-o", base_name_extended + ".bam", base_name_extended + ".sam"])
        subprocess.run(["samtools", "sort", base_name_extended + ".bam", base_name_extended + ".sorted"])
        mpileup_out = open(base_name_extended + ".mpileup", "w")
        mpileup_ps = subprocess.run(["samtools", "mpileup", "-d", "5000", "-A", "-Q", "0", base_name_extended + ".sorted.bam"], stdout=mpileup_out)
        subprocess.run(["ivar", "consensus", "-p", base_name_extended, "-n", "'N'", "-m", "5", "-t", "0"], stdin=mpileup_ps)
        # remove leading and trailing Ns
        # subprocess.run(["sed", "-r", "'/^>/! s/n+$|^n+//g'", base_name_extended + ".fa"])
        new_record = next(SeqIO.parse(base_name_extended + ".fa", "fasta"))
        new_seq_str = str(new_record.seq)        
        if extend_direction == ExtendDirection.BOTH:
            return new_record
        elif extend_direction == ExtendDirection.LEFT:
            begin_original_contig = new_seq_str.index(str(contig_record.seq))
            end_original_contig = begin_original_contig + contig_length
            left_record = SeqRecord(
                new_record.seq[:end_original_contig],
                id=contig_record.name,
                name=contig_record.name,
                description="",
            )
            return left_record
        elif extend_direction == ExtendDirection.RIGHT:
            begin_original_contig = new_seq_str.index(str(contig_record.seq))
            end_original_contig = begin_original_contig + contig_length
            right_record = SeqRecord(
                new_record.seq[begin_original_contig:],
                id=contig_record.name,
                name=contig_record.name,
                description="",
            )
            return right_record

def extend_contig(contig_record, extend_direction):
    contig_length = len(contig_record.seq)
    if contig_length >= 4000:
        truncate_seq = 4000
    else:
        truncate_seq = contig_length
    if extend_direction == ExtendDirection.LEFT:
        left_tag = SeqRecord(
            contig_record.seq[:truncate_seq],
            id=contig_record.name + "_left",
            name=contig_record.name + "_left",
            description="",
        )
        contig_left_extended = SeqIO.parse(extend_contig(left_tag, ExtendDirection.LEFT))
        extended_contig_str = str(contig_left_extended.seq) + str(contig_record.seq[truncate_seq:])
    elif extend_direction == ExtendDirection.RIGHT:
        right_tag = SeqRecord(
            contig_record.seq[truncate_seq:],
            id=contig_record.name + "_right",
            name=contig_record.name + "_right",
            description="",
        )
        contig_right_extended = extend_contig(right_tag, ExtendDirection.RIGHT)
        extended_contig_str = str(contig_record.seq[:4000]) + str(contig_right_extended.seq)



    # else:
    #     both_tag = SeqRecord(
    #         contig_record.seq,
    #         id=contig_record.name + "_both",
    #         name=contig_record.name + "_both",
    #         description="",
    #     )
    #     contig_both_extended = extend_contig(both_tag, ExtendDirection.BOTH)