import pandas as pd
import pysam
import numpy as np
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, os
import subprocess
import itertools
import logging
from enum import Enum


# 5 command line args are required - scaffolds_file, sample_id, fastq_1, fastq_2, cpus
scaffolds_file = sys.argv[1]
meta_id = sys.argv[2]
fastq_1 = sys.argv[3]
fastq_2 = sys.argv[4]
cpus = sys.argv[5]

if ".gz" in fastq_1:
    subprocess.run(["gunzip", "-f", fastq_1])
    fastq_1 = fastq_1[:fastq_1.index(".gz")]
if ".gz" in fastq_2:
    subprocess.run(["gunzip", "-f", fastq_2])
    fastq_2 = fastq_2[:fastq_2.index(".gz")]

# Curated list of HCMV genomes and the regions extracted from their published 
# annotation.  Region info may or may not be trustworthy. 
gb_genomes = "/genome_identification/cmv/cmv_genomes.gb"
gb_genomes_regions = "/genome_identification/cmv/cmv_genomes_regions.json"

logging.basicConfig(filename=meta_id + '_build_reference.log', level=logging.DEBUG, filemode = 'w')

# Pull genome regions into a dictionary for easy lookup
regions_dict = json.load(open(gb_genomes_regions,))

# Pull any genome genbank record by name
def get_genbank_record(name):
    for gb_record in SeqIO.parse(open(gb_genomes,"r"), "genbank"):
        if gb_record.name in name:
            return gb_record
    return None

# Gather information from a standard Merlin reference to use as a
# defining, consistent coordinate system for large regions.
merlin_ref_name = "NC_006273"

merlin_record = get_genbank_record(merlin_ref_name)
merlin_regions = regions_dict[merlin_ref_name]
merlin_ul_start = merlin_regions["regions"]["UL"][0]
merlin_ul_end = merlin_regions["regions"]["UL"][1]
merlin_us_start = merlin_regions["regions"]["US"][0]
merlin_us_end = merlin_regions["regions"]["US"][1]
merlin_pretrs1_start = merlin_regions["regions"]["PRETRS1"][0]
merlin_pretrs1_end = merlin_regions["regions"]["PRETRS1"][1]

# Create fasta files of the flanking 300bp regions for the "unique regions" 
merlin_ul_start_file_name = merlin_ref_name + "_ul_start.fasta"
merlin_ul_end_file_name = merlin_ref_name + "_ul_end.fasta"
merlin_us_start_file_name = merlin_ref_name + "_us_start.fasta"
merlin_us_end_file_name = merlin_ref_name + "_us_end.fasta"
merlin_pretrs1_start_file_name = merlin_ref_name + "_pretrs1_start.fasta"
merlin_pretrs1_end_file_name = merlin_ref_name + "_pretrs1_end.fasta"
merlin_complete_file_name = merlin_ref_name + "_complete.fasta"

with open(merlin_complete_file_name, "w") as merlin_complete_outfile:
    SeqIO.write(merlin_record, merlin_complete_outfile, 'fasta')

with open(merlin_ul_end_file_name, "w") as ul_end_outfile:
    ul_end_record = merlin_record[merlin_ul_end-300:merlin_ul_end]
    SeqIO.write(ul_end_record, ul_end_outfile, 'fasta')
with open(merlin_ul_start_file_name, "w") as ul_start_outfile:
    ul_start_record = merlin_record[merlin_ul_start:merlin_ul_start + 300]
    SeqIO.write(ul_start_record, ul_start_outfile, 'fasta')    

with open(merlin_us_end_file_name, "w") as us_end_outfile:
    us_end_record = merlin_record[merlin_us_end-300:merlin_us_end]
    SeqIO.write(us_end_record, us_end_outfile, 'fasta')
with open(merlin_us_start_file_name, "w") as us_start_outfile:
    us_start_record = merlin_record[merlin_us_start:merlin_us_start + 300]
    SeqIO.write(us_start_record, us_start_outfile, 'fasta')

with open(merlin_pretrs1_end_file_name, "w") as pretrs1_end_outfile:
    pretrs1_end_record = merlin_record[merlin_pretrs1_end:merlin_pretrs1_end + 300]
    SeqIO.write(pretrs1_end_record, pretrs1_end_outfile, 'fasta')
with open(merlin_pretrs1_start_file_name, "w") as pretrs1_start_outfile:
    pretrs1_start_record = merlin_record[merlin_pretrs1_start-300:merlin_pretrs1_start]
    SeqIO.write(pretrs1_start_record, pretrs1_start_outfile, 'fasta')    

# return the longest of the top query hits
def top_hit(file_name):
    try:
        df = pd.read_csv(file_name,sep='\t',header=None)
        return df.iloc[df[11].idxmax()][1]
    except pd.errors.EmptyDataError:
        return None
        logging.debug("Failed to get top hit from " + file_name + ".")

# write region file
def write_region(region_name, sequence):
    with open(region_name + ".fasta", 'w') as out:
        out.write(">" + region_name + "\n" + sequence)

# Get the top hit for the hyper variable region of UL and top hit for the 'a region' and TRS1
hyper_hit_file = sys.argv[1] + ".hyper.txt"
hyper_hit = top_hit(hyper_hit_file)
hyper_hit_regions = regions_dict[hyper_hit[:hyper_hit.index(".")]]
hyper_hit_record = get_genbank_record(hyper_hit)
hyper_complete_file_name = hyper_hit + "_complete.fasta"
with open(hyper_complete_file_name, "w") as hyper_complete_outfile:
    SeqIO.write(hyper_hit_record, hyper_complete_outfile, 'fasta')

# Extract UL_end from hyper_hit to use for IRL extraction in the case of UL substitution
hyper_hit_ul_end = hyper_hit_regions["regions"]["UL"][1]
hyper_hit_ul_start = hyper_hit_regions["regions"]["UL"][0]
hyper_hit_ul_end_file_name = hyper_hit_record.name + "_ul_end.fasta"
hyper_hit_ul_start_file_name = hyper_hit_record.name + "_ul_start.fasta"
with open(hyper_hit_ul_end_file_name, "w") as hyper_ul_end_outfile:
    hyper_hit_ul_end_record = hyper_hit_record[hyper_hit_ul_end-300:hyper_hit_ul_end]
    SeqIO.write(hyper_hit_ul_end_record, hyper_ul_end_outfile, 'fasta')
with open(hyper_hit_ul_start_file_name, "w") as hyper_ul_start_outfile:
    hyper_hit_ul_start_record = hyper_hit_record[hyper_hit_ul_start:hyper_hit_ul_start+300]
    SeqIO.write(hyper_hit_ul_start_record, hyper_ul_start_outfile, 'fasta')

irl_irs_hit_file = sys.argv[1] + ".irl_irs.txt"
irl_irs_hit = top_hit(irl_irs_hit_file)
irl_irs_hit_regions = regions_dict[irl_irs_hit[:irl_irs_hit.index(".")]]
irl_irs_hit_record = get_genbank_record(irl_irs_hit)
irl_irs_start = irl_irs_hit_regions["regions"]["IRL_IRS"][0]
irl_irs_end = irl_irs_hit_regions["regions"]["IRL_IRS"][1]
irl_irs_only_record = str(irl_irs_hit_record.seq[irl_irs_start:irl_irs_end])
with open(meta_id + "_irl_irs_hit.fasta", 'w') as irl_irs_out:
    irl_irs_out.write(">" + irl_irs_hit_record.name + "\n" + irl_irs_only_record)

irl_irs_complete_file_name = irl_irs_hit + "_complete.fasta"
with open(irl_irs_complete_file_name, "w") as irl_irs_complete_outfile:
    SeqIO.write(irl_irs_hit_record, irl_irs_complete_outfile, 'fasta')

trl_hit_file = sys.argv[1] + ".trl.txt"
trl_hit = top_hit(trl_hit_file)
if not trl_hit:
    trl_hit = irl_irs_hit
trl_hit_regions = regions_dict[trl_hit[:trl_hit.index(".")]]
trl_hit_record = get_genbank_record(trl_hit)
trl_start = trl_hit_regions["regions"]["TRL"][0]
trl_end = trl_hit_regions["regions"]["TRL"][1]
trl_only_record = str(trl_hit_record.seq[trl_start:trl_end])
with open(meta_id + "_trl_hit.fasta", 'w') as trl_out:
    trl_out.write(">" + trl_hit_record.name + "\n" + trl_only_record)  

trs_hit_file = sys.argv[1] + ".trs.txt"
trs_hit = top_hit(trs_hit_file)
if not trs_hit:
    trs_hit = irl_irs_hit
trs_hit_regions = regions_dict[trs_hit[:trs_hit.index(".")]]
trs_hit_record = get_genbank_record(trs_hit)
trs_start = trs_hit_regions["regions"]["TRS"][0]
trs_end = trs_hit_regions["regions"]["TRS"][1]
trs_only_record = str(trs_hit_record.seq[trs_start:trs_end])
with open(meta_id + "_trs_hit.fasta", 'w') as trs_out:
    trs_out.write(">" + trs_hit_record.name + "\n" + trs_only_record)  

a_hit_file = sys.argv[1] + ".a.txt"
a_hit = top_hit(a_hit_file)
if not a_hit:
    a_hit = irl_irs_hit
a_hit_regions = regions_dict[a_hit[:a_hit.index(".")]]
a_hit_record = get_genbank_record(a_hit)
a_start = a_hit_regions["regions"]["A"][0]
a_end = a_hit_regions["regions"]["A"][1]
a_only_record = str(a_hit_record.seq[a_start:a_end])
with open(meta_id + "_a_hit.fasta", 'w') as a_out:
    a_out.write(">" + a_hit_record.name + "\n" + a_only_record)  

trs1_hit_file = sys.argv[1] + ".trs1.txt"
trs1_hit = top_hit(trs1_hit_file)    
trs1_record = get_genbank_record(trs1_hit)
trs1_regions = regions_dict[trs1_hit[:trs1_hit.index(".")]]  
trs1_start = trs1_regions["regions"]["TRS1"][0]
trs1_end = trs1_regions["regions"]["TRS1"][1]
trs1_only_record = str(trs1_record.seq[trs1_start:trs1_end])

with open(meta_id + "_trs1_hit.fasta", 'w') as out:
    out.write(">" + trs1_record.name + "\n" + trs1_only_record)
trs1_start_tag_file_name = trs1_record.name + "_trs1_start_tag.fasta"
with open(trs1_start_tag_file_name, "w") as trs1_start_tag_outfile:
    trs1_start_tag_outfile.write(">" + trs1_record.name + "_trs1_start_tag\n" + trs1_only_record[0:300])
trs1_end_tag_file_name = trs1_record.name + "_trs1_end_tag.fasta"
with open(trs1_end_tag_file_name, "w") as trs1_end_tag_outfile:
    trs1_end_tag_outfile.write(">" + trs1_record.name + "_trs1_end_tag\n" + trs1_only_record[len(trs1_only_record)-300:])

pretrs1_hit_file = sys.argv[1] + ".pretrs1.txt"
pretrs1_hit = top_hit(pretrs1_hit_file)
if not pretrs1_hit:
    pretrs1_hit = "NC_006273.2"
pretrs1_record = get_genbank_record(pretrs1_hit)
pretrs1_regions = regions_dict[pretrs1_hit[:pretrs1_hit.index(".")]]  
pretrs1_start = pretrs1_regions["regions"]["PRETRS1"][0]
pretrs1_end = pretrs1_regions["regions"]["PRETRS1"][1]
pretrs1_only_record = str(pretrs1_record.seq[pretrs1_start:pretrs1_end])
with open(meta_id + "_pretrs1_hit.fasta", 'w') as pretrs1_out:
    pretrs1_out.write(">" + pretrs1_record.name + "_pretrs1\n" + pretrs1_only_record) 

def find_top_combined_hit(trl_hit_file, trs_hit_file, irl_irs_hit_file):
    try:
        df_1 = pd.read_csv(trl_hit_file,sep='\t',header=None)
        df_2 = pd.read_csv(trs_hit_file,sep='\t',header=None)
        df_3 = df_1.merge(df_2[[1, 11]], left_on=1,
                    right_on=1)
        if df_3.shape[0] > 0:
            best_hit = df_3.iloc[df_3[['11_x', '11_y']].sum(axis = 1).idxmax()][1]
        else:
            raise pd.errors.EmptyDataError
    except pd.errors.EmptyDataError:
        best_hit = top_hit(irl_irs_hit_file)
        if not best_hit:
            logging.debug("Failed to get any hit from TRL, TRS, IRL_IRS. Probably not worth proceeding with this sample.")
    return best_hit

class ExtendDirection(Enum):
    BOTH = 0
    LEFT = 1
    RIGHT = 2

def extend_contig_map(contig_record, extend_direction, fastq_1, fastq_2, ignore_back_branch):
    base_name = meta_id + "_" + contig_record.name
    fasta_name = base_name + ".fasta"
    base_name_extended = base_name + "_extended"
    extended_fasta = base_name_extended + ".fasta"
    #interleaved_fastq = meta_id +"_interleaved.fastq"
    contig_length = len(contig_record.seq)
    SeqIO.write(contig_record, fasta_name, "fasta")
    
    if extend_direction == ExtendDirection.LEFT:
        extend_flag = "extendleft=300"
    elif extend_direction == ExtendDirection.RIGHT:
        extend_flag = "extendright=300"
    elif extend_direction == ExtendDirection.BOTH:
        extend_flag = "extendleft=300 extendright=300"

    if ignore_back_branch:
        ibb = "ibb=t"
    else:
        ibb = "ibb=f"
 
    logging.debug("Extending with tadpole.sh.")
    subprocess.run(["tadpole.sh", "in=" + fasta_name, "extra=" + fastq_1 + "," + fastq_2, "out=" + extended_fasta, extend_flag, "overwrite=true" , ibb, "mode=extend"])
    extended_fasta_exists = os.path.isfile("./" + extended_fasta)
    if not extended_fasta_exists:
        logging.debug("Tadpole did not generate an extension fasta.  This may be a complete contig.")
        return contig_record
    new_record = next(SeqIO.parse(extended_fasta, "fasta"))
    new_seq_str = str(new_record.seq)        
    if extend_direction == ExtendDirection.LEFT:
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
    elif extend_direction == ExtendDirection.BOTH:
        return new_record

def extend_contig(contig_record, extend_direction, fastq_1, fastq_2, trim_contig_length, ignore_back_branch):
    contig_length = len(contig_record.seq)
    if contig_length >= trim_contig_length:
        truncate_seq = trim_contig_length
    else:
        truncate_seq = contig_length
    logging.debug("Attempting to extend contig " + contig_record.name + " length=" + str(contig_length) + " trim_contig_length=" + str(trim_contig_length))
    if extend_direction == ExtendDirection.LEFT:
        left_tag = SeqRecord(
            contig_record.seq[:truncate_seq],
            id=contig_record.name + "_left",
            name=contig_record.name + "_left",
            description="",
        )
        contig_extended = extend_contig_map(left_tag, ExtendDirection.LEFT, fastq_1, fastq_2, ignore_back_branch)        
    elif extend_direction == ExtendDirection.RIGHT:
        right_tag = SeqRecord(
            contig_record.seq[contig_length-truncate_seq:],
            id=contig_record.name + "_right",
            name=contig_record.name + "_right",
            description="",
        )
        contig_extended = extend_contig_map(right_tag, ExtendDirection.RIGHT, fastq_1, fastq_2, ignore_back_branch)
    elif extend_direction == ExtendDirection.BOTH:
        contig_extended = extend_contig_map(contig_record, ExtendDirection.BOTH, fastq_1, fastq_2, ignore_back_branch)
    logging.debug("Completed extension. New contig length=" + str(len(contig_extended.seq)))
    return contig_extended

def get_contig_by_tag_fasta(scaffolds_file, tag_fasta_file_name, return_name_and_coordinates_only):
    logging.debug("Looking for tag in contig.")
    for record in list(SeqIO.parse(scaffolds_file, "fasta")):
        base_name = meta_id + tag_fasta_file_name[:tag_fasta_file_name.index(".fasta")] + record.name
        db_fasta_name = base_name + ".fasta"
        SeqIO.write(record, db_fasta_name, 'fasta')
        subprocess.run(["makeblastdb", "-in", db_fasta_name, "-out", base_name + "_id", "-parse_seqids", "-dbtype", "nucl"])            
        ul_start_out = open(base_name + "_tag.txt", "w")
        subprocess.run(["blastn", "-db", base_name + "_id", "-query", tag_fasta_file_name, "-outfmt", "6"], stdout=ul_start_out)
        try:
            df_tag_start = pd.read_csv(base_name + "_tag.txt",sep='\t',header=None)
            tag_start = df_tag_start.iloc[df_tag_start[11].idxmax()][8]
            # if contig is in reverse orientation of merlin_ul_start
            if df_tag_start.iloc[df_tag_start[11].idxmax()][8] > df_tag_start.iloc[df_tag_start[11].idxmax()][9]:
                logging.debug("Found reverse complemented tag_start.")
                if return_name_and_coordinates_only:
                    return {"name": record.name, "reverse": True, "start": df_tag_start.iloc[df_tag_start[11].idxmax()][8], "end": df_tag_start.iloc[df_tag_start[11].idxmax()][9]}
                else:                
                    return SeqRecord(
                        record.seq.reverse_complement(),
                        id=record.name,
                        name=record.name,
                        description="",
                    )
            else:
                logging.debug("Found tag_start.")
                if return_name_and_coordinates_only:
                    return {"name": record.name, "reverse": False, "start": df_tag_start.iloc[df_tag_start[11].idxmax()][8], "end": df_tag_start.iloc[df_tag_start[11].idxmax()][9]}
                else:  
                    return record
        except pd.errors.EmptyDataError:
            continue
    logging.debug("Did not find tag_start.")       
    return None

#python3 /genome_identification/cmv/extract_build_reference.py CMV_S1_0M_best_ref.scaffolds.fa CMV_S1_0M_best_ref CMV_S1_0M_best_ref.preprocessed_1.fastq CMV_S1_0M_best_ref.preprocessed_2.fastq
def run_minimap2(target_sequence_file, query_sequence_file):

    paf_filename = meta_id + ".paf"

    # Run minimap2 to generate PAF output
    subprocess.run(["minimap2", "-x", "asm5", target_sequence_file, query_sequence_file, "--secondary=no", "-t", "1", "-o", paf_filename], check=True)

    if os.path.exists(paf_filename) and os.path.getsize(paf_filename) > 0:
        return paf_filename
    else:
        return None

# This finds exact match on ends
def find_max_overlap(consensus, contig):
    max_overlap = 0

    for i in range(1, min(len(consensus), len(contig)) + 1):
        if consensus[-i:] == contig[:i]:
            max_overlap = i

    return max_overlap

def generate_consensus_sequence(bam_file_path, reference_fasta_path, output_fasta_path):
    # Open BAM file and reference FASTA file
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    reference = pysam.FastaFile(reference_fasta_path)

    # Create a dictionary to store the consensus base for each position
    consensus_dict = {}

    # Iterate over each position in the reference genome
    for i in range(reference.get_reference_length(reference.references[0])):
        # Collect the bases at the current position from all reads in the BAM file
        bases_at_position = [
            read.query_sequence[pos] 
            for read in bam_file.fetch(reference.references[0], i, i + 1)
            for pos, ref_pos in read.get_aligned_pairs(matches_only=True)
            if ref_pos == i
        ]
        # If there are bases at the current position, calculate the most frequent base as the consensus base
        if bases_at_position:
            consensus_base = max(set(bases_at_position), key=bases_at_position.count)            
        else:
            # If no bases are present, use the reference base as the consensus base
            consensus_base = reference.fetch(reference.references[0], i, i + 1)

        # Store the consensus base in the dictionary
        consensus_dict[i] = consensus_base

    # Close the BAM file and reference FASTA file
    bam_file.close()
    reference.close()

    # Write the consensus sequence to a FASTA file
    with open(output_fasta_path, "w") as output_fasta:
        output_fasta.write(f">{reference.references[0]}\n")
        output_fasta.write("".join(consensus_dict.values()))

    print(f"Consensus sequence generated and saved to {output_fasta_path}")

def overlap_consensus(fasta_file, start_sequence_name):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    logging.debug("Attempting to create overlap consensus.")
    if start_sequence_name not in sequences:
        logging.debug(f"Sequence with name {start_sequence_name} not found in the input FASTA file.")
        return ""

    consensus_sequence = str(sequences[start_sequence_name].seq)
    logging.debug(f"Start sequence: {sequences[start_sequence_name].name}")
    remaining_sequences = list(sequences.keys())
    remaining_sequences.remove(start_sequence_name)
    logging.debug(f"Remaining sequences: {str(remaining_sequences)}")

    last_remaining_count = len(remaining_sequences)
    current_remaining_count = last_remaining_count - 1
    while len(remaining_sequences) > 0 and current_remaining_count < last_remaining_count:
        last_remaining_count = len(remaining_sequences)
        for seq_name in remaining_sequences:
            current_sequence = str(sequences[seq_name].seq)
            overlap = find_max_overlap(consensus_sequence, current_sequence)
            if overlap >= 10:
                consensus_sequence = consensus_sequence[:len(consensus_sequence)-overlap] + current_sequence
                remaining_sequences.remove(seq_name)
            else:
                current_sequence_reverse = str(sequences[seq_name].seq.reverse_complement())
                overlap = find_max_overlap(consensus_sequence, current_sequence_reverse)
                if overlap >= 10:
                    consensus_sequence = consensus_sequence[:len(consensus_sequence)-overlap] + current_sequence_reverse
                    remaining_sequences.remove(seq_name)

        current_remaining_count = len(remaining_sequences)

    logging.debug(f"Overlap complete. Consensus sequence length: {str(len(consensus_sequence))} Remaining sequences: {str(remaining_sequences)}")
    overlap_records = []
    overlap_records.append(SeqRecord(
            Seq(consensus_sequence),
            id="1_Overlap_consensus",
            name="1_Overlap_consensus",
            description="",            
        ))
    
    for name in remaining_sequences:
        print(str(remaining_sequences))
        overlap_records.append(SeqRecord(
            sequences[name].seq,
            id=sequences[name].id,
            name=sequences[name].name,
            description="",            
        ))        

    return overlap_records

def parse_paf(paf_filename):
    with open(paf_filename, 'r') as f:
        lines = f.readlines()
    if lines:
        fields = lines[0].strip().split('\t')
        #qstart, qend, strand, target, tlength, tstart, tend = map(int, fields[2:9])
        print(str(fields))
        return {"strand": fields[4], "coords": (fields[2], fields[3], fields[7], fields[8])}
    return None

def map_and_de_novo_assemble_region(region, region_ref_fasta, perform_de_novo):
    logging.debug("Bowtie2 index and mapping for best " + region + " hit.")
    subprocess.run(["mkdir", meta_id + "_" + region + "_bowtie2/"])
    bowtie2_ref = meta_id + "_" + region + "_bowtie2/" + meta_id + "_" + region
    subprocess.run(["bowtie2-build", region_ref_fasta, bowtie2_ref])
    subprocess.run(["bowtie2", "-x",  bowtie2_ref, "-1", fastq_1, "-2", fastq_2, "-S", meta_id + "_" + region + ".sam", "-p", cpus])
    bam_out = open(meta_id + "_" + region + ".bam", "w")
    subprocess.run(["samtools", "view", "-h", "-Sb", meta_id + "_" + region + ".sam"], stdout=bam_out)
    F12_bam_out = open(meta_id + "_F12_" + region + ".bam", "w")
    subprocess.run(["samtools", "view", "-F", "12", "-b", meta_id + "_" + region + ".bam"], stdout=F12_bam_out)
    F8_bam_out = open(meta_id + "_F8_" + region + ".bam", "w")
    subprocess.run(["samtools", "view", "-f", "4", "-F", "8", "-b", meta_id + "_" + region + ".bam"], stdout=F8_bam_out)
    F4_bam_out = open(meta_id + "_F4_" + region + ".bam", "w")
    subprocess.run(["samtools", "view", "-f", "8", "-F", "4", "-b", meta_id + "_" + region + ".bam"], stdout=F4_bam_out)
    subprocess.run(["samtools", "merge", "-f", meta_id + "_full_" + region + "_reads.bam", meta_id + "_F12_" + region + ".bam", meta_id + "_F8_" + region + ".bam", meta_id + "_F4_" + region + ".bam"])
    subprocess.run(["bam2fastq", "-f", "-o", meta_id + "_" + region + "_filtered_reads#.fastq", meta_id + "_full_" + region + "_reads.bam"])
    if not perform_de_novo:
        bam_file_path = meta_id + "_full_" + region + "_reads.bam"
        sorted_bam_file_path = meta_id + "_full_" + region + "_reads_sorted.bam"
        reference_fasta_path = region_ref_fasta
        output_fasta_path = meta_id + "_" + region + "_scaffold.fasta"
        subprocess.run(["samtools", "sort", "-o", sorted_bam_file_path, bam_file_path])
        subprocess.run(["samtools", "index", sorted_bam_file_path])
        generate_consensus_sequence(sorted_bam_file_path, reference_fasta_path, output_fasta_path)
        return meta_id + "_" + region + "_scaffold.fasta"
    else:
        logging.debug("De novo assembly of extracted reads mapped to best " + region + " hit.")
        subprocess.run(["unicycler", "--threads", cpus, "--mode", "bold", "-1", meta_id + "_" + region + "_filtered_reads_1.fastq", "-2", meta_id + "_" + region + "_filtered_reads_2.fastq", "--out", "./"])
        subprocess.run(["cp", "assembly.fasta", meta_id + "_" + region + "_scaffold.fasta"])
        return meta_id + "_" + region + "_scaffold.fasta"

def trim_for_trs1_sequence(us_contig, trs_contig):
    logging.debug("Attempting to trim at US_end and TRS_start and insert best TRS1 hit.")
    
    trs1_scaffold_fasta = map_and_de_novo_assemble_region("trs1", meta_id + "_trs1_hit.fasta", False)
    records = SeqIO.parse(trs1_scaffold_fasta, "fasta")

    num_records = sum(1 for record in records)
    if num_records > 1:
        logging.debug("De novo assembly of extracted trs1 reads generated more than 1 contig. Scaffold_builder with best trs1 hit reference.")
        subprocess.run(["python2.7", "/scaffold_builder.py", "-q", trs1_scaffold_fasta, "-r", meta_id + "_trs1_hit.fasta", "-p", meta_id + "_trs1_final"])
        subprocess.run(["mv", meta_id + "_trs1_final_Scaffold.fasta", meta_id + "_best_trs1.fasta"])
    else:
        logging.debug("De novo assembly of extracted trs1 reads generated 1 contig.")
        subprocess.run(["mv", trs1_scaffold_fasta, meta_id + "_best_trs1.fasta"])
    
    new_trs1_record = next(SeqIO.parse(meta_id + "_best_trs1.fasta","fasta"))
    extended_trs1_record = extend_contig(new_trs1_record, ExtendDirection.BOTH, fastq_1, fastq_2, 0, True)
    extended_trs1_record.id = meta_id + "_extended_trs1"
    extended_trs1_record.name = meta_id + "_extended_trs1"
    with open(meta_id + "_extended_trs1.fasta", "w") as trs1_extended_out:
        SeqIO.write(extended_trs1_record, trs1_extended_out, 'fasta')    
    
    subprocess.run(["makeblastdb", "-in", meta_id + "_extended_trs1.fasta", "-out", meta_id + "_trs1_hit_id", "-parse_seqids", "-dbtype", "nucl"])

    # Only search for US_end at len(us_contig) - 10000 so that IRS1 hit doesn't confuse matters.
    # Just remove that possibility temporarily. 
    with open(meta_id + "_US.fasta", "w") as us_out:
        us_out.write(">" + meta_id + "_US\n" + us_contig[len(us_contig)-10000:])
    with open(meta_id + "_TRS.fasta", "w") as trs_out:
        trs_out.write(">" + meta_id + "_TRS\n" + trs_contig)
    us_results_out = open(meta_id + "_us.txt", "w")
    subprocess.run(["blastn", "-db", meta_id + "_trs1_hit_id", "-query", meta_id + "_US.fasta", "-outfmt", "6"], stdout=us_results_out)
    trs_results_out = open(meta_id + "_trs.txt", "w")
    subprocess.run(["blastn", "-db", meta_id + "_trs1_hit_id", "-query", meta_id + "_TRS.fasta", "-outfmt", "6"], stdout=trs_results_out) 
    try:
        logging.debug("Finding US_end...")
        df = pd.read_csv(meta_id + "_us.txt",sep='\t',header=None)
        if df.shape[0] > 0:
            adjusted_trs1_start = (len(us_contig)-10000) - 1 + df.iloc[df[11].idxmax()][6]
            us = us_contig[:adjusted_trs1_start]
            logging.debug("Found US_end: " + str(adjusted_trs1_start))
        logging.debug("Finding TRS_start...")
        df = pd.read_csv(meta_id + "_trs.txt",sep='\t',header=None)
        logging.debug("Found TRS_start.")
        if df.shape[0] > 0:
            logging.debug("Found TRS_start: " + str(df.iloc[df[11].idxmax()][7]))
            
            # start_trs1_info = get_contig_by_tag_fasta(meta_id + "_extended_trs1.fasta", trs1_start_tag_file_name, True)
            # end_trs1_info = get_contig_by_tag_fasta(meta_id + "_extended_trs1.fasta", trs1_end_tag_file_name, True)

            # start_trs1 = start_trs1_info["start"]
            # end_trs1 = end_trs1_info["end"]

            # if start_trs1 > end_trs1:
            #     trs1_insert = str(extended_trs1_record.seq[end_trs1:start_trs1].reverse_complement())
            # else:
            #     trs1_insert = str(extended_trs1_record.seq[start_trs1:end_trs1])

            full_contig_with_trs1_inserted = us + str(extended_trs1_record.seq) + trs_contig[df.iloc[df[11].idxmax()][7]:]
        logging.debug("Full contig with TRS1 inserted length: " + str(len(full_contig_with_trs1_inserted)))
        return {"final_contig": full_contig_with_trs1_inserted, "trs1_start": adjusted_trs1_start}                 
    except pd.errors.EmptyDataError:
        return None
        logging.debug("Failed to trim and insert best TRS1 hit.")
        
def extend_short_us_end(new_contigs, contig_name):
    logging.debug("Attempting to extend contig to best possible US.")
    updated_contigs = []
    for record in new_contigs:
        if record.name == contig_name:
            contig = record
        else:
            updated_contigs.append(record)
    base_name = meta_id + "_search_us_end"
    fasta_name = base_name + ".fasta"
    blast_out_name = base_name + ".txt"
    with open(meta_id + "_search_us_end.fasta", "w") as search_out:
        search_out.write(">" + meta_id + "_search_us_end.fasta\n" + str(contig.seq))
    search_out = open(blast_out_name, "w")    
    # subprocess.run(["blastn", "-db", "/genome_identification/cmv/cmv_ussub_id", "-query", fasta_name, "-outfmt", "6"], stdout=search_out)
    with open(meta_id + "_hyper_ussub.fasta", "w") as search_subject_out:
        search_subject_out.write(">" + meta_id + "_hyper_ussub.fasta\n" + str(hyper_hit_record.seq[hyper_hit_regions["regions"]["USSUB"][0]:hyper_hit_regions["regions"]["USSUB"][1]]))    
    subprocess.run(["makeblastdb", "-in", meta_id + "_hyper_ussub.fasta", "-out", meta_id + "_hyper_us_hit_id", "-parse_seqids", "-dbtype", "nucl"])
    subprocess.run(["blastn", "-db", meta_id + "_hyper_us_hit_id", "-query", fasta_name, "-outfmt", "6"], stdout=search_out)
    try:
        df = pd.read_csv(blast_out_name ,sep='\t', header=None)
        logging.debug("Found best US hit.")
        # subject = get_genbank_record(df.iloc[df[11].idxmax()][1])
        # subject_regions = regions_dict[subject.name]
        subject = hyper_hit_record
        subject_regions = hyper_hit_regions

        subject_seq = subject.seq[subject_regions["regions"]["USSUB"][0] + df.iloc[df[11].idxmax()][9]:subject_regions["regions"]["TRS1"][0]+2700]
        logging.debug("Got extension. length:" + str(len(subject_seq)) + " start: " + str(subject_regions["regions"]["USSUB"][0] + df.iloc[df[11].idxmax()][9]) + " end: " +str(subject_regions["regions"]["TRS1"][0]+2000) )
        if df.iloc[df[11].idxmax()][6] > df.iloc[df[11].idxmax()][7]:
            oriented_contig = contig.seq[df.iloc[df[11].idxmax()][7]:].reverse_complement()
        else:
            oriented_contig = contig.seq[:df.iloc[df[11].idxmax()][7]]            
        updated_contigs.append(SeqRecord(
            Seq(str(oriented_contig) + str(subject_seq)),
            id=contig.name + "_us_extended",
            name=contig.name + "_us_extended",
            description="",            
        ))
        return updated_contigs
    except pd.errors.EmptyDataError:
        logging.debug("Failed to extend contig to best possible US.")
        None

def remove_contigs_that_map_to_sequence(contig_records, sequence):
    logging.debug("Mapping contigs to reference sequence for removal.")
    with open(meta_id + "_remove_reads_reference.fasta", "w") as remove_ref_out:
        remove_ref_out.write(">" + meta_id + "_remove_reads_reference\n" + sequence + "\n") 
    subprocess.run(["makeblastdb", "-in", meta_id + "_remove_reads_reference.fasta", "-out", meta_id + "_remove_reads_id", "-parse_seqids", "-dbtype", "nucl"])
    with open("remove_contigs.fasta", "w") as remove_outfile:
        SeqIO.write(contig_records, remove_outfile, 'fasta')
    search_remove_out = open("remove_contigs.txt", "w")
    subprocess.run(["blastn", "-db", meta_id + "_remove_reads_id", "-query", "remove_contigs.fasta", "-outfmt", "6"], stdout=search_remove_out)
    try:
        df = pd.read_csv("remove_contigs.txt" ,sep='\t', header=None)
        contig_names_to_remove = df[0].unique().tolist()
        contig_names = [str(x) for x in contig_names_to_remove]
        logging.debug("Contig names to be removed: " + str(contig_names_to_remove))
        new_contigs = []
        for record in contig_records:
            logging.debug("Checking to remove record.id="+record.id)
            if record.id in contig_names:
                logging.debug("Removing record.id="+record.id)
            else:
                new_contigs.append(record)
        logging.debug("Completed removing mapped contigs.  Contigs: " + str(new_contigs))
        return new_contigs
    except pd.errors.EmptyDataError:
        logging.debug("Failed to find any contigs to remove.")
        return contig_records    

def remove_contigs_replace_region(replace_region_start, replace_region_end, contigs, reference_complete_file_name, reference_record, merlin_start_tag_file, merlin_end_tag_file):
    logging.debug("Attempting to insert region into either replace_region_start or replace_region_end.")
    reference_tag_start = get_contig_by_tag_fasta(reference_complete_file_name, merlin_start_tag_file, True)
    reference_tag_end = get_contig_by_tag_fasta(reference_complete_file_name, merlin_end_tag_file, True)
    no_hit_contigs = []
    if replace_region_start:
        logging.debug("Found replace_region_start. Start: " + str(replace_region_start)) 
        # Get UL_start contig and extend with hyper variable hit UL
        for unfinished_record in contigs:
            if unfinished_record.name == replace_region_start["name"]:
                if replace_region_start["reverse"]:
                    trimmed_start_seq = unfinished_record.seq[replace_region_start["start"]:].reverse_complement()
                else:
                    trimmed_start_seq = unfinished_record.seq[:replace_region_start["start"]]
                logging.debug("Creating new region start contig. Trimmed contig length: " + str(len(trimmed_start_seq)) + " Region length: " + 
                              str(len(reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]])))
                new_ul_contig_str = str(trimmed_start_seq) + str(reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]])
                new_contig = SeqRecord(
                                            Seq(new_ul_contig_str),
                                            id=replace_region_start["name"] + "_insert",
                                            name=replace_region_start["name"] + "_insert",
                                            description="",
                                        )
            else:
                no_hit_contigs.append(unfinished_record)
    elif replace_region_end:                
        for unfinished_record in contigs:
            if unfinished_record.name == replace_region_end["name"]:
                if replace_region_end["reverse"]:
                    trimmed_start_seq = unfinished_record.seq[:replace_region_end["end"]].reverse_complement()
                else:
                    trimmed_start_seq = unfinished_record.seq[replace_region_end["end"]:]
                new_ul_contig_str =  str(reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]]) + str(trimmed_start_seq)
                logging.debug("Creating new region end contig. Complete Length: " + str(len(new_ul_contig_str)))
                new_contig = SeqRecord(
                                            Seq(new_ul_contig_str),
                                            id=replace_region_end["name"] + "_insert",
                                            name=replace_region_end["name"] + "_insert",
                                            description="",
                                        )    
            else:
                no_hit_contigs.append(unfinished_record)
    else:
        logging.debug("Neither start or end found in contigs. Adding complete region as new contig.")
        new_contig = SeqRecord(
                                    reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]],
                                    id=reference_record.id + "_insert",
                                    name=reference_record.name + "_insert",
                                    description="",
                                )
        no_hit_contigs = contigs       
    updated_contigs = remove_contigs_that_map_to_sequence(no_hit_contigs, str(reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]]))
    updated_contigs.append(new_contig)
    return updated_contigs
 
def replace_region(replace_region_start, replace_region_end, contigs, reference_complete_file_name, reference_record, merlin_start_tag_file, merlin_end_tag_file):
    logging.debug(f"Attempting to insert region between {replace_region_start['name']} and {replace_region_end['name']}.")
    reference_tag_start = get_contig_by_tag_fasta(reference_complete_file_name, merlin_start_tag_file, True)
    reference_tag_end = get_contig_by_tag_fasta(reference_complete_file_name, merlin_end_tag_file, True)
    ul_contigs = []
    if replace_region_start and replace_region_end:
        logging.debug("Found both ends of UL. Start: " + str(replace_region_start) + " End: " + str(replace_region_end)) 
        # Get UL contig and extend with hyper variable hit UL
        for unfinished_record in contigs:
            if unfinished_record.name == replace_region_start["name"]:
                if replace_region_start["reverse"]:
                    trimmed_start_seq = unfinished_record.seq[replace_region_start["start"]:].reverse_complement()
                else:
                    trimmed_start_seq = unfinished_record.seq[:replace_region_start["start"]]
                logging.debug("Creating new region. Trimmed contig length: " + str(len(trimmed_start_seq)) + " Region length: " + 
                              str(len(reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]])))
                new_ul_contig_str = str(trimmed_start_seq) + str(reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]])
                #print(new_ul_contig_str)
            else:
                ul_contigs.append(unfinished_record)
        new_ul_contigs = []
        for remaining_contig in ul_contigs:
            if remaining_contig.name == replace_region_end["name"]:
                if replace_region_end["reverse"]:
                    trimmed_start_seq = remaining_contig.seq[:replace_region_end["end"]].reverse_complement()
                else:
                    trimmed_start_seq = remaining_contig.seq[replace_region_end["end"]:]
                new_ul_contig_str = new_ul_contig_str + str(trimmed_start_seq)
                logging.debug("Creating new contig. Complete Length: " + str(len(new_ul_contig_str)))
                new_contig = SeqRecord(
                                Seq(new_ul_contig_str),
                                id=replace_region_start["name"] + "_insert",
                                name=replace_region_start["name"] + "_insert",
                                description="",
                            )   
            else:
                new_ul_contigs.append(remaining_contig)
        updated_contigs = remove_contigs_that_map_to_sequence(new_ul_contigs, str(reference_record.seq[reference_tag_start["start"]:reference_tag_end["end"]]))
        updated_contigs.append(new_contig)                
        return updated_contigs
    else:
        logging.debug("Failed to get both ends of region. Start: " + (str(replace_region_start) if replace_region_start else "Failed") + " End: " + (str(replace_region_end) if replace_region_end else "Failed"))

def join_contigs_by_irl_irs(fasta_file, start_sequence_name, overlap_end_ul_contig, overlap_start_us_contig):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    if start_sequence_name not in sequences:
        print(f"Sequence with name {start_sequence_name} not found in the input FASTA file.")
        return ""
    
    start_name = meta_id + "_start_contig_for_irl_irs_join"
    start_name_fasta = start_name + ".fasta"
    with open(start_name_fasta, "w") as start_outfile:
        SeqIO.write(sequences[start_sequence_name], start_outfile, 'fasta')

    top_irl_irs_hit = find_top_combined_hit(trl_hit_file, trs_hit_file, irl_irs_hit_file)
    top_irl_irs_record = get_genbank_record(top_irl_irs_hit)
    top_irl_irs_regions = regions_dict[top_irl_irs_hit[:top_irl_irs_hit.index(".")]]
    irl_irs_start = top_irl_irs_regions["regions"]["IRL_IRS"][0]
    irl_irs_end = top_irl_irs_regions["regions"]["IRL_IRS"][1]
    irl_irs = str(top_irl_irs_record.seq[irl_irs_start:irl_irs_end])

    top_irl_irs_fasta = meta_id + "_" + top_irl_irs_hit + "_top_irl_irs.fasta"
    with open(top_irl_irs_fasta, "w") as top_irl_irs_outfile:
        top_irl_irs_outfile.write(">" + top_irl_irs_hit + "_irl_irs\n" + irl_irs + "\n")

    irl_irs_scaffold_fasta = map_and_de_novo_assemble_region("irl_irs", top_irl_irs_fasta, True)
    records = SeqIO.parse(irl_irs_scaffold_fasta, "fasta")

    num_records = sum(1 for record in records)
    if num_records > 1:
        logging.debug("De novo assembly of extracted irl_irs reads generated more than 1 contig. Scaffold_builder with best irl_irs hit reference.")
        subprocess.run(["python2.7", "/scaffold_builder.py", "-q", irl_irs_scaffold_fasta, "-r", top_irl_irs_fasta, "-p", meta_id + "_isl_irs_final"])
        subprocess.run(["mv", meta_id + "_isl_irs_final_Scaffold.fasta", meta_id + "_best_irl_irs.fasta"])
    else:
        logging.debug("De novo assembly of extracted irl_irs reads generated 1 contig.")
        subprocess.run(["mv", irl_irs_scaffold_fasta, meta_id + "_best_irl_irs.fasta"])
    new_irl_irs_record = next(SeqIO.parse(meta_id + "_best_irl_irs.fasta","fasta"))
    extended_irl_irs_record = extend_contig(new_irl_irs_record, ExtendDirection.BOTH, fastq_1, fastq_2, 0, True)
    extended_irl_irs_record.id = meta_id + "_extended_irl_irs"
    extended_irl_irs_record.name = meta_id + "_extended_irl_irs"

    new_records = []
    new_records.append(sequences[start_sequence_name])
    start_seq_length = len(sequences[start_sequence_name].seq)
    new_records.append(extended_irl_irs_record)
    remaining_sequences = list(sequences.keys())
    remaining_sequences.remove(start_sequence_name)
    second_sequence_name = remaining_sequences[0]
    second_seq_length = len(sequences[second_sequence_name].seq)
    new_records.append(sequences[second_sequence_name])
    with open(meta_id + "_extended_plus_irl_irs.fasta", "w") as extended_irl_irs_outfile:
        SeqIO.write(new_records, extended_irl_irs_outfile, 'fasta')
    overlap_records = overlap_consensus(meta_id + "_extended_plus_irl_irs.fasta", start_sequence_name)

    # if the contigs don't overlap completely/perfectly
    if len(overlap_records) > 1:
        for record in overlap_records:
            # Extra contigs are probably host
            if len(record.seq) > 228000:
                logging.debug("Join contigs by irl_irs generated a multiple contigs with one contig > 228000.")
                return overlap_records
        # join contigs by using merlin coordinates
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        ul_end_record = sequences[overlap_end_ul_contig["name"]]
        us_start_record = sequences[overlap_start_us_contig["name"]]
        joined_contig = SeqRecord(
                                Seq(str(ul_end_record.seq[:overlap_end_ul_contig["end"]]) + irl_irs + str(us_start_record.seq[overlap_start_us_contig["start"]:])),
                                id="irl_irs_joined",
                                name="irl_irs_joined",
                                description="",
                            )
        logging.debug("Join contigs by irl_irs generated a single contig after appending at merlin coordinates length: " + str(len(joined_contig.seq)))
        return [joined_contig]
    
    logging.debug("Join contigs by irl_irs generated a single contig after overlap.")
    return overlap_records

def join_contigs_by_pretrs1(fasta_file, start_sequence_name, contig_1, contig_2):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    if start_sequence_name not in sequences:
        print(f"Sequence with name {start_sequence_name} not found in the input FASTA file.")
        return ""
    
    start_name = meta_id + "_start_contig_for_pretrs1_join"
    start_name_fasta = start_name + ".fasta"
    with open(start_name_fasta, "w") as start_outfile:
        SeqIO.write(sequences[start_sequence_name], start_outfile, 'fasta')

    top_pretrs1_fasta = meta_id + "_pretrs1_hit.fasta"
    pretrs1_scaffold_fasta = map_and_de_novo_assemble_region("pretrs1", top_pretrs1_fasta, True)
    records = SeqIO.parse(pretrs1_scaffold_fasta, "fasta")

    num_records = sum(1 for record in records)
    if num_records > 1:
        logging.debug("De novo assembly of extracted pretrs1 reads generated more than 1 contig. Scaffold_builder with best pretrs1 hit reference.")
        subprocess.run(["python2.7", "/scaffold_builder.py", "-q", pretrs1_scaffold_fasta, "-r", top_pretrs1_fasta, "-p", meta_id + "_pretrs1_final"])
        subprocess.run(["mv", meta_id + "_pretrs1_final_Scaffold.fasta", meta_id + "_best_pretrs1.fasta"])
    else:
        logging.debug("De novo assembly of extracted pretrs1 reads generated 1 contig.")
        subprocess.run(["mv", pretrs1_scaffold_fasta, meta_id + "_best_pretrs1.fasta"])
    new_pretrs1_record = next(SeqIO.parse(meta_id + "_best_pretrs1.fasta","fasta"))
    extended_pretrs1_record = extend_contig(new_pretrs1_record, ExtendDirection.BOTH, fastq_1, fastq_2, 0, True)
    extended_pretrs1_record.id = meta_id + "_extended_pretrs1"
    extended_pretrs1_record.name = meta_id + "_extended_pretrs1"

    new_records = []
    new_records.append(sequences[start_sequence_name])
    start_seq_length = len(sequences[start_sequence_name].seq)
    new_records.append(extended_pretrs1_record)
    remaining_sequences = list(sequences.keys())
    remaining_sequences.remove(start_sequence_name)
    second_sequence_name = remaining_sequences[0]
    second_seq_length = len(sequences[second_sequence_name].seq)
    new_records.append(sequences[second_sequence_name])
    with open(meta_id + "_extended_plus_pretrs1.fasta", "w") as extended_pretrs1_outfile:
        SeqIO.write(new_records, extended_pretrs1_outfile, 'fasta')
    overlap_records = overlap_consensus(meta_id + "_extended_plus_pretrs1.fasta", start_sequence_name)

    # if the contigs don't overlap completely/perfectly
    if len(overlap_records) > 1:
        for record in overlap_records:
            # Extra contigs are probably host
            if len(record.seq) > 232000:
                logging.debug("Join contigs by pretrs1 generated a multiple contigs with one contig > 232000: " + str(len(record.seq)))
                return overlap_records
        # join contigs by using merlin coordinates
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        start_record = sequences[contig_1["name"]]
        end_record = sequences[contig_2["name"]]
        joined_contig = SeqRecord(
                                Seq(str(start_record.seq) + pretrs1_only_record + str(end_record.seq)),
                                id="pretrs1_joined",
                                name="pretrs1_joined",
                                description="",
                            )
        logging.debug("Join contigs by pretrs1 generated a single contig after appending at merlin coordinates length: " + str(len(joined_contig.seq)))
        return [joined_contig]
    
    logging.debug("Join contigs by pretrs1 generated a single contig after overlap.")
    return overlap_records

def find_TATTTA(sequence):
    pattern = "TATTTA"
    index = sequence.find(pattern)

    if index != -1:
        logging.debug(f"TATTTA found at index {index}.")
        return index
    else:
        logging.debug("TATTTA not found.")
        return None


# Remove any old scaffold files
subprocess.run(["rm *scaffold.fasta"], shell=True)

# Begin scaffolding, joining if necessary, extracting IRL, IRS, and setting TRL and TRS
logging.debug("Initial scaffold file: more than one contig detected. Extend and join contigs.")
# Extend all contigs
new_contigs = []
logging.debug("Extending contigs.")
for record in list(SeqIO.parse(scaffolds_file, "fasta")):
    extended_contig = extend_contig(record, ExtendDirection.BOTH, fastq_1, fastq_2, 0, True)
    new_contigs.append(extended_contig)

with open(meta_id + "_new_scaffold.fasta", "w") as new_scaffold_outfile:
    SeqIO.write(new_contigs, new_scaffold_outfile, 'fasta')
start_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_start_file_name, True)
end_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_end_file_name, True)

# If missing start_UL_contig or end_UL_contig, replace the UL region and remove any mapping contigs.
# The contigs that map will be replaced during Mugsy/make_reference step.
if not start_UL_contig or not end_UL_contig:
    new_contigs = remove_contigs_replace_region(start_UL_contig, end_UL_contig, new_contigs, hyper_complete_file_name, hyper_hit_record, merlin_ul_start_file_name, merlin_ul_end_file_name)
    with open(meta_id + "_new_scaffold.fasta", "w") as new_scaffold_outfile:
        SeqIO.write(new_contigs, new_scaffold_outfile, 'fasta')   
    start_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_start_file_name, True)
    end_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_end_file_name, True)     
# If start_UL_contig is not the same as end_UL_contig, replace the UL region with the 
# best hypervariable reference hit.
elif start_UL_contig["name"] != end_UL_contig["name"]:   
    new_contigs = replace_region(start_UL_contig, end_UL_contig, new_contigs, hyper_complete_file_name, hyper_hit_record, merlin_ul_start_file_name, merlin_ul_end_file_name)
    with open(meta_id + "_new_scaffold.fasta", "w") as new_scaffold_outfile:
        SeqIO.write(new_contigs, new_scaffold_outfile, 'fasta')   
    start_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_start_file_name, True)
    end_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_end_file_name, True)     

end_US_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", trs1_start_tag_file_name, True)
replace_end_US_contig = None
# If unable to find end_US contig then extend
if end_US_contig == None:
    start_us_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_us_start_file_name, True)
    new_contigs = extend_short_us_end(new_contigs, start_us_contig["name"])
    with open(meta_id + "_new_scaffold.fasta", "w") as new_scaffold_outfile:
        SeqIO.write(new_contigs, new_scaffold_outfile, 'fasta')
    end_US_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_us_end_file_name, True)
    start_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_start_file_name, True)
    end_UL_contig = get_contig_by_tag_fasta(meta_id + "_new_scaffold.fasta", merlin_ul_end_file_name, True)    


extension_contigs = []
for contig in new_contigs:
    if contig.name == start_UL_contig["name"] and contig.name == end_US_contig["name"]:
        logging.debug("Found single contig for trimming.")
        if end_US_contig["reverse"]:
            logging.debug("Contig is reverse orientation.")
            trimmed_end_seq = contig.seq[end_US_contig["end"]-300:].reverse_complement()
            end_seq = contig.seq[:end_US_contig["end"]-300].reverse_complement()
            logging.debug("trimmed_end_seq length: " + str(len(trimmed_end_seq)) + " end_seq length: " + str(len(end_seq)))
            trimmed_start_seq = trimmed_end_seq[len(contig.seq) - start_UL_contig["end"] - 300:]
            beginning_seq = trimmed_end_seq[:len(contig.seq) - start_UL_contig["end"] - 300]    
        else:
            logging.debug("Contig is forward orientation.")
            trimmed_end_seq = contig.seq[:end_US_contig["end"]+300]
            end_seq = contig.seq[end_US_contig["end"]+300:]
            logging.debug("trimmed_end_seq length: " + str(len(trimmed_end_seq)) + " end_seq length: " + str(len(end_seq)))
            trimmed_start_seq = trimmed_end_seq[start_UL_contig["end"]+300:]
            beginning_seq = trimmed_end_seq[:start_UL_contig["end"]+300]    
        logging.debug("Writing trimmed single contig length: " + str(len(trimmed_start_seq)))
        extension_contigs.append(SeqRecord(
                                    trimmed_start_seq,
                                    id=contig.name,
                                    name=contig.name,
                                    description="",
                                ))
    elif contig.name == start_UL_contig["name"]:
        if start_UL_contig["reverse"]:
            trimmed_start_seq = contig.seq[:start_UL_contig["end"]-300].reverse_complement()
            beginning_seq = contig.seq[start_UL_contig["end"]-300:].reverse_complement()
        else:
            trimmed_start_seq = contig.seq[start_UL_contig["end"]+300:]
            beginning_seq = contig.seq[:start_UL_contig["end"]+300]
        logging.debug("Writing trimmed UL start contig length: " + str(len(trimmed_start_seq)))
        extension_contigs.append(SeqRecord(
                                    trimmed_start_seq,
                                    id=contig.name,
                                    name=contig.name,
                                    description="",
                                ))
        print("Trimmed UL length:" + str(len(trimmed_start_seq)))
    elif contig.name == end_US_contig["name"]:
        if replace_end_US_contig:
            logging.debug("Replacing US_end contig with PRETRS1 extended contig.")
            logging.debug("Writing US_end contig with PRETRS1 extension length: " + str(len(replace_end_US_contig)))
            extension_contigs.append(replace_end_US_contig)
        else:
            if end_US_contig["reverse"]:
                trimmed_end_seq = contig.seq[end_US_contig["end"]-300:].reverse_complement()
                end_seq = contig.seq[:end_US_contig["end"]-300].reverse_complement()
            else:
                trimmed_end_seq = contig.seq[:end_US_contig["end"]+300]
                end_seq = contig.seq[end_US_contig["end"]+300:]
            logging.debug("Writing US_end contig with length: " + str(len(trimmed_end_seq)))
            extension_contigs.append(SeqRecord(
                                        trimmed_end_seq,
                                        id=contig.name,
                                        name=contig.name,
                                        description="",
                                    ))
            print("Trimmed US length:" + str(len(trimmed_end_seq)))
    else:
        extension_contigs.append(contig)

with open(meta_id + "_extension_scaffold.fasta", "w") as extension_scaffold_outfile:
    SeqIO.write(extension_contigs, extension_scaffold_outfile, 'fasta')
overlap_records = overlap_consensus(meta_id + "_extension_scaffold.fasta", start_UL_contig["name"])
with open(meta_id + "_overlap.fasta", "w") as overlap_outfile:
    SeqIO.write(overlap_records, overlap_outfile, 'fasta')

# Check to see if US region failed to complete after overlap
overlap_start_us_contig = get_contig_by_tag_fasta(meta_id + "_overlap.fasta", merlin_us_start_file_name, True)
overlap_end_us_contig = get_contig_by_tag_fasta(meta_id + "_overlap.fasta", merlin_us_end_file_name, True)

# Missing either the start or end of US region in contigs
if not overlap_start_us_contig or not overlap_end_us_contig:
    logging.debug("Either US_start or US_end is missing in contigs.  Attempting to insert best US region.")
    # insert best US region
    overlap_records = remove_contigs_replace_region(overlap_start_us_contig, overlap_end_us_contig, overlap_records, hyper_complete_file_name, hyper_hit_record, merlin_us_start_file_name, merlin_us_end_file_name)
    with open(meta_id + "_overlap.fasta", "w") as overlap_outfile:
        SeqIO.write(overlap_records, overlap_outfile, 'fasta')

# Both the start and end of US region are in the contigs
elif overlap_start_us_contig["name"] != overlap_end_us_contig["name"]:
    logging.debug("US_start and US_end are on different contigs.  Attempting to insert best US region.")
    # insert best US region and join with surrounding contigs
    overlap_records = replace_region(overlap_start_us_contig, overlap_end_us_contig, overlap_records, hyper_complete_file_name, hyper_hit_record, merlin_us_start_file_name, merlin_us_end_file_name)
    with open(meta_id + "_overlap.fasta", "w") as overlap_outfile:
        SeqIO.write(overlap_records, overlap_outfile, 'fasta')

logging.debug("overlap_records size: " + str(len(overlap_records)))

# Check to see if UL and US are on different contigs and, if so, join with irl_irs.
overlap_end_ul_contig = get_contig_by_tag_fasta(meta_id + "_overlap.fasta", merlin_ul_end_file_name, True)
overlap_start_us_contig = get_contig_by_tag_fasta(meta_id + "_overlap.fasta", merlin_us_start_file_name, True)

if overlap_end_ul_contig["name"] != overlap_start_us_contig["name"]:
    logging.debug("Attempting to join contigs with best irl_irs hit region.")
    overlap_records = join_contigs_by_irl_irs(meta_id + "_overlap.fasta", overlap_end_ul_contig["name"], overlap_end_ul_contig, overlap_start_us_contig)
    logging.debug("irl_irs joined contigs size: " + str(len(overlap_records)))

if len(overlap_records) == 1:
    result = str(overlap_records[0].seq)
else:
    for record in overlap_records:
        # Check remaining 
        if len(record.seq) > 228000:
            result = str(record.seq)
            break
    if not result:
        print("Failed to join all contigs. overlap_records size: " + str(len(overlap_records)))
        sys.exit(0)
    

consensus_record = SeqRecord(
    Seq(result),
    id=meta_id + "_consensus_scaffold",
    name=meta_id + "_consensus_scaffold",
    description="",
)
with open(meta_id + "_initial_irl_irs_scaffold.fasta", "w") as initial_scaffold_out:
    SeqIO.write(consensus_record, initial_scaffold_out, 'fasta')

logging.debug("Attempting to split out PRETRS1, de novo assemble, and re-join the region.")
# Split contig at merlins PRETRS1 coordinates.
pretrs1_records = []
pretrs1_initial_end_contig = get_contig_by_tag_fasta(meta_id + "_initial_irl_irs_scaffold.fasta", merlin_pretrs1_end_file_name, True)
pretrs1_contig_2 = SeqRecord(
                            consensus_record.seq[pretrs1_initial_end_contig["start"]:],
                            id=consensus_record.id + "_2",
                            name=consensus_record.id + "_2",
                            description=""
                        )
pretrs1_initial_start_contig = get_contig_by_tag_fasta(meta_id + "_initial_irl_irs_scaffold.fasta", merlin_pretrs1_start_file_name, True)
pretrs1_contig_1 = SeqRecord(
                    consensus_record.seq[:pretrs1_initial_start_contig["end"]],
                    id=consensus_record.id + "_1",
                    name=consensus_record.id + "_1",
                    description=""
                )
pretrs1_records.append(pretrs1_contig_1)
pretrs1_records.append(pretrs1_contig_2)
with open(meta_id + "_pretrs1_contigs.fasta", "w") as pretrs1_outfile:
    SeqIO.write(pretrs1_records, pretrs1_outfile, 'fasta')

pretrs1_start_contig = {"name": consensus_record.id + "_1", "start": pretrs1_initial_start_contig["start"], "end": pretrs1_initial_start_contig["end"]-1}
pretrs1_end_tag_end = pretrs1_initial_start_contig["end"] - pretrs1_initial_start_contig["start"] -1
pretrs1_end_contig = {"name": consensus_record.id + "_2", "start": 0, "end": pretrs1_end_tag_end}

logging.debug("Attempting to join contigs with best PRETRS1 hit region.")
joined_records = join_contigs_by_pretrs1(meta_id + "_pretrs1_contigs.fasta", pretrs1_start_contig["name"], pretrs1_start_contig, pretrs1_end_contig)
logging.debug("USSUB joined contigs size: " + str(len(joined_records)))
if len(joined_records) == 1:
    consensus_record = SeqRecord(
        joined_records[0].seq,
        id=meta_id + "_consensus_pretrs1_scaffold",
        name=meta_id + "_consensus_pretrs1_scaffold",
        description="",
    )
    with open(meta_id + "_initial_scaffold.fasta", "w") as initial_joined_scaffold_out:
        SeqIO.write(consensus_record, initial_joined_scaffold_out, 'fasta')
else:
    print("Exiting joined_records size: " + str(len(joined_records)))

#     logging.debug("Mapping to initial_scaffold.fasta")
#     subprocess.run(["mkdir", meta_id + "_bowtie2"])
#     subprocess.run(["bowtie2-build", meta_id + "_initial_scaffold.fasta", meta_id + "_bowtie2/" + meta_id])
#     sp_sam = subprocess.run(["bowtie2", "-x", meta_id + "_bowtie2/" + meta_id, "-1", fastq_1, "-2", fastq_2, "--threads", cpus, "--no-unal", "--local", "--very-sensitive-local", "--seed", "1"], stdout=subprocess.PIPE)
#     sp_sorted_bam = subprocess.run(["samtools", "sort", "-@", cpus, "-o", meta_id + "_initial_mapped_sorted.bam"], input=sp_sam.stdout)
#     sp_mpileup = subprocess.run(["samtools", "mpileup", "-d", "5000", "-A", "-Q", "0", meta_id + "_initial_mapped_sorted.bam"], stdout=subprocess.PIPE)
#     subprocess.run(["ivar", "consensus", "-p", meta_id + "_initial_consensus", "-n", "'N'", "-m", "5", "-t", "0.75"], input=sp_mpileup.stdout)
#     logging.debug("Consensus created from mapping to initial_scaffold.fasta.")
    
#     if os.path.isfile(meta_id + "_initial_consensus.fa") == True:
#         initial_consensus_record = next(SeqIO.parse(meta_id + "_initial_consensus.fa","fasta"))
#         initial_consensus_str = str(initial_consensus_record.seq)
#         # Check for Ns in initial consensus.  If there, run GapFiller        
#         if "N" in initial_consensus_str:
#             logging.debug("Found Ns in initial_scaffold mapping. Using picard CollectInsertSizeMetrics to prepare for GapFiller.")
#             subprocess.run(["picard", "CollectInsertSizeMetrics", "I=" + meta_id + "_initial_mapped_sorted.bam", "O=" + meta_id 
#                             + "_insert_metrics.txt", "H=" + meta_id + "_insert_histogram.pdf", "M=0.5"])
#             with open(meta_id + "_insert_metrics.txt", 'r') as file:
#                 lines = file.readlines()
#                 if len(lines) >= 8:
#                     metrics = lines[7].strip().split("\t")
#                     insert_size = str(round(float(metrics[5]), 2))
#                     std_dev = str(round(float(metrics[6])/float(metrics[5]), 2))
#                     print("CollectInsertSizeMetrics insert_size: " + insert_size + " std_dev: " + std_dev)
#                     with open("gflib.txt","w") as gf_lib_out:
#                         gf_lib_out.write("lib1 bwa " + fastq_1 + " " + fastq_2 + " " + insert_size + " " + std_dev + " FR")
#                     logging.debug("GapFiller with gflib.txt: lib1 bwa " + fastq_1 + " " + fastq_2 + " " + insert_size + " " + std_dev + " FR")
#                     subprocess.run(["perl", "/genome_identification/cmv/GapFiller", "-q", "./", "-l", "gflib.txt", "-s", meta_id + "_initial_consensus.fa"])
#                     new_consensus_record = next(SeqIO.parse("./standard_output/standard_output.gapfilled.final.fa","fasta"))
#                     new_consensus_str = str(new_consensus_record.seq)
#                     logging.debug("Writing to initial_scaffold_consensus fasta.")
#                     with open(meta_id + "_initial_mapped_consensus.fasta", "w") as new_consensus_out:
#                         new_consensus_out.write(">" + meta_id + "_initial_mapped_consensus\n" + new_consensus_str + "\n")
#                     #subprocess.run(["cp", "./standard_output/standard_output.gapfilled.final.fa", meta_id + "_initial_mapped_consensus.fasta"])
#         else:
#             logging.debug("No Ns in initial_scaffold mapping. Writing to initial_scaffold_consensus.fasta.")
#             with open(meta_id + "_initial_mapped_consensus.fasta", "w") as new_consensus_out:
#                 new_consensus_out.write(">" + meta_id + "_initial_mapped_consensus\n" + initial_consensus_str + "\n")            
#         mapped_consensus_record = next(SeqIO.parse(meta_id + "_initial_mapped_consensus.fasta", "fasta"))

# mapped_consensus_str = str(mapped_consensus_record.seq)

logging.debug("Attempting to find UL_end and US_start in initial_scaffold_consensus.fasta.")
subprocess.run(["makeblastdb", "-in", meta_id + "_initial_scaffold.fasta", "-out", consensus_record.name + "_id", "-parse_seqids", "-dbtype", "nucl"])            
ul_end_out = open(consensus_record.name + "_ul_end.txt", "w")
subprocess.run(["blastn", "-db", consensus_record.name + "_id", "-query", merlin_ul_end_file_name, "-outfmt", "6"], stdout=ul_end_out)
us_start_out = open(consensus_record.name + "_us_start.txt", "w")
subprocess.run(["blastn", "-db", consensus_record.name + "_id", "-query", merlin_us_start_file_name, "-outfmt", "6"], stdout=us_start_out)       
logging.debug("Attempting to find 'a region' in consensus.")
a_out = open(consensus_record.name + "_a.txt", "w")
subprocess.run(["blastn", "-db", "/genome_identification/cmv/cmv_a_id", 
                "-query", meta_id + "_initial_scaffold.fasta", "-outfmt", "6"], stdout=a_out)
try:
    df_ul_end = pd.read_csv(consensus_record.name + "_ul_end.txt",sep='\t',header=None)
    ul_end = df_ul_end.iloc[df_ul_end[11].idxmax()][9]
    # if hyper_hit_trl_offset > 0:
    #     logging.debug("Extending merlin_ul_end by hyper_hit_trl_offset: " + str(hyper_hit_trl_offset))
    #     ul_end = ul_end + hyper_hit_trl_offset
    df_us_start = pd.read_csv(consensus_record.name + "_us_start.txt",sep='\t',header=None)
    us_start = df_us_start.iloc[df_us_start[11].idxmax()][8]       
    df_a = pd.read_csv(consensus_record.name + "_a.txt",sep='\t',header=None)
    scaffold_a_start = df_a.iloc[df_a[11].idxmax()][6]
    scaffold_a_stop = df_a.iloc[df_a[11].idxmax()][7]
    trs_sub_seq = consensus_record.seq[scaffold_a_start:us_start]
    trs = str(trs_sub_seq.reverse_complement())            
    trl_sub_seq = consensus_record.seq[ul_end+1:scaffold_a_stop+1]
    trl = str(trl_sub_seq.reverse_complement())
    
    print("TRL length: " + str(len(trl)) + " UL_end: " + str(ul_end) + " a_end: " + str(scaffold_a_stop))
    print("TRS length: " + str(len(trs)))

    with open(meta_id + "_beginning.fasta", "w") as beginning_outfile:
        beginning_outfile.write(">" + meta_id + "_beginning\n" + str(beginning_seq))
    with open(meta_id + "_end.fasta", "w") as end_outfile:
        end_outfile.write(">" + meta_id + "_end\n" + str(end_seq))
    with open(meta_id + "_trl.fasta", "w") as trl_outfile:
        trl_outfile.write(">" + meta_id + "_trl\n" + trl)
    with open(meta_id + "_trs.fasta", "w") as trs_outfile:
        trs_outfile.write(">" + meta_id + "_trs\n" + trs)                                 
    begin_paf_file = run_minimap2(meta_id + "_beginning.fasta", meta_id + "_trl.fasta")

    if not begin_paf_file:
        logging.debug("TRL paf file was not generated.  Attempting to extract TRL from hyper hit and append.")
        TATTTA_index = find_TATTTA(str(beginning_seq))
        if TATTTA_index:
            no_trl_complete = str(beginning_seq[TATTTA_index:]) + result
            no_trs = trl + no_trl_complete
        else:  
            no_trl_complete = str(beginning_seq[len(beginning_seq)-601:]) + result
            trl_search_out = open(meta_id + "_trl_best_hit.txt", "w")
            subprocess.run(["blastn", "-db", "/genome_identification/cmv/cmv_trl_id", "-query", meta_id + "_trl.fasta", "-outfmt", "6"], stdout=trl_search_out)
            try:
                df = pd.read_csv(meta_id + "_trl_best_hit.txt" ,sep='\t', header=None)
                subject = get_genbank_record(df.iloc[df[11].idxmax()][1])
                subject_regions = regions_dict[subject.name]
                with open(meta_id + "_best_trl_" + subject.name + ".fasta", "w") as best_trl_out:
                    SeqIO.write(subject, best_trl_out, 'fasta')
                #hyper_start_UL_contig = get_contig_by_tag_fasta(hyper_complete_file_name, merlin_ul_start_file_name, True)
                #hyper_trl = hyper_hit_record.seq[:hyper_start_UL_contig["start"]]
                best_trl_start_UL_contig = get_contig_by_tag_fasta(meta_id + "_best_trl_" + subject.name + ".fasta", merlin_ul_start_file_name, True)
                best_trl = subject.seq[:best_trl_start_UL_contig["start"]]                
                no_trs = str(best_trl) + no_trl_complete
            except pd.errors.EmptyDataError:
                    logging.debug("Could not extract trl from best hit trl db.")
                    sys.exit()                 
    else:
        begin_overlap_info = parse_paf(begin_paf_file)
        bqstart, bqend, btstart, btend = tuple(map(int, begin_overlap_info["coords"]))
        new_begin = trl[:bqstart] + str(beginning_seq[btstart:])
        no_trs = new_begin + result
    trs1_insert_info = trim_for_trs1_sequence(no_trs, trs)
    final_consensus = trs1_insert_info["final_contig"]
    trs1_start = trs1_insert_info["trs1_start"]
    logging.debug("Attempting to write complete final reference length: " + str(len(final_consensus)))
    final_consensus_record = SeqRecord(
        Seq(final_consensus),
        id=meta_id + "_final_consensus_scaffold",
        name=meta_id + "_final_consensus_scaffold",
        description="",
    )    
    with open(meta_id + "_final_scaffold.fasta", "w") as final_out:
        SeqIO.write(final_consensus_record, final_out, 'fasta')
        logging.debug("Final reference file generated.")
except pd.errors.EmptyDataError:
        logging.debug("Could not extract irl_irs from consensus.")
        sys.exit()

scaffolds_file_fasta = meta_id + "_final_scaffold.fasta"

subprocess.run(["makeblastdb", "-in", scaffolds_file_fasta, "-out", final_consensus_record.name + "_id", "-parse_seqids", "-dbtype", "nucl"])            
ul_start_out = open(final_consensus_record.name + "_ul_start.txt", "w")
subprocess.run(["blastn", "-db", final_consensus_record.name + "_id", "-query", merlin_ul_start_file_name, "-outfmt", "6"], stdout=ul_start_out)
ul_end_out = open(final_consensus_record.name + "_ul_end.txt", "w")
subprocess.run(["blastn", "-db", final_consensus_record.name + "_id", "-query", merlin_ul_end_file_name, "-outfmt", "6"], stdout=ul_end_out)
us_start_out = open(final_consensus_record.name + "_us_start.txt", "w")
subprocess.run(["blastn", "-db", final_consensus_record.name + "_id", "-query", merlin_us_start_file_name, "-outfmt", "6"], stdout=us_start_out)  
try:
    df_ul_start = pd.read_csv(final_consensus_record.name + "_ul_start.txt",sep='\t',header=None)
    ul_start = df_ul_start.iloc[df_ul_start[11].idxmax()][8]
    df_ul_end = pd.read_csv(final_consensus_record.name + "_ul_end.txt",sep='\t',header=None)
    ul_end = df_ul_end.iloc[df_ul_end[11].idxmax()][9]
    df_us_start = pd.read_csv(final_consensus_record.name + "_us_start.txt",sep='\t',header=None)

    trs1_start
    final_trl = str(final_consensus_record.seq[:ul_start])
    final_ul = str(final_consensus_record.seq[ul_start:ul_end])
    final_irl_irs = str(final_consensus_record.seq[ul_end:us_start])
    final_us = str(final_consensus_record.seq[us_start:trs1_start])
    final_trs = str(final_consensus_record.seq[trs1_start:])

    region_name = "region_TRL_" + meta_id
    write_region(region_name, final_trl)
    region_name = "region_UL_" + meta_id
    write_region(region_name, final_ul)
    region_name = "region_IRL_IRS_" + meta_id
    write_region(region_name, final_irl_irs)
    region_name = "region_US_" + meta_id
    write_region(region_name, final_us)
    region_name = "region_TRS_" + meta_id
    write_region(region_name, final_trs)

    # this will only be used for the region_order (it doesn't matter what coordinates are)
    with open("CMV-" + meta_id + ".json", 'w') as jsonfile:
        jsonfile.write(json.dumps(hyper_hit_regions))

    with open(meta_id + ".gb", "w") as gb_file:
        gb_file.write("This is intended to be empty.")

except pd.errors.EmptyDataError:
    logging.debug("Failed to find and extract region files.")  
sys.exit()

