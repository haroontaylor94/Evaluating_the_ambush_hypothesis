from Bio import SeqIO
import numpy as np


def extract_coding_sequences(input_file, file_type): 
# Code adapted from http://wilke.openwetware.org/Parsing_Genbank_files_with_Biopython.html 
# Read in embl file
    genome = SeqIO.read(input_file, file_type)
    genome_CDS = ''
    for feature in genome.features: if feature.type == 'CDS':
        seq = feature.extract(genome.seq)
        # Verify it is an actual protein coding sequence
        if check_for_protein_sequence(seq): genome_CDS += seq
    return genome_CDS


def check_for_protein_sequence(input_seq): 
    if len(input_seq) % 3 != 0:
        # Check if sequence is a multiple of 3
        protein_seq = False
    elif input_seq[1:3] != 'TG':
        # Check if sequence starts with NTG
        protein_seq = False
    elif input_seq[-3:] not in ['TAG', 'TAA', 'TGA']:
        # Check if sequence ends with a stop codon
        protein_seq = False 
    else:
        protein_seq = True 
    return protein_seq


def calculate_OSC_percentage(input_seq, frame_shift, transl_table): 
    frame_shifts = 0
    if transl_table == 11:
        stop_codon_list = ['TAG', 'TAA', 'TGA']
    elif transl_table == 4: 
        stop_codon_list = ['TAG', 'TAA']
    else:
        stop_codon_list = 'NA'
    for i in range(frame_shift, len(input_seq), 3): codon = input_seq[i:i+3]
        if codon in stop_codon_list:
            frame_shifts += 1
    OSC_percentage = (frame_shifts * 100)/float(len(input_seq)) 
    return OSC_percentage


def calculate_simulation_means_and_standard_deviations(input_seq, transl_table):
    plus_one_OSCs_list = []
    plus_two_OSCs_list = []
    total_OSCs_list = []
    input_seq = str(input_seq)
    seq_list = [input_seq[i:i+3] for i in range(0, len(input_seq), 3)] 
    for i in xrange(100):
        random.shuffle(seq_list)
        random_CDS = ''.join(seq_list)
        plus_one_OSCs = calculate_OSC_percentage(random_CDS, 1, transl_table)
        plus_two_OSCs = calculate_OSC_percentage(random_CDS, 2, transl_table)
        total_OSCs = plus_one_OSCs + plus_two_OSCs
        plus_one_OSCs_list.append(plus_one_OSCs) 
        plus_two_OSCs_list.append(plus_two_OSCs) total_OSCs_list.append(total_OSCs)
    plus_one_mean = np.mean(plus_one_OSCs_list)
    plus_one_stdev = np.std(plus_one_OSCs_list)
    plus_two_mean = np.mean(plus_two_OSCs_list)
    plus_two_stdev = np.std(plus_two_OSCs_list)
    total_mean = np.mean(total_OSCs_list)
    total_stdev = np.std(total_OSCs_list)
    return plus_one_mean, plus_one_stdev, plus_two_mean, plus_two_stdev, total_mean, total_stdev


def calculate_GC_content(input_seq):
    total_G = input_seq.count('G')
    total_C = input_seq.count('C')
    GC_content = (total_G + total_C)/float(len(input_seq)) 
    return GC_content


def calculate_GC3_content(input_seq): 
    GC3_count = 0
    for i in range(2, len(input_seq), 3): base = input_seq[i]
        if base in ['G', 'C']:
        GC3_count += 1
    GC3_content = GC3_count/float(len(input_seq)/3) 
    return GC3_content


def calculate_intervals_between_OSCs(input_seq, transl_table): 
    if transl_table == 11:
        stop_codon_list = ['TAG', 'TAA', 'TGA'] 
    elif transl_table == 4:
        stop_codon_list = ['TAG', 'TAA'] 
    else:
        stop_codon_list = 'NA' 
    OSC_gap_list = []
    OSC_index = 0
    for i in range(1, len(input_seq)):
        codon = input_seq[i:i+3]
        if codon in stop_codon_list and i % 3 != 0:
            OSC_gap = i - OSC_index 
            OSC_gap_list.append(OSC_gap) 
            OSC_index = i
        # Reset indexing at end of each gene
        elif codon in stop_codon_list and i % 3 == 0: 
            OSC_index = i + 2
    gap_mean = np.mean(OSC_gap_list)
    gap_stdev = np.std(OSC_gap_list) 
    return gap_mean, gap_stdev
 
