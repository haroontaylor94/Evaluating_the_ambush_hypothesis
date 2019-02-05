import os
from OSC_analysis import *

output_dir = '../Results'
os.mkdir(output_dir)

for directory in ['../filtered_genomes_table_11', '../bacterial_genomes_table_4']:
    translation_table = int(directory.split('_')[3])
    output_file_path = output_dir + '/Table_' + str(translation_table) + '_results.csv'
    with open(output_file_path, 'w') as outfile: 
        header = 
'File_name,Plus_1_OSC_percentage,Plus_1_OSC_simulation_mean,' \ 'Plus_1_OSC_simulation_stdev,Plus_2_OSC_percentage,Plus_2_OSC_simulation_mean,' \ 'Plus_2_OSC_simulation_stdev,Total_OSC_percentage,Total_OSC_simulation_mean,' \ 'Total_OSC_simulation_stdev,GC_content,GC3_content,OSC_gap_mean,OSC_gap_stdev,' \ 'Sequence_length\n' 
        outfile.write(header)
        for file in os.listdir(directory):
            input_file_path = directory + file
            CDS = extract_coding_sequences(input_file_path) 
            plus1_OSCs = calculate_OSC_percentage(CDS, 1, translation_table)
            plus2_OSCs = calculate_OSC_percentage(CDS, 2, translation_table)
            total_OSCs = plus1_OSCs + plus2_OSCs
            plus1_mean, plus1_stdev, plus2_mean, plus2_stdev, total_mean, 
total_stdev = \
                calculate_simulation_means_and_standard_deviations(CDS, translation_table)
            GC_content = calculate_GC_content(CDS)
            GC3_content = calculate_GC3_content(CDS)
            gap_mean, gap_stdev = calculate_intervals_between_OSCs(CDS, translation_table)
            seq_length = len(CDS)
            results_list = [plus1_OSCs, plus1_mean, plus1_stdev, plus2_OSCs, plus2_mean, plus2_stdev, total_OSCs, total_mean, total_stdev, GC_content, GC3_content, gap_mean, gap_stdev, seq_length]
            line = ','.join([file[:-5]] + [str(x) for x in results_list] + '\n')
            outfile.write(line)