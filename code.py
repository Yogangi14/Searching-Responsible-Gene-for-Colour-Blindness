#!/usr/bin/env python
# coding: utf-8

# In[39]:


from IPython.display import display
import ipywidgets as widgets
from collections import defaultdict
import warnings, os
import numpy as np
import pandas as pd
from sympy import binomial as dist
warnings.filterwarnings('ignore')
from tqdm import tqdm
import re


# In[40]:


def load_fasta(fasta_file_path):
    with open(fasta_file_path, 'r') as f:
        lines = f.readlines()
    # Ignore the header and concatenate the sequence lines
    reference = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return reference

# Load the reference sequence from your local file
fasta_file_path = "C:/Users/yogan/Downloads/chrX_bwt/chrX_bwt/chrX.fa"  
reference_sequence = load_fasta(fasta_file_path)


# In[41]:


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# In[42]:


# Define the gene start sequence
gene_start_sequence = "ATGGCCCAGCAGTGGAGCCTC" 

# Function to search for the gene start sequence in the reference
def find_gene_start_positions(reference, gene_start_sequence):
    positions = []
    search_len = len(gene_start_sequence)
    
    for i in range(len(reference) - search_len + 1):
        if reference[i:i + search_len] == gene_start_sequence:
            positions.append(i)
    
    return positions

# Find gene start positions in the reference sequence
gene_start_positions = find_gene_start_positions(reference_sequence, gene_start_sequence)

# Print the positions of the gene start sequence
print(f"Gene start positions: {gene_start_positions}")


# In[12]:


# Define the test read
test_read = "GAGGACAGCACCCAGTCCAGCATCTTCACCTACACCAACAGCAACTCCACCAGAGGTGAGCCAGCAGGCCCGTGGAGGCTGGGTGGCTGCACTGGGGGCCA"

# Function to check if a read aligns at a specific position in the reference
def check_read_alignment(reference, read, start_position):
    # Extract the subsequence from the reference starting at the given position
    ref_subseq = reference[start_position:start_position + len(read)]
    return ref_subseq == read

# Check if the test read aligns at position 149249814
alignment_result_test = check_read_alignment(reference_sequence, test_read, 149249814)

# Output the result of the alignment check
if alignment_result_test:
    print(f"The test read aligns at position 149249814.")
else:
    print(f"The test read does not align at position 149249814.")


# In[13]:


# Define the exon ranges for red and green genes
red_exons = [(149249757, 149249868), (149256127, 149256423), (149258412, 149258580), 
             (149260048, 149260213), (149261768, 149262007), (149264290, 149264400)]

green_exons = [(149288166, 149288277), (149293258, 149293554), (149295542, 149295710), 
               (149297178, 149297343), (149298898, 149299137), (149301420, 149301530)]

# Function to check if a position falls within the exons
def is_in_exon(position, exons):
    for start, end in exons:
        if start <= position < end:
            return True
    return False

# Function to count testreads mapping to red and green gene exons
def count_reads_in_exons(reads, reference, red_exons, green_exons):
    red_count = 0
    green_count = 0

    for read in reads:
        # Check if the read aligns to the reference at any of the red or green exon positions
        for start in range(len(reference) - len(read) + 1):
            if check_read_alignment(reference, read, start):
                if is_in_exon(start, red_exons):
                    red_count += 1
                elif is_in_exon(start, green_exons):
                    green_count += 1

    return red_count, green_count

# Example usage (replace 'reads' with your actual reads list)
testreads = [test_read]  # Example with the test read, replace with your full list of testreads

# Count testreads mapping to the red and green gene exons
red_count_test, green_count_test = count_reads_in_exons(testreads, reference_sequence, red_exons, green_exons)

print(f"Red gene read count: {red_count_test}")
print(f"Green gene read count: {green_count_test}")
red_counts = [256, 57, 42.5, 100, 203.5, 285]
green_counts = [256, 193.5, 132.5, 133.5, 334.5, 285 ]


def load_reads(reads_file):
    with open(reads_file, 'r') as f:
        reads = [line.strip().replace('N', 'A') for line in f.readlines()]
    return reads
reads = np.loadtxt('C:/Users/yogan/Downloads/chrX_bwt/chrX_bwt/reads', dtype=str)
red_count, green_count = count_reads_in_exons(reads, reference_sequence, red_exons, green_exons)
print (red_count, green_count)



# In[38]:


import gc
red_counts_p = [57, 42.5, 100, 203.5]
green_count_p = [193.5, 132.5, 133.5, 334.5 ]
# Force garbage collection
collected = gc.collect()
print(f"Garbage collected: {collected} objects")


# In[45]:


def compute_probabilities(red_counts_p: np.ndarray, green_count_p: np.ndarray) -> List[float]:
    """Compute probabilities for each row of configuration_red using 
    probability_red = configuration_red * red_counts_p[i] / (red_counts_p[i] + green_count_p[i])
    where i ranges from 0 to 3.
    """
    
    # Define 4x4 matrix configuration for red
    configuration_red = [
        [0.5, 0.5, 0.5, 0.5],
        [1.0, 1.0, 0.0, 0.0],
        [0.33, 0.33, 1.0, 1.0],
        [0.33, 0.33, 0.33, 1.0]
    ]
    
    probabilities = []
    
    # Calculate probability for each row in configuration_red
    for red_prob_row in configuration_red:
        config_probability = 1.0  # Start with multiplicative identity
        
        # Apply the formula to each exon position i (0 to 3)
        for i in range(4):  # Index from 0 to 3 for each element in the row
            # Avoid division by zero by ensuring both counts are non-zero
            if red_counts_p[i] + green_count_p[i] > 0:
                prob_component = red_prob_row[i] * (red_counts_p[i] / (red_counts_p[i] + green_count_p[i]))
            else:
                prob_component = 0.0  # If counts are zero, set probability component to zero
            
            # Multiply to accumulate product of probabilities across the configuration row
            config_probability *= prob_component
        
        probabilities.append(config_probability)
    
    return probabilities
    

red_counts_p = [57, 42.5, 100, 203.5]
green_count_p = [193.5, 132.5, 133.5, 334.5 ]
probability_of_configration = compute_probabilities(red_counts_p, green_count_p)
print(probability_of_configration)


# In[ ]:





# In[36]:


## loading and prepping the data
def load_and_prepare_data(file_paths: List[str], delta: int) -> Tuple[np.ndarray, np.ndarray]:
    """Loads and preprocesses the last column and mapping data files."""
    
    # Load and process the last column data
    content_chunks = []
    with open(file_paths[0], 'r') as file:
        while True:
            chunk = file.read(1024)
            if not chunk:
                break
            content_chunks.append(chunk.replace('\n', ''))
    last_column = np.array(list(''.join(content_chunks)), dtype='S1')
    
    # Load mapping data, skipping rows based on delta
    mapping_data = pd.read_csv(file_paths[2], skiprows=lambda x: x % delta != 0, header=None, dtype=np.uint32)
    index_mapping = mapping_data.squeeze().to_numpy()
    
    return last_column, index_mapping

#processing data with setting bounds and milestones
def adjust_bounds(character, upper_bound, lower_bound, bwt_column, milestones, counts, char_map, delta):
    """Adjusts the upper and lower bounds for a character."""
    
    def get_next_bound(start, end, step, char):
        """Finds the next occurrence of char in the specified range."""
        for idx in range(start, end + step, step):
            if bwt_column[idx] == char:
                rank = rank_query(bwt_column, idx, milestones[char_map[char]], delta)
                return select_query(char, counts, rank)
        return None

    new_upper = get_next_bound(int(upper_bound), int(lower_bound), 1, character)
    new_lower = get_next_bound(int(lower_bound), int(upper_bound), -1, character)
    
    return new_upper, new_lower


def build_milestone_matrix(data_vector, delta):
    """Creates a milestone matrix to store character counts up to each milestone."""
    total_milestones = (len(data_vector) - 1) // delta
    nucleotide_tags = [b'A', b'C', b'G', b'T']
    milestones = np.zeros((4, total_milestones + 1), dtype=np.uint32)
    
    progress = tqdm(desc='Building Milestones', total=4 * (total_milestones + 1))
    
    masks = [data_vector == nucleotide for nucleotide in nucleotide_tags]
    milestones[nucleotide_tags.index(data_vector[0])][0] = 1
    
    for j, mask in enumerate(masks):
        for i in range(1, total_milestones + 1):
            milestones[j][i] = np.sum(mask[delta * (i - 1) + 1:delta * i + 1]) + milestones[j][i - 1]
            progress.update(1)
    
    progress.close()
    return milestones, milestones[:, -1]

# determining rank
def rank_query(data, index, milestone, delta):
    """Calculates the rank of a character up to a specific index using milestones."""
    nearest_milestone = index // delta
    rank_from_milestone = milestone[nearest_milestone]
    start_position = nearest_milestone * delta
    additional_count = np.sum(data[start_position + 1:index + 1] == data[index])
    return rank_from_milestone + additional_count

# selecting the esential
def select_query(char, cumulative_counts, rank):
    """Returns the index of the specified rank for a character."""
    offsets = {
        b'A': 0,
        b'C': cumulative_counts[0],
        b'G': cumulative_counts[0] + cumulative_counts[1],
        b'T': sum(cumulative_counts[:3])
    }
    return offsets[char] + rank - 1

# locating the index
def locate_index(mapping_data, delta, initial_index, last_column, milestones, char_map, cumulative_counts, offset):
    """Finds the mapped index given the initial position and the milestone structure."""
    steps = 0
    
    while initial_index % delta != 0:
        char = last_column[initial_index]
        rank_value = rank_query(last_column, initial_index, milestones[char_map[char]], delta)
        initial_index = select_query(char, cumulative_counts, rank_value)
        steps += 1

    return mapping_data[initial_index // delta] + steps - offset

#partially scanning the data
def partial_scan(mappings, partial_read, bwt_column, milestones, counts, char_tags, delta, offset):
    """Performs a partial scan on the BWT column to locate a sequence."""
    
    lower, upper = 0, 0
    
    # Ensure the last character is valid, or default to 'A'
    last_char = partial_read[-1] if partial_read[-1] in char_tags else b'A'
    
    if last_char != b'A':
        lower = np.sum(counts[:char_tags[last_char]]) - 1
        upper = lower + 1 - counts[char_tags[last_char]]
    else:
        lower = counts[char_tags[last_char]] - 1
        upper = lower + 1 - counts[char_tags[last_char]]
    
    for i in range(1, len(partial_read)):
        # Ensure character validity within each iteration
        current_char = partial_read[-i-1] if partial_read[-i-1] in char_tags else b'A'
        upper, lower = adjust_bounds(current_char, upper, lower, bwt_column, milestones, counts, char_tags, delta)
        
        if upper is None or lower is None:
            return [None]

    # Final result handling based on the adjusted bounds
    if upper == lower:
        return [locate_index(mappings, delta, upper, bwt_column, milestones, char_tags, counts, offset)]
    return [locate_index(mappings, delta, idx, bwt_column, milestones, char_tags, counts, offset) for idx in range(upper, lower + 1)]

def locate_reads(read_sequence, last_column, milestones, cumulative_counts, char_tags, delta, mappings, is_complement=False):
    """Locates the position of a read sequence in the BWT column."""
    
    if not is_complement:
        read_vector = np.array(list(read_sequence), dtype='S1')
        read_vector = read_vector[read_vector != b'\n']
        read_vector[read_vector == b'N'] = b'A'
    else:
        read_vector = read_sequence

    indices = partial_scan(mappings, read_vector, last_column, milestones, cumulative_counts, char_tags, delta, 0)
    if np.sum(np.array(indices) != None) > 0:
        return set(indices), read_vector

    if len(read_sequence) > 49:
        read_batches, offsets = [], []
        for i in range(len(read_sequence) // 33 - 1):
            read_batches.append(read_vector[i * 33:(i + 1) * 33])
            offsets.append(i * 33)
        read_batches.append(read_vector[(i + 1) * 33:])
        offsets.append((i + 1) * 33)

        initial_batch_indices = set(partial_scan(mappings, read_batches[0], last_column, milestones, cumulative_counts, char_tags, delta, 0))
        batch_indices = [set(partial_scan(mappings, read_batches[i], last_column, milestones, cumulative_counts, char_tags, delta, offsets[i])) for i in range(1, len(read_batches))]
        batch_indices.append(initial_batch_indices)
        
        intersected_indices = set.intersection(*[indices - {None} for indices in batch_indices])
        if intersected_indices:
            return intersected_indices, read_vector

    if not is_complement:
        comp_vector = np.copy(read_vector)
        comp_vector[read_vector == b'A'] = b'T'
        comp_vector[read_vector == b'T'] = b'A'
        comp_vector[read_vector == b'C'] = b'G'
        comp_vector[read_vector == b'G'] = b'C'
        comp_vector = np.flip(comp_vector)
        
        indices, reverse_vector = locate_reads(comp_vector, last_column, milestones, cumulative_counts, char_tags, delta, mappings, True)
        if len(indices - {None}) > 0:
            return indices, reverse_vector
    
    return {None}, read_vector

## findding exon counts
def find_exon_counts(exon_ranges, read_positions):
    """Counts occurrences of reads in specified exon ranges."""
    exon_count = [0] * len(exon_ranges)
    
    for idx, (start, end) in enumerate(exon_ranges):
        for positions in read_positions:
            valid_positions = []
            for pos in positions:
                try:
                    # Ensure numeric conversion; decode if in bytes format
                    if isinstance(pos, bytes):
                        pos = pos.decode()
                    valid_positions.append(int(pos))
                except (ValueError, AttributeError):
                    # Skip non-numeric or malformed entries
                    continue
            
            # Check valid positions within the exon range
            if any(start <= pos <= end for pos in valid_positions):
                exon_count[idx] += 1
                
    return exon_count

## main function



# In[ ]:


if __name__ == "__main__":
    # Define file paths and processing parameters
    file_paths = ["C:/Users/yogan/Downloads/chrX_bwt/chrX_bwt/chrX_last_col.txt", 'C:/Users/yogan/Downloads/chrX_bwt/chrX_bwt/reads', "C:/Users/yogan/Downloads/chrX_bwt/chrX_bwt/chrX_map.txt"]
    delta = 13

    print("\nLoading and Preprocessing Data\n")
    last_column_data, index_map = load_and_prepare_data(file_paths, delta)

    milestone_matrix, nucleotide_counts = build_milestone_matrix(last_column_data, delta)
    char_mapping = {b'A': 0, b'C': 1, b'G': 2, b'T': 3}

    # Define exon ranges for analysis
    red_exons = [
        (149249757, 149249868), 
        (149256127, 149256423), 
        (149258412, 149258580), 
        (149260048, 149260213), 
        (149261768, 149262007), 
        (149264290, 149264400)
    ]
     green_exons = [
        (149288166, 149288277), 
        (149293258, 149293554), 
        (149295542, 149295710), 
        (149297178, 149297343), 
        (149298898, 149299137), 
        (149301420, 149301530)
    ]

    # Process reads and calculate exon counts
    interest_indices = set(range(149249600, 149301650))
    found_positions = locate_reads(file_paths[1], last_column_data, milestone_matrix, nucleotide_counts, char_mapping, delta, index_map)
    
    print("Red Exon Counts: ", red_counts)
    print("Green Exon Counts: ", green_counts)

