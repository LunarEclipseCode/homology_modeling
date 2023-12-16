import requests, signal, re
import os, subprocess, tempfile
from io import StringIO
from Bio import SeqIO, Align, PDB
from Bio.Blast import NCBIXML
from concurrent.futures import ThreadPoolExecutor, as_completed
from pard.grantham import grantham

# protein = "AF_AFB9K865F1"
protein = "AF_AFA5IEB3F1"

# given the protein id, return the sequence, and first line of the fasta file
# first line of fasta file to know whether the protein has multiple identical chains
def fetch_sequence_from_fasta(comp_seq):
    url = f'https://www.rcsb.org/fasta/entry/{comp_seq}/display'
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = StringIO(response.text)

        for record in SeqIO.parse(fasta_data, "fasta"):
            return str(record.seq), record.description
    else:
        print(f"Error: Unable to fetch data from {url}")
        return None, None
    
query_sequence, query_fasta = fetch_sequence_from_fasta(protein)

# given a input sequence, returns the list of proteins that matches the best
def blast_search_local(sequence):
    blastp_path = r'./blast-2.15.0+/bin/blastp.exe'
    db_path = r'./blast-2.15.0+/blast/db/pdbaa'

    # cannot directly pass in sequence, need to store in a temporary file first
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_seq_file:
        temp_seq_file.write(sequence)
        temp_seq_file_path = temp_seq_file.name

    cmd = f'"{blastp_path}" -db {db_path} -query {temp_seq_file_path} -outfmt 5'

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        blast_records = NCBIXML.read(StringIO(result.stdout.decode('utf-8')))
        return blast_records

    finally:
        # delete that temporary file
        if temp_seq_file_path:
            os.remove(temp_seq_file_path)

# blast_search_local doesn't directly give the protein ids, it's hidden in the blast_records output
# post-processing function to extract the pdb_ids from blast_records
def get_pdb_ids(blast_record):
    pdb_ids = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            pdb_ids.append(alignment.accession)
    return pdb_ids

# Run local BLAST search
blast_records = blast_search_local(query_sequence)

# In all_pdb_ids, the protein ids are from best match to least match order
all_pdb_ids = get_pdb_ids(blast_records)
all_pdb_ids = [pdb_id.split('_')[0] for pdb_id in all_pdb_ids]
all_pdb_ids = all_pdb_ids[:100]
# print(all_pdb_ids)

# Store the sequence of all the protein ids in this array
all_protein = []

# when running multi-threading task like downloading files in parallel using concurrent package
# Ctrl + C doesn't work to terminate code, need explicit signal handler for keyboard interrupt like Ctrl + C
def handle_sigint(signum, frame):
    print("Received SIGINT. Stopping execution.")
    os._exit(1)

# When we have list of proteins that matches best with input sequence, get the sequence
# and also the pdb file of all the proteins.
def fetch_and_save_pdb_files(pdb_ids):
    os.makedirs("pdb_files", exist_ok=True)
    def fetch_sequence_with_id(pdb_id):
        template_sequence, template_fasta = fetch_sequence_from_fasta(pdb_id)
        pdb_file_path = os.path.join("pdb_files", f"{pdb_id}.pdb")
        
        if os.path.exists(pdb_file_path):
            return pdb_id, template_sequence, template_fasta
        
        url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            with open(pdb_file_path, "w") as pdb_file:
                pdb_file.write(response.text)
        else:
            pdb_id = None
            
        return pdb_id, template_sequence, template_fasta

    # download multiple files in parallel
    with ThreadPoolExecutor() as executor:
        future_to_pdb_id = {executor.submit(fetch_sequence_with_id, pdb_id): pdb_id for pdb_id in pdb_ids}

        for future in as_completed(future_to_pdb_id):
            pdb_id = future_to_pdb_id[future]
            try:
                result = future.result()
                if None not in result:
                    all_protein.append(result)
            except Exception as e:
                print(f"Error fetching data for {pdb_id}: {e}")
        
            except KeyboardInterrupt:
                print("Received KeyboardInterrupt. Shutting down gracefully.")
                for future in future_to_pdb_id:
                    future.cancel()     
                executor.shutdown(wait=True)

fetch_and_save_pdb_files(all_pdb_ids)   

# Now do sequence alignment between input sequence and all the protein sequences
# The alignment will be stored as the 4th element in the tuple
def perform_alignment(protein_tuple):
    aligner = Align.PairwiseAligner()
    protein_id, sequence, fasta = protein_tuple
    # index 0 because that's the highest scoring alignment
    # Actually, there are a lot of highest scoring alignment with same score, to choose between these takes long
    # time even with multi-threading just because alignment array length is really long. So, we will do some post-processing later.
    alignment = aligner.align(query_sequence, sequence)[0]
    return (protein_id, sequence, fasta, alignment)

# multi-threading again to speed up the alignment process
with ThreadPoolExecutor() as executor:
    futures = {executor.submit(perform_alignment, protein): protein for protein in all_protein}

    aligned_proteins = []
    for future in as_completed(futures):
        protein_tuple = futures[future]
        try:
            result = future.result()
            aligned_proteins.append(result)
        except Exception as e:
            print(f"Error processing {protein_tuple}: {e}")
        except KeyboardInterrupt:
                print("Received KeyboardInterrupt. Shutting down gracefully.")
                for future in futures:
                    future.cancel()     
                executor.shutdown(wait=True)

# When running things in parallel, some tasks will finish before others, and the order in all_pdb_ids
# will not be maintained anymore. Use sorting to put the proteins back into best match to least match order
sorted_protein = sorted(aligned_proteins, key=lambda x: all_pdb_ids.index(x[0]) if x[0] in all_pdb_ids else float('inf'))

# At this point the sorted_protein data structure has got quite deep, so let's recap.
# sorted_protein[0] gives the tuple with structure (protein_id, sequence, first line of fasta file, sequence alignment)
# Now, to see sequnece alignment, you would use sorted_protein[0][3]
# Now, inside this, sorted_protein[0][3][0] gives aligned input sequence, sorted_protein[0][3][0] -> aligned protein sequence
# print(sorted_protein[0][3][0])
# print(sorted_protein[0][3][1])

# Now by trying different sequence alignment, some observations:
# In this scenario.
# Aligned Input: M------------FRGVGTAIVTPFKN...
# Aligned Query: MGSDKIHHHHHHMFRGVGTAIVTPFKN...

# However, we would like the following alignment 
# Aligned Input: ------------MFRGVGTAIVTPFKN...
# Aligned Query: MGSDKIHHHHHHMFRGVGTAIVTPFKN...

# The second alignment has multiple advantages. Firstly, as we are starting with gaps in input sequence, that means
# we can just ignore those amino acids in the template protein sequence. Not to mention, even after 'swapping' the gaps
# and the 'M' , our 'M' still matches with the template protein sequence in the new index

# Now, it happens to be the scoring parameters Align.PairwiseAligner() uses doesn't like alignment starting with gap in input
# sequence. The scoring function of PairwiseAligner is match_score = 2, mismatch_score = 0, gap opening penalty = 0, 
# gap extension penalty = -0.5. I tried altering these parameters little bit, but they change the alignment output quite
# drastically. So, I used option (2) post-processing

# The sequence alignment from Pairwise Aligner is great but it has some less optimal choices. One example is shown above.
# Other common scenario is:
# AKMT-LVSYPATLD-YNV
# AKMTK-VSYP-TL-KYNV
# Notice there are some "one-off gaps" in the alignment. For example:
# -L       D-
# K-       -K
# However, it's easier for us to handle mismatch amino acids using Grantham distance than these gaps. So, we want
# to merge these type of gaps and get the following alignment:
# AKMTLVSYPATLDYNV
# AKMTKVSYP-TLKYNV
# As an extension to this problem, we might need to merge more general gaps too. Example:
# --RN      should be merged into    RN
# KM--                               KM

# Same goes for the first problem too. We might have an alignment like below:
# MKL--------------FRGVGTAIVTPFKN...
# MKLGSDKIHHHHHHMKLFRGVGTAIVTPFKN...
# and we should turn it into:
# --------------MKLFRGVGTAIVTPFKN...
# MKLGSDKIHHHHHHMKLFRGVGTAIVTPFKN...

# These two post-processing are done for all the sequence alignments using the process_alignment 
# # and merge_gaps function.
# print(query_sequence)

def merge_gaps(input1, input2):
    output1 = ""
    output2 = ""
    i = 0
    
    while i < len(input1):
        if input1[i] == '-' and input2[i] != '-':
            track = 1
            for j in range (i+1, len(input1), 1):
                if input1[j] == '-' and input2[j] != '-':
                    track += 1
                else:
                    break
                
            if input1[i + track: i+2*track].isalpha() and input2[i + track: i+2*track] == '-' * track:
                output1 += input1[i + track: i+2*track]
                output2 += input2[i:i+track]
                i += 2*track
            else:
                output1 += input1[i]
                output2 += input2[i]
                i += 1
        
        elif input1[i] != '-' and input2[i] == '-':
            track = 1
            for j in range (i+1, len(input1), 1):
                if input1[j] != '-' and input2[j] == '-':
                    track += 1
                else:
                    break
                
            if input1[i + track: i+2*track]== '-' * track and input2[i + track: i+2*track].isalpha():
                output1 += input1[i:i+track]
                output2 += input2[i + track: i+2*track]
                i += 2*track
            else:
                output1 += input1[i]
                output2 += input2[i]
                i += 1
        else:
            output1 += input1[i]
            output2 += input2[i]
            i += 1
            
    return output1, output2

def merge_gaps_v2(query, template):
    output_query = ""
    output_template = ""
    i = 0
    while i < len(query):
        if query[i] == '-' and template[i] != '-':
            track = 1
            track_end = -1
            for j in range (i+1, len(query), 1):
                if query[j] == template[j]:
                    track += 1
                elif query[j] != '-' and template[j]== '-':
                    track += 1
                    track_end = 1
                    break
                else:
                    break
                
            if track != 0 and track_end == 1 and query[i+2:i+track] == template[i+1:i+track-1]:
                output_query += query[i+1:i+track]
                output_template += template[i:i+track -1]
                i += track
            else:
                output_query += query[i]
                output_template += template[i]
                i += 1
                
        elif query[i] != '-' and template[i] == '-':
            track = 1
            track_end = -1
            for j in range (i+1, len(query), 1):
                if query[j] == template[j]:
                    track += 1
                elif query[j] == '-' and template[j] != '-':
                    track += 1
                    track_end = 1
                    break
                else:
                    break

            if track != 0 and track_end == 1 and template[i+2:i+track] == query[i+1:i+track-1]:
                output_query += query[i:i+track -1]
                output_template += template[i+1:i+track]
                i += track
            else:
                output_query += query[i]
                output_template += template[i]
                i += 1
                
        else:
            output_query += query[i]
            output_template += template[i]
            i += 1
        
    return output_query, output_template   

def process_alignment(alignment):
    aligned_query = alignment[0]
    aligned_template = alignment[1]
    
    # swapping gaps and amino acids section
    start_index = aligned_query.find("-")
    end_index = None
    
    if start_index != -1:
        for i in range(start_index + 1, len(aligned_query), 1):
            if aligned_query[i] == '-':
                end_index = i
            else:
                break
    
    if start_index is not None and end_index is not None:
        swapped_string = aligned_query[start_index:end_index + 1] + aligned_query[:start_index] + aligned_query[end_index + 1:]
        
        move_len = len(aligned_query[:start_index])
        index1 = end_index - start_index + 1    # start_check
        index2 = index1 + move_len              # end_check
        
        if swapped_string[index1:index2] == aligned_template[index1:index2]:
            aligned_query = swapped_string
            
    # gap merging adjustment 
    aligned_query, aligned_template = merge_gaps(aligned_query, aligned_template)
    aligned_query, aligned_template = merge_gaps_v2 (aligned_query, aligned_template)

    return [aligned_query, aligned_template]
    
    
for k in range (len(sorted_protein)):
    sorted_protein[k] = (sorted_protein[k][0], sorted_protein[k][1], sorted_protein[k][2], 
                         process_alignment(sorted_protein[k][3]))
    
# print(sorted_protein[0][3][0])
# print(sorted_protein[0][3][1])

combine_pdb = []

def find_longest_exact_match(tuples_array):
    longest_match_length, longest_true_index = 0, 0
    matching_tuple, match_id = None, None
    stop_index, true_index = 0, 0
    initial_run = True

    while true_index < len(query_sequence) -1:
        if not initial_run:
            true_index = combine_pdb[-1][3]
        longest_match_length = 0
        longest_true_index = 0
        for protein_id, sequence, fasta_file, (aligned_input, aligned_template) in tuples_array:
            match_length = 0
            
            if initial_run == True:
                for i, char in enumerate(aligned_input):
                    if char.isalpha():
                        start_index = i
                        break
            else:
                count = -1
                for i, char in enumerate(aligned_input):
                    if char.isalpha():
                        count += 1
                        if count == true_index:
                            break
            
            start_index = i
            for i in range(start_index, len(aligned_input)):
                if aligned_input[i] == aligned_template[i]:
                    match_length += 1
                elif aligned_input[i] != aligned_template[i] and aligned_input[i].isalpha() and aligned_template[i].isalpha():
                    if grantham(aligned_input[i], aligned_template[i]) <= 150:
                        match_length += 1
                    else:
                        break
                # elif aligned_input[i] == '-' and aligned_template[i] != '-':
                #     match_length += 1
                else:
                    break 
            
            # Update the result if a longer match is found
            if match_length > longest_match_length:
                true_index = - 1
                for m in range(len(aligned_input)):
                    if aligned_input[m].isalpha():
                        true_index += 1
                    if m == i:
                        break
                    
                if true_index > longest_true_index:
                    longest_match_length = match_length
                    longest_true_index = true_index
                    match_id = protein_id
                    stop_index = i
                    match_start_index = start_index
                    
                    backtrack = 0
                    if not initial_run:
                        backtrack_end = combine_pdb[-1][1]
                        protein_tuple = next((tup for tup in sorted_protein if tup[0] == combine_pdb[-1][0]), None)
                        convert_index = -1
                        for i, char in enumerate(protein_tuple[3][0]):
                            if char.isalpha():
                                convert_index += 1
                                if i == backtrack_end:
                                    break
                        
                        current_index = -1
                        for i, char in enumerate(aligned_input):
                            if char.isalpha():
                                current_index += 1
                                if current_index == convert_index:
                                    break
                            
                        
                        for w in range(start_index -1, i - 1, -1):
                            if aligned_input[w] == aligned_template[w]:
                                backtrack += 1
                            elif aligned_input[w] == '-' and aligned_template[w] != '-':
                                 backtrack += 1
                            else:
                                break
   
        combine_pdb.append((match_id, match_start_index, stop_index, longest_true_index, backtrack)) 
        initial_run = False                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    return combine_pdb

protein_list = find_longest_exact_match(sorted_protein)
print(protein_list)

unique_protein_ids = set(item[0] for item in protein_list)

# Iterate through unique protein_ids and find corresponding tuple in the first array
for protein_id in unique_protein_ids:
    matching_tuple = next((tup for tup in sorted_protein if tup[0] == protein_id), None)
    if matching_tuple is not None:
        aligned_input, aligned_template = matching_tuple[3]
        print(f"For protein_id {protein_id}:")
        print("Aligned Input:", aligned_input)
        print("Aligned Template:", aligned_template)
        print()
    else:
        print(f"No matching tuple found for protein_id {protein_id}")

    
# For next step, among all the sequence alignments, find the one that exactly matches with the input sequence upto the highest index,
# let's say protein 4K1Y exactly matches upto index 50 of input sequence. That means all other proteins matches with input some index before that.
# Now, look for proteins that matches exactly with our input from index 51 onwards. Let's say protein 5NEM exactly matches with upto index 91.
# Now, backtrack and check upto which point 5NEM's match continues.Let's say, we find from index 43 to 91 it's exact match. 
# So, take index 0 to 47 of 4K1Y, then index 48 to "91" and keep on going until the end of sequence. Then merge all these amino acid
# chunks into one pdb file.

# I put 91 in quotation above, because when you find the next match and backtrack, it probably will be less than 91.
# Also, we can relax the exact match condition, when we have the grantham distance thing fully implemented.
# I haven't yet start doing that part, so if you want to do that, go ahead...

output_pdb_path = f"{protein}.pdb"

aa_code = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
           'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

swapped_code = {value: key for key, value in aa_code.items()}


def update_ter_line(input_lines):
    # Find the index of the line before TER
    lines = input_lines
    for i, line in enumerate(lines):
        if line.startswith("ATOM"):
            # Extract information from the ATOM line
            atom_serial_number = int(line[6:11])
            amino_acid = line[17:20]
            residue_number = int(line[22:26])
            chain_id = line[21]

            # Update TER line
            if lines[i + 1].startswith("TER"):
                lines[i +1] = f"TER    {atom_serial_number + 1:4}      {amino_acid} {chain_id} {residue_number:3}\n"
    return lines

def ret_atom_count():
    count_map = {one_letter: 0 for one_letter in swapped_code.values()}
    file_path = 'aa_format_pdb.txt'

    with open(file_path, 'r') as file:
        for line in file:
            codes = line.split()
            for code in codes:
                one_letter = swapped_code.get(code, None)
                if one_letter is not None:
                    count_map[one_letter] += 1

    count_map = {key: value - 1 for key, value in count_map.items()}
    return count_map

def extract_chain(input_pdb_path, chain_id, output_pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('template', input_pdb_path)

    io = PDB.PDBIO()
    io.set_structure(structure[0][chain_id])
    io.save(output_pdb_path)
    with open(output_pdb_path, 'r') as infile:
        lines = infile.readlines()

    # Keep only ATOM lines and TER lines
    filtered_lines = [line for line in lines if line.startswith('ATOM') or
                      line.startswith('TER') or line.startswith('END')]
    filtered_lines = update_ter_line(filtered_lines)

    with open(output_pdb_path, 'w') as outfile:
        outfile.writelines(filtered_lines)
        

def extract_aa_sequence(pdb_file_path):
    parser = PDB.PDBParser()
    structure = parser.get_structure('protein', pdb_file_path)
    amino_acid_sequence = ''

    for model in structure:
        for chain in model:
            total_residues = len(list(chain.get_residues()))

            for i, residue in enumerate(chain):
                if PDB.is_aa(residue):
                    resname = residue.get_resname()
                    amino_acid_sequence += swapped_code[resname]

                    count = 0
                    for atom in residue:
                        count += 1
                        
                    if i == total_residues - 1:
                        count -= 1

                    if count != aa_atom_count[swapped_code[resname]]:
                        print("Warning!", residue)

    return amino_acid_sequence

def aa_processing(current, new, pdb_lines):
    hydrophobic = ['CYS', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
    special_case = ['CYS', 'GLY', 'PRO']
    polar_uncharged = ['SER', 'THR', 'ASN', 'GLN']
    polar_charged = ['ARG']
    
    if new == swapped_code['GLY']:
        pdb_lines = pdb_lines[:4]
        pdb_lines = [line.replace(line.split()[3], 'GLY') for line in pdb_lines]

    return pdb_lines
        

ideal_entries = {}

def ideal_aa_entries():
    current_amino_acid = None
    with open('aa_format_pdb.txt', 'r') as file:
        for line in file:
            if line.startswith('='):
                continue
            elif line.startswith('ATOM'):
                # Extract the amino acid name from the line
                amino_acid_name = line.split()[3]

                if current_amino_acid != amino_acid_name:
                    current_amino_acid = amino_acid_name
                    ideal_entries[current_amino_acid] = []

                ideal_entries[current_amino_acid].append(line)


def fix_number_in_lines(input_content):
    modified_lines = []

    atom_serial_number = 1
    residue_sequence_number = 0
    previous_amino_acid = None

    for line in input_content:
        if line.startswith('ATOM') and 'OXT' not in line:
            # Check if the residue number in the current line is different from the previous line
            current_amino_acid = line[17:20]
            if current_amino_acid != previous_amino_acid:
                residue_sequence_number += 1
                previous_amino_acid = current_amino_acid
                track = aa_atom_count[swapped_code[current_amino_acid]]
            else:
                track -= 1
                if track == 0:
                    residue_sequence_number += 1
                    track = aa_atom_count[swapped_code[current_amino_acid]]

            line = line[:6] + f'{atom_serial_number: >5}' + \
                line[11:22] + f'{residue_sequence_number: >4}' + line[26:]
            atom_serial_number += 1

            modified_line = line.rstrip() + '\n'
            modified_lines.append(modified_line)

    modified_lines.extend(input_content[-3:])
    line = input_content[-3]
    modified_lines[-3] = line[:6] + f'{atom_serial_number: >5}' + line[11:22] + f'{residue_sequence_number: >4}' + line[26:]
   
    return modified_lines
    

def create_pdb(aligned_sequences, begin_index, end_index, fasta_aa):
    selected_lines = []
    aligned_query_seq = aligned_sequences[0]
    aligned_template_seq = aligned_sequences[1]
    
    with open('temp.pdb', 'r') as infile:
        lines = infile.readlines()
        line_index = 0

    collision = []
    pdb_aa = extract_aa_sequence("temp.pdb")
    
    aligner = Align.PairwiseAligner()
    aligner.extend_gap_score = 0.1
    alignment = aligner.align(pdb_aa, fasta_aa)[0]
    print("Aligned PDB Sequence:", alignment[0])
    print("Aligned FASTA Sequence:", alignment[1])

    for i in range(begin_index):
        if aligned_template_seq != '-' and alignment[0][i] != '-':
            line_index += aa_atom_count[aligned_template_seq[i]]

    # print(line_index)
    for i in range(begin_index, end_index+1, 1):
        if aligned_query_seq[i] == aligned_template_seq[i]:
            increment = aa_atom_count[aligned_template_seq[i]]
            atom_records = lines[line_index:line_index + increment]     
            residues = set()
            for line in atom_records:
                residue_name = line[17:20].strip()
                residues.add(residue_name)
            
            # print(residues)
            
            
            if len(residues) > 1:
                ideal_entry = ideal_entries[aa_code[aligned_template_seq[i]]]
                
                expected_entries = []
                specific_lines = [element for element in atom_records if aa_code[aligned_template_seq[i]] in element]
                current_entries = []
                for line in ideal_entry:
                    atom_type = line[12:16].strip()
                    expected_entries.append(atom_type)
                
                # print(expected_entries)
                # print(current_entries)
            #     # Get the missing entries
            #     missing_entries = set(expected_entries) - set(amino_acids)
            #     increment -= len(missing_entries)
            #     # print(missing_entries)
            #     # for missing_entry in missing_entries:
            #     #     for record in ideal_entry:
            #     #         if f' {missing_entry} ' in record:
            #     #             lines_to_add = [record.rstrip()]
            #     #             lines[line_index:line_index] = lines_to_add
            #     #             line_index += len(lines_to_add)
            #     print(amino_acids)
            selected_lines.extend(atom_records)
            line_index += increment
        
        if aligned_query_seq[i] != aligned_template_seq[i]:
            if aligned_query_seq[i].isalpha() and aligned_template_seq[i].isalpha():
                process_lines = []
                increment = aa_atom_count[aligned_template_seq[i]]
                process_lines = lines[line_index:line_index + increment]
                amino_acids = list(set([re.search(r'ATOM\s+\d+\s+\w+\s+(\w+)', record).group(1) for record in process_lines]))
                if len(amino_acids) > 1:
                    amino_acid = aligned_template_seq[i]
                    ideal_entry = ideal_entries.get(amino_acid, [])
                    expected_entries = [re.search(r'ATOM\s+\d+\s+(\w+)\s+', record).group(1) for record in ideal_entry]

                    # Get the missing entries
                    missing_entries = set(expected_entries) - set(amino_acids)
                    increment -= len(missing_entries)
                else:
                    line_index += aa_atom_count[aligned_template_seq[i]]
                # for missing_entry in missing_entries:
                #     for record in ideal_entry:
                #         if f' {missing_entry} ' in record:
                #             lines_to_add = [record.rstrip()]
                #             lines[line_index:line_index] = lines_to_add
                #             line_index += len(lines_to_add)
                # print(amino_acids)
                # print(amino_acids)

                selected_lines.extend(aa_processing(aligned_template_seq[i], aligned_query_seq[i], process_lines))
            
            
            if aligned_query_seq[i] == '-' and aligned_template_seq[i].isalpha():
                line_index += aa_atom_count[aligned_template_seq[i]]
    #       
                
    # while i < len(aligned_query_seq) and i < len(aligned_template_seq) - 1:
    #     if aligned_query_seq[i] != aligned_template_seq[i]:
    #         if any(i in tup for tup in collision):
    #             tuple = next(tup for tup in collision if i in tup)
    #             line_update = aa_atom_count[tuple[0]]
    #             selected_lines.extend(
    #                 lines[line_index:line_index + line_update])
    #             collision_pdbfix.append((tuple[0], tuple[1], line_index))
    #             line_index += aa_atom_count[tuple[1]]
    #             i += 2

    #         else:
    #             line_index += aa_atom_count[aligned_template_seq[i]]
    #             i += 1

    #     if aligned_query_seq[i] == aligned_template_seq[i]:
    #         increment = aa_atom_count[aligned_template_seq[i]]
    #         selected_lines.extend(lines[line_index:line_index + increment])
    #         line_index += increment
    #         i += 1

    selected_lines.extend(lines[-3:])
    
    return selected_lines

aa_atom_count = ret_atom_count()
print(aa_atom_count)

ideal_aa_entries()

for i in range(len(protein_list)):
    protein_tuple = next((tup for tup in sorted_protein if tup[0] == protein_list[i][0]), None)
    template_pdb_path = os.path.join(os.getcwd(), f"pdb_files", f"{protein_list[i][0]}.pdb")
    
    if 'chains' in protein_tuple[2].lower():
        extract_chain(template_pdb_path, 'A', "temp.pdb")
        content = create_pdb(protein_tuple[3], protein_list[i][1], protein_list[i][2], protein_tuple[1])
    
    content = fix_number_in_lines(content)
    
    with open('selected_output.pdb', 'w') as outfile:
        outfile.writelines(content)

    with open('selected_output.pdb', 'r') as infile:
        testing = infile.readlines()
        
    testing = update_ter_line(testing)
    
    os.makedirs("trial_models", exist_ok=True)
    file_path = os.path.join("trial_models", f"{protein}_predict.pdb")
    with open(file_path, 'w') as outfile:
         outfile.writelines(testing)

    os.remove('temp.pdb')
    os.remove('selected_output.pdb')