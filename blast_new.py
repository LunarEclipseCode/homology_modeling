import requests
import os, subprocess, tempfile
from io import StringIO
from Bio import SeqIO, Align, PDB
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import PDB
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import signal

protein = "AF_AFB9K865F1"

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

def blast_search_local(sequence):
    blastp_path = r'./blast-2.15.0+/bin/blastp.exe'
    db_path = r'./blast-2.15.0+/blast/db/pdbaa'

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_seq_file:
        temp_seq_file.write(sequence)
        temp_seq_file_path = temp_seq_file.name

    cmd = f'"{blastp_path}" -db {db_path} -query {temp_seq_file_path} -outfmt 5'

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        blast_records = NCBIXML.read(StringIO(result.stdout.decode('utf-8')))
        return blast_records

    finally:
        if temp_seq_file_path:
            os.remove(temp_seq_file_path)

def get_pdb_ids(blast_record):
    pdb_ids = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            pdb_ids.append(alignment.accession)
    return pdb_ids

# Run local BLAST search
blast_records = blast_search_local(query_sequence)
all_pdb_ids = get_pdb_ids(blast_records)
all_pdb_ids = [pdb_id.split('_')[0] for pdb_id in all_pdb_ids]
print(all_pdb_ids)

all_sequence = []
def handle_sigint(signum, frame):
    print("Received SIGINT. Stopping execution.")
    os._exit(1)

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

    with ThreadPoolExecutor() as executor:
        future_to_pdb_id = {executor.submit(fetch_sequence_with_id, pdb_id): pdb_id for pdb_id in pdb_ids}

        for future in as_completed(future_to_pdb_id):
            pdb_id = future_to_pdb_id[future]
            try:
                result = future.result()
                if None not in result:
                    all_sequence.append(result)
            except Exception as e:
                print(f"Error fetching data for {pdb_id}: {e}")
        
            except KeyboardInterrupt:
                print("Received KeyboardInterrupt. Shutting down gracefully.")
                for future in future_to_pdb_id:
                    future.cancel()     
                executor.shutdown(wait=True)

fetch_and_save_pdb_files(all_pdb_ids)   
sorted_sequence = sorted(all_sequence, key=lambda x: all_pdb_ids.index(x[0]) if x[0] in all_pdb_ids else float('inf'))

# template_sequence = sorted_sequence[0][1]
# # sequence alignment
# aligner = Align.PairwiseAligner()
# alignments = aligner.align(query_sequence, template_sequence)

# max_val = alignments[0].score
# max_val_count = sum(1 for alignment in alignments if alignment.score == max_val)

# print(max_val_count)

# best_alignment = max(alignments, key=lambda x: x.score)

# print(np.max(alignment_scores))

# for i in range(len(alignments)):
#     print(alignments[i][0])
#     print(alignments[i][1])

# for align in alignments:
#     print(align.score)

# best_align1 = alignments[0]
# if len(alignments) > 1:
#     best_align2 = alignments[1]

# aligned_query_seq = best_alignment[0]
# aligned_template_seq = best_alignment[1]
# print("Aligned Query Sequence:", aligned_query_seq)   #best_align1[0]
# print("Aligned Template Sequence:", aligned_template_seq)     #best_align[1]

# # Generate 3D model
# output_pdb_path = f"{protein}.pdb"
# template_pdb_path = f"{template_id}.pdb"


# def update_ter_line(input_lines):
#     # Find the index of the line before TER
#     lines = input_lines
#     for i, line in enumerate(lines):
#         if line.startswith("ATOM"):
#             # Extract information from the ATOM line
#             atom_serial_number = int(line[6:11])
#             amino_acid = line[17:20]
#             residue_number = int(line[22:26])
#             chain_id = line[21]

#             # Update TER line
#             if lines[i + 1].startswith("TER"):
#                 lines[i +
#                       1] = f"TER    {atom_serial_number + 1:4}      {amino_acid} {chain_id} {residue_number:3}\n"

#     return lines

# # Case 1: exact match but template has multiple chains, query has one chain


# def extract_chain(input_pdb_path, chain_id, output_pdb_path):
#     parser = PDB.PDBParser(QUIET=True)
#     structure = parser.get_structure('template', input_pdb_path)

#     io = PDB.PDBIO()
#     io.set_structure(structure[0][chain_id])
#     io.save(output_pdb_path)
#     with open(output_pdb_path, 'r') as infile:
#         lines = infile.readlines()

#     # Keep only ATOM lines and TER lines
#     filtered_lines = [line for line in lines if line.startswith('ATOM') or
#                       line.startswith('TER') or line.startswith('END')]
#     filtered_lines = update_ter_line(filtered_lines)

#     with open(output_pdb_path, 'w') as outfile:
#         outfile.writelines(filtered_lines)


# if (best_align1[0] == best_align1[1]) and 'chains' in template_fasta.lower():
#     extract_chain(template_pdb_path, 'A', output_pdb_path)

# # Case 2: almost matches with one chain from the template protein
# # # five cases:
# # # exact match of letters
# # # letters don't match
# # # gap in query but not in template
# # # gap in template but not in query
# # # gap in both

# # # letters don't match, need to estimate atom position


# def choose_alignment(align1, align2):
#     len1 = len(align1[0])
#     len2 = len(align2[0])
#     index1 = 0
#     index2 = 0
#     for char1, char2 in zip(align1[0], align1[1]):
#         if char1 != char2:
#             index1 += 1
#         else:
#             break

#     for char3, char4 in zip(align2[0], align2[1]):
#         if char3 != char4:
#             index2 += 1
#         else:
#             break

#     index = max(index1, index2) + 3
#     if align1[0][index:len1] == align2[0][index:len1] and align1[1][index:len1] == align2[1][index:len1]:
#         if align2[0][0] == '-' and align1[0][0] != '-':
#             return align2
#         else:
#             return align1


# best_align = best_align1

# if len(alignments) > 1:
#     best_align = choose_alignment(best_align1, best_align2)

# aligned_query_seq = best_align[0]
# aligned_template_seq = best_align[1]


# aa_code = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
#            'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

# swapped_code = {value: key for key, value in aa_code.items()}


# def ret_atom_count():
#     count_map = {one_letter: 0 for one_letter in swapped_code.values()}
#     file_path = 'aa_format_pdb.txt'

#     with open(file_path, 'r') as file:
#         for line in file:
#             codes = line.split()
#             for code in codes:
#                 one_letter = swapped_code.get(code, None)
#                 if one_letter is not None:
#                     count_map[one_letter] += 1

#     count_map = {key: value - 1 for key, value in count_map.items()}
#     return count_map


# aa_atom_count = ret_atom_count()
# collision_pdbfix = []

# print(aa_atom_count)


# def create_pdb():
#     selected_lines = []
#     with open('temp.pdb', 'r') as infile:
#         lines = infile.readlines()
#         line_index = 0

#     collision = []
#     i = 0
#     while i < len(aligned_query_seq) - 1:
#         if aligned_query_seq[i] == '-' and aligned_template_seq[i+1] == '-':
#             collide = (aligned_query_seq[i+1], aligned_template_seq[i], i+1, i)
#             collision.append(collide)
#         if aligned_query_seq[i+1] == '-' and aligned_template_seq[i] == '-':
#             collide = (aligned_query_seq[i], aligned_template_seq[i+1], i, i+1)
#             collision.append(collide)
#         i += 1

#     i = 0

#     while i < len(aligned_query_seq) and i < len(aligned_template_seq) - 1:
#         if aligned_query_seq[i] != aligned_template_seq[i]:
#             if any(i in tup for tup in collision):
#                 tuple = next(tup for tup in collision if i in tup)
#                 line_update = aa_atom_count[tuple[0]]
#                 selected_lines.extend(
#                     lines[line_index:line_index + line_update])
#                 collision_pdbfix.append((tuple[0], tuple[1], line_index))
#                 line_index += aa_atom_count[tuple[1]]
#                 i += 2

#             else:
#                 line_index += aa_atom_count[aligned_template_seq[i]]
#                 i += 1

#         if aligned_query_seq[i] == aligned_template_seq[i]:
#             increment = aa_atom_count[aligned_template_seq[i]]
#             selected_lines.extend(lines[line_index:line_index + increment])
#             line_index += increment
#             i += 1

#     selected_lines.extend(lines[-3:])
#     return selected_lines


# def fix_number_in_lines(input_content):
#     modified_lines = []

#     atom_serial_number = 1
#     residue_sequence_number = 0
#     previous_amino_acid = None

#     for line in input_content:
#         if line.startswith('ATOM') and 'OXT' not in line:
#             # Check if the residue number in the current line is different from the previous line
#             current_amino_acid = line[17:20]
#             if current_amino_acid != previous_amino_acid:
#                 residue_sequence_number += 1
#                 previous_amino_acid = current_amino_acid
#                 track = aa_atom_count[swapped_code[current_amino_acid]]
#             else:
#                 track -= 1
#                 if track == 0:
#                     residue_sequence_number += 1
#                     track = aa_atom_count[swapped_code[current_amino_acid]]

#             line = line[:6] + f'{atom_serial_number: >5}' + \
#                 line[11:22] + f'{residue_sequence_number: >4}' + line[26:]
#             atom_serial_number += 1

#             modified_line = line.rstrip() + '\n'
#             modified_lines.append(modified_line)
#         # else:
#         #     modified_lines.append(line)
#     modified_lines.extend(input_content[-2:])

#     return modified_lines


# def update_residue_name(pdb_param, atom_range, new_residue_name):
#     pdb_lines = pdb_param
#     for i, line in enumerate(pdb_lines):
#         if line.startswith("ATOM"):
#             atom_serial_number = int(line[6:11].strip())
#             if atom_range[0] <= atom_serial_number <= atom_range[1]:
#                 pdb_lines[i] = line[:17] + new_residue_name + line[20:]
#     return pdb_lines


# if (aligned_query_seq != aligned_template_seq) and ('chains' in template_fasta.lower()):
#     extract_chain(template_pdb_path, 'A', f"temp.pdb")
#     content = create_pdb()

#     for i in range(len(collision_pdbfix)):
#         atom_start = collision_pdbfix[i][2] + 1
#         atom_end = collision_pdbfix[i][2] + \
#             aa_atom_count[collision_pdbfix[i][0]]
#         content = update_residue_name(
#             content, (atom_start, atom_end), aa_code[collision_pdbfix[i][0]])

#     content = fix_number_in_lines(content)

#     with open('selected_output.pdb', 'w') as outfile:
#         outfile.writelines(content)

#     with open('selected_output.pdb', 'r') as infile:
#         testing = infile.readlines()

#     testing = update_ter_line(testing)
#     filename = f"{protein}_predict.pdb"
#     with open(filename, 'w') as outfile:
#         outfile.writelines(testing)

#     os.remove('temp.pdb')
#     os.remove('selected_output.pdb')

#     print(best_align[0])
#     print(best_align[1])


