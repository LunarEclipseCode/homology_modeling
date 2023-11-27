import requests, os
from io import StringIO
from Bio import SeqIO, Align, PDB
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio import PDB
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, PDBIO

protein = "AF_AFB5EZH0F1"

def fetch_sequence_from_fasta(comp_seq):
    url = f'https://www.rcsb.org/fasta/entry/{comp_seq}/display'
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = StringIO(response.text)

        # Parse the FASTA file
        for record in SeqIO.parse(fasta_data, "fasta"):
            return str(record.seq), record.description
    else:
        print(f"Error: Unable to fetch data from {url}")
        return None
     
def blast_search(sequence, database='pdb', num_alignments=1):
    result_handle = NCBIWWW.qblast("blastp", database, sequence, alignments=num_alignments)
    blast_records = NCBIXML.parse(result_handle)
    return blast_records

def get_pdb_id(blast_record):
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            return alignment.accession  # Assuming you want the top hit

def read_mapping_file():
    mapping = {}
    with open('mapping.txt', 'r') as f:
        for line in f:
            query, pdb_id = line.strip().split(' -> ')
            mapping[query] = pdb_id
    return mapping

query_sequence, query_fasta = fetch_sequence_from_fasta(protein)

# if found most similar protein once before, use that from mapping
mapping = read_mapping_file()

if protein in mapping:
    top_hit_pdb_id = mapping[protein]
# elif query_sequence in mapping:
#     top_hit_pdb_id = mapping[query_sequence]
else:
    blast_records = blast_search(query_sequence)
    top_hit_pdb_id = get_pdb_id(next(blast_records))
    with open('mapping.txt', 'a') as f:
        f.write(f"{protein} -> {top_hit_pdb_id}\n")
        # f.write(f"{query_sequence} -> {top_hit_pdb_id}\n")

# download the template protein PDB file
template_id = top_hit_pdb_id[:4]
template_sequence, template_fasta = fetch_sequence_from_fasta(template_id)
url = f"https://files.rcsb.org/view/{template_id}.pdb"

response = requests.get(url)
if response.status_code == 200:
    pdb_file_path = f"{template_id}.pdb"
    with open(pdb_file_path, "w") as pdb_file:
        pdb_file.write(response.text)
    print(f"File saved successfully as {pdb_file_path}")
else:
    print(f"Error: Unable to fetch data from the URL. Status code: {response.status_code}")
    
# sequence alignment
aligner = Align.PairwiseAligner()
alignments = aligner.align(query_sequence, template_sequence)

best_align1 =  alignments[0]
best_align2 = alignments[1]

# print("Aligned Query Sequence:", aligned_query_seq)   #best_align1[0]
# print("Aligned Template Sequence:", aligned_template_seq)     #best_align[1]

# Generate 3D model
output_pdb_path = f"{protein}.pdb"
template_pdb_path = f"{template_id}.pdb"

# Case 1: exact match but template has multiple chains, query has one chain
def extract_chain(input_pdb_path, chain_id, output_pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('template', input_pdb_path)

    io = PDB.PDBIO()
    io.set_structure(structure[0][chain_id])
    io.save(output_pdb_path)
    with open(output_pdb_path, 'r') as infile:
        lines = infile.readlines()

    # Keep only ATOM lines and TER lines
    filtered_lines = [line for line in lines if line.startswith('ATOM')]

    with open(output_pdb_path, 'w') as outfile:
        outfile.writelines(filtered_lines)
    
if (best_align1[0] == best_align1[1]) and 'chains' in template_fasta.lower():
    extract_chain(template_pdb_path, 'A', output_pdb_path)

# Case 2: almost matches with one chain from the template protein
# # five cases:
# # exact match of letters
# # letters don't match
# # gap in query but not in template
# # gap in template but not in query
# # gap in both

# # letters don't match, need to estimate atom position

def choose_alignment(align1, align2):
    len1 = len(align1[0])
    len2 = len(align2[0])
    index1 = 0
    index2 = 0
    for char1, char2 in zip(align1[0], align1[1]):
        if char1 != char2:
            index1 += 1
        else:
            break
        
    for char3, char4 in zip(align2[0], align2[1]):
        if char3 != char4:
            index2 += 1
        else:
            break
    
    index = max(index1, index2)
    
    if align1[0][index:len1] == align2[0][index:len1] and align1[1][index:len1] == align2[1][index:len1]:
        if align1[0]
        
    
    
if (aligned_query_seq != aligned_template_seq) and ('chains' in template_fasta.lower()):
#     extract_chain(template_pdb_path, 'A', f"temp.pdb")
#     create_pdb(query_sequence, f"temp.pdb", aligned_query_seq, aligned_template_seq, f"{protein}.pdb")