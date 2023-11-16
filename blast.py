from Bio import SeqIO
from Bio import pairwise2
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import seq3
from Bio.PDB import PDBParser, PPBuilder
import requests
import io

print("Step 1")

input_sequence = "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG"

print("Step 2")
result_handle = NCBIWWW.qblast("blastp", "pdb", input_sequence)
blast_records = NCBIXML.read(result_handle)
print(blast_records)
hit = blast_records.alignments[0]
hit_id = hit.title.split('|')[3]

print("Step 3")
pdb_file = f"{hit_id}.pdb"
url = f"https://files.rcsb.org/download/{pdb_file}"
response = requests.get(url)
with open(pdb_file, 'wb') as f:
    f.write(response.content)

print("Step 4")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

pdb_sequence = ''
for model in structure:
    for chain in model:
        pdb_sequence += PPBuilder().build_peptides(chain)[0].get_sequence()
print(pdb_sequence)

print("Step 5")
input_seq_record = SeqIO.read(io.StringIO(
    f">{hit_id}\n{input_sequence}\n"), "fasta")
pdb_seq_record = SeqIO.read(io.StringIO(
    f">{hit_id}\n{pdb_sequence}\n"), "fasta")

alignment = pairwise2.align.globalxx(
    input_seq_record.seq, pdb_seq_record.seq, one_alignment_only=True)[0]
input_aligned, pdb_aligned, score, start, end = alignment
print(input_aligned)
print(pdb_aligned)
print(score)
print(start)
print(end)
