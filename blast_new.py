import requests, signal
import os, subprocess, tempfile
from io import StringIO
from Bio import SeqIO, Align, PDB
from Bio.Blast import NCBIXML
from concurrent.futures import ThreadPoolExecutor, as_completed

protein = "AF_AFB9K865F1"

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
print(all_pdb_ids)

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
print(sorted_protein[0][3][0])
print(sorted_protein[0][3][1])

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
# and merge_gaps function.

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
    return [aligned_query, aligned_template]
    
    
for k in range (len(sorted_protein)):
    sorted_protein[k] = (sorted_protein[k][0], sorted_protein[k][1], sorted_protein[k][2], 
                         process_alignment(sorted_protein[k][3]))
    
print(sorted_protein[0][3][0])
print(sorted_protein[0][3][1])
    
# For next step, among all the sequence alignments, find the one that exactly matches with the input sequence upto the highest index,
# let's say protein 4K1Y exactly matches upto index 50 of input sequence. That means all other proteins matches with input some index before that.
# Now, look for proteins that matches exactly with our input from index 51 onwards. Let's say protein 5NEM exactly matches with upto index 91.
# Now, backtrack and check upto which point 5NEM's match continues.Let's say, we find from index 43 to 91 it's exact match. 
# So, take index 0 to 47 of 4K1Y, then index 48 to "91" and keep on going until the end of sequence. Then merge all these amino acid
# chunks into one pdb file.

# I put 91 in quotation above, because when you find the next match and backtrack, it probably will be less than 91.
# Also, we can relax the exact match condition, when we have the grantham distance thing fully implemented.
# I haven't yet start doing that part, so if you want to do that, go ahead...
