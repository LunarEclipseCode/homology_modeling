from Bio import PDB
from Bio.PDB import Superimposer

def extract_all_atoms(structure):
    return [atom for atom in structure.get_atoms() if atom.id == 'CA']

def calculate_rmsd(pdb1, pdb2):
    # Load PDB structures
    parser = PDB.PDBParser(QUIET=True)
    structure1 = parser.get_structure('structure1', pdb1)
    structure2 = parser.get_structure('structure2', pdb2)

    # Extract all heavy atoms
    atoms1 = extract_all_atoms(structure1)
    atoms2 = extract_all_atoms(structure2)

    # Check if both structures have the same number of atoms
    if len(atoms1) != len(atoms2):
        raise ValueError("The structures have different numbers of atoms.")

    # Superimpose structures
    superimposer = Superimposer()
    superimposer.set_atoms(atoms1, atoms2)
    superimposer.apply(structure2.get_atoms())

    # Calculate RMSD
    rmsd = superimposer.rms

    return rmsd


if __name__ == "__main__":
    # Replace 'path_to_pdb1.pdb' and 'path_to_pdb2.pdb' with the actual paths to your PDB files
    pdb_file1 = "output.pdb"
    pdb_file2 = "AF_AFA0A023GPI8F1.pdb"

    try:
        rmsd_value = calculate_rmsd(pdb_file1, pdb_file2)
        print(f"RMSD between {pdb_file1} and {pdb_file2}: {rmsd_value:.3f} Å")
    except ValueError as e:
        print(f"Error: {e}")
