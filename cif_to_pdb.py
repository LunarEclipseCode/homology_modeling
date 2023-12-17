from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser, PDBIO
import os


def convert_cif_to_pdb(cif_file, pdb_file):
    # Parse CIF file
    mmcif_dict = MMCIF2Dict(cif_file)
    structure_id = mmcif_dict["_entry.id"]
    parser = MMCIFParser()
    structure = parser.get_structure(structure_id, cif_file)

    # Write PDB file
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(pdb_file)


if __name__ == "__main__":
    cif_file = os.path.join(os.getcwd(), f"generated_models", f"AF-P30820-F1-model_v4.cif")
    pdb_file = os.path.join(os.getcwd(), f"generated_models", f"AF_AFP30820F1_alphafold.pdb")
    convert_cif_to_pdb(cif_file, pdb_file)

