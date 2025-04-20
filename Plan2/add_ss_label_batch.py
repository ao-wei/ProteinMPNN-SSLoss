import argparse
import os
import torch
import subprocess
from Bio.PDB.DSSP import make_dssp_dict

aa_map = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
    'X': 'GLY'
}

mainchain_atoms = ['N', 'CA', 'C', 'O']

def get_chain_id_from_filename(filename):
    name = os.path.splitext(os.path.basename(filename))[0]
    if "_" not in name:
        raise ValueError(f"{filename}")
    return name.split("_")[1]

def pt_to_pdb(pt_path, pdb_path):
    data = torch.load(pt_path)
    seq = data['seq']
    xyz = data['xyz']
    mask = data['mask']
    with open(pdb_path, 'w') as f:
        atom_index = 1
        for res_idx, aa in enumerate(seq):
            resname = aa_map.get(aa, 'GLY')
            for atom_i, atom_name in enumerate(mainchain_atoms):
                if not mask[res_idx, atom_i]:
                    continue
                x, y, z = xyz[res_idx, atom_i].tolist()
                line = (
                    f"ATOM  {atom_index:5d} {atom_name:>4s} {resname:>3s} A"
                    f"{res_idx+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                    f"  1.00  0.00           {atom_name[0]}\n"
                )
                f.write(line)
                atom_index += 1
        f.write("END\n")

def run_dssp(pdb_path, dssp_path):
    subprocess.run(['mkdssp', pdb_path, dssp_path], check=True)

def extract_q3_from_dssp(dssp_file, chain_id, seq_len):
    dssp_dict = make_dssp_dict(dssp_file)[0]
    q3_map = {'H': 0, 'E': 1}
    labels = []
    for i in range(1, seq_len + 1):
        ss = dssp_dict.get((chain_id, i), ('', ' '))[1]
        labels.append(q3_map.get(ss, 2))
    return torch.tensor(labels, dtype=torch.long)

def process_single_pt(pt_path, pdb_dir, dssp_dir):
    base = os.path.splitext(os.path.basename(pt_path))[0]
    chain_id = get_chain_id_from_filename(pt_path)

    tmp_pdb = os.path.join(pdb_dir, f"{base}.pdb")
    tmp_dssp = os.path.join(dssp_dir, f"{base}.dssp")

    try:
        pt_to_pdb(pt_path, tmp_pdb)
        run_dssp(tmp_pdb, tmp_dssp)

        data = torch.load(pt_path)
        seq_len = len(data['seq'])
        ss_label = extract_q3_from_dssp(tmp_dssp, chain_id, seq_len)

        assert len(ss_label) == len(data['seq']), "Length mismatch between ss_label and seq"

        data['ss_label'] = ss_label
        torch.save(data, pt_path)
        print(f"OK: {base}")
    except Exception as e:
        print(f"Error: {base} -> {e}")
    finally:
        if os.path.exists(tmp_pdb): os.remove(tmp_pdb)
        if os.path.exists(tmp_dssp): os.remove(tmp_dssp)

def batch_process(root_dir, pdb_dir, dssp_dir):
    os.makedirs(pdb_dir, exist_ok=True)
    os.makedirs(dssp_dir, exist_ok=True)

    # pt_files = [f for f in os.listdir(pt_dir) if f.endswith(".pt") and "_" in f]
    # for f in sorted(pt_files):
    #     pt_path = os.path.join(pt_dir, f)
    #     process_single_pt(pt_path, pdb_dir, dssp_dir)

    for subdir, _, files in os.walk(root_dir):
        pt_files = [f for f in files if f.endswith(".pt") and "_" in f]
        for f in sorted(pt_files):
            pt_path = os.path.join(subdir, f)
            process_single_pt(pt_path, pdb_dir, dssp_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch process .pt files to add ss_label")
    parser.add_argument("--pt_dir", type=str, required=True, help= "input .pt file directory")
    parser.add_argument("--pdb_dir", type=str, required=True, help=".pdb output file directory")
    parser.add_argument("--dssp_dir", type=str, required=True, help=".dssp file directory")
    args = parser.parse_args()

    batch_process(args.pt_dir, args.pdb_dir, args.dssp_dir)