import os
import subprocess
import torch
import numpy as np
import glob


ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'

def find_netsurfp_txt(out_prefix):
    """
    Find the .netsurfp.txt file in the output directory.
    """
    search_path = os.path.join(out_prefix, "**", "*.netsurfp.txt")
    matches = glob.glob(search_path, recursive=True)
    if len(matches) == 0:
        raise FileNotFoundError(f"[NetSurfP] No .netsurfp.txt found under: {out_prefix}")
    if len(matches) > 1:
        print(f"[NetSurfP] Multiple .netsurfp.txt found, returning first.")
    print(f"[DEBUG] Found .netsurfp.txt: {matches[0]}")
    return matches[0]

def write_batch_fasta(S_sample, save_dir="./netsurfp_tmp", batch_index=0):
    """
    Write a batch of sequences to a FASTA file.
    """
    os.makedirs(save_dir, exist_ok=True)
    fasta_paths = []

    for i, s in enumerate(S_sample):
        seq = "".join([ALPHABET[idx] for idx in s.cpu().tolist()])
        fasta_path = os.path.join(save_dir, f"sample_{batch_index}_{i}.fasta")
        with open(fasta_path, 'w') as f:
            f.write(f">sample_{batch_index}_{i}\n{seq}\n")
        fasta_paths.append(fasta_path)

    return fasta_paths

def run_netsurfp(fasta_path):
    """
    Run NetSurfP on a given FASTA file and return the path to the output file.
    """
    base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    out_dir = os.path.dirname(fasta_path)
    out_prefix = os.path.join(out_dir, base_name)

    script_path = "/home/aowei/protein_ss_evaluation/tools/NetSurfP-3.0_standalone/nsp3.py"
    model_path = "/home/aowei/protein_ss_evaluation/tools/NetSurfP-3.0_standalone/models/nsp3.pth"

    cmd = f"python3 {script_path} -m {model_path} -i {fasta_path} -o {out_prefix} > {out_prefix}.log 2>&1"
    subprocess.run(cmd, shell=True)

    # for ext in [".json", ".csv", ".fasta", ".log"]:
    #     for file in glob.glob(os.path.join(out_prefix, f"**/*{ext}"), recursive=True):
    #         try:
    #             os.remove(file)
    #         except Exception as e:
    #             print(f"[NetSurfP][WARN] Failed to delete {file}: {e}")

    netsurfp_dir = os.path.join(out_prefix, "01", f"0000_{base_name}")
    netsurfp_file = os.path.join(netsurfp_dir, f"0000_{base_name}.netsurfp.txt")

    print(f"[DEBUG] NetSurfP output file: {netsurfp_file}")
    # ss3_path = find_netsurfp_txt(out_prefix)
    return netsurfp_file

def parse_netsurfp_q3(netsurfp_path):
    """
    Parse the NetSurfP output file to extract the Q3 secondary structure predictions.
    """
    # ss_pred = []
    probs = []
    with open(netsurfp_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) >= 10:
                prob_H = float(parts[7])
                prob_E = float(parts[8])
                prob_C = float(parts[9])
                probs.append([prob_H, prob_E, prob_C])
                # q3_index = int(np.argmax([prob_H, prob_E, prob_C]))
                # ss_pred.append(q3_index)
    return torch.tensor(probs, dtype=torch.float32)
    # return torch.tensor(ss_pred, dtype=torch.long)

def run_netsurfp_and_parse(fasta_path):
    """
    Run NetSurfP on a given FASTA file and parse the output to get secondary structure predictions.
    """
    ss3_path = run_netsurfp(fasta_path)
    ss_pred = parse_netsurfp_q3(ss3_path)
    return ss_pred

