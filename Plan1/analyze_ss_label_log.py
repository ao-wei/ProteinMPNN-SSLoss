import argparse
from collections import defaultdict

def analyze_log(log_path, save_failed=False):
    total_lines = 0
    success_count = 0
    mkdssp_errors = []
    empty_errors = []
    other_errors = []

    with open(log_path, 'r') as f:
        lines = f.readlines()

    total_lines = len(lines)

    for i, line in enumerate(lines):
        line = line.strip()

        if line.startswith("OK:"):
            success_count += 1

        elif "mkdssp" in line and line.startswith("Error:"):
            entry = line.split(":")[1].strip().split("->")[0].strip()
            mkdssp_errors.append(entry)

        elif "Empty file" in line:
            prev_line = lines[i-1].strip()
            if prev_line.startswith("Error:"):
                entry = prev_line.split(":")[1].strip().split("->")[0].strip()
                empty_errors.append(entry)

        elif line.startswith("Error:"):
            entry = line.split(":")[1].strip().split("->")[0].strip()
            other_errors.append(entry)

    print(f"\n=== Log Summary ===")
    print(f"Total lines in log:      {total_lines}")
    print(f"Success count (OK):      {success_count}")
    print(f"DSSP failures (mkdssp):  {len(mkdssp_errors)}")
    print(f"Empty PDB files:        {len(empty_errors)}")
    print(f"Other errors:            {len(other_errors)}\n")

    print(f"Example mkdssp errors:   {mkdssp_errors[:5]}")
    print(f"Example empty errors:    {empty_errors[:5]}")
    print(f"Example other errors:    {other_errors[:5]}")

    if save_failed:
        with open("failed_mkdssp.txt", "w") as f1:
            f1.writelines([f"{x}\n" for x in mkdssp_errors])
        with open("failed_empty.txt", "w") as f2:
            f2.writelines([f"{x}\n" for x in empty_errors])
        with open("failed_other.txt", "w") as f3:
            f3.writelines([f"{x}\n" for x in other_errors])
        print("\nSaved failed sample lists: failed_mkdssp.txt, failed_empty.txt, failed_other.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=str, default="log.txt", help="Path to log.txt file")
    parser.add_argument("--save_failed", action="store_true", help="Save failed sample names to file")
    args = parser.parse_args()

    analyze_log(args.log, save_failed=args.save_failed)
