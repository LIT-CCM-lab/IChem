import subprocess
from pathlib import Path
import numpy as np
import csv
import os


BASE_DIRECTORY = Path.cwd()

ICHEM_COMMAND = BASE_DIRECTORY / "../../../build/IChem"  


# GRIM Score calculation
def calculate_grsc(n_lig, n_center, n_prot, sum_cl, rmsd, diff_i):
    grsc = (
        0.5006
        + 0.0151 * n_lig
        + 0.0039 * n_center
        + 0.0143 * n_prot
        + 0.2098 * sum_cl
        - 0.0720 * rmsd
        - 0.0003 * diff_i
    )
    return round(grsc, 4)

#  Calculate DiffI
def calculate_diffi(r_lig, r_center, r_prot, c_lig, c_center, c_prot):
    ref_total = r_lig + r_center + r_prot
    comp_total = c_lig + c_center + c_prot
    return abs(ref_total - comp_total)

# Generate Grifp_res.csv using the IChem command
def run_grim_ichem_command(ref_file, comp_file, output_file="Grifp_res.csv"):
    command = f"{ICHEM_COMMAND} grim {ref_file} {comp_file}"
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"IChem command failed: {result.stderr}")
        print(f"Generated {output_file} successfully.")
        return output_file
    except Exception as e:
        raise RuntimeError(f"Error running IChem command: {e}")

# Parse the file and compute DiffI and GrSc
def process_file(file_path):
    results = []
    with open(file_path, mode="r") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            # Parse input values
            n_lig = int(row["LIG"])
            n_center = int(row["CENTER"])
            n_prot = int(row["PROT"])
            sum_cl = float(row["SumCl"])
            rmsd = float(row["RMSD"])
            r_lig = int(row["RLig"])
            r_center = int(row["RCent"])
            r_prot = int(row["RProt"])
            c_lig = int(row["CLig"])
            c_center = int(row["CCent"])
            c_prot = int(row["CProt"])

            # Calculate DiffI
            diff_i = calculate_diffi(r_lig, r_center, r_prot, c_lig, c_center, c_prot)

            # Calculate GrSc
            grsc = calculate_grsc(n_lig, n_center, n_prot, sum_cl, rmsd, diff_i)

            # Append results for the row
            results.append({
                "NCli": row["NCli"],
                "Ref": row["Ref"],
                "Comp": row["Comp"],
                "DiffI": diff_i,
                "GrSc": grsc,
                "ExpectedGrSc": float(row["GrSc"]),
            })
    return results

# # Calculate RMSD Alignment from matched pseudoatom coordinates
# def calculate_rmsdal(ref_coords, comp_coords):
#     ref_coords = np.array(ref_coords)
#     comp_coords = np.array(comp_coords)

#     # squared differences
#     squared_diffs = np.sum((ref_coords - comp_coords) ** 2, axis=1)

#     # mean squared difference and RMSD
#     mean_squared_diff = np.mean(squared_diffs)
#     rmsdal = np.sqrt(mean_squared_diff)
#     return round(rmsdal, 4) # Round to 4


## Compare calculated and expected GrSc
def compare_results(results):
    test_summary = []
    for result in results:
        calculated = result["GrSc"]
        expected = result["ExpectedGrSc"]
        status = "Passed" if abs(calculated - expected) < 1e-4 else "Failed"
        
        test_summary.append({
            "NCli": result["NCli"],
            "Ref": result["Ref"],
            "Comp": result["Comp"],
            "DiffI": result["DiffI"],
            "CalculatedGrSc": calculated,
            "ExpectedGrSc": expected,
            "Status": status
        })
    
    return test_summary


# Display results
def display_results(test_summary):
    print("\nTest Summary\n")
    total_tests = len(test_summary)
    passed_tests = [t for t in test_summary if t["Status"] == "Passed"]
    failed_tests = [t for t in test_summary if t["Status"] == "Failed"]

    print(f"Total tests: {total_tests}")
    print(f"Passed: {len(passed_tests)}")
    print(f"Failed: {len(failed_tests)}\n")

    for test in test_summary:
        print(
            f"NCli: {test['NCli']}, Ref: {test['Ref']}, Comp: {test['Comp']}, "
            f"DiffI: {test['DiffI']}, CalculatedGrSc: {test['CalculatedGrSc']}, "
            f"ExpectedGrSc: {test['ExpectedGrSc']}, Status: {test['Status']}"
        )


if __name__ == "__main__":
    ref_file = "grim-test-data/RefInts.mol2"
    comp_file = "grim-test-data/CompInts.mol2"

    try:
        # Generate Grifp_res.csv
        csv_file_path = run_grim_ichem_command(ref_file, comp_file)

        # Process and compare results
        results = process_file(csv_file_path)
        test_summary = compare_results(results)
        display_results(test_summary)
    finally:
        # Delete the generated file
        if Path(csv_file_path).exists():
            os.remove(csv_file_path)
            # print(f"Deleted file: {csv_file_path}")