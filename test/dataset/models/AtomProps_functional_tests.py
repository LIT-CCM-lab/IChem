import subprocess
from pathlib import Path
import re

BASE_DIRECTORY = Path.cwd()
ATOMPROPS_TEST_DIRECTORY = BASE_DIRECTORY / "atomprops-test-data"
ICHEM_COMMAND = BASE_DIRECTORY / "../../../build/IChem"

# Run AtomProps command for a single file
def run_atomprops_command(mol_file):
    command = f"{ICHEM_COMMAND} AtomProps {mol_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"IChem command failed for {mol_file.name}: {result.stderr}")
    return result.stdout


def count_molecules_in_output(output, mol_file):
    # get full path
    full_path = str(mol_file.resolve())
    # We are looking for the line "Ligand..." ; the start of the molecule
    pattern = rf"Molecule \d+ with name \S+ from file: {re.escape(full_path)}"
    matches = re.findall(pattern, output)
    print(f"Tell me you found a match {matches}")
    return len(matches)

# Process a single file test in atomprops-test-data directory
def process_single_file_test(mol_file):
    try:
        output = run_atomprops_command(mol_file)
        if not output.strip():  # Output must not be empty
            return {
                "file": mol_file.name,
                "status": "Failed",
                "message": "No output produced for this file",
                "molecule_count": 0
            }
    except Exception as e:
        return {
            "file": mol_file.name,
            "status": "Error",
            "message": f"Command execution failed: {e}",
            "molecule_count": 0
        }

    # Count the molecules 
    # print(output)
    molecule_count = count_molecules_in_output(output, mol_file)

    # At least one molecule so the test passes
    if molecule_count > 0:
        return {
            "file": mol_file.name,
            "status": "Passed",
            "message": f"Processed {molecule_count} molecule(s)",
            "molecule_count": molecule_count
        }
    else:
        return {
            "file": mol_file.name,
            "status": "Failed",
            "message": "No molecules found in output",
            "molecule_count": 0
        }

# Run AtomProps tests for all files in atomprops-test-data
def run_all_atomprops_tests():
    if not ATOMPROPS_TEST_DIRECTORY.exists() or not ATOMPROPS_TEST_DIRECTORY.is_dir():
        print(f"Error: AtomProps test directory {ATOMPROPS_TEST_DIRECTORY} does not exist or is not a directory")
        return []

    mol_files = list(ATOMPROPS_TEST_DIRECTORY.glob("*.mol2"))
    if not mol_files:
        print("No MOL2 files found in atomprops-test-data")
        return [{
            "status": "Failed",
            "message": "No MOL2 files found in atomprops-test-data"
        }]

    results_summary = []
    total_molecules = 0  # Molecules processed

    for mol_file in mol_files:
        print(f"Processing test case: {mol_file.name}")
        test_result = process_single_file_test(mol_file)
        total_molecules += test_result.get("molecule_count", 0)
        results_summary.append(test_result)

    # Summary
    passed_tests = [res for res in results_summary if res.get("status") == "Passed"]
    failed_tests = [res for res in results_summary if res.get("status") in ["Failed", "Error"]]

    print("\n\t\t\t############## Test Summary ##############\n")
    print(f"Total tests: {len(results_summary)}")
    print(f"Total molecules processed: {total_molecules}")
    print(f"\t=> Passed: {len(passed_tests)}")
    print(f"\t=> Failed: {len(failed_tests)}\n\n")

    for res in results_summary:
        print(f"{res['status']}: {res['file']} - {res['message']} (Molecules: {res['molecule_count']})")

    return results_summary

if __name__ == "__main__":
    run_all_atomprops_tests()
