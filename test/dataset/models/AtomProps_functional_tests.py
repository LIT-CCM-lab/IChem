import subprocess
from pathlib import Path
import re

BASE_DIRECTORY = Path.cwd()
ATOMPROPS_TEST_DIRECTORY = BASE_DIRECTORY / "atomprops-test-data"
ICHEM_COMMAND = BASE_DIRECTORY / "../../../build/IChem"

# Expected properties for each molecule
EXPECTED_PROPERTIES = {
    "acetamide.mol2": {
        "atoms": [
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": 0},
            {"name": "N", "type": "N.am", "props": "", "charge": 0},
            {"name": "H1", "type": "H", "props": "", "charge": 0},
            {"name": "H2", "type": "H", "props": "", "charge": 0},
            {"name": "C3", "type": "C.3", "props": "Hydrophobic|Apolar", "charge": 0},
            {"name": "H31", "type": "H", "props": "", "charge": 0},
            {"name": "H32", "type": "H", "props": "", "charge": 0},
            {"name": "H33", "type": "H", "props": "", "charge": 0},
            {"name": "O4", "type": "O.2", "props": "Acceptor", "charge": 0}
        ]
    },
    "acetic_acid.mol2": {
        "atoms": [
            {"name": "C", "type": "C.3", "props": "Hydrophobic|Apolar", "charge": 0},
            {"name": "H1", "type": "H", "props": "", "charge": 0},
            {"name": "H2", "type": "H", "props": "", "charge": 0},
            {"name": "C3", "type": "C.2", "props": "Apolar", "charge": 0},
            {"name": "O31", "type": "O.3", "props": "Acceptor", "charge": 0},
            {"name": "H31", "type": "H", "props": "", "charge": 0},
            {"name": "O32", "type": "O.2", "props": "Acceptor", "charge": 0},
            {"name": "H4", "type": "H", "props": "", "charge": 0}
        ]
    },
    "acetone.mol2": {
        "atoms": [
            {"name": "C", "type": "C.3", "props": "Hydrophobic|Apolar", "charge": 0},
            {"name": "H1", "type": "H", "props": "", "charge": 0},
            {"name": "C2", "type": "C.2", "props": "Apolar", "charge": 0},
            {"name": "C2", "type": "C.3", "props": "Hydrophobic|Apolar", "charge": 0},
            {"name": "H21", "type": "H", "props": "", "charge": 0},
            {"name": "H22", "type": "H", "props": "", "charge": 0},
            {"name": "H23", "type": "H", "props": "", "charge": 0},
            {"name": "H3", "type": "H", "props": "", "charge": 0},
            {"name": "H4", "type": "H", "props": "", "charge": 0},
            {"name": "O", "type": "O.2", "props": "Acceptor", "charge": 0}
        ]
    },
    "benzene.mol2": {
        "atoms": [
            {"name": "C", "type": "C.ar", "props": "Hydrophobic|Apolar|Aromatic", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "C", "type": "C.ar", "props": "Hydrophobic|Apolar|Aromatic", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "C", "type": "C.ar", "props": "Hydrophobic|Apolar|Aromatic", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "C", "type": "C.ar", "props": "Hydrophobic|Apolar|Aromatic", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "C", "type": "C.ar", "props": "Hydrophobic|Apolar|Aromatic", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "C", "type": "C.ar", "props": "Hydrophobic|Apolar|Aromatic", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0}
        ]
    },
    "furan.mol2": {
        "atoms": [
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": 0},
            {"name": "O", "type": "O.3", "props": "Acceptor", "charge": 0},
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0},
            {"name": "H", "type": "H", "props": "", "charge": 0}
        ]
    },
    "methane.mol2": {
        "atoms": [
            {"name": "C", "type": "C.3", "props": "Hydrophobic|Apolar", "charge": 0},
            {"name": "H2", "type": "H", "props": "", "charge": 0},
            {"name": "H3", "type": "H", "props": "", "charge": 0},
            {"name": "H1", "type": "H", "props": "", "charge": 0},
            {"name": "H4", "type": "H", "props": "", "charge": 0}
        ]
    },
    "methanesulfonamide.mol2": {
        "atoms": [
            {"name": "C", "type": "C.3", "props": "Hydrophobic|Apolar", "charge": 0},
            {"name": "H1", "type": "H", "props": "", "charge": 0},
            {"name": "H3", "type": "H", "props": "", "charge": 0},
            {"name": "S2", "type": "S.o2", "props": "", "charge": 0},
            {"name": "N2", "type": "N.am", "props": "", "charge": 0},
            {"name": "H21", "type": "H", "props": "", "charge": 0},
            {"name": "H22", "type": "H", "props": "", "charge": 0},
            {"name": "H4", "type": "H", "props": "", "charge": 0},
            {"name": "O", "type": "O.2", "props": "Acceptor", "charge": 0},
            {"name": "O", "type": "O.2", "props": "Acceptor", "charge": 0}
        ]
    },
    "methanol.mol2": {
        "atoms": [
            {"name": "C", "type": "C.3", "props": "Hydrophobic|Apolar", "charge": 0},
            {"name": "O1", "type": "O.3", "props": "Acceptor", "charge": 0},
            {"name": "H1", "type": "H", "props": "", "charge": 0},
            {"name": "H2", "type": "H", "props": "", "charge": 0},
            {"name": "H3", "type": "H", "props": "", "charge": 0},
            {"name": "H4", "type": "H", "props": "", "charge": 0}
        ]
    },
    "methylammonium.mol2": {
        "atoms": [
            {"name": "C", "type": "C.3", "props": "Apolar", "charge": 0},
            {"name": "H1", "type": "H", "props": "", "charge": 0},
            {"name": "H2", "type": "H", "props": "", "charge": 0},
            {"name": "N3", "type": "N.4", "props": "Cation", "charge": 1},
            {"name": "H31", "type": "H", "props": "", "charge": 0},
            {"name": "H32", "type": "H", "props": "", "charge": 0},
            {"name": "H33", "type": "H", "props": "", "charge": 0},
            {"name": "H4", "type": "H", "props": "", "charge": 0}
        ]
    },
    "pyrrole.mol2": {
        "atoms": [
             {"name": "C", "type": "C.2", "props": "Apolar", "charge": -0.1500},
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": -0.1500},
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": -0.3016},
            {"name": "N", "type": "N.pl3", "props": "", "charge": 0.0332},
            {"name": "C", "type": "C.2", "props": "Apolar", "charge": -0.3016},
            {"name": "H", "type": "H", "props": "", "charge": 0.1500},
            {"name": "H", "type": "H", "props": "", "charge": 0.1500},
            {"name": "H", "type": "H", "props": "", "charge": 0.2700},
            {"name": "H", "type": "H", "props": "", "charge": 0.1500}
        ]
    }   
}

# Run AtomProps command for a single file
def run_atomprops_command(mol_file):
    command = f"{ICHEM_COMMAND} AtomProps {mol_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"IChem command failed for {mol_file.name}: {result.stderr}")
    return result.stdout

#Parses the AtomProps output to extract atom details
def parse_atomprops_output(output):
    # Capture the result from the command AtomProps 
    atom_pattern = r"ATOM\s+(\S+)\s+\|\d+\s+\|\S+\s+\d+-\S+\s*\n\s*type:\s+(\S+)\s*\n\s*props\s*:\s*(.*?)\s*\n\s*charge:\s+([-+]?\d*\.\d+|\d+)"
    matches = re.findall(atom_pattern, output)
    parsed_atoms = []
    for match in matches:
        atom_name, atom_type, atom_props, atom_charge = match
        parsed_atoms.append({
            "name": atom_name,
            "type": atom_type.strip(),
            "props": atom_props.strip(),
            "charge": float(atom_charge)  # Convert to float for precision
        })
    print(f"Parsed atoms from output:\n{parsed_atoms}")
    return parsed_atoms

# Verify parsed output against expectations
# Compares the parsed atom properties with the expected properties
def verify_atom_properties(parsed_atoms, expected_atoms):
    for expected_atom in expected_atoms:
        # Find all matching atoms by name
        print(f"Expected atom: {expected_atom}\n")
        # print(f"Parsed atoms: {parsed_atoms}")
        matching_atoms = [atom for atom in parsed_atoms if atom["name"] == expected_atom["name"]]
        # print(f"Matching atoms for {expected_atom['name']}: {matching_atoms}")
        if not matching_atoms:
            return False, f"Missing atom {expected_atom['name']}"

        # matching atom: satisfies all properties
        for actual_atom in matching_atoms:
            if all(actual_atom[key] == expected_atom[key] for key in ["type", "props", "charge"]):
                break
        else:  # Do not match
            return False, (
                f"Atom {expected_atom['name']} does not match: "
                f"expected {expected_atom}, got {matching_atoms}"
            )

    return True, "All properties match"

# Process a single file test in atomprops-test-data directory and compares the output with expected values
def process_single_file_test(mol_file):
    try:
        output = run_atomprops_command(mol_file)
        if not output.strip():
            return {
                "file": mol_file.name,
                "status": "Failed",
                "message": "No output produced for this file"
            }
    except Exception as e:
        return {
            "file": mol_file.name,
            "status": "Error",
            "message": f"Command execution failed: {e}"
        }

    # Parse the output
    parsed_atoms = parse_atomprops_output(output)

    # Verify against expected properties
    expected_data = EXPECTED_PROPERTIES.get(mol_file.name)
    if not expected_data:
        return {
            "file": mol_file.name,
            "status": "Error",
            "message": "No expected data defined for this molecule"
        }

    is_valid, message = verify_atom_properties(parsed_atoms, expected_data["atoms"])

    return {
        "file": mol_file.name,
        "status": "Passed" if is_valid else "Failed",
        "message": message
    }

# Run AtomProps tests for all files in atomprops-test-data
def run_all_atomprops_tests():
    if not ATOMPROPS_TEST_DIRECTORY.exists() or not ATOMPROPS_TEST_DIRECTORY.is_dir():
        print(f"Error: AtomProps test directory {ATOMPROPS_TEST_DIRECTORY} does not exist or is not a directory")
        return []

    mol_files = list(ATOMPROPS_TEST_DIRECTORY.glob("*.mol2"))
    if not mol_files:
        print("No MOL2 files found in atomprops-test-data")
        return [{"status": "Failed", "message": "No MOL2 files found in atomprops-test-data"}]

    results_summary = []

    for mol_file in mol_files:
        print(f"Processing test case: {mol_file.name}")
        test_result = process_single_file_test(mol_file)
        results_summary.append(test_result)

    # Summary
    passed_tests = [res for res in results_summary if res.get("status") == "Passed"]
    failed_tests = [res for res in results_summary if res.get("status") in ["Failed", "Error"]]

    print("\n\t\t\t############## Test Summary ##############\n")
    print(f"Total tests: {len(results_summary)}")
    print(f"\t=> Passed: {len(passed_tests)}")
    print(f"\t=> Failed: {len(failed_tests)}\n\n")

    for res in results_summary:
        print(f"{res['status']}: {res['file']} - {res['message']}")

    return results_summary

if __name__ == "__main__":
    run_all_atomprops_tests()