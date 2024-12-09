import subprocess
import re
from pathlib import Path

BASE_DIRECTORY = Path.cwd()

ICHEM_COMMAND = BASE_DIRECTORY / "../../../build/IChem"  

BIT_POSITION_MAP = {
    0: "hydrophobic",
    1: "aromatic_face_to_face",
    2: "aromatic_edge_to_face",
    3: "h_bond_protein_donor",
    4: "h_bond_ligand_donor",
    5: "ionic_protein_charged_plus",
    6: "ionic_ligand_charged_plus",
    7: "pication",
    8: "metal"
}

# Find mol2 files in interaction subdir 
def find_mol2_files(base_directory):
    test_files = {}
    for subdirectory in base_directory.iterdir():
        if subdirectory.is_dir():
            mol2_files = sorted([f for f in subdirectory.glob("*.mol2")])
            files = [(mol2_files[i], mol2_files[i + 1]) for i in range(0, len(mol2_files) - 1, 2)]
            test_files[subdirectory.name] = files
    return test_files


# Run IChem command
def run_chem_command(file_a, file_b, interaction_type):
    # print(f"In run_chem_command(), displaying file_a: {file_a} and file_b: {file_b} and interaction_type: {interaction_type}")
    if interaction_type == "metal" or interaction_type == "pication":
        # print("I am in metal")
        command = f"{ICHEM_COMMAND} --extended IFP  {file_a} {file_b}"
    # elif interaction_type == "pication":
    #     command = f"{ICHEM_COMMAND} --extended IFP {file_a} {file_b}"
    else:
        command = f"{ICHEM_COMMAND} IFP {file_a} {file_b}"

    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"IChem command failed: {result.stderr}")
    # print(f"Command output for {file_a} and {file_b} with interaction {interaction_type}:\n{result.stdout}")
    return result.stdout


# Parse result for either 7-bit or 9-bit string
def parse_result(output, extended=False):
    # print(f"Output defined when entering parse_result: {output}")
    bit_length_string = 9 if extended else 7
    # rf for raw f-strings
    match = re.search(rf'[01]{{{bit_length_string}}}', output)  # Match either 7 or 9 bits based on --extended option
    if match:
        return match.group(0)
    else:
        print(f"Failed to parse {bit_length_string}-bit string in output:\n{output}")
    return None


# Identify interactions based on bit positions
def detect_interactions(result_string):
    detected_interactions = []
    for position, interaction_type in BIT_POSITION_MAP.items():
        if position < len(result_string) and result_string[position] == "1":
            detected_interactions.append(interaction_type)
    # print(f"FINAL INTERACTIONS {detected_interactions}")
    return detected_interactions # Return an array of interactions: exp ['hydrophobic', 'aromatic_edge_to_face']


# Log results for each case
def process_test_results(file_a, file_b, result_string, interaction_type):
    if result_string is None:
        # print(f"No valid bit string for files {file_a.name} and {file_b.name}")
        return {
            "interaction_type": interaction_type,
            "file_a": file_a.name,
            "file_b": file_b.name,
            "result_string": None,
            "interactions_found": [],
            "status": "Failed",
            "message": f"No interactions found for the test case {interaction_type} with file_a: ({file_a.name} and file_b:{file_b.name})"
        }
    
    # print(f"The result string before the interactions get detected {result_string}")
    detected_interactions = detect_interactions(result_string)
    # print(f"Interactions detected {detected_interactions}") 
    status = "Passed" if detected_interactions else "Failed" # Status defined
    message = f"Interactions found for the test case {interaction_type}: {detected_interactions} with the bit string: {result_string}" if detected_interactions else f"No interactions found for the test case {interaction_type} with file_a:({file_a.name} and  file_b:{file_b.name})"
    
    # print(f"{status} - Result for files {file_a.name} and {file_b.name}: {result_string} - {message}")
    
    return {
        "interaction_type": interaction_type,
        "file_a": file_a.name,
        "file_b": file_b.name,
        "result_string": result_string,
        "interactions_found": detected_interactions,
        "status": status,
        "message": message
    }

# Run all the tests
def run_all_tests():
    results_summary = []
    test_files = find_mol2_files(BASE_DIRECTORY)
    
    for test_dir, file_pairs in test_files.items():
        for file_a, file_b in file_pairs:
            interaction_type = test_dir  # Interaction type: name of the file containing the test
            extended = interaction_type in ["metal", "pication"] # Extend parsing for those two cases 
            # print(f"Displaying extended value {extended}")
            
            try:
                output = run_chem_command(file_a, file_b, interaction_type)
                # print(f"Output result in run_all_tests {output}")
                result_string = parse_result(output, extended=extended)
                test_result = process_test_results(file_a, file_b, result_string, interaction_type)
                results_summary.append(test_result)
                
            except Exception as e:
                error_message = f"Error running test with files {file_a.name} and {file_b.name}: {e}"
                print(error_message)
                results_summary.append({
                    "interaction_type": interaction_type,
                    "file_a": file_a.name,
                    "file_b": file_b.name,
                    "error": str(e),
                    "status": "Failed",
                    "message": error_message
                })

    # Summary
    passed_tests = [res for res in results_summary if res["status"] == "Passed"]
    failed_tests = [res for res in results_summary if res["status"] == "Failed"]

    print("\n\t\t\t############## Test Summary ##############\n")
    print(f"Total tests: {len(results_summary)}\n")
    print(f"\t=> Passed: {len(passed_tests)}")
    print(f"\t=> Failed: {len(failed_tests)}\n\n")
    print("Details:\n")
    
    for res in results_summary:
        # print(results_summary)
        print(f"\t=> {res['status']} - {res['message']}")

    print("\n")
    return results_summary

# Run all tests
if __name__ == "__main__":
    run_all_tests()
    
    
# if __name__ == "__main__":
    
#     results_summary = []
#     test_files = find_mol2_files(BASE_DIRECTORY)
#     for test_dir, file_pairs in test_files.items():
#         for file_a, file_b in file_pairs:
#             interaction_type = test_dir 
#             print(interaction_type)
#             # extended = interaction_type in ["metal", "pication"]
            
#             output = run_chem_command(file_a, file_b, interaction_type)
#             # print(output)
#             result_string = parse_result(output, extended=False)
#             # print(result_string) # The bit string 
#             test_result = process_test_results(file_a, file_b, result_string)
#             # print(test_result)
#             # results_summary.append(test_result)
            
            # print(results_summary)
            # passed_tests = [res for res in results_summary if res["status"] == "Passed"]
            # # print(f" Passed: {passed_tests}")
            # failed_tests = [res for res in results_summary if res["status"] == "Failed"]
            # print(f" Passed: {failed_tests}")

    
