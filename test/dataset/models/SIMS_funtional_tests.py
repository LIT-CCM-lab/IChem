import subprocess
import re
from pathlib import Path

BASE_DIRECTORY = Path.cwd()
SIMS_DIRECTORY = BASE_DIRECTORY / "sims-test-data"  # Directory containing functional test cases
ICHEM_COMMAND = BASE_DIRECTORY / "../../../build/IChem"

# Read the bitstring from a file
# Should be in the format of: namefile\tbitstring
def read_bit_string(file_path):
    with open(file_path, "r") as f:
        line = f.readline().strip()  # Read the first line
        parts = line.split("\t")  # Split into namefile and bitstring
        if len(parts) != 2:
            raise ValueError(f"Invalid file format in {file_path}. Expected 'namefile\\tbitstring'")
        return parts[1]  # Return the bitstring

# Parse result.txt for the similarity values
# Parse result.txt for the format 'IF1\tIF2\tSimilarity
def parse_result_file(result_file):
    results = []
    if not result_file.exists():
        print(f"Result file {result_file} not found")
        return results
    
    with open(result_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                print(f"Invalid line in {result_file}: {line}")
                continue
            results.append((parts[0], parts[1], float(parts[2])))  # IF1, IF2, Similarity
    return results

# Calculate Tanimoto similarity between two binary fingerprints
def tanimoto_similarity(ifp1, ifp2):
    if len(ifp1) != len(ifp2):
        raise ValueError("Bitstrings must have the same length")

    intersection = sum(1 for a, b in zip(ifp1, ifp2) if a == b == "1")
    union = sum(1 for a, b in zip(ifp1, ifp2) if a == "1" or b == "1")

    return intersection / union if union > 0 else 0

# Run the similarity command
# Generate result.txt
def run_similarity_command(file_a, file_b, result_file):
    command = f"{ICHEM_COMMAND} sims {file_a} {file_b} >> {result_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"IChem command failed: {result.stderr}")


# Process a single test directory and compare results
# Process a single test directory to compare computed and command-generated similarities.
def process_test_directory(case_dir):
    if not case_dir.exists() or not case_dir.is_dir():
        print(f"Test directory {case_dir} does not exist")
        return []

    result_file = case_dir / "result.txt"
    result_file.unlink(missing_ok=True)  # Remove previous result file if it exists
    
    file_pairs = list(case_dir.glob("*.txt"))
    if len(file_pairs) != 2:
        print(f"Test directory {case_dir.name} must contain exactly 2 .txt files")
        return []

    file_a, file_b = file_pairs

    # Extract bitstrings from file_a and file_b
    try:
        bitstring_a = read_bit_string(file_a)
        bitstring_b = read_bit_string(file_b)
    except Exception as e:
        print(f"Error reading bitstrings in {case_dir.name}: {e}")
        return [{
            "test_case": case_dir.name,
            "status": "Error",
            "message": f"Failed to read bitstrings: {e}"
        }]

    # Run the similarity command
    try:
        run_similarity_command(file_a, file_b, result_file)
    except Exception as e:
        print(f"Error running similarity command for {case_dir.name}: {e}")
        return [{
            "test_case": case_dir.name,
            "status": "Error",
            "message": f"Command execution failed: {e}"
        }]

    # Parse result.txt
    if not result_file.exists() or result_file.stat().st_size == 0:
        return [{
            "test_case": case_dir.name,
            "status": "Failed",
            "message": f"No result written in {result_file}."
        }]
    
    parsed_results = parse_result_file(result_file)
    validation_results = []

    # Compare each parsed result with the computed Tanimoto similarity
    for if1, if2, similarity in parsed_results:
        try:
            computed_similarity = tanimoto_similarity(bitstring_a, bitstring_b)
            status = "Passed" if abs(similarity - computed_similarity) < 1e-6 else "Failed"
            validation_results.append({
                "if1": if1,
                "if2": if2,
                "parsed_similarity": similarity,
                "computed_similarity": computed_similarity,
                "status": status
            })
        except Exception as e:
            validation_results.append({
                "if1": if1,
                "if2": if2,
                "parsed_similarity": similarity,
                "computed_similarity": None,
                "status": "Error",
                "message": str(e)
            })

    return validation_results

# Run similarity tests for all directories dynamically detected in sims-test-data
def run_all_tests():
    if not SIMS_DIRECTORY.exists() or not SIMS_DIRECTORY.is_dir():
        print(f"Error: Sims test directory {SIMS_DIRECTORY} does not exist or is not a directory")
        return []

    results_summary = []
    test_directories = [d for d in SIMS_DIRECTORY.iterdir() if d.is_dir()]
    
    if not test_directories:
        print(f"No test directories found in {SIMS_DIRECTORY}.")
        return results_summary

    for case_dir in test_directories:
        print(f"Processing test case: {case_dir.name}")
        try:
            test_results = process_test_directory(case_dir)
            results_summary.extend(test_results)
        except Exception as e:
            print(f"Error processing {case_dir.name}: {e}")
            results_summary.append({
                "test_case": case_dir.name,
                "status": "Error",
                "message": str(e)
            })
    
    # Summary
    passed_tests = [res for res in results_summary if res.get("status") == "Passed"]
    failed_tests = [res for res in results_summary if res.get("status") == "Failed"]

    print("\n\t\t\t############## Test Summary ##############\n")
    print(f"Total tests: {len(results_summary)}\n")
    print(f"\t=> Passed: {len(passed_tests)}")
    print(f"\t=> Failed: {len(failed_tests)}\n\n")
    print("Details:\n")
    
    for res in results_summary:
        print(f"\t=> {res['status']} - Parsed: {res['parsed_similarity']}, Computed: {res['computed_similarity']}")

    return results_summary

if __name__ == "__main__":
    run_all_tests()