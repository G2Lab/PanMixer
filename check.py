import os
import glob
import sys
from constants import BASE_PATH

error_keywords = ["error", "failed", "exception", "traceback", "fatal", "abort", "unable", "denied", "core dumped"]

def find_most_recent_directory(parent_dir, name_of_job = None):
    """Finds the most recently created directory in the given parent directory."""
    subdirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]

    if name_of_job:
        subdirs = [d for d in subdirs if name_of_job in d]
    if not subdirs:
        return None
    return max(subdirs, key=os.path.getctime)

def check_for_errors(directory):
    """Checks for error files in the directory and prints relevant outputs."""
    error_files = sorted(glob.glob(os.path.join(directory, "error_*")))

    directory_job_name = directory.split("/")[-1].split("_")[2:]
    directory_job_name = "_".join(directory_job_name)
    print(f"Checking for errors in recently run {directory_job_name}...")
    
    for error_file in error_files:
        with open(error_file, 'r') as ef:
            lines = ef.readlines()
        
        for line in lines:
            if any(keyword in line.lower() for keyword in error_keywords):
                print(f"Error found! In file: {os.path.basename(error_file)}")
                
                # Find the corresponding out_* file
                print("------ Error ------")
                print("".join(lines))
                print("------ Output ------")
                out_file = error_file.replace("error_", "out_", 1)
                if os.path.exists(out_file):
                    with open(out_file, 'r') as of:
                        print(of.read())

                return  # Stop after finding the first error
    print("No errors found.")
    
if __name__ == "__main__":
    parent_directory = f"{BASE_PATH}/slurm"  # Change this to the directory containing job folders
    
    if sys.argv[1:]:
        name_of_job = sys.argv[1]
    else:
        name_of_job = None
    recent_directory = find_most_recent_directory(parent_directory, name_of_job = name_of_job)
    
    if recent_directory:
        check_for_errors(recent_directory)
    else:
        print("No directories found.")
