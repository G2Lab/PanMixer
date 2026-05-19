import numpy as np
import datetime
import os
import subprocess

from constants import BASE_PATH, PYTHON_ENV

ROOT_TOOL_PATH = f"{BASE_PATH}/tools"

RESTRICT_RESOURCES = False

NUMBER_CHROMSOMES = 22

# Global dependency: when set, all submitted jobs will depend on this job ID
_pending_dependency = None

def set_dependency(job_id):
    """Set a dependency for the next batch of job submissions."""
    global _pending_dependency
    _pending_dependency = job_id

def clear_dependency():
    """Clear the current dependency."""
    global _pending_dependency
    _pending_dependency = None

def _submit_sbatch(script_path):
    """Submit a sbatch script and return the job ID."""
    global _pending_dependency
    cmd = ["sbatch", "--parsable"]
    if _pending_dependency is not None:
        cmd.append(f"--dependency=afterok:{_pending_dependency}")
    cmd.append(script_path)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"sbatch failed for {script_path}: {result.stderr.strip()}")
    job_id = result.stdout.strip()
    print(f"Submitted {script_path} -> Job ID: {job_id}")
    return job_id

month_number_to_name = {
    1: "Jan",
    2: "Feb",
    3: "Mar",
    4: "Apr",
    5: "May",
    6: "Jun",
    7: "Jul",
    8: "Aug",
    9: "Sep",
    10: "Oct",
    11: "Nov",
    12: "Dec"
}

def get_nyc_time_in_month_day_time():
    nyc_time = datetime.datetime.now() - datetime.timedelta(hours=4)
    month = month_number_to_name[int(nyc_time.strftime("%m"))]
    time = nyc_time.strftime(f"%d_%H.%M.%S")
    return month + time

def create_slurm_header(tool_name, path, memory = "64g", cpus = "1", num_tasks = "1"):
    if RESTRICT_RESOURCES:
        return f"""#!/bin/bash
#SBATCH --job-name={tool_name}
#SBATCH --time=100:00:00
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --out={path}/out_%j_%a.log
#SBATCH --error={path}/error_%j_%a.log
#SBATCH --array=0-{str(int(num_tasks) - 1)}
#SBATCH --nodelist=ne1dc6-001,ne1dc6-002  # Specify nodes
"""
    else:
        return f"""#!/bin/bash
#SBATCH --job-name={tool_name}
#SBATCH --time=100:00:00
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --out={path}/out_%j_%a.log
#SBATCH --error={path}/error_%j_%a.log
#SBATCH --array=0-{str(int(num_tasks) - 1)}
"""
    
def create_slurm_header_multichromosome(tool_name, path, memory = "64g", cpus = "1", num_tasks = "1", max_concurrent = None):
    array_spec = f"0-{str(int(num_tasks) - 1)}"
    if max_concurrent:
        array_spec += f"%{max_concurrent}"
    return f"""#!/bin/bash
#SBATCH --job-name={tool_name}
#SBATCH --time=24:00:00
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --out={path}/out_%j_%a.log
#SBATCH --error={path}/error_%j_%a.log
#SBATCH --array={array_spec}

ROW_ID=$((SLURM_ARRAY_TASK_ID / {NUMBER_CHROMSOMES}))
CHR=$(( (SLURM_ARRAY_TASK_ID % {NUMBER_CHROMSOMES}) + 1 ))
"""

def create_slurm_header_by_chromosome(tool_name, path, memory = "64g", cpus = "1"):
    return f"""#!/bin/bash
#SBATCH --job-name={tool_name}
#SBATCH --time=100:00:00
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --out={path}/out_%j_%a.log
#SBATCH --error={path}/error_%j_%a.log
#SBATCH --array=1-{str(NUMBER_CHROMSOMES)}
CHR=$((SLURM_ARRAY_TASK_ID))
"""

def create_slurm_header_by_rowid(tool_name, path, memory = "64g", cpus = "1", num_tasks = "1"):
    return f"""#!/bin/bash
#SBATCH --job-name={tool_name}
#SBATCH --time=100:00:00
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --out={path}/out_%j_%a.log
#SBATCH --error={path}/error_%j_%a.log
#SBATCH --array=0-{str(int(num_tasks) - 1)}
ROW_ID=$((SLURM_ARRAY_TASK_ID))
"""

PYTHON = f"{PYTHON_ENV}/bin/python3"

def _write_env_setup(file):
    """Write environment setup commands to a slurm script."""
    file.write(f"export PATH={PYTHON_ENV}/bin:$PATH\n")
    file.write(f"cd {BASE_PATH}\n")
    file.write(f"export PYTHONPATH=$PYTHONPATH:{BASE_PATH}\n")

def launch_job(tool_name, args, memory = "64g", cpus = "1", num_tasks = "1"):
    os.makedirs(f"{BASE_PATH}/slurm", exist_ok=True)
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header(tool_name, path, memory, cpus, num_tasks)
    args_to_string = " ".join([str(arg) for arg in args])
    args_to_string += " $SLURM_ARRAY_TASK_ID"
    with open(path + "/slurm.sh", "w") as file:
        file.write(header)
        _write_env_setup(file)
        file.write(f"{PYTHON} tools/{tool_name}.py " + args_to_string)
    job_id = _submit_sbatch(f"{path}/slurm.sh")
    return job_id

def launch_job_multichromosome(tool_name, args, memory = "64g", cpus = "1", num_tasks = "1", max_concurrent = None):
    os.makedirs(f"{BASE_PATH}/slurm", exist_ok=True)
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_multichromosome(tool_name, path, memory, cpus, num_tasks, max_concurrent)
    args_to_string = " ".join([str(arg) for arg in args])
    args_to_string += " $CHR $ROW_ID"
    with open(path + "/slurm.sh", "w") as file:
        file.write(header)
        _write_env_setup(file)
        file.write(f"{PYTHON} tools/{tool_name}.py " + args_to_string)
    job_id = _submit_sbatch(f"{path}/slurm.sh")
    return job_id

def launch_job_with_custom_command_multichromosome(tool_name, id_counter, command, memory = "64g", cpus = "1", num_tasks = "1"):
    os.makedirs(f"{BASE_PATH}/slurm", exist_ok=True)
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_multichromosome(tool_name, path, memory, cpus, num_tasks)
    with open(path + f"/slurm_{id_counter}.sh", "w") as file:
        file.write(header)
        _write_env_setup(file)
        file.write(command)
    job_id = _submit_sbatch(f"{path}/slurm_{id_counter}.sh")
    return job_id

def launch_job_with_custom_command_by_chromosome(tool_name, id_counter, command, memory = "64g", cpus = "1", num_tasks = "1"):
    os.makedirs(f"{BASE_PATH}/slurm", exist_ok=True)
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_by_chromosome(tool_name, path, memory, cpus)
    with open(path + f"/slurm_{id_counter}.sh", "w") as file:
        file.write(header)
        _write_env_setup(file)
        file.write(f"for ROW_ID in $(seq 0 {str(int(num_tasks) - 1)}); do\n")
        file.write(command)
        file.write("done\n")
    job_id = _submit_sbatch(f"{path}/slurm_{id_counter}.sh")
    return job_id

def launch_job_with_custom_command_by_rowid(tool_name, id_counter, command, memory = "64g", cpus = "1", num_tasks = "1"):
    os.makedirs(f"{BASE_PATH}/slurm", exist_ok=True)
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_by_rowid(tool_name, path, memory, cpus, num_tasks)
    with open(path + f"/slurm_{id_counter}.sh", "w") as file:
        file.write(header)
        _write_env_setup(file)
        file.write(command)
    job_id = _submit_sbatch(f"{path}/slurm_{id_counter}.sh")
    return job_id



