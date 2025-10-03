import numpy as np
import datetime
import os

from constants import BASE_PATH, PYTHON_ENV

ROOT_TOOL_PATH = f"{BASE_PATH}/tools"

RESTRICT_RESOURCES = False

NUMBER_CHROMSOMES = 22

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
    
def create_slurm_header_multichromosome(tool_name, path, memory = "64g", cpus = "1", num_tasks = "1"):
    return f"""#!/bin/bash
#SBATCH --job-name={tool_name}
#SBATCH --time=1:00:00
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --out={path}/out_%j_%a.log
#SBATCH --error={path}/error_%j_%a.log
#SBATCH --array=0-{str(int(num_tasks) - 1)}

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

def launch_job(tool_name, args, memory = "64g", cpus = "1", num_tasks = "1"):
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header(tool_name, path, memory, cpus, num_tasks)
    args_to_string = " ".join([str(arg) for arg in args])
    args_to_string += " $SLURM_ARRAY_TASK_ID"
    with open(path + "/slurm.sh", "w") as file:
        file.write(header)
        file.write("if [[ $(hostname) == ne1* ]]; then\n")
        file.write("    module load Python/3.10.8-GCCcore-12.2.0\n")
        file.write(f"    source {PYTHON_ENV}/bin/activate\n")
        file.write("    module load bcftools\n")
        file.write("else\n")
        file.write("    module load jupyter3\n")
        file.write("fi\n")
        file.write(f"cd {BASE_PATH}\n")
        file.write(f"export PYTHONPATH=$PYTHONPATH:{BASE_PATH}\n")
        file.write(f"python3 tools/{tool_name}.py " + args_to_string)
    os.system(f"sbatch {path}/slurm.sh")
    return path

def launch_job_multichromosome(tool_name, args, memory = "64g", cpus = "1", num_tasks = "1"):
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_multichromosome(tool_name, path, memory, cpus, num_tasks)
    args_to_string = " ".join([str(arg) for arg in args])
    args_to_string += " $CHR $ROW_ID"
    with open(path + "/slurm.sh", "w") as file:
        file.write(header)
        file.write("if [[ $(hostname) == ne1* ]]; then\n")
        file.write("    module load Python/3.10.8-GCCcore-12.2.0\n")
        file.write(f"    source {PYTHON_ENV}/bin/activate\n")
        file.write("    module load bcftools\n")
        file.write("else\n")
        file.write("    module load jupyter3\n")
        file.write("fi\n")
        file.write(f"cd {BASE_PATH}\n")
        file.write(f"export PYTHONPATH=$PYTHONPATH:{BASE_PATH}\n")
        file.write(f"python3 tools/{tool_name}.py " + args_to_string)
    os.system(f"sbatch {path}/slurm.sh")
    return path

def launch_job_with_custom_command_multichromosome(tool_name, id_counter, command, memory = "64g", cpus = "1", num_tasks = "1"):
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_multichromosome(tool_name, path, memory, cpus, num_tasks)
    with open(path + f"/slurm_{id_counter}.sh", "w") as file:
        file.write(header)
        file.write("if [[ $(hostname) == ne1* ]]; then\n")
        file.write("    module load Python/3.10.8-GCCcore-12.2.0\n")
        file.write(f"    source {PYTHON_ENV}/bin/activate\n")
        file.write("    module load bcftools\n")
        file.write("else\n")
        file.write("    module load jupyter3\n")
        file.write("fi\n")
        file.write(f"cd {BASE_PATH}\n")
        file.write(f"export PYTHONPATH=$PYTHONPATH:{BASE_PATH}\n")
        file.write(command)
    os.system(f"sbatch {path}/slurm_{id_counter}.sh")

def launch_job_with_custom_command_by_chromosome(tool_name, id_counter, command, memory = "64g", cpus = "1", num_tasks = "1"):
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_by_chromosome(tool_name, path, memory, cpus)
    with open(path + f"/slurm_{id_counter}.sh", "w") as file:
        file.write(header)
        file.write("if [[ $(hostname) == ne1* ]]; then\n")
        file.write("    module load Python/3.10.8-GCCcore-12.2.0\n")
        file.write(f"    source {PYTHON_ENV}/bin/activate\n")
        file.write("    module load bcftools\n")
        file.write("else\n")
        file.write("    module load jupyter3\n")
        file.write("fi\n")
        file.write(f"cd {BASE_PATH}\n")
        file.write(f"export PYTHONPATH=$PYTHONPATH:{BASE_PATH}\n")
        file.write(f"for ROW_ID in $(seq 0 {str(int(num_tasks) - 1)}); do\n")
        file.write(command)
        file.write("done\n")
    os.system(f"sbatch {path}/slurm_{id_counter}.sh")

def launch_job_with_custom_command_by_rowid(tool_name, id_counter, command, memory = "64g", cpus = "1", num_tasks = "1"):
    path = f"{BASE_PATH}/slurm/{get_nyc_time_in_month_day_time()}_{tool_name}"
    os.mkdir(path)
    header = create_slurm_header_by_rowid(tool_name, path, memory, cpus, num_tasks)
    with open(path + f"/slurm_{id_counter}.sh", "w") as file:
        file.write(header)
        file.write("if [[ $(hostname) == ne1* ]]; then\n")
        file.write("    module load Python/3.10.8-GCCcore-12.2.0\n")
        file.write(f"    source {PYTHON_ENV}/bin/activate\n")
        file.write("    module load bcftools\n")
        file.write("else\n")
        file.write("    module load jupyter3\n")
        file.write("fi\n")
        file.write(f"cd {BASE_PATH}\n")
        file.write(f"export PYTHONPATH=$PYTHONPATH:{BASE_PATH}\n")
        file.write(command)
    os.system(f"sbatch {path}/slurm_{id_counter}.sh")



