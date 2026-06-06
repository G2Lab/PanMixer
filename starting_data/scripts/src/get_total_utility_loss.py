import json
import numpy as np
import time
from pathlib import Path

from paths import BASE_PATH

STARTING_DATA_PATH = Path(BASE_PATH) / "starting_data"


def load_subjects():
    subjects_path = STARTING_DATA_PATH / "subjects.txt"
    if not subjects_path.exists():
        subjects_path = STARTING_DATA_PATH / "subjects_files" / "subjects.txt"
    return np.loadtxt(subjects_path, dtype=str)


def compute_total_utility_loss():
    subjects = load_subjects()
    t_total = time.time()
    true_max_utility_loss = {}
    for i, subject in enumerate(subjects):
        t_subject = time.time()
        true_max_utility_loss[subject] = 0
        for chrom in range(1, 23):
            utility_loss = np.load(STARTING_DATA_PATH / f"chr{chrom}/subjects/{i}/utility_loss.npy")
            true_max_utility_loss[subject] += 2 * np.sum(utility_loss)
        print(f"[TIMING] Subject {subject} ({i+1}/{len(subjects)}) completed in {time.time() - t_subject:.2f}s")

    with open(STARTING_DATA_PATH / "utility_loss.json", "w") as f:
        json.dump(true_max_utility_loss, f)

    print(f"[TIMING] Total utility loss computation completed in {time.time() - t_total:.2f}s")


if __name__ == "__main__":
    compute_total_utility_loss()
