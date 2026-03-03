import sys
from pathlib import Path
import numpy as np

if len(sys.argv) != 3:
    print("Compare two aurora runs and print differences.\nUsage:")
    print("diff.py [AURORA_OUTPUT_FILE1] [AURORA_OUTPUT_FILE2]")
    sys.exit(1)


def load_regions(file_path):
    regions = {}

    with open(file_path, "r") as f:
        for line in f:
            amb_count, target, tstart, tend, query, qstart, qend, direction, confidence, join_id, region = line.strip().split()

            region = int(region)
            join_id = int(join_id)
            tstart = int(tstart)
            tend = int(tend)

            if region not in regions:
                regions[region] = []

            regions[region].append((tstart, tend, line.strip()))

    return regions


def get_edit_path(region1_lines, region2_lines):
    cost_matrix = np.zeros((len(region1_lines) + 1, len(region2_lines) + 1))

    for i in range(len(region1_lines)):
        cost_matrix[i + 1, 0] = cost_matrix[i, 0] - 1
    
    for j in range(len(region2_lines)):
        cost_matrix[0, j + 1] = cost_matrix[0, j] - 1

    for i in range(len(region1_lines)):
        for j in range(len(region2_lines)):
            cost_matrix[i + 1, j + 1] = max(
                cost_matrix[i, j] - (region1_lines[i][2] != region2_lines[j][2]),
                cost_matrix[i, j + 1] - 1,
                cost_matrix[i + 1, j] - 1
            )
    
    i = len(region1_lines)
    j = len(region2_lines)

    final_trace = []

    while i > 0 or j > 0:
        best_val = float("-INFINITY")
        best_j = 0
        best_i = 0
        is_match = False
        if i > 0 and j > 0:
            is_match = (region1_lines[i - 1][2] == region2_lines[j - 1][2])
            cost = cost_matrix[i - 1, j - 1] - (not is_match)
            best_val = cost
            best_i = -1
            best_j = -1

        if i > 0:
            cost = cost_matrix[i - 1, j] - 1
            if cost > best_val:
                best_val = cost
                best_i = -1
                best_j = 0
        
        if j > 0:
            cost = cost_matrix[i, j - 1] - 1
            if cost > best_val:
                best_val = cost
                best_i = 0
                best_j = -1
        
        final_trace.append((
            is_match,
            region1_lines[i - 1] if best_i != 0 else None,
            region2_lines[j - 1] if best_j != 0 else None
        ))
        i += best_i
        j += best_j
    
    final_trace.reverse()
    return final_trace


def print_per_region(file1_name, file2_name, regions1, regions2):
    all_regions = sorted(regions1.keys() | regions2.keys())

    file_names = [file1_name, file2_name]
    max_file_name_len = max(len(f) for f in file_names)
    file_names = ["  " + (" " * (max_file_name_len - len(f))) + f for f in file_names]

    no_mismatches_found_yet = True

    for region_idx in all_regions:
        try:
            r1 = regions1[region_idx]
            r2 = regions2[region_idx]
        except KeyError:
            print(f"Region {region_idx} not found in both files! Skipping....")

        trace = get_edit_path(r1, r2)

        prior_mismatches = []
        for is_match, line1, line2 in trace:
            if not is_match:
                if no_mismatches_found_yet:
                    print(f"Region: {region_idx}")
                    no_mismatches_found_yet = False
                if line1 is not None:
                    prior_mismatches.append((*line1, 0))
                if line2 is not None:
                    prior_mismatches.append((*line2, 1))
            elif len(prior_mismatches) > 0:
                bound_start = min(l[0] for l in prior_mismatches)
                bound_end = max(l[1] for l in prior_mismatches)
                print(f"Target from {bound_start} to {bound_end} does not match. Differences:")
                for l in prior_mismatches:
                    print(f"{file_names[l[3]]}: {l[2]}")
                prior_mismatches.clear()
        
        if len(prior_mismatches) > 0:
            bound_start = min(l[0] for l in prior_mismatches)
            bound_end = max(l[1] for l in prior_mismatches)
            print(f"Target from {bound_start} to {bound_end} does not match. Differences:")
            for l in prior_mismatches:
                print(f"{file_names[l[3]]}: {l[2]}")
            prior_mismatches.clear()
        

aurora_output_file1 = Path(sys.argv[1])
aurora_output_file2 = Path(sys.argv[2])

regions1 = load_regions(aurora_output_file1)
regions2 = load_regions(aurora_output_file2)

print_per_region(aurora_output_file1.name, aurora_output_file2.name, regions1, regions2)
