import sys
from pathlib import Path

if len(sys.argv) not in [2, 3]:
    print("Count inversions for a given aurora run.\nUsage:")
    print("inversions.py [AURORA_AMB_OUTPUT_FILE] [OPTIONAL_TSV_SAVE_PATH]")
    sys.exit(1)


prior_region = -1
join_ids_to_prior_direction = {}
inversions_per_region = {}

aurora_output_file = Path(sys.argv[1])

with open(aurora_output_file, "r") as f:
    for line in f:
        amb_count, target, tstart, tend, query, qstart, qend, direction, confidence, join_id, region = line.strip().split()

        region = int(region)
        join_id = int(join_id)
        # Just take first direction if there is multiple...
        direction = direction.split(",")[0]

        if region != prior_region:
            join_ids_to_prior_direction.clear()
            inversions_per_region[region] = 0
            prior_region = region
        
        if join_id in join_ids_to_prior_direction and join_ids_to_prior_direction[join_id] != direction:
            inversions_per_region[region] = inversions_per_region.get(region, 0) + 1
        
        join_ids_to_prior_direction[join_id] = direction


def write_results(file, inversions_per_region):
    file.write("Region\tInversion Count\n")
    file.write(f"ALL\t{sum(inversions_per_region.values())}\n")

    for region in sorted(inversions_per_region.keys()):
        file.write(f"{region}\t{inversions_per_region[region]}\n")

if len(sys.argv) == 3:
    if Path(sys.argv[2]).resolve() == aurora_output_file:
        print("Error, save path matches aurora file path, not saving...")
    else:
        with open(sys.argv[2], "w") as f:
            write_results(f, inversions_per_region)

write_results(sys.stdout, inversions_per_region)
