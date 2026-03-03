import sys
from pathlib import Path

if len(sys.argv) not in [2, 3]:
    print("Collect and print out or save family statistics for a aurora run.\nUsage:")
    print("family_stats.py [AURORA_AMB_OUTPUT_FILE] [OPTIONAL_TSV_SAVE_PATH]")
    sys.exit(1)


prior_region = -1
already_seen_joins = set()
family_counts = {}
family_coverage = {}

aurora_output_file = Path(sys.argv[1])

with open(aurora_output_file, "r") as f:
    for line in f:
        amb_count, target, tstart, tend, query, qstart, qend, direction, confidence, join_id, region = line.strip().split()

        region = int(region)
        join_id = int(join_id)
        query_families = query.split(",")
        # Just take first direction if there is multiple...
        direction = direction.split(",")[0]

        if region != prior_region:
            already_seen_joins.clear()
            prior_region = region
        
        if join_id not in already_seen_joins:
            for family in query_families:
                family_counts[family] = family_counts.get(family, 0) + 1
            already_seen_joins.add(join_id)

        for family in query_families:
            family_coverage[family] = family_coverage.get(family, 0) + abs(int(tend) - int(tstart))


def write_results(file, family_counts, family_coverage):
    file.write("Family\tCount\tCoverage\n")

    for family in sorted(family_counts, key=lambda k: (family_counts[k], family_coverage[k]), reverse=True):
        file.write(f"{family}\t{family_counts[family]}\t{family_coverage[family]}\n")

if len(sys.argv) == 3:
    if Path(sys.argv[2]).resolve() == aurora_output_file:
        print("Error, save path matches aurora file path, not saving...")
    else:
        with open(sys.argv[2], "w") as f:
            write_results(f, family_counts, family_coverage)

write_results(sys.stdout, family_counts, family_coverage)