import sys
from pathlib import Path

if len(sys.argv) != 2:
    print("Extract history trace statistics from debug visualization output.")
    print("Usage: stats.py VISUALIZATION_DIRECTORY") 
    sys.exit(1)

worst_ranks = []
worst_scores = []

for dir in Path(sys.argv[1]).iterdir():
    if not dir.is_dir():
        continue

    worst_rank = 0
    worst_rel_threshold = 0.0

    with open(dir / "final_trace_stats.csv") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            vals = [v.strip() for v in line.split(",")]
            worst_rank = max(worst_rank, int(vals[2]))
            worst_rel_threshold = min(worst_rel_threshold, float(vals[3]))

    print(f"Region {dir.name}: rank: {worst_rank}, threshold: {worst_rel_threshold}")
    worst_ranks.append(worst_rank)
    worst_scores.append(worst_rel_threshold)

print(f"Worst Rank Found: {max(worst_ranks)}")
print(f"Worst Score Found: {min(worst_scores)}")

