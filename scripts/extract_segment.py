import sys
from pathlib import Path

if len(sys.argv) != 4:
    print("Extract a Segment from a CAF file.\nUsage:")
    print("extract_segment.py [CAF_FILE_PATH] [TARGET_START] [TARGET_END]")
    print("Target start and end are the start and end positions on the target sequence.")
    sys.exit(1)

caf_file = Path(sys.argv[1]).resolve()
start = int(sys.argv[2])
end = int(sys.argv[3])


with open(caf_file, "r") as cf:
    with open(caf_file.parent / (caf_file.stem + f"_{start}_{end}.caf"), "w") as cf_new:
        for line in cf:
            if line.strip() == "":
                continue
            al_start, al_end = [int(v) for v in line.strip().split(",")[5:7]]
            if al_start <= end and al_end >= start:
                cf_new.write(line)


            
