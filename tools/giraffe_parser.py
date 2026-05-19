import sys

import sys

def read_stats_file(stats_filepath):
    results = {}
    with open(stats_filepath, "r") as f:
        # 1) Total alignments: 11755181
        line = f.readline().strip()
        results["total_alignments"] = int(line.split()[2])

        # 2) Total primary: 11755181
        line = f.readline().strip()
        results["total_primary"] = int(line.split(":", 1)[1].strip().split()[0])

        # 3) Total secondary: 0
        line = f.readline().strip()
        results["total_secondary"] = int(line.split(":", 1)[1].strip().split()[0])

        # 4) Total aligned: 11521845
        line = f.readline().strip()
        results["total_aligned"] = int(line.split(":", 1)[1].strip().split()[0])

        # 5) Total perfect: 8949888
        line = f.readline().strip()
        results["total_perfect"] = int(line.split(":", 1)[1].strip().split()[0])

        # 6) Total gapless (softclips allowed): 11039914
        line = f.readline().strip()
        results["total_gapless_softclips_allowed"] = int(line.split(":", 1)[1].strip().split()[0])

        # 7) Total paired: 0
        line = f.readline().strip()
        results["total_paired"] = int(line.split(":", 1)[1].strip().split()[0])

        # 8) Total properly paired: 0
        line = f.readline().strip()
        results["total_properly_paired"] = int(line.split(":", 1)[1].strip().split()[0])

        # 9) Alignment score: mean 154.327, median 160, stdev 17.9166, max 160 (…)
        line = f.readline().strip()
        rhs = line.split(":", 1)[1].split("(", 1)[0]  # drop the "(... reads)" part
        parts = [p.strip() for p in rhs.split(",")]
        
        results["alignment_score_mean"]   = float(parts[0].split()[1])
        results["alignment_score_median"] = float(parts[1].split()[1])
        results["alignment_score_stdev"]  = float(parts[2].split()[1])
        results["alignment_score_max"]    = float(parts[3].split()[1])

        results["alignment_score_max_160_reads"] = int(line.split(" ")[-2][1:])

        # 10) Mapping quality: mean 49.6319, median 60, stdev 21.4919, max 60 (…)
        line = f.readline().strip()
        rhs = line.split(":", 1)[1].split("(", 1)[0]
        parts = [p.strip() for p in rhs.split(",")]
        results["mapping_quality_mean"]   = float(parts[0].split()[1])
        results["mapping_quality_median"] = float(parts[1].split()[1])
        results["mapping_quality_stdev"]  = float(parts[2].split()[1])
        results["mapping_quality_max"]    = float(parts[3].split()[1])

        results["mapping_quality_max_60_reads"] = int(line.split(" ")[-2][1:])


        # 11) Insertions: 1075032 bp in 295719 read events
        line = f.readline().strip()
        rhs = line.split(":", 1)[1].strip()
        results["insertions_bp"] = int(rhs.split()[0])
        results["insertions_events"] = int(rhs.rsplit(" in ", 1)[1].split()[0])

        # 12) Deletions: 884781 bp in 296886 read events
        line = f.readline().strip()
        rhs = line.split(":", 1)[1].strip()
        results["deletions_bp"] = int(rhs.split()[0])
        results["deletions_events"] = int(rhs.rsplit(" in ", 1)[1].split()[0])

        # 13) Substitutions: 7064082 bp in 7064082 read events
        line = f.readline().strip()
        rhs = line.split(":", 1)[1].strip()
        results["substitutions_bp"] = int(rhs.split()[0])
        results["substitutions_events"] = int(rhs.rsplit(" in ", 1)[1].split()[0])

        # 14) Matches: 1697100535 bp (144.37 bp/alignment)
        line = f.readline().strip()
        rhs = line.split(":", 1)[1].strip()
        results["matches_bp"] = int(rhs.split()[0])
        # value inside parentheses before space, e.g. "144.37 bp/alignment"
        paren = rhs.split("(", 1)[1].split(")", 1)[0]
        results["matches_bp_per_alignment"] = float(paren.split()[0])

        # 15) Softclips: 23037101 bp in 541609 read events
        line = f.readline().strip()
        rhs = line.split(":", 1)[1].strip()
        results["softclips_bp"] = int(rhs.split()[0])
        results["softclips_events"] = int(rhs.rsplit(" in ", 1)[1].split()[0])

        # 16) Total time: 1025.24 seconds
        line = f.readline().strip()
        results["total_time_seconds"] = float(line.split(":", 1)[1].strip().split()[0])

        # 17) Speed: 11465.8 reads/second
        line = f.readline().strip()
        results["speed_reads_per_second"] = float(line.split(":", 1)[1].strip().split()[0])

    return results