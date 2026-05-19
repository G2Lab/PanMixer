import numpy as np
import gzip

cache_path = "./cache"

def get_variant_info(vcf_file):
    variant_level_info = {}
    total = 0
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7]
            info_parts = info.split(";")
            level = -1
            for part in info_parts:
                if part.startswith("LV="):
                    level = int(part.split("=")[1])
                    break
            
            if len(alt) > 1:
                variant_type = "MULTIALLELIC"
            elif len(ref) == 1 and len(alt[0]) == 1 and (alt[0] != "." and ref != "."):
                variant_type = "SNP"
            elif len(ref) > len(alt[0]) and alt[0] != ".":
                variant_type = "DEL"
            elif len(alt[0]) > len(ref) and ref != ".":
                variant_type = "INS"
            elif alt[0] == "." and ref != ".":
                variant_type = "DEL"
            elif ref == "." and alt[0] != ".":
                variant_type = "INS"
            else:
                variant_type = "OTHER"

            if level not in variant_level_info:
                variant_level_info[level] = {}
            
            if variant_type not in variant_level_info[level]:
                variant_level_info[level][variant_type] = 0
            
            variant_level_info[level][variant_type] += 1
            total += 1

    return variant_level_info, total
            
import sys
import json
if __name__ == "__main__":
    vcf_file = sys.argv[1]
    variant_level_info, total = get_variant_info(vcf_file)
    json.dump(variant_level_info, open(f"{cache_path}/variant_level_info.json", "w"), indent=4)
    
    print(f"Total variants: {total}")
    for level, types in variant_level_info.items():
        print(f"Level {level}:")
        for variant_type, count in types.items():
            print(f"  {variant_type}: {count}")

def get_block_info():
    for i in range(1,23):
        
    