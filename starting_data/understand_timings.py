import pandas as pd

df = pd.read_csv("timings.csv")

def parse_time_to_s(s):
    parts = s.split(":")
    parts = [float(p) for p in parts]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    else:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]

def parse_memory(s):
    # Handles K/M/G suffixes — adjust if your data differs
    unit = s[-1]
    value = float(s[:-1])
    multiplier = {"K": 1, "M": 1024, "G": 1024**2}.get(unit, 1)
    return value * multiplier

# Apply to columns (replace with your actual column names)
df["elapsed_s"] = df["Elapsed"].apply(parse_time_to_s)
df["cpu_s"] = df["TotalCPU"].apply(parse_time_to_s)
df["mem_kb"] = df["MaxRSS"].apply(parse_memory)


# jobs that are graph preperation: get_pmi_weights
# jobs that are obfuscation: optimizer, stacker



# get me max optimizer and stacker time

df_prep = df[df["JobName"] == "get_pmi_weights"]
df_optimizer = df[df["JobName"] == "optimizer"]
df_stacker = df[df["JobName"] == "stacker"]

def to_hours(s):
    return s / (60 * 60)
def to_gb(mem_kb):
    return mem_kb / (1000 * 1000)
print("graph preperation:")
print("cpu time:", to_hours(df_prep["cpu_s"].sum()), "hr", to_gb(df_prep["mem_kb"].max()), "GB")

print("obfuscation (per subject/capacity, aggregated across 22 chromosomes):")
print("optimizer  — avg elapsed:", to_hours(df_optimizer["elapsed_s"].mean()), "hr",
      "max elapsed:", to_hours(df_optimizer["elapsed_s"].max()), "hr",
      "max RSS:", to_gb(df_optimizer["mem_kb"].max()), "GB")
print("stacker    — avg elapsed:", to_hours(df_stacker["elapsed_s"].mean()), "hr",
      "max elapsed:", to_hours(df_stacker["elapsed_s"].max()), "hr",
      "max RSS:", to_gb(df_stacker["mem_kb"].max()), "GB")
print("total obfuscation (optimizer + stacker) avg elapsed:",
      to_hours(df_optimizer["elapsed_s"].mean() + df_stacker["elapsed_s"].mean()), "hr")