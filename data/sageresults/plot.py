# /// script
# dependencies = [
#   "polars",
#   "matplotlib",
# ]
# ///

import json
import numpy as np
import matplotlib.pyplot as plt

data = json.load(open("query_results/results.json"))
labels = json.load(open("ubb_elution_groups.json"))
id_to_label = {x["id"]: x["human_id"] for x in labels}

num_elution_groups = len(data)
fig, ax = plt.subplots(ncols=num_elution_groups, nrows=2, figsize=(10, 7))

# Top is the MS1 and bottom is the MS2

# Getting the data for ms1 is like so:
# x["result"]["ms1_stats"]["retention_time_miliseconds"]
# x["result"]["ms1_stats"]["transition_intensities"][y]

# Getting the data for ms2 is like so:
# x["result"]["ms2_stat"]["retention_time_miliseconds"]
# x["result"]["ms2_stats"]["transition_intensities"][y]


def ms_to_mins(x):
    return (x / 1_000) / 60


def infinite_color_cycle():
    colors = ["#308AAD", "#C8102E", "#96D8D8", "#B2EE79"]
    while True:
        for c in colors:
            yield c


for i, x in enumerate(data):
    ms1_rts = np.array(x["result"]["ms1_stats"]["retention_time_miliseconds"])
    ms1_rts = ms_to_mins(ms1_rts)
    ms2_rts = np.array(x["result"]["ms2_stats"]["retention_time_miliseconds"])
    ms2_rts = ms_to_mins(ms2_rts)

    colors = infinite_color_cycle()
    sorted_ms1_keys = sorted(x["result"]["ms1_stats"]["transition_intensities"].keys())
    sorted_ms2_keys = sorted(x["result"]["ms2_stats"]["transition_intensities"].keys())

    for j in sorted_ms1_keys:
        j = x["result"]["ms1_stats"]["transition_intensities"][j]
        ax[0, i].plot(ms1_rts, j, color=next(colors))

    for j in sorted_ms2_keys:
        j = x["result"]["ms2_stats"]["transition_intensities"][j]
        ax[1, i].plot(ms2_rts, j, color=next(colors))

    # Label axes ...
    ax[0, i].set_xlabel("Retention Time (min)")
    ax[1, i].set_xlabel("Retention Time (min)")
    ax[0, i].set_ylabel("Intensity")
    ax[1, i].set_ylabel("Intensity")

    # Label titles
    human_lab = id_to_label[x["elution_group"]["id"]]
    ax[0, i].set_title(human_lab + " MS1")
    ax[1, i].set_title(human_lab + " MS2")

plt.tight_layout()
plt.savefig("ubb_peptide_plot.png")
