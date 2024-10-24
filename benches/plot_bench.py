# /// script
# dependencies = [
#   "altair",
#   "polars",
#   "vl-convert-python",
# ]
# ///

from pathlib import Path
import polars as pl
import altair as alt
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("benchmark_file")
args = parser.parse_args()

data = (
    pl.read_json(args.benchmark_file)
    .explode("results")
    .with_columns(
        bench=pl.col("results").struct.field("name"),
        context=pl.col("results").struct.field("context"),
        iterations=pl.col("results").struct.field("iterations"),
    )
    .explode("iterations")
    .with_columns(time_seconds=pl.col("iterations").struct.field("duration_seconds"))
)
print(data)

for context, sdf in data.group_by("context"):
    context = context[0]
    print(context)
    # Sample to max 100 points on each bench
    sdf = sdf.filter(pl.int_range(pl.len()).shuffle().over("bench") < 100)
    alt.Chart(sdf).mark_point().encode(
        x="bench",
        y=alt.Y("time_seconds", scale=alt.Scale(type="log")),
        color="bench",
    ).properties(
        width=600,
        height=600,
        title=f"{context}",
    ).save(f"{context}.png")
