
import pandas as pd

df = pd.DataFrame({"clones_file": snakemake.input.clones, "gex_data": snakemake.input.gex})
df["gex_data_type"] = "10x_h5"

with open(snakemake.output.tsv, mode="w") as f:
    df.to_csv(f, sep="\t", index=False)