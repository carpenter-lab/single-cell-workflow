import pandas as pd
import tcrdist
from tcrdist.repertoire import TCRrep
from tcrdist.public import TCRpublic

df = pd.read_csv(snakemake.input["csv"])


tr = TCRrep(cell_df = df, organism = "human", chains = ["beta"])

tp = TCRpublic(tcrrep = tr, output_html_name = "quasi_public_clones.html")

tp.query_str = 'K_neighbors > 2'
public = tp.report()


public['quasi_public_df'].to_csv(snakemake.output["quasi_public"])
public['nn_summary'].to_csv(snakemake.output["nn_summary"])
public['clone_df'].to_csv(snakemake.output["clone_df"])