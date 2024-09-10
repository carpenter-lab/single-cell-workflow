# Single-cell RNA Sequencing Workflow

This repository contains all the scripts to analyze single-cell RNA sequencing data after running
[Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest).
Each script is prepared for use in a [Snakemake](https://snakemake.github.io) workflow.


## Setup
There are two enviornments that need to be set up prior to running the workflow.

### R Enviornment
The R enviornment is managed using [renv](https://rstudio.github.io/renv/).

After installing R, run the following code in the terminal:

```bash
Rscript -e "renv::restore()"
```

This will install all the required packages into a *project-specific* directory.
This will *not* impact your global packages.

## Python
There are two options for Python, [venv](https://docs.python.org/3/library/venv.html) or [Conda](https://docs.conda.io)/[Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).
There is some difference between these tools, but if you need a starting point, choose Mamba.

### Mamba (or Conda)

To set up the enviornment just run the following code in the terminal:

```bash
conda env create -f enviornment.yml
```

### venv (TODO)

## Configuring the Workflow

### Creating PEP

To manage samples the workflow uses [PEPkit](https://pep.databio.org).
First edit the file in `pep/sample_table.csv` to include your samples.
Then edit `pep/config.yml` to include the paths to the sample inputs.
These paths should include some wildcard patterns, see example in repository, to automatically create file paths.

### Configuring the Workflow

The configuration for the workflow is complex, but is *technically* optional.
All rules can be included in another Snakefile and bypass the pre-configured workflow to include other rules ([see here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows))

#### YAML header

You can use YTE to run Python code in the YAML configuration file.
Some examples of this can be seen later.

```yaml
__use_yte__: true                   # requried if using YTE
__definitions__:
  - from itertools import product   # required if importing Python packages
```

#### Preprocessing

```yaml
preprocessing:
  normalisation:
    method: "sctransform"         # Normalisation Method (LogNormalize or sctransform)
    vars_to_regress: null         # Variables to regress during normalisation, null for no regression
```

#### Subclustering

```yaml
subcluster:
  all_data_key: "all_data"         # Key to use for the filename containing all cells. Can be anything.
  subclusters:
    - name: t_cells                # Name to use for the subset. Can be anything.
      params: null                 # Set to null for collecting only TCR+ cells
    - name: macs
      params:
        idents_use: "predicted.ann_finest_level"                           # To do other subsets, set this to the metadata column
        idents_keep: [ "Alveolar Mφ CCL3+", "Non-classical monocytes",     # And include a list of the identities to subset
                       "Alveolar macrophages", "Monocyte-derived Mφ",
                       "Alveolar Mφ MT-positive", "Classical monocytes" ]
```

#### Differential Expression

```yaml
differential_expression:
  test: "MAST"                     # DE Test to use. Choose any of: wilcox, wilcox_limma, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2
  params:
    subsets: [ "t_cells", "macs", "all_data" ]  # Subsets to use for DE testing
    assays: [ "ADT", "SCT" ]                    # Assays to use for DE testing
    group_by: [ "predicted.ann_finest_level", "predicted.ann_level_2", "tcr_ident" ] # Metadata columns to use for DE testing
```

#### Plotting

The plotting configuration is the most complex as this forms the basis of every plot.
A condensed format is included below.

```yaml
plotting:
  dot_plot|umap_plot:
    subsets: [ "t_cells", "macs", "all_data" ]
    assays: [ "ADT", "SCT" ]
    reductions: [ "harmony", "pca" ] # only for UMAP
    features: # only for Dot Plot
      ADT: [...]
      SCT: [...]
```

The group by parameter is different between the two plot functions.
For dot plots the parameter is simelar to other options such as in [Differential Expression](#differential-expression).
```yaml
plotting:
  dot_plot:
    group_by: [ "predicted.ann_finest_level", "predicted.ann_level_2" ]
```

However, for UMAP plots this becomes more complex.
There are two `group_by` parameters, `group_by` and `group_by_yte_only`.
Both have the same structure, however `group_by_yte_only` will be parsed by YTE.
`group` is the metadata column to **color** the UMAP plot by.
`split` is the metadata column to **split** by, i.e. create separate plots.
This can be a list, and when `null` there is no split.

```yaml
plotting:
  umap_plot:
    group_by_yte_only:
      ?for assay in ["SCT"]:
        - group: ?f"{assay}_{{reduction}}_clusters"   # Two curly braces to prevent replacement now
          split: [ null, "condition" ]
    group_by:
      - group: "clonal_expansion_type"
        split: "condition"
      - group: "predicted.ann_finest_level"
        split: [ null, "predicted.ann_level_2" ]
      - group: "tcr_ident"
        split: null
```

## Running the Workflow

There are a few ways to run the workflow.
You can show what will be run with the `--dry-run` option.
The simplest is below.

``` bash
snakemake --cores NUMBER_OF_CORES
```

To use the included `run_workflow.py` script run either of these:

```bash
# for a dry run
DRYRUN=1 ./run_workflow.py 

# setup upload to Google Drive does not run workflow
# DO NOT SHARE THE "token.pickle" file. 
SETUP_GDRIVE=1 ./run_workflow.py 

# Run workflow with 10 cores, creates report.html, and uploads report to Google Drive
WORKFLOW_CORES=10 ./run_workflow.py
```

This cannot be emphasized enough. 
**DO NOT SHARE THE `token.pickle` FILE!** 
This file can allow **anyone** the access to upload to your Google Drive.