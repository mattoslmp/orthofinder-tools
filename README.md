# OrthoFinder Tools

## Idea

* Calculate the most common gene name of each orthogroup by majority vote: `annotate_orthogroups`
* Create publication-ready comparative-genomics plots from OrthoFinder results: `orthofinder_plots`

## Setup

```shell
pip install orthofinder-tools
```

For the Databiomics fork directly from GitHub:

```shell
pip install "git+https://github.com/mattoslmp/orthofinder-tools.git@master"
```

## Usage

### annotate_orthogroups

#### Prerequisites

Your FASTA sequences must have some description, e.g.:

```text
>gnl|extdb|STRAIN-XY_000001 DNA-directed RNA polymerase subunit beta' [Pediococcus stilesii]
MIDVNKFESMQIGLASPDKIRMWSYGEVKKPETINYRTLKPEKDGLFDERIFGPTKDYECACGKYKRIRY
...
```

From this protein, `DNA-directed RNA polymerase subunit beta'` will be extracted.

#### Command line usage

```shell
annotate_orthogroups --help

annotate_orthogroups \
    --orthogroups_tsv /path/to/N0_or_Orthogroups.tsv \
    --hog True \
    --fasta_dir /path/to/fastas \
    --file_endings faa \
    --out outfile.tsv \
    --simple True
```

If `--simple=False` resulting tsv looks like this:

| HOG | Best Gene Name | Gene Name Occurrences |
| --- | --- | --- |
| N0.HOG0000000 | amino acid ABC transporter | `{JSON}` |
| N0.HOG0000001 | IS30 family transposase | `{JSON}` |
| N0.HOG0000002 | IS5/IS1182 family transposase | `{JSON}` |

The JSON is a dictionary with key='gene name' -> value=occurrence, for example:

```json5
{
  'Integrase core domain protein': 47,
  'hypothetical protein': 15,
  'IS30 family transposase': 126
}
```

If `--simple=True` the resulting TSV contains the orthogroup and its majority-vote gene name.

#### Usage as Python class

```python
from orthofinder_tools import OrthogroupToGeneName

PATH_TO_ORTHOFINDER_FASTAS = '/path/to/OrthoFinder/fastas'
CURRENT_FOLDER = 'Results_Mon00'

otg = OrthogroupToGeneName(
    fasta_dir=PATH_TO_ORTHOFINDER_FASTAS,
    file_endings='faa',
)
otg.load_hog(
    hog_tsv=f'{PATH_TO_ORTHOFINDER_FASTAS}/OrthoFinder/{CURRENT_FOLDER}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
)
```

`otg.majority_dict` is a Python dictionary with orthogroup IDs as keys and majority-vote gene names as values.

### orthofinder_plots

This script started as a port of [roary_plots](https://github.com/sanger-pathogens/Roary/tree/master/contrib/roary_plots) by Marco Galardini and now includes publication-ready Databiomics extensions.

The legacy command remains valid:

```shell
orthofinder_plots \
  --tree data/SpeciesTree_rooted.txt \
  --orthogroups_tsv data/Orthogroups.tsv \
  --out output
```

For manuscript-ready exports, generate vector and high-resolution raster files together:

```shell
orthofinder_plots \
  --tree data/SpeciesTree_rooted.txt \
  --orthogroups_tsv data/Orthogroups.tsv \
  --out output \
  --formats svg,pdf,png \
  --dpi 600 \
  --permutations 100 \
  --max_patterns 20
```

For hierarchical orthogroups (`N0.tsv`), add `--hog True`.

#### Generated figures

The plotting workflow now produces:

* `pangenome_frequency` — orthogroup frequency across genomes.
* `pangenome_matrix` — species tree aligned to the orthogroup presence/absence matrix.
* `pangenome_pie` — publication-ready donut summarising core, soft-core, shell and cloud orthogroups. The historic filename is retained for compatibility.
* `pangenome_accumulation` — pan-genome and core-genome accumulation curves across random genome orders, with 95% empirical intervals.
* `orthogroup_patterns` — the most frequent exact presence/absence patterns across genomes.
* `orthofinder_plot_summary.json` — machine-readable counts, thresholds, output formats and figure manifest.

Supported figure formats are `svg`, `pdf`, `png` and `tiff`. SVG/PDF are recommended for vector artwork; PNG/TIFF default to 600 dpi for publication-quality raster output.

The pangenome categories follow the existing package convention: core >=99% of genomes, soft-core >=95% and <99%, shell >=15% and <95%, and cloud <15%. Integer thresholds are calculated with `ceil`, and the exact values used are recorded in the JSON summary.

#### Usage as Python class

```python
from orthofinder_tools import create_plots

summary = create_plots(
    tree='/path/to/SpeciesTree_rooted.txt',
    orthogroups_tsv='/path/to/Orthogroups.tsv',
    out='/path/to/output/folder',
    formats='svg,pdf,png',
    dpi=600,
    permutations=100,
    max_patterns=20,
    no_labels=False,
)

print(summary)
```

### Notes for reproducible figures

* Species names in `SpeciesTree_rooted.txt` must match the columns in `Orthogroups.tsv`.
* The accumulation plot uses a deterministic random seed internally, so the same input and parameters reproduce the same curves.
* The plotting layer does not modify OrthoFinder results; it only reads the tree and orthogroup tables and writes derived figures/summary files.
