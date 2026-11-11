# PHANES

**Proteome-wide Human Atlas of iNterpretable fEatureS**

PHANES is an interactive web application that provides biologically interpretable descriptors for the human proteome. It is meant as a general-purpose feature resource for exploratory analysis, comparative studies and machine-learning workflows, not as a predictor of one specific phenomenon.

**Web server:** https://tools.tartaglialab.com/phanes/

---

## Features

For every human protein, PHANES provides:

| Group | Examples |
|---|---|
| **Physico-chemical** (sequence) | hydrophobicity, charge, aggregation, disorder, nucleic-acid binding, secondary-structure and membrane propensities, at residue level and as protein-level summaries |
| **Structural** (AlphaFoldDB models) | size, compactness, solvent accessibility, secondary-structure content, pLDDT confidence |
| **Graph** (residue interaction networks) | degree, centrality, hubness, clustering, path length, entropy |
| **Thermostability** | interaction-energy–based stability descriptors |
| **Expression** | tissue abundance and tissue specificity (τ) |

In the web interface you can:

- query proteins by UniProt accession or gene name
- select and filter features, and view distributions, correlations, PCA/UMAP and per-residue profiles
- download feature tables (CSV) and per-residue profiles
- **upload your own PDB structures** (up to 10) and compute the same descriptors on demand. Expression features are not computed for uploads. Results are e-mailed as a download link.

---

## Repository layout

```
app.py                     Dash application (entry point, `app:server` for gunicorn)
api_v6.py                  analysis & plotting functions
data.py                    loaders for the precomputed feature tables
descriptions.py            feature / analysis descriptions shown in the UI
assets/                    CSS and JS served by Dash
content/                   About, Documentation, Tutorial, Feature description pages
code/                      feature calculator for uploaded structures
  FeatureCalculatorPcAfGraphs.py
  catgranuleFunctions.py
  src/                     physico-chemical scales, trained models, helper files
thermometer/               thermostability pipeline (NAMD minimisation + R scripts)
combine_thermo_features.py merges thermostability into the job results
package_results.py         packages job results into the downloadable zip
dash.txt / catg.txt / thermo.txt   conda environment specs (linux-64)
```

---

## Installation

The full app (including on-demand feature calculation) needs **Linux** and **three conda environments**:

| Env | Used for | Spec file |
|---|---|---|
| `profet_db_dash` | the web app | `dash.txt` + `dash_pip_requirements.txt` |
| `catgranule2` | physico-chemical / structural / graph features | `catg.txt` |
| `thermo` | thermostability (Thermometer, R) | `thermo.txt` |

```bash
git clone https://github.com/tartaglialabIIT/phanes_pub.git
cd phanes_pub

conda create -n profet_db_dash --file dash.txt -y
conda activate profet_db_dash
pip install -r dash_pip_requirements.txt

conda create -n catgranule2 --file catg.txt -y
conda create -n thermo      --file thermo.txt -y
```

---

## Running

**Development**

```bash
conda activate profet_db_dash
python app.py
```

Then open <http://localhost:8080/phanes/>.

## Citation

If you use PHANES, please cite:

```
https://doi.org/10.64898/2026.10.14.756639
```

## License

Code: see [LICENSE](LICENSE). Data: CC-BY-NC 4.0.
