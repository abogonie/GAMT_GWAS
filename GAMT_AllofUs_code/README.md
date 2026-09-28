# GAMT deficiency carrier frequency in the All of Us Research Program

Code for: Bogoniewski, A.A.; Lipshutz, G.S. Carrier Frequency and Predicted Birth Prevalence of Guanidinoacetate N-Methyltransferase Deficiency in 535,606 Participants from the All of Us Research Program.

Python 3.10.7; pandas 2.3.3; NumPy 2.2.6; SciPy 1.15.3; matplotlib 3.11.1.

## Files

| Script | Runs in | Produces |
|---|---|---|
| 01_allofus_workbench_analysis.py | All of Us Researcher Workbench (Controlled Tier, CDR v9) | data/GAMT_summary_export_safe.csv, data/GAMT_variants_export_safe.csv, Table 3 counts |
| 02_compare_gnomad.py | Local | results/gnomad_vs_allofus_summary.csv, gnomad_common_variant_comparison.csv, gnomad_qualifying_variants.csv |
| 03_pool_and_test.py | Local | results/pooled_allofus_gnomad.csv (Table 4 p values and pooled estimates) |
| 04–06_figure*.py | Local | figures/Figure1–3 (PNG and PDF) |

## Data

- `data/GAMT_summary_export_safe.csv` and `data/GAMT_variants_export_safe.csv`: aggregate All of Us results exported under the All of Us Data and Statistics Dissemination Policy (counts of 1–19 masked). Individual-level All of Us data are available only to registered researchers through the Researcher Workbench.
- `data/gnomAD_v4.1.1_ENSG00000130005_*.csv`: GAMT variants exported from the gnomAD v4.1.1 browser (https://gnomad.broadinstitute.org).

## Run

Step 01 must be run inside the All of Us Researcher Workbench. Steps 02–06 run locally from this folder:

```
pip install -r requirements.txt
python 02_compare_gnomad.py
python 03_pool_and_test.py
python 04_figure1_carrier_frequency.py
python 05_figure2_top_variants.py
python 06_figure3_variant_positions.py
```
