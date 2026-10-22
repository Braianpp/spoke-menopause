# SPOKE menopause and brain aging

Analysis code and processed results accompanying **Accelerating the discovery of biological pathways linking menopause to brain aging using biomedical knowledge graphs**.

This release covers the ten research questions (Q1–Q10), with the data behind all 40 submission supplementary tables. It includes the historical analysis sources, a portable command-line interface, frozen enrichment responses, the complete Q6 random-pair record, and plotting code.

## Start here: reproduce saved results offline

Python 3.10 or later is required. Python 3.12 was used to test this release.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-tested.txt
python -m pip install -e .
spoke-menopause validate
spoke-menopause export-tables
spoke-menopause figures
python -m pytest -q
```

The validation command checks reported counts, replays the Q6 empirical P values and BH correction, recomputes Q5 Fisher tests, and checks saved Q1–Q3 enrichment responses. CSV export creates 62 blocks, with an index giving table numbers, titles, question mappings and notes. Plotting creates 14 numerical panels as PNG and PDF. They reproduce the selected numerical results; the final manuscript composites also contain manual layout and Cytoscape work.

These commands do not contact SPOKE or Enrichr. Internet access is needed to install dependencies. `requirements-tested.txt` records this release's tested environment, not the software versions of every historical analysis.

## Repository contents

| Location | Contents |
|---|---|
| `src/spoke_menopause/` | Tested CLI for offline checks, table exports, numerical panels, graph DWPC, PPR, path extraction and explicit API queries |
| `analysis/q01/` … `analysis/q10/` | Selected historical scripts and notebook cells, organized by manuscript question |
| `analysis/shared/` | Original DWPC/metapath definitions and shared code |
| `analysis/figure_sources/` | Historical figure scripts, including Cytoscape exports and later layout revisions |
| `data/supplementary_tables.json` | Publication-facing content of the final 40-table workbook, with October 22 Figure 3B flags |
| `data/enrichment/` | Frozen Q1, full-PD Q2 and Q3 gene inputs and Enrichr responses |
| `data/q06/` | All 1,500 sampled disease pairs, their gene counts, real paths and schemas |
| `queries/q07_q08.json` | Parameterized queries used for the refreshed Q7–Q8 checks |
| `docs/` | Analysis map, source provenance, graph hashes, verification results and limitations |

**Historical sources are reference material, not an independently executable workflow.** Notebook cells depend on earlier notebook state; legacy scripts use historical directory layouts. Private paths were replaced by descriptive `${SPOKE_*}` placeholders. Those placeholders are not automatically expanded. Use the tested CLI for portable operations. See [the analysis guide](docs/ANALYSIS_GUIDE.md) before adapting historical code.

## Recompute from a SPOKE graph

The full SPOKE graphs and database credentials are not distributed here. Obtain authorized access to the relevant graph separately. Exact file hashes and graph properties are in [GRAPH_INPUTS.json](docs/GRAPH_INPUTS.json). Graph files are trusted Python pickles and must never be loaded from untrusted sources.

For the typed DWPC graph:

```bash
spoke-menopause dwpc --graph /path/to/V1_homosapiens_1N_9E.gpickle \
  --source Disease:DOID:10787 --target Disease:DOID:1307 \
  --max-length 3 --out results/q1_metapaths.csv
```

For the human-focused PPR graph, use that graph's native node identifiers (they differ from the DWPC graph):

```bash
spoke-menopause pagerank --graph /path/to/spoke_subnetwork_human_focused.gpickle \
  --seeds DOID:1307 --out results/dementia_pagerank.csv
```

Check that identifiers are present before running. PPR uses alpha 0.85, at most 500 iterations and tolerance 1e-20. Multiple seeds receive equal personalization weights. The command returns every node; question-specific filtering and ranking are documented in the analysis guide and historical sources. DWPC counts walks, which may revisit nodes. It uses binary predicate-specific adjacency and damping −0.4 at both ends of each edge, downweighting high-degree nodes. Large-graph runs require substantial memory and time; no full-graph runtime benchmark is claimed for the new CLI.

`extract-genes` takes a JSON list of ordered predicate lists, a path cap, source and target. Capped paths depend on graph neighbor order. Q1's historical extractor counts parallel matching edges separately; select `--parallel-edges` for that behavior. Other extractions generally collapse matching parallel edges per neighbor. Follow the exact question's schema selection and cap.

The portable full Q6 replay is available as:

```bash
python analysis/q06/rerun_matched_null.py --graph /path/to/V1_homosapiens_1N_9E.gpickle --out results/q06
```

It uses the saved schemas and sampled pairs, checks observed counts against Table 25, and applies the same extraction to observed and null pairs. This full-graph command has not been rerun as part of release testing; its extraction algorithm is tested against ordered DFS on small graphs, and the saved full-graph output is independently replayed by `validate`.

## Explicit online operations

These commands contact external services only when invoked:

```bash
spoke-menopause enrichr --genes data/enrichment/q01/genes.txt \
  --libraries GO_Biological_Process_2023 --out results/new_enrichment
spoke-menopause neo4j --query-file queries/q07_q08.json --out results/new_queries
```

Set `NEO4J_URI`, `NEO4J_USER`, `NEO4J_PASSWORD`, and optionally `NEO4J_DATABASE` in your shell before the Neo4j command. Do not commit credentials. Live Enrichr/library and development-SPOKE updates can change results. The archived responses and supplementary tables are the record of the reported analysis.

## Scope and reproducibility limits

All 40 tables are included as processed results, but an exact fresh rerun of every historical analysis cannot be guaranteed. The development graph lacked an immutable release identifier; Q1 random-pair identities and Q7 random-compound identities/seed were not saved. Some scripts depend on notebook state or cached PPR exports. The complete original notebooks are excluded because they mix credentials, exploratory work and superseded analyses. Selected sanitized source cells retain provenance hashes.

The portable panels use the current numerical tables. In particular, Figure 3B remains the shared-27 enrichment panel, while the full-PD analysis is in Table 11. The Figure 5A reconstruction locates generic dementia by **DOID:1307**; an older graphic placed its dementia marker at ALS-PDC1's score. See [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the complete limitations and figure-version notes.

## Citation and reuse

The analysis code is released under the MIT license; see [LICENSE](LICENSE).

Use the versioned repository release when citing the code. `CITATION.cff` provides software citation metadata. A Zenodo DOI can be added after the archive is actually published. No DOI has been claimed in this source package.

This repository does not grant rights to the underlying SPOKE database or its upstream sources. Data access and upstream terms remain separate from the license for the analysis code.
