# Reproducibility scope and version notes

## What is verified in this release

- Every one of the 76,893 processed data rows matches the final October 22 submission workbook. Table 11 correctly marks the retained Figure 3B as shared-set enrichment only.
- The offline validator checks Q1–Q10 counts, Q2/Q9 DWPC sums, Q5 Fisher tests, Q6 observed and null gene counts, empirical P values and BH adjustment, and saved Q1–Q3 Enrichr overlaps and adjusted P values.
- Q7 counts are independently recovered from distinct nonblank disease identifiers for the five refreshed compounds. Diethylstilbestrol's historical score is retained but lacks a refreshed route export.
- The portable DWPC implementation is compared with the original notebook's full matrix products on a multigraph with parallel edges and a self-loop. Repeated-node walks, schema enumeration, path caps, and PPR stationary probabilities are tested.
- The portable plot command renders 14 numerical panels. These are table-based reconstructions, not claims of identical manuscript typography or manual composite layout.

## What requires external inputs

SPOKE access is required to repeat the graph extraction. Two locally preserved graph files have recorded SHA-256 hashes, but the live development Neo4j database has no immutable snapshot identifier. The graph files are not included. Full graph runs were not repeated during code packaging.

PPR filtering differs across questions. Q4's gene ranking mapped numeric PPR node identifiers back to the DWPC graph's typed `Gene:` identifiers. The October rank check reproduced all 14 overlap genes and 42 reported ranks within the available candidate mapping; the full historical 19,507-gene filter was not reconstructed. Table 16 records the reported overlap and ranks. The portable generic PPR command computes node scores; it does not silently substitute a different Q4 candidate universe.

The historical Q5 workflow includes multisource PPR, BFS restriction, degree-matched sampling, overlap tests and threshold sensitivity. Its selected notebook cells are preserved, but a new end-to-end rerun requires the original graph and notebook context. The saved results are in Tables 21–23. Q5 Fisher tests are replayed offline; the full degree-matched permutation experiment is not rerun by that command.

## Randomization records

Q1 saved the null-yield distribution and gene-level summaries but did not retain all random-pair identities in the supplied records. Q7 saved 1,000 random compound scores but not compound identities or a sampling seed. Their exact historical draws cannot be reconstructed from this package. A new random run is a new analysis, not an exact replay.

Q6 retains all 1,500 sampled ordered pairs, per-pair gene counts and real paths. The same ordered extraction is applied to real and null pairs, with a cap of 50 paths per selected schema. Sampling was with replacement among predicate-compatible disease pairs; it was not degree matched. SD in Table 25 is the population SD (`ddof=0`). The portable Q6 command reuses the saved pair identities and selected schemas. Neighbor order still matters if a graph is rebuilt.

## Enrichment

Frozen Enrichr responses are supplied for the updated Q1 55-gene set, Q2 full 165-gene PD set, and Q3 15-gene route set. Other reported enrichment results are preserved in the full tables, with the historical source code where recovered. Enrichr's reported adjusted P values apply to returned overlapping terms. The full-library BH values in some tables are separately labelled sensitivity calculations; they are not substituted for the reported values.

Different libraries and gene universes should not be treated as one pooled statistical test. The Q3 display combines selected biological themes across libraries, with library identity encoded by color; adjusted P values remain specific to each enrichment run. Term-name annotation in Q1/Q2 is a documented descriptive grouping, not a new enrichment test.

## Figure versions

Historical plotting sources disclose the original rendering logic but can refer to earlier selections, private paths, intermediate files, or an active Cytoscape session. They are not executed by the release validation or portable figure command.

- **Figure 2B:** 34 network genes are those on the selected gene-regulation schemas. This is a mechanistic subset of the Q1 extraction, not an arbitrary selection of 34 labels from 55 genes. Typed edge data are in Table 8.
- **Figure 3B/C:** The retained B panel shows shared-27 enrichment. Full-PD mitophagy enrichment is reported in Table 11. C includes the six PD-only mitophagy genes alongside the 27 shared genes. The archived October 8 full-layout script also contains an alternative PD enrichment display; it is a historical revision, not the retained B panel.
- **Figure 4A:** The current symptom set excludes compassion fatigue and zero-scoring sparse hair. The axis is log10(PageRank + 1e-10), not a P value.
- **Figure 5A:** The historical graphic's dementia marker used the score for ALS-PDC1 (`DOID:0111246`). The portable reconstruction uses generic dementia (`DOID:1307`). It plots a neutral background histogram rather than reusing the historical category label as a disease identifier. The corresponding Table 19 scores are unchanged.
- **Figure 5B:** This is a metapath/DWPC analysis, not enrichment of the 14 Q4 PPR-overlap genes. Those overlap enrichment results are in Table 17. The focused SPP1 network is a separate typed relationship extraction (Table 20).
- **Figure 6/Q7–Q8:** Use current verified target and route tables; earlier network images may contain superseded annotations. Q7 CGPS measures annotated graph connectivity, not clinical cancer risk.
- **Figure 8:** The reconstructed bars use DWPC and downweight hubs. Historical images labelled “Embeddings” or claiming hub upweighting should not be treated as the current method.

No manuscript document, reviewer comment history, passwords, raw graph pickle, or full exploratory notebook is part of this release.
