# Analysis map

All table numbers below refer to the submission workbook. The notes in each table are part of the methods. `export-tables` preserves them in its index.

| Question | Supplementary tables | Source |
|---|---|---|
| Methods | 1–3 | General design, identifiers and enrichment rules |
| Q1 | 4–8 | Premature menopause–dementia DWPC, 55 intermediate genes and regulation network |
| Q2 | 9–12 | Three neurodegenerative targets, shared intermediates and separate PD enrichment |
| Q3 | 13–15 | Symptom PPR, complete symptom–gene–disease routes and enrichment |
| Q4 | 16–20 | Bone–brain top-200 gene overlap, DWPC, PPR and focused SPP1 network |
| Q5 | 21–23 | APO, menopause and cognition PPR sets, overlaps and enrichment |
| Q6 | 24–27 | Delirium metapaths, observed gene sets and matched extraction null |
| Q7 | 28–31 | Compound target/cancer connectivity and mechanistic networks |
| Q8 | 32–36 | Breast-cancer therapies: PPR, typed dementia routes, gene overlap and enrichment |
| Q9 | 37 | Estradiol/estrone DWPC by route category |
| Q10 | 38–40 | MHT ingredients, queried pharmacogenes and typed relationships |

## Recorded selection rules

### Q1: DWPC and intermediate genes

- **Source or seed nodes:** Premature menopause
- **Target nodes:** Dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** All enumerated schemas; gene intermediates from selected schemas. Maximum length: 3.
- **Selection and parameters:** 346 schemas; 11 active; up to 50 paths/schema for 55 genes; damping −0.4
- **Statistical comparison:** 1,500 compatible random disease pairs; seed 42; empirical P with +1 correction and BH
- **Interpretation:** Matrix DWPC counts walks; gene lists use capped path enumeration.

### Q1: Figure 2B network

- **Source or seed nodes:** Premature menopause
- **Target nodes:** Dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** Four gene-regulation schemas: GPuG, KGuG, GPdG, KGdG with DaG flanks. Maximum length: 3.
- **Selection and parameters:** Up to 50 paths/schema; 34 displayed genes
- **Statistical comparison:** No separate null comparison reported for the displayed network
- **Interpretation:**

### Q2: DWPC and shared genes

- **Source or seed nodes:** Premature menopause
- **Target nodes:** Parkinson’s disease, vascular dementia and frontotemporal dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** All enumerated schemas; gene intermediates from selected schemas. Maximum length: 3.
- **Selection and parameters:** Top five schemas per target; up to 200 paths/schema; 242 union genes, 27 shared; damping −0.4
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Figure 3C shows all 27 shared genes plus the six PD-only contributors to GO:1903146 from separate enrichment of all 165 saved PD intermediates.

### Q2: Disease proximity

- **Source or seed nodes:** Premature menopause
- **Target nodes:** 2,422 connected diseases
- **Knowledge graph:** SPOKE Neo4j queries
- **Relationships and maximum path length:** Disease–ASSOCIATES_DaG–Gene–any edge–Gene–ASSOCIATES_DaG–Disease. Maximum length: 3.
- **Selection and parameters:** Count distinct second-gene intermediates per target disease
- **Statistical comparison:** No null comparison reported
- **Interpretation:** This analysis is a shared-gene count, not PageRank.

### Q3: Symptom proximity

- **Source or seed nodes:** Dementia
- **Target nodes:** 1,891 symptom nodes
- **Knowledge graph:** Human-focused SPOKE graph
- **Relationships and maximum path length:** Personalized PageRank over graph edges. Maximum length: Global.
- **Selection and parameters:** alpha 0.85; max_iter 500; tol 1e-20; rank among symptoms
- **Statistical comparison:** Rank/N only; no degree-controlled null
- **Interpretation:**

### Q3: Direct symptom routes

- **Source or seed nodes:** Insomnia; paresthesia
- **Target nodes:** Dementia-related disease nodes in saved export
- **Knowledge graph:** SPOKE Neo4j queries
- **Relationships and maximum path length:** Symptom–ASSOCIATES_GaS–Gene–ASSOCIATES_DaG–Disease. Maximum length: 2.
- **Selection and parameters:** All routes in the network export: 15 genes, 36 routes and 50 unique edges
- **Statistical comparison:** No null comparison reported
- **Interpretation:** All genes in the route export are included.

### Q4: Disease and gene proximity

- **Source or seed nodes:** Premature menopause; osteoporosis; dementia
- **Target nodes:** Disease rankings; condition-specific gene rankings
- **Knowledge graph:** Human-focused SPOKE graph
- **Relationships and maximum path length:** Personalized PageRank over graph edges. Maximum length: Global.
- **Selection and parameters:** alpha 0.85; max_iter 500; tol 1e-20; top 200 genes per condition
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Gene rankings used the 19,507-gene identifier universe of the DWPC subgraph; the condition-specific PageRank scores came from the human-focused graph.

### Q4: Osteoporosis–dementia DWPC

- **Source or seed nodes:** Osteoporosis
- **Target nodes:** Dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** All enumerated schemas. Maximum length: 3.
- **Selection and parameters:** 346 schemas, 11 active; damping −0.4
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Independent of the 14-gene enrichment; active schemas contain no Pathway node.

### Q4: Focused SPP1 network

- **Source or seed nodes:** SPP1; osteopontin
- **Target nodes:** Premature menopause; osteoporosis; selected dementia-related partners
- **Knowledge graph:** SPOKE Neo4j queries
- **Relationships and maximum path length:** Typed SPOKE relationships in exported edge table. Maximum length: Query-specific.
- **Selection and parameters:** Selected partners include IL6, IL1B, TNF, TREM2 and APOE
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Focused export; not an unbiased centrality screen or representation of all 14 overlap genes.

### Q5: Multiseed PageRank and overlap

- **Source or seed nodes:** 7 APO, 16 menopause-symptom and 4 cognitive seeds
- **Target nodes:** Gene nodes
- **Knowledge graph:** Human-focused SPOKE graph
- **Relationships and maximum path length:** Personalized PageRank; candidate genes within two BFS hops. Maximum length: Global PPR; 2-hop candidates.
- **Selection and parameters:** Uniform seed weights; alpha 0.85; tol 1e-10; max_iter 500; top 2% rounded up
- **Statistical comparison:** 1,000 degree-matched permutations; 20 log1p-degree bins; 18,753 connected-gene universe; seed 42
- **Interpretation:** The complete seed sets, including mapped nodes with degree zero, are listed in Supplementary Table 2.

### Q6: DWPC and path-gene overlap

- **Source or seed nodes:** Premature menopause
- **Target nodes:** Subacute delirium; dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** All schemas; genes from selected schemas. Maximum length: 3.
- **Selection and parameters:** Top five schemas; up to 200 paths/schema: 68 and 37 genes; damping −0.4
- **Statistical comparison:** No matched null for this set comparison
- **Interpretation:** The top-five-schema sets contain 68 premature-menopause–delirium genes and 37 premature-menopause–dementia genes; the latter is a subset of the Q1 55-gene set.

### Q6: Delirium–dementia intermediate genes

- **Source or seed nodes:** Subacute delirium
- **Target nodes:** Dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** 11 selected schemas. Maximum length: 3.
- **Selection and parameters:** Up to 50 paths/schema; 36 genes
- **Statistical comparison:** 1,500 predicate-compatible random disease pairs; seed 42; identical ordered extraction for real/null; no expansion cap
- **Interpretation:** The null conditions on the 11 observed schemas and does not match disease-node degree.

### Q7: Cancer Gene Proximity Score

- **Source or seed nodes:** Fezolinetant, elinzanetant, tamoxifen, diethylstilbestrol, letrozole and goserelin
- **Target nodes:** 1,884 cancer disease nodes
- **Knowledge graph:** SPOKE Neo4j queries
- **Relationships and maximum path length:** Compound–BINDS_CbP–Protein–ENCODES_GeP–Gene–cancer; cancer edges DaG, GeiD, GmpD or GmnD. Maximum length: 3.
- **Selection and parameters:** Count distinct cancer diseases; binding and encoding retain direction, gene–disease edges traversed either way
- **Statistical comparison:** 1,000 eligible random compounds with ≥1 binding edge; 868 nonzero
- **Interpretation:** The random-compound identities and random seed are unavailable. CGPS measures graph connectivity and does not estimate clinical risk.

### Q7: Neurokinin-antagonist network exports

- **Source or seed nodes:** Fezolinetant; elinzanetant
- **Target nodes:** Dementia and breast-cancer pathway WP4262
- **Knowledge graph:** SPOKE Neo4j queries
- **Relationships and maximum path length:** Binding, encoding and typed mechanistic and protein-interaction relationships. Maximum length: Query-specific.
- **Selection and parameters:** Selected mechanistic and PPI neighborhoods
- **Statistical comparison:** No null comparison reported
- **Interpretation:** The network comprises selected neighborhoods; it is not an exhaustive target screen.

### Q8: Compound PageRank

- **Source or seed nodes:** Dementia
- **Target nodes:** Tamoxifen; leuprolide; letrozole; goserelin
- **Knowledge graph:** Human-focused SPOKE graph
- **Relationships and maximum path length:** Personalized PageRank over graph edges. Maximum length: Global.
- **Selection and parameters:** alpha 0.85; max_iter 500; tol 1e-20; ranks among 2,432,054 compounds
- **Statistical comparison:** Descriptive compound-rank background
- **Interpretation:**

### Q8: Drug–gene–dementia paths

- **Source or seed nodes:** Tamoxifen, letrozole, goserelin and leuprolide
- **Target nodes:** Dementia
- **Knowledge graph:** SPOKE Neo4j queries
- **Relationships and maximum path length:** Compound–[CdG/CuG]–Gene–[GPdG/KGdG/KGuG/GPuG]–Gene–DaG–Disease. Maximum length: 3.
- **Selection and parameters:** Identical query on 1 October 2026: tamoxifen 927, letrozole 174, goserelin 237 unique genes; leuprolide zero routes
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Genes are deduplicated over both intermediate positions. A zero route count applies only to the specified schema.

### Q8: Longer letrozole schemas

- **Source or seed nodes:** Letrozole
- **Target nodes:** Dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** All enumerated schemas. Maximum length: 4.
- **Selection and parameters:** 7,511 schemas through length four; damping −0.4
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Lengths three and four are reported separately.

### Q9: MHT DWPC

- **Source or seed nodes:** Estradiol and estrone
- **Target nodes:** Alzheimer’s disease, dementia and vascular dementia
- **Knowledge graph:** SPOKE subgraph for DWPC
- **Relationships and maximum path length:** All schemas; route class assigned by saved notebook rule. Maximum length: 3.
- **Selection and parameters:** 3,006 compound–target–schema rows; damping −0.4
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Estrone represents one component of conjugated estrogens, not the mixture. Connectivity does not establish a protective effect.

### Q10: Pharmacogene query

- **Source or seed nodes:** Seven estrogen/SERM ingredients
- **Target nodes:** 14 queried pharmacogenes
- **Knowledge graph:** SPOKE Neo4j queries
- **Relationships and maximum path length:** Binding/interacting protein links through ENCODES plus gene up/downregulation. Maximum length: Direct gene or encoded-protein link.
- **Selection and parameters:** Completed seven-ingredient run; progestins were not included
- **Statistical comparison:** No null comparison reported
- **Interpretation:** Output matrix has eight genes and six connected ingredients; all 14 queried genes are tabulated.

## Source conventions

Notebook source filenames use zero-based cell positions, not notebook execution counts. Sanitization removes outputs, credentials, private locations and IPython shell/magic statements; it also normalizes Python formatting. See `SOURCE_MAP.csv` for hashes before and after extraction.

Historical sources may require cached CSVs, graph objects, functions from other cells and Python packages beyond the portable CLI. They are included to disclose the analysis logic, not advertised as turnkey scripts. The Q6 portable replay is an explicitly adapted entry point. New CLI implementations live under `src/` and are tested separately.

The original full graph computations were not rerun during release packaging. Offline validation checks agreement among frozen inputs, tables and reported statistics. It does not validate the biological conclusions independently.
