# Phase positive control — discovery fit and F1 freeze (protocol step 3)

Date: 2026-10-27. Executed on the native analysis host (SRV01) in the project's
statistical runtime (R 4.4.3, edgeR 4.4.0, limma 3.62.1, statmod 1.5.2; Python 3.10.22,
NumPy 2.2.6), the same versions as the v1 native run. **No HECA phase model has been
fitted.** Receipt: [`phase_positive_control_f1_discovery_freeze.json`](../receipts/phase_positive_control_f1_discovery_freeze.json);
published copies of the frozen products: [`reports/phase_positive_control_f1/`](phase_positive_control_f1/).

## What was fitted

The 35-donor discovery matrix (22 v1 proliferative + 13 secretory donors; 26,109 stable
genes) with the primary design frozen at the feasibility step: intercept, `SEC_minus_PRO`,
`ENDO_minus_non_ENDO`, centered age, eight capture-batch contrasts, adenomyosis
(rank 13, residual df 22). The estimator is the v1 configuration unchanged
(`filterByExpr` min.count 10 / min.total.count 15 / large.n 10 / min.prop 0.7, TMM, voom span
0.5, least-squares fit, robust empirical Bayes, BH over the full expression-pass family),
with cycle phase as the grouping variable, as the protocol specifies. The driver
([`phase_positive_control_fit_discovery.R`](../scripts/phase_positive_control_fit_discovery.R))
reproduces the v1 library step by step; the F0 freeze was re-verified before fitting
(every data and configuration artifact unchanged; the one documentation artifact that
differs on the analysis host, `docs/REPLICATION_EVALUATION.md`, holds pre-revision bytes
there and matches the F0 hash on this published checkout).

## Results frozen at F1

| Quantity | Value |
|---|---|
| Expression-pass family (tested, all estimable) | 17,341 genes (CPM cutoff 0.32 in ≥12.1 donors; median library 31.2 M) |
| **BH q < 0.05** | **7,367 genes** (3,611 higher in secretory, 3,756 lower); 5,182 at q < 0.01 |
| P1, discovery component | threshold ≥100: **met**. P1 itself needs the HECA component, not yet fitted |
| Ranking universe | 11,650 of 17,341 expression-pass genes carry a frozen-network edge; 122,355 induced edges of 159,312; 5,691 isolated genes excluded, as the v1 rule requires |
| G1 propagation | converged in 28 iterations (L1 change 5.3 × 10⁻¹¹) |
| Fit warnings | none |

**B0 top-20 (expression ranking, all 20 at q < 0.05; so are all of the top 50):** BCAT1,
CFD, ISOC1 (−), EVA1C, MAOA, PDLIM1, APOD, ZBTB16, TRIB2 (−), SCARA5, ENOX2, MTHFD2, MYOCD,
CRYAB, ZEB1, MASP1, IL1R1, FMO1 (−), DKK1, ALDH1A3. Signs are secretory minus
proliferative. Several of these are long-established mid-secretory stromal genes in the
cited literature (for example MAOA, APOD, DKK1 and CRYAB); the list is reported as frozen,
not curated.

**B1 top-20 (network degree):** H4C6, TP53, MAPT, ANLN, UBC, BRD4, ACTB, H3C1, BRCA1, EGFR,
EP300, UBB, H2AC17, HDAC1, H4C9, ACTG1, H4C14, H4C12, H4C8, H4C2 — 10 of 20 at q < 0.05.
**G1 top-20 (network propagation):** TRIM28, TP53, MAPT, ANLN, NUDCD1, GNB1, UBC, GNG2,
EGFR, CUL3, H4C6, CTNNB1, GNAI1, ACTB, HDAC1, BRCA1, BRD4, CALM3, EP300, YWHAB — 9 of 20 at
q < 0.05. As in v1, both network rankings are dominated by high-degree genes rather than
by the strongest phase effects; whether that transfers is exactly what the HECA stage
measures.

## Marker panel (S1), discovery component

| Gene | Expected | In family | log2 FC | BH q | Sign | q < 0.05 |
|---|---|---|---:|---:|---|---|
| IGFBP1 | + | **no** | — | — | — | — |
| PRL | + | yes | +2.23 | 6.5 × 10⁻³ | agrees | yes |
| DKK1 | + | yes | +6.37 | 4.1 × 10⁻¹¹ | agrees | yes |
| LEFTY2 | + | yes | +1.89 | 0.43 | agrees | **no** |
| PLCL1 | + | yes | +4.75 | 1.0 × 10⁻¹⁰ | agrees | yes |
| ESR1 | − | yes | −1.51 | 2.1 × 10⁻⁶ | agrees | yes |
| PGR | − | yes | −1.25 | 7.3 × 10⁻⁹ | agrees | yes |

Six of six fitted genes have the expected sign; five of seven meet the prespecified
discovery criterion (expected sign **and** q < 0.05). The threshold is six, so **S1 fails
on the discovery component.** The two shortfalls have an identifiable cause, recorded here
without changing the rule: IGFBP1 was removed by the expression filter — its counts are
essentially confined to the two late-secretory controls (sample104, sample131) and a few
mid-secretory donors, far fewer than the ≥12.1 donors the group-size-based filter requires
in a pooled secretory arm of 13 — and LEFTY2, likewise a late-secretory gene, is positive
but imprecise. Both are late-secretory markers in an arm that is mostly early and mid
secretory (4 early, 7 mid, 2 late after the clinical-flag exclusion). Per the protocol, S1
is explanatory, not decisive; its thresholds are not revisited after effects.

## What is frozen, and what remains

Frozen at F1: the full `SEC_minus_PRO` table (17,341 rows), the expression filter, the
normalization factors, the design, the B0/B1/G1 full rankings with top-20/top-50 lists, the
saved discovery signs, the ranking universe and induced graph, and the marker readout.
Their hashes and the commit that holds the code are in the receipt. Two published files
had the analysis host's project root replaced by the repository placeholder after the
freeze (the ranking provenance records absolute paths by v1 design, and the receipt's
software record names the interpreter); the receipt records the hashes before and after
that substitution, following the repository's published-history convention.

Not done, by protocol order: the HECA `SEC_minus_PRO` fit (step 5) and the P2/S1-HECA/S2
evaluations, which may now proceed because F1 exists. Nothing here is disease evidence;
a secretory-phase gene is not an endometriosis gene, and the B0 list is a pipeline
positive-control input, not a target list.
