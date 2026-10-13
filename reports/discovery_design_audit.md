# Discovery metadata-only design audit

These are candidate designs, not frozen models. No expression values or disease effects were read.

| Candidate set | Terms in addition to intercept and disease | n | Rank / columns | Residual df | Structural gate |
|---|---|---:|---:|---:|---|
| author_adjudicated_phase22 | age, capture | 22 | 11 / 11 | 11 | True |
| author_adjudicated_phase22 | age, capture, adenomyosis | 22 | 12 / 12 | 10 | True |
| author_adjudicated_phase22 | age, sequencing | 22 | 4 / 4 | 18 | True |
| author_adjudicated_phase22 | age, sequencing, adenomyosis | 22 | 5 / 5 | 17 | True |
| author_adjudicated_phase22 | age, capture, sequencing | 22 | 11 / 12 | 11 | False |
| author_adjudicated_phase22 | age | 22 | 3 / 3 | 19 | True |
| strict_histology_concordance20 | age, capture | 20 | 11 / 11 | 9 | True |
| strict_histology_concordance20 | age, capture, adenomyosis | 20 | 12 / 12 | 8 | True |
| strict_histology_concordance20 | age, sequencing | 20 | 4 / 4 | 16 | True |
| strict_histology_concordance20 | age, sequencing, adenomyosis | 20 | 5 / 5 | 15 | True |
| strict_histology_concordance20 | age, capture, sequencing | 20 | 11 / 12 | 9 | False |
| strict_histology_concordance20 | age | 20 | 3 / 3 | 17 | True |
| author_phase_adenomyosis_free18 | age, capture | 18 | 11 / 11 | 7 | True |
| author_phase_adenomyosis_free18 | age, capture, adenomyosis | 18 | 11 / 12 | 7 | False |
| author_phase_adenomyosis_free18 | age, sequencing | 18 | 4 / 4 | 14 | True |
| author_phase_adenomyosis_free18 | age, sequencing, adenomyosis | 18 | 4 / 5 | 14 | False |
| author_phase_adenomyosis_free18 | age, capture, sequencing | 18 | 11 / 12 | 7 | False |
| author_phase_adenomyosis_free18 | age | 18 | 3 / 3 | 15 | True |
| strict_phase_adenomyosis_free16 | age, capture | 16 | 11 / 11 | 5 | False |
| strict_phase_adenomyosis_free16 | age, capture, adenomyosis | 16 | 11 / 12 | 5 | False |
| strict_phase_adenomyosis_free16 | age, sequencing | 16 | 4 / 4 | 12 | False |
| strict_phase_adenomyosis_free16 | age, sequencing, adenomyosis | 16 | 4 / 5 | 12 | False |
| strict_phase_adenomyosis_free16 | age, capture, sequencing | 16 | 11 / 12 | 5 | False |
| strict_phase_adenomyosis_free16 | age | 16 | 3 / 3 | 13 | False |

The author-adjudicated phase scenario contains12 ENDO and10 non-ENDO donors. The strict histology-concordance scenario contains11 and9. All22 meet>=50 cells under the externally derived fibroblast definition, which remains subject to our annotation crosswalk. Under author-adjudicated phase, both disease groups coexist in capture batches12,14,16 and18; under strict histology concordance they coexist only in14 and16. The remaining batches contribute nuisance-parameter/variance information; they do not create within-batch disease comparisons. A full-rank design is not evidence that confounding is absent.

Do not drop capture-batch adjustment simply because the smaller model gives more residual degrees of freedom. Sequencing batch is nested in capture batch; adding both produces redundant columns. Adenomyosis-adjusted and adenomyosis-free analyses answer different questions. Under author-adjudicated phase, excluding adenomyosis leaves8 cases and10 controls; structural results are shown above. Under strict histology concordance, it leaves only7 cases and9 controls and fails the fixed8-per-group minimum, even if a numerical regression could run. When adenomyosis is constant after exclusion, its redundant model column must be omitted, with this rule stated before effects.

Machine-readable inputs, cohort IDs, covariate columns, source hashes, structural checks and batch tables are in discovery_design_audit.json. The script checks a known perfectly confounded artificial metadata copy as a negative control; it does not simulate biological results.

YELLOW: understanding not yet assessed. Explain why an identifiable adjusted disease coefficient can still rely on few mixed batches, and why this is neither external validation nor a guarantee of power.
