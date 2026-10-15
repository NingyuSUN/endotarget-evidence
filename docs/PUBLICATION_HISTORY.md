# Source and execution provenance

The published repository retains the original scientific results and the evidence
needed to distinguish recorded analyses from later documentation changes. This page
records every change made to the published history and what each change means for
the verification records.

## Commit identifiers

The published history has been rewritten twice. Both rewrites kept the commit
sequence, authors, committers, timestamps and scientific content. The complete
original history is kept by the maintainer as a verified local Git bundle.

- **14 October 2026.** The descriptions of 45 commits were standardized in English.
  Only commit messages and parent identifiers changed; every Git tree was preserved.
- **15 October 2026.** Content not intended for publication was removed from every
  commit (see below). Git trees therefore changed, and every commit identifier
  changed again.

The [commit map](../receipts/publication_commit_map.tsv) links each identifier in
historical execution records to its corresponding published commit. The
`original_tree` column is the Git tree of the original, unpublished version. A
mapped commit is the published counterpart of that version, not evidence that a new
analysis was run under its new identifier. Older published commit URLs no longer
resolve. Use the mapped identifier to inspect a historical file at
`https://github.com/NingyuSUN/endotarget-evidence/blob/<published_commit>/<path>`.

## 15 October 2026 history revision

Applied identically to every commit:

- Internal planning, progress, working-state and server-access documents, and three
  host permission probe receipts, were removed. Earlier versions of the research-brief
  PDF were removed; the current PDF was rendered again from the current brief.
- The analysis host name, account name and project root were replaced with the
  neutral placeholders `SRV01`, `analyst` and `/home/analyst/te` in every file,
  file name and commit message.
- Links to removed documents were unwrapped, and two figure files were renamed to
  `scripts/render_summary_figures_v1.py` and `receipts/summary_figures_v1.json`.

**Consequence for hash verification.** SHA-256 values in receipts were not edited.
For files whose published bytes changed in this revision — including scripts and
method documents referenced by the scientific freezes — the recorded hashes refer to
the original bytes and cannot be checked against this public history. The original
bytes are retained in the maintainer's local bundle. The
[revision verification record](../receipts/publication_history_revision_2026-10-15.json)
lists, per commit, which files changed, and confirms that each frozen file differs
from its original only by the placeholder substitution and that no other file was
altered. Scientific numbers, candidate lists and figures were not changed.

## Current documentation

The repository homepages, method card, reproduction guide, evidence index and
research brief present the completed study.

- Learning-status annotations were rewritten as neutral interpretation notes in
  editable documentation. The nine documents recorded by the scientific freezes or
  the M0 gate decision, and the generated reports, keep their original wording.
- [results_summary_en.md](../reports/results_summary_en.md) adds an English summary
  of the generated bilingual replication report, including two interpretation notes
  derived from the reported numbers. It contains no new analysis.
- The generated replication reports list a removed working-state note as a build
  input; their receipts keep its recorded hash.

Some frozen documents include status statements and relative links appropriate to
their original version; they are not current execution instructions. The
[documentation index](README.md) links to the current research and reproduction
entry points.

## Verification scope

Historical receipts retain their original source identifiers and hashes. In
particular, `local_delivery_verification.json` and the prior PDF receipt describe
the files present at their recorded time; they do not verify the edited current
homepages or research brief. The current PDF has a separate
[rendering record](../receipts/publication_brief_pdf.json).

The original `verify_local_delivery.py` also checks original Git ancestry and file
hashes. Neither is preserved in the published history, so running that checker on
a published clone is not a supported validation. Its original checks have not been
relaxed or represented as passing. Reproducing the historical check requires the
original Git lineage and its complete artifact set, which the maintainer retains
locally.

None of these revisions refit statistical models, altered candidate lists,
regenerated scientific figures or added biological evidence. The scientific
limitations remain as reported.
