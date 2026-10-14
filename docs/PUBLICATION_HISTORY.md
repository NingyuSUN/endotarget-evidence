# Source and execution provenance

The published repository retains the original scientific results and the evidence
needed to distinguish recorded analyses from later documentation changes.

## Commit identifiers

On 14 October 2026, the descriptions of 45 existing commits were standardized in
English. Only commit messages and the resulting parent identifiers changed.
Every corresponding Git tree, author, committer and timestamp was preserved.
The complete original history was backed up locally before the rewrite.

The [commit map](../receipts/publication_commit_map.tsv) links each identifier in
historical execution records to its corresponding published commit. The
`identical_tree` column is the shared Git tree identifier: it covers file contents,
paths and modes at that historical version. A mapped commit is a source-equivalent
version, not evidence that a new analysis was run under its new identifier.

Old commit URLs may no longer resolve. Use the mapped identifier to inspect a
historical file at `https://github.com/NingyuSUN/endotarget-evidence/blob/<published_commit>/<path>`.

## Current documentation

The repository homepages, method card, reproduction guide, evidence index and
research brief present the completed study. Obsolete planning, progress, server
access documents were removed from the current file tree.
Their original versions remain in the mapped Git history and the local backup.

The original protocol, pre-effect decisions, statistical implementation and
evaluation contracts remain unchanged. So do the three primary scientific
freeze receipts and their 95 distinct referenced artifacts. Some frozen documents
include status statements and relative links appropriate to their original
version; those should be read at the corresponding historical commit. They are
not current execution instructions. The [documentation index](README.md) links
to the current research and reproduction entry points.

## Verification scope

Historical receipts retain their original source identifiers and hashes. In
particular, `local_delivery_verification.json` and the prior PDF receipt describe
the files present at their recorded time; they do not verify the edited current
homepages or research brief. The current documentation has a separate
[verification record](../receipts/publication_documentation_check.json), and the
current PDF has a [separate rendering record](../receipts/publication_brief_pdf.json).

The original `verify_local_delivery.py` also checks original Git ancestry.
Rewriting commit messages preserves historical file trees but does not preserve
that ancestry relationship to the new branch. Consequently, running that checker
on the publication branch is not a supported way to validate this documentation
revision. Its original checks have not been relaxed or represented as passing.
Reproducing the historical check requires the original Git lineage and its
complete artifact set, not just the published source clone. The maintainer's
verified Git bundle and archive ref preserve that lineage locally.

This revision did not refit statistical models, alter candidate lists, regenerate
scientific figures, or add biological evidence. It changes presentation and source
navigation; the scientific limitations remain as reported.
