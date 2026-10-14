# Repository working rules

- Use public data only. Keep donor-aware study design, source provenance and scientific uncertainty explicit.
- Read `docs/METHOD_CARD.md`, `docs/LOCAL_EXECUTION.md` and `docs/PUBLICATION_HISTORY.md` before changes. Frozen historical protocols describe their recorded versions, not a pending execution queue.
- Preserve the three scientific freeze receipts and all referenced artifact bytes. Never silently rerun or overwrite frozen analyses, outputs or checkpoints.
- Local execution uses `/mnt/g/projects/te`; historical scripts retain `/home/analyst/te` through the documented runtime mount. Do not alter unrelated projects or host security settings.
- Plan measured disk and memory requirements before large downloads. Do not commit raw data, credentials, environments or caches.
- Use only `gpt-6-astra` and `gpt-6.1-sol` for delegated work, with no other-model fallback.
- Write concise English commit messages describing technical changes. Preserve partial or failed status when applicable. Commit only intentional changes; push or rewrite history only with explicit user authorization.
- Keep project documentation focused on methods, results, limitations and reproducibility. Keep internal planning and progress notes out of published project documentation.
