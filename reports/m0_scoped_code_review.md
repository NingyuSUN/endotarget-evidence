# M0 scoped rereview
Reviewed fix commit: 0290d997e73e9e5ec6e8e4f1be4d116d4b8ccd63
Previous head: 6c8f975ebbb35bdd23506cdf39888af3333c5ebf
Scope: the two previous findings and regressions introduced by their fix. No tests rerun, downloads, code changes, subagents, scientific analysis or host mutations.

## Finding 1: NOT ADDRESSED completely (P2)
The original blocking body-read issue is addressed by SIGALRM, and saved socketpair red/green evidence supports this. However, the claimed overall deadline can still be swallowed during address connection retries.

Location: scripts/m0_metadata.py:108-109.
The signal handler raises built-in TimeoutError, which belongs to the OSError hierarchy. Actual SRV01 Python 3.10 /usr/lib/python3.10/socket.py:824-841 loops through resolved addresses and catches socket.error/OSError around socket.connect. If the one-shot alarm fires during that operation, the exception can be caught and another address attempted. The timer has already expired, leaving subsequent connection and body reads unprotected by the overall alarm. This is a static control-flow finding; no adversarial rerun was performed.

The slow-opener regression replaces OpenerDirector.open with sleep and therefore does not exercise this standard-library catch. Minimal correction: have the deadline raise a dedicated exception outside OSError, for example a RuntimeError subclass, so network retry code cannot consume it; adapt the deadline tests and add one narrow catch/retry regression. Preserve the existing Unix main-thread limitation.

## Finding 2: ADDRESSED
tests/test_m0_metadata.py:11-14 defines absolute project-local SCRATCH and creates it in setUpModule; all temporary directories use it. Saved old clean-snapshot log demonstrates six FileNotFoundError errors, and fixed clean-snapshot log demonstrates 24 passing tests. The actual script/test and four log hashes match receipts/m0_metadata_review_evidence.json.

## Other scoped observations
No other load-bearing regression identified in the fix diff. Scientific donor rules, thresholds and readiness are unchanged.
The uncommitted Docker diagnostic receipt and runtime/plan clarification accurately distinguish a password-required noninteractive sudo probe from privileged execution or host modification. No container success is inferred.

Quality verdict: changes still requested for the one deadline exception-path defect.
Compliance verdict: scoped work and evidence preserved; recorded tests are valid for their covered paths but do not prove the untested connection-retry deadline.
Scientific verdict: unchanged OPEN/NOT_VERIFIED.
Understanding label: YELLOW, not assessed; socket inactivity timeout, absolute deadline, and exceptions swallowed by retries are distinct issues.

Only this ignored review artifact was written. No tracked change or review commit was made; parent changes were left untouched.
