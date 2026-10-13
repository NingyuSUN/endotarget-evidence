# Task 1 independent-review fixes

VERIFIED: two review findings fixed in one scoped pass on SRV01 under /home/analyst/te. Scientific donor definitions, thresholds, audit outputs and readiness were not changed. No public metadata or expression input was downloaded during this fix.

## Overall deadline
The prior pre-read elapsed-time check did not interrupt a single HTTPResponse.read call. A real socketpair-backed HTTPResponse test supplied 30 body bytes at 10 ms intervals. With TOTAL_TIMEOUT patched to 50 ms, the old implementation returned from its first read only after 0.305326 seconds, exceeding the asserted 200 ms upper bound. A second test demonstrated that a slow opener could complete after 300 ms without raising TimeoutError.

The fetch now arms a Unix ITIMER_REAL/SIGALRM deadline before opener.open, keeps it active through redirects/connection/body reads and streaming fsync, and restores the prior signal handler afterward. The raising signal handler interrupts blocking open/read operations. Standalone stream_to_part also gets a deadline; fetch_one supplies its covering deadline to avoid nesting timers. Existing active alarms are refused instead of overwritten. A timed-out stream removes only its newly created partial file. Receipts now name the setting overall_deadline_seconds.

This is an explicit Unix main-thread CLI implementation, matching SRV01. It is not a promise of portability to Windows or invocation from worker threads. The old receipts preserve their original configured stream_deadline_seconds field; the earlier claim that this was a strictly enforced overall bound was too strong and is superseded by this review fix.

## Clean test execution
The old suite used TemporaryDirectory(dir="work") without creating work. A git archive of the committed HEAD was extracted into work/m0_review_clean_red, with no ignored work directory. The documented command actually failed with six FileNotFoundError errors.

setUpModule now creates an absolute project-local SCRATCH parent derived from the test file's location. All TemporaryDirectory calls use that parent. A second clean tracked-code snapshot (HEAD archive plus the two pending fix files only) was extracted into work/m0_review_clean_green. Its work directory did not exist before the run. The exact same command succeeded and created its own project-local scratch directory.

## Commands and verified results
All operations occurred under /home/analyst/te. No external environment or model was invoked.

- python3 -m unittest discover -s tests -v
  - Before fix with new deadline regressions: 24 tests, two assertion failures, exit 1.
  - After fix: 24 tests, all passed, exit 0.
- Clean snapshots prepared by Python subprocess calling git archive --format=tar HEAD, then tarfile extraction under work/. Green snapshot overlaid only scripts/m0_metadata.py and tests/test_m0_metadata.py.
- subprocess.run(["python3","-m","unittest","discover","-s","tests","-v"], cwd=clean_snapshot)
  - Old clean snapshot: 22 tests, six errors, exit 1.
  - Fixed clean snapshot: 24 tests, all passed, exit 0.

Both deadline regressions configure a 50 ms timer and assert completion in less than 200 ms; the trickling-response case additionally verifies partial-file cleanup. The socketpair is local interprocess communication, not an HTTP connection to a research source.

Logs, exact code hashes, command metadata and full outputs are saved in receipts/m0_metadata_review_evidence.json. Original logs are in logs/m0_metadata_{review_red,review_green,clean_red,clean_green}.log.

## Interpretation and remaining scope
These results establish the two engineering fixes on SRV01. They do not alter or verify biological eligibility, donor mapping or raw counts. Stage 0 remains OPEN and final eligibility NOT_VERIFIED. No additional user decision is needed. YELLOW: understanding not yet assessed; a per-socket inactivity timeout and an overall elapsed-time deadline guard different failure modes.


## Scoped rereview: prevent deadline swallowing by address retries
VERIFIED on SRV01 Python 3.10.12: socket.create_connection catches socket.error (OSError) around socket.connect. Built-in TimeoutError is an OSError subclass, so the one-shot alarm could be caught and the next address tried after the timer had fired.

The minimal correction is DeadlineExceeded(RuntimeError), deliberately outside the OSError hierarchy. The signal handler raises this dedicated exception; existing streaming/opening regression assertions now expect it. Alarm mechanics and the accepted scratch-directory fix are unchanged.

The new regression invokes the actual stdlib socket.create_connection retry loop. Only getaddrinfo and the endpoint socket are replaced: two documentation-range addresses are supplied, the first connect blocks, and a real 50 ms ITIMER_REAL alarm fires. Before the correction, the test observed both address attempts and no escaping exception (25 tests, one failure). After the correction, DeadlineExceeded escapes on the first address within the 200 ms assertion bound; all 25 tests pass.

Exact command for both red and green runs, from /home/analyst/te:
    python3 -m unittest discover -s tests -v
Red log: logs/m0_metadata_alarm_red.log, exit 1.
Green log: logs/m0_metadata_alarm_green.log, exit 0.
Full outputs, implementation hashes, stdlib socket.py hash and command metadata: receipts/m0_metadata_alarm_evidence.json.

No public source fetch or scientific analysis occurred. Original metadata, donor counts and scientific readiness remain unchanged. The Unix main-thread restriction still applies. No new user decision is needed.
