# AWS deployment guide — documented, not executed

Current execution remains SRV01:/home/analyst/te as requested. No AWS account resources,
image registry, bucket, instance, IAM role or billable job were created. This is an optional
portability guide, not evidence of cloud deployment experience or a completed cloud release.

## Target and prerequisites

Use one Linux x86_64 EC2 host for a first replay of the audited donor-level workflow.
Single-host replay is sufficient for these small donor matrices; an AWS Batch port is outside
this first guide. A proposed initial capacity is 8 vCPU,16 GiB RAM and80 GiB attached storage,
to be revised from the measured receipts before purchase. This is not an AWS performance
measurement or price quotation. Reserve at least20 GiB free; cap numerical threads at1
per model and bootstrap workers at4. Keep new raw acquisition outside the replay.

Prerequisites for an actual attempt:
- Explicitly authorize the AWS account, region, spend limit and shutdown plan. None is assumed.
- Finish and verify the SRV01 container build/run first if the intended claim is container portability.
  The current Docker profile is unexecuted and no validated image digest exists yet.
- Export a versioned source commit and the complete referenced local artifact closure. Git alone
  omits raw/count/work artifacts. Include the files hashed by the three freezes and the immutable
  discovery products; also preserve all original output paths needed by provenance receipts.
- Store the transfer bundle privately with a SHA256 manifest. Preserve source access/licensing
  information. Exclude credentials, local environments, unrelated projects and shell histories.
- Retain the absolute project layout /home/analyst/te on the replay host because archived
  provenance and some frozen interfaces reference that path. Another layout needs an explicit
  portability change and a new manifest; do not silently rewrite frozen files.

The AWS documentation describes installing Docker on an EC2 host and granting the configured
account access before verifying the daemon. These are future administrator actions, not commands
executed on SRV01: [AWS container/EC2 setup](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/create-container-image.html).

## Replay sequence

1. Select the EC2 host and private storage according to the approved budget. Permit only the
   necessary administration access; use a scoped instance role for the artifact location.
   Record region, AMI ID, architecture, vCPU/RAM/storage and launch/termination timestamps.
2. Restore the source and artifact bundle into the preserved project layout. Recreate the
   statistical runtime from requirements/statistical-runtime.explicit.txt using a verified
   package manager; do not copy a machine-dependent conda prefix blindly. Use the pinned
   Nextflow25.10.2 and Java21 runtime. Record actual package versions and the new machine.
3. Verify the supplied checkpoints before any computation:

       python scripts/freeze_analysis.py --verify receipts/statistical_pre_effect_freeze.json
       python scripts/freeze_analysis.py --verify receipts/discovery_rank_freeze.json
       python scripts/freeze_analysis.py --verify receipts/replication_implementation_freeze.json

   Run these from the restored project with the required Python available. Failure to find a
   referenced ignored file means the bundle is incomplete. Failure of a hash is a stop, not
   permission to update the freeze.
4. For a native portability check, source scripts/statistical_runtime_env.sh and run a NEW
   external output tree; preserve the SRV01 result:

       nextflow -C replication.config -log logs/aws_heca_replay.log run replication.nf -profile srv01 -params-file config/replication_run.json --outdir work/aws_heca_replay -with-trace logs/aws_heca_replay_trace.tsv

   The profile name srv01 selects the native locked environment; it does not establish host
   identity. The receipt must explicitly say EC2 when that is the actual host.
   This command is a proposed replay and has NOT been executed on AWS.
5. For container portability, use the separately verified image digest and an audited runtime
   overlay with the correct host UID/GID and project/task mounts. Verify read access to frozen
   inputs and write access only to task outputs. Nextflow stages inputs through filesystem
   paths, so host-path/symlink visibility matters. See [Nextflow containers](https://docs.seqera.io/nextflow/container).
   Do not claim the existing untested docker profile is ready merely because it parses.

## Acceptance and shutdown

A successful cloud replay requires actual exit0 for both Nextflow stages, real host/image/runtime
identities, input hashes matching the SRV01 checkpoint, unchanged selected gene IDs/directions,
matching tested-family sizes and agreement counts, and effect/statistic differences within a
declared numerical tolerance. Timestamps and host-specific receipt paths will differ; compare
the scientific tables and explain any numerical differences. A proposed tolerance of1e-8
absolute for finite effect/statistic values must be declared before the cloud comparison.

Preserve execution trace, console/error logs, output hashes, measured resources, any deviations,
and the exact command. A passing replay demonstrates engineering portability, not independent
biological validation: it reruns the same data. Public cohort replication is a separate result.

After the authorized attempt, copy the final evidence bundle to the approved location, then
stop/terminate the instance and review persistent volume, object, snapshot and registry storage.
Stopping compute alone does not demonstrate that every storage charge has ended. Verify the
account's actual resource inventory; this guide does not delete resources.

## Current limits

VERIFIED: native SRV01 workflow and its measured receipts.
DOCUMENTED: this EC2 replay and acceptance plan.
NOT VERIFIED: image build, actual container execution, AWS execution, cloud cost or performance.
