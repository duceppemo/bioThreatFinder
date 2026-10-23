<!-- One concern per PR. Say what changed and why, not what the diff shows. -->

## What this changes, and why

## Evidence

<!-- If this rests on a measurement, give it. If it rests on judgement, say so — the
     distinction is tracked in docs/decisions.md and is more useful than a confident
     claim either way. -->

## Checklist

- [ ] `python3 tests/test_pipeline.py` passes
- [ ] A bug fix comes with the check that would have caught it
- [ ] A new check was **mutation-tested**: I broke the thing it catches and watched it fail
- [ ] The check cannot fire on a correct installation (a guard that cries wolf gets switched off)
- [ ] No absolute path anywhere, including comments and docs (`-k roots` covers this)
- [ ] Pins changed in `modules/local/*.nf` were followed by `scripts/sync_container_specs.py`
- [ ] `docs/decisions.md` has an entry if the reason is not obvious from the diff
