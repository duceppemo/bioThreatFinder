# bioThreatFinder

Metagenomic screening and in-silico confirmation of high-consequence pathogens, across
human, animal and plant health.

> **Research software, not a validated diagnostic assay.** Sequence detection does not
> establish a viable organism, an active infection, or exposure. No threshold here is a
> measured operating point, and **nothing this produces should be used to decide that a
> sample is safe.** Read
> [DISCLAIMER.md](https://github.com/duceppemo/bioThreatFinder/blob/main/DISCLAIMER.md).

## The problem this is shaped around

Stage-1 metagenomic classifiers are built for *community profiling*: "what is the
composition of this sample". They are good at it. They are not built for "is organism X
present", which carries a much harsher false-positive cost — a handful of k-mers shared
with a common environmental relative is enough to put a high-consequence name in a Kraken2
report. Run a classifier over a soil or sputum metagenome and you will see plenty of
alarming names. Essentially all of them are noise.

So the question is split in two, and the second half is where the tool earns its keep:

| Stage | Question | Posture |
|---|---|---|
| 1 | What *might* be here? | Deliberately permissive — any single arm raises a candidate, because missing a lead costs more than carrying a false one forward |
| 2 | Does the alignment evidence hold up? | Strict — this is what decides. Nothing above `screened_only` is reported on stage-1 evidence alone |

Stage 2 maps each candidate's reads against an index holding the **target together with its
near neighbours**, so a read that fits a harmless relative better is scored to the relative
rather than to the agent. A marker-gene tier gates the top tier wherever the chromosome is
not the actionable signal.

[Methods](Methods) covers all of it, including the three-valued marker logic and why
"absent" and "untested" must not be collapsed.

## Where to go

- **Installing it** → [Installation](Installation), then [Databases](Databases)
- **Running it** → [Usage](Usage), then [Outputs](Outputs)
- **Deciding whether to trust a result** → [Validation](Validation). Read this before
  drawing a conclusion from a run
- **Something looks wrong** → [FAQ](FAQ)
- **Working on it** → [Development](Development)

## The thing worth knowing before anything else

Almost every bug found in this project so far shares one shape: *the pipeline produced an
empty or degraded result that was indistinguishable from a clean negative.* An arm that ran
and reported nothing. A reference scored against a mask map that did not contain it. A
sketch that silently omitted sixteen genomes, including every one that a divergent
filovirus would have matched. A marker whose database row was missing, so the rule
evaluated UNKNOWN and the finding settled one tier lower with no message.

That is why there are so many checks, why each distinguishes *no* from *could not tell*,
and why a negative from this pipeline is weaker than it looks. [Validation](Validation) is
specific about which of those gaps are still open.

## Provenance

- Software: [github.com/duceppemo/bioThreatFinder](https://github.com/duceppemo/bioThreatFinder)
- Databases: [doi:10.5281/zenodo.23602998](https://doi.org/10.5281/zenodo.23602998)
- The engineering record, with the evidence behind each choice and what would reverse it:
  [`docs/decisions.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/decisions.md)
