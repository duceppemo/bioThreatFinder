"""Specimen composition per sample, for labelling a NEIGHBOUR present or absent.

The manifests' truth_present lists WATCHLIST taxa only; almost every neighbour that
matters here (L. innocua, S. sonnei, C. deneoformans, Citrobacter, Y. enterocolitica)
is a near_neighbour reference and appears in no manifest column. So composition is
written out from what each specimen actually IS.

THREE CLASSES, NOT TWO. A sibling reference that is absent still recruits the reads of
a present relative across most of its length, so lumping it with distant absentees
would hide the only case the metric has to get right.
"""
D6300 = {"Pseudomonas aeruginosa", "Escherichia coli", "Salmonella enterica",
         "Limosilactobacillus fermentum", "Lactobacillus fermentum",
         "Enterococcus faecalis", "Staphylococcus aureus", "Listeria monocytogenes",
         "Bacillus subtilis", "Saccharomyces cerevisiae", "Cryptococcus neoformans"}

# benchmarks/scenarios/human_core.yaml `background:` plus the spiked target and host.
SIM_BACKGROUND = {"Escherichia coli", "Neisseria lactamica", "Haemophilus haemolyticus",
                  "Listeria innocua", "Bacillus thuringiensis", "Homo sapiens"}

# Established by D18/defa5c7 as REAL trace contamination, confirmed by Kraken2.
ISOLATE_CONTAM = {"bacillus_subtilis": {"Listeria monocytogenes"},
                  "saccharomyces_cerevisiae": {"Listeria monocytogenes"}}

# Same-species-complex pairs: absent, but genomically inseparable from a present one.
SIBLING_OF = {
    "Shigella": "Escherichia coli", "Shigella dysenteriae": "Escherichia coli",
    "Shigella flexneri": "Escherichia coli", "Shigella sonnei": "Escherichia coli",
    "Shigella boydii": "Escherichia coli",
    "Cryptococcus deneoformans": "Cryptococcus neoformans",
    "Cryptococcus gattii": "Cryptococcus neoformans",
    "Listeria innocua": "Listeria monocytogenes",
    "Bacillus thuringiensis": "Bacillus cereus",
    "Yersinia pseudotuberculosis": "Yersinia pestis",
}

def present_set(run: str, sample: str) -> set:
    if run.startswith("zymo"):
        if "_isolate_" in sample:
            org = sample.split("_isolate_")[1].rsplit("_illumina", 1)[0]
            name = org.replace("_", " ").title()
            # Title() breaks binomials like "Escherichia Coli"; fix the epithet.
            parts = name.split()
            name = " ".join([parts[0]] + [p.lower() for p in parts[1:]])
            out = {name} | ISOLATE_CONTAM.get(org, set())
            if name == "Lactobacillus fermentum":
                out.add("Limosilactobacillus fermentum")
            return out
        return set(D6300)                      # zymo_even / zymo_log
    return set(SIM_BACKGROUND) | {"Yersinia pestis"}   # lod_v4 / d12fine

def label(run: str, sample: str, nb_name: str) -> str:
    if not nb_name:
        return "unknown"
    pres = present_set(run, sample)
    if nb_name in pres:
        return "PRESENT"
    sib = SIBLING_OF.get(nb_name)
    if sib and sib in pres:
        return "ABSENT_SIBLING"
    return "ABSENT_DISTANT"
