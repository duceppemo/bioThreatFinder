#!/usr/bin/env nextflow
/*
 * bioThreatFinder
 *
 * Two-stage metagenomic detection of high-consequence pathogens.
 *
 *   Stage 1  broad k-mer / ANI screen  -> candidate taxa (permissive)
 *   Stage 2  targeted competitive alignment -> evidence (strict)
 *
 * The split exists because stage-1 classifiers are built for community
 * profiling, not for single-organism calls. A handful of k-mers shared with an
 * environmental relative is enough to put a high-consequence taxon in a Kraken2
 * report. That is a lead, not a detection. Stage 2 is what decides.
 */

nextflow.enable.dsl = 2

include { FASTP; CHOPPER; SEQKIT_STATS      } from './modules/local/preprocess'
include { MINIMAP2_HOST_DEPLETE             } from './modules/local/preprocess'
include { RESOLVE_WATCHLIST                 } from './modules/local/detect'
include { KRAKEN2; BRACKEN; SYLPH_PROFILE   } from './modules/local/detect'
include { KRAKENUNIQ; METABULI              } from './modules/local/detect'
include { PARSE_SCREEN                      } from './modules/local/detect'
include { PAVIAN_EXPORT; PAVIAN_REPORT      } from './modules/local/detect'
include { MINIMAP2_COMPETITIVE              } from './modules/local/detect'
include { COVERAGE_STATS                    } from './modules/local/detect'
include { MARKER_ALIGN_NT; MARKER_SEARCH_AA } from './modules/local/detect'
include { MARKER_SEARCH_HMM; MARKER_EVIDENCE } from './modules/local/detect'
include { SCORE_EVIDENCE; MAKE_REPORT       } from './modules/local/detect'

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

def helpMessage() {
    log.info """
    bioThreatFinder v${workflow.manifest.version}

    Usage:
      nextflow run . -profile conda,workstation \\
          --input samplesheet.csv \\
          --pack human \\
          --outdir results

    Required:
      --input            Samplesheet CSV: sample,platform,fastq_1,fastq_2,role
                         platform = illumina | ont
                         role     = sample | control   (controls drive contamination subtraction)
      --pack             Domain pack: human | animal | plant   [${params.pack}]

    Key options:
      --outdir           Output directory                      [${params.outdir}]
      --db_root          Shared read-only DBs (taxdump, sylph, nt)  [${params.db_root}]
      --btf_db           Project-owned DBs on fast NVMe             [${params.btf_db}]

    Stage-1 classifier arms (any-classifier OR; consensus recorded, not required):
      --skip_kraken2     Disable the PRIMARY arm. Kraken2 with
                         --report-minimizer-data supplies distinct-minimizer counts
                         and a duplicity ratio, which is what separates reads spread
                         across a genome from reads piled on a conserved locus.
      --skip_krakenuniq false
                         OPT-IN second opinion (default off - note that OMITTING this
                         flag leaves the arm OFF; `--krakenuniq` is not a parameter and
                         Nextflow accepts unknown params silently). Kraken2's minimizer
                         columns are literally KrakenUniq's HyperLogLog code at the
                         same default precision, counting minimizers rather than
                         k-mers, and the only published comparison finds comparable
                         F1 once the threshold is shifted. KrakenUniq has had no
                         commit since Dec 2023. Prefer the prebuilt index over a
                         build if you enable this.
      --skip_metabuli    The DNA+amino-acid arm, DEFAULT TRUE (opt in with
                         `--skip_metabuli false`, not by omitting the flag).
                         Independent benchmarks rank its precision poorly and
                         several upstream bugs are open. Measured here it changed
                         no decision on any scenario run, including the divergent_*
                         pair meant to favour it (D15); its untested case is an
                         agent absent from the index. Treat its calls as leads.
      --skip_sylph       Disable ANI/containment.
      --min_distinct_kmers / --max_kmer_dup   stage-1 k-mer quality flags
      --high_unclassified 0.40
                         Unclassified-read fraction above which a sample is
                         flagged and bannered. UNCALIBRATED, and matrix- and
                         database-dependent. It is the only signal that
                         survives a divergent agent being confidently
                         mislabelled as a known taxon rather than missed.
      --skip_host_depletion
      --skip_markers     Disable marker gating. Taxa whose calls depend on
                         plasmid/toxin loci will then be promoted on alignment
                         evidence alone, which is NOT defensible for several.
      Stage-2 threshold OVERRIDES. packs/<pack>/pack.yaml is the source of truth;
      these only override it, and null means "use the pack value". An unknown key is
      rejected rather than ignored.
      --min_unique_reads    on-target reads      [${params.min_unique_reads ?: 'pack.yaml'}]
      --min_breadth         reference breadth    [${params.min_breadth ?: 'pack.yaml'}]
      --max_depth_cv        depth CV ceiling     [${params.max_depth_cv ?: 'pack.yaml'}]
      --min_read_identity   identity floor       [${params.min_read_identity ?: 'pack.yaml (per platform AND group)'}]

    Profiles: conda, mamba, docker, singularity, workstation, test
    """.stripIndent()
}

/* Pack YAML carries ${db_root} placeholders so a pack is portable between
 * machines without editing it. Resolve them here, once. */
def resolvePlaceholders(Object v, Map ctx) {
    if (v instanceof String) {
        def out = v
        ctx.each { k, val -> out = out.replace('${' + k + '}', val.toString()) }
        return out
    }
    if (v instanceof Map)  return v.collectEntries { k, val -> [k, resolvePlaceholders(val, ctx)] }
    if (v instanceof List) return v.collect { resolvePlaceholders(it, ctx) }
    return v
}

def loadPack(String packName) {
    def packDir = file("${params.packs_dir}/${packName}")
    if (!packDir.exists()) {
        exit 1, "Pack '${packName}' not found at ${packDir}. Available: " +
                file(params.packs_dir).list().sort().join(', ')
    }
    def yml = packDir.resolve('pack.yaml')
    if (!yml.exists()) exit 1, "Pack '${packName}' has no pack.yaml (${yml})"
    def parsed = new org.yaml.snakeyaml.Yaml().load(yml.text)
    def resolved = resolvePlaceholders(
        parsed,
        [db_root: params.db_root, btf_db: params.btf_db, projectDir: projectDir])
    resolved.dir = packDir
    return resolved
}

/* Pads an optional per-sample channel out to cover every sample.
 *
 * `base` carries (meta, placeholder) for all samples; `produced` carries
 * (meta, file) only for samples where the arm actually ran. A remainder join keeps
 * the missing ones and substitutes the placeholder, so a skipped or failed
 * classifier degrades that one arm rather than dropping the sample from the screen
 * entirely - which would look like a clean negative. */
def fillMissing(base, produced) {
    return base
        .join(produced, remainder: true)
        .map { meta, placeholder, actual -> tuple(meta, actual ?: placeholder) }
}

/* Guards against the worst failure mode this pipeline has.
 *
 * Classifier output is joined to packs/<pack>/watchlist.resolved.tsv, which is keyed
 * by NCBI taxid. Metabuli publishes GTDB-based indexes as well, and GTDB taxids are
 * NOT NCBI taxids. Point the pipeline at a GTDB index and every watchlist lookup
 * misses - producing not an error but an empty result, indistinguishable from a
 * clean negative across the entire screen.
 *
 * A silent false negative on every watchlist organism is worse than a crash, so this
 * is refused rather than warned about. The build scripts record the namespace in
 * BTF_BUILD_MANIFEST.txt precisely so it can be checked here. */
// A PACK CAN BE MASKED AND STILL HAVE REFERENCES THAT WERE NEVER EXAMINED.
//
// The check below this used to be all-or-nothing: the masked/ directory exists, or
// it does not. It cannot see the state the human pack was actually in for twelve
// days (D34) - 371 references, 359 mask BEDs, with every file of the Filoviridae
// group entry unmasked because the entry was added to the pack after the masking
// pass. Six of those had contigs in no BED anywhere, and coverage_stats reads every
// BED into one contig->intervals map, so a contig in no BED is simply absent from it
// - which is exactly what a contig with no masked intervals looks like.
//
// mask_pack_refs.py now reconciles this and exits non-zero, but that only helps
// someone who re-runs the masker. A run uses whatever is on disk, so the run says so
// too. WARN rather than exit: a partially masked pack is usable and the finding is
// still worth having, which is not true of the GTDB-namespace case below, where the
// arm contributes literally nothing.
def warnIfPackPartlyMasked(refsManifest, maskedDir, packName) {
    if (!refsManifest || !refsManifest.exists() || !maskedDir || !maskedDir.exists()) {
        return
    }
    def lines = refsManifest.readLines()
    if (lines.size() < 2) { return }
    def header = lines[0].split('\t') as List
    def iPath = header.indexOf('path')
    def iName = header.indexOf('name')
    if (iPath < 0) { return }

    // Basenames, never full paths: bed_path is recorded absolute, so a pack mounted
    // somewhere other than where it was masked would otherwise look wholly unmasked
    // - a false alarm in a guard teaches people to stop reading it.
    def haveBeds = (maskedDir.listFiles() ?: [])
        .findAll { it.name.endsWith('.mask.bed') }
        .collect { it.name } as Set

    // Counted over rows that HAVE a path, which is what mask_pack_refs.py masks and
    // therefore the only denominator the two can agree on. The manifest carries a
    // few rows with no path (a watchlist entry whose genome could not be fetched);
    // they are a different gap, reported as `no_reference` at scoring time.
    def considered = 0
    def missing = []
    lines.drop(1).each { line ->
        def f = line.split('\t') as List
        if (f.size() <= iPath) { return }
        def pth = f[iPath]?.trim()
        if (!pth) { return }
        // `considered++` is rejected by Nextflow's strict script syntax.
        considered += 1
        def stem = pth.split('/')[-1].replaceAll(/\.fna(\.gz)?$/, '')
        if (!haveBeds.contains("${stem}.mask.bed".toString())) {
            missing << (iName >= 0 && f.size() > iName ? f[iName] : stem)
        }
    }
    if (!missing) { return }

    // toUnique(), NOT unique(): the no-arg List.unique() deduplicates IN PLACE, so
    // building `shown` used to shrink `missing` before ${missing.size()} was evaluated -
    // and these names come from the manifest's `name` column, where twelve Filoviridae
    // references share one name. The twelve-unmasked-references state this guard exists to
    // report would have printed as "1 of 375".
    def uniq = missing.toUnique()
    def shown = uniq.take(6).join(', ')
    log.warn "${missing.size()} of ${considered} references in the ${packName} " +
             "pack have NO masking result (${shown}${uniq.size() > 6 ? ', ...' : ''}). " +
             "A reference that was never examined is NOT a clean one: coverage_stats " +
             "loads every BED into one contig map, so a contig that appears in no BED " +
             "is indistinguishable there from a contig with nothing to mask. Findings " +
             "on these references carry a full-length breadth denominator and no " +
             "reads_in_masked_regions. Fix with:\n" +
             "  python3 scripts/mask_pack_refs.py --pack ${packName} --dry-run   " +
             "(read-only; lists them and exits 1)"
}

// A SKETCH THAT EXISTS IS NOT A SKETCH OF THIS PACK.
//
// `file(sy_db_path, checkIfExists: true)` asks the only question Nextflow can ask about a
// database - is it there - and the sylph sketch spent 22 days answering it correctly while
// being wrong. It was written 2026-09-29 over 359 genomes; the pack reached 375, and the
// sixteen references added after it (all twelve Filoviridae, D42's second L. monocytogenes
// assembly, three Listeria neighbours) were in no sketch. The arm ran, the database loaded,
// and those genomes produced no rows - so screen_summary.tsv recorded sylph as `ok` for
// every run in that window. A reference nothing ever sketched reports exactly like a
// reference nothing was found on.
//
// D45 established the predicate for the competitive bundles - a bundle is rebuilt when its
// INPUTS changed, not merely when it is absent - and the lesson did not travel here.
//
// A RUN IS WHERE THE OPERATOR HAS TO BE TOLD, so this is deliberately not a schedule: the
// only consumer of the sketch is a run, and a nightly job would report staleness to a log
// on the days nobody was going to use the database anyway.
//
// Two signals, which is D45's pattern because it is the same problem:
//
//   the recorded genome set  EXACT, and the reason <db>.genomes.txt exists at all. The
//                            head node has no sylph - it lives in the btf-sylph env, used
//                            inside the process - so a run cannot read the genome list out
//                            of the database. The builder records it from `sylph inspect`
//                            at build time and the run compares against that, which sees a
//                            reference added, removed or swapped with every mtime settled.
//                            That is the state that actually shipped, and mtime is blind
//                            to it.
//   mtime                    catches content rewritten under a name the record knows.
//
// A sketch with no record still gets the mtime check, and the absence is reported only when
// something else is already wrong - a guard that nags on every run of a correct pack is a
// guard that gets ignored. `scripts/build_sylph_db.sh --stamp` writes the record for a
// database that predates it.
//
// WARN rather than exit, like the masking guard: the arm is still worth having for the
// genomes it does cover, and the other arms are unaffected.
def warnIfSketchStale(refsManifest, sylphDb, packName) {
    if (!refsManifest || !refsManifest.exists() || !sylphDb || !sylphDb.exists()) {
        return
    }
    def sketchMtime = sylphDb.lastModified()
    if (!sketchMtime) { return }

    def lines = refsManifest.readLines()
    if (lines.size() < 2) { return }
    def header = lines[0].split('\t') as List
    def iPath = header.indexOf('path')
    def iName = header.indexOf('name')
    if (iPath < 0) { return }

    // The genome set the builder read back out of the database, if it recorded one.
    // Basenames, because that is what sylph reports and what parse_screen.py reads a taxid
    // out of - and comparing absolute paths would call a pack that moved volumes stale in
    // its entirety.
    def sidecar = file("${sylphDb}.genomes.txt")
    def recorded = null
    if (sidecar.exists()) {
        def names = sidecar.readLines()
            .collect { it.trim() }
            .findAll { it && !it.startsWith('#') }
        if (names) { recorded = names as Set }
    }

    // Rows WITH a path, the same denominator the masking guard uses: a watchlist entry
    // whose genome could not be fetched is a different gap, reported as `no_reference`.
    def considered = 0
    def newer = []
    def notSketched = []
    def havePath = [] as Set
    lines.drop(1).each { line ->
        def f = line.split('\t') as List
        if (f.size() <= iPath) { return }
        def pth = f[iPath]?.trim()
        if (!pth) { return }
        considered += 1
        def base = pth.split('/')[-1]
        havePath << base
        def label = (iName >= 0 && f.size() > iName ? f[iName] : base)
        if (recorded != null && !recorded.contains(base)) {
            notSketched << label
        }
        def ref = file(pth)
        if (!ref.exists()) { return }
        if (ref.lastModified() > sketchMtime) {
            newer << label
        }
    }
    // Recorded but no longer in the pack. Not silent like the case above, but wrong in the
    // other direction: sylph reports filenames, so a hit on one is screen support for a
    // reference stage 2 has nothing to map against.
    def notInPack = recorded == null ? [] : (recorded - havePath) as List

    // The count the builder recorded, when it built this sketch. Only consulted when there
    // is no recorded genome set, which is strictly better: a count can only ever say HOW
    // MANY are missing, never which. Every pack's sketch shares one BTF_BUILD_MANIFEST.txt
    // in ${btf_db}/sylph/, so `db_file=` has to be checked: without it, a second pack's
    // build would make this guard compare the human sketch against the animal pack's
    // genome count and report a shortfall that is not there.
    def shortBy = 0
    def buildManifest = file("${sylphDb.parent}/BTF_BUILD_MANIFEST.txt")
    if (recorded == null && buildManifest.exists()) {
        def ml = buildManifest.readLines()
        def dbLine = ml.find { it.startsWith('db_file=') }
        def gLine = ml.find { it.startsWith('genomes=') }
        if (dbLine && gLine &&
            dbLine.split('=', 2)[1].trim() == sylphDb.toString() &&
            gLine.split('=', 2)[1].trim().isInteger()) {
            def sketched = gLine.split('=', 2)[1].trim() as Integer
            if (sketched < considered) { shortBy = considered - sketched }
        }
    }

    if (!newer && !shortBy && !notSketched && !notInPack) { return }

    def detail = []
    if (notSketched) {
        def uniq = notSketched.toUnique()
        detail << "${notSketched.size()} reference(s) are NOT in the sketch " +
                  "(${uniq.take(6).join(', ')}${uniq.size() > 6 ? ', ...' : ''})"
    }
    if (notInPack) {
        detail << "${notInPack.size()} sketched genome(s) are no longer in the pack " +
                  "(${notInPack.take(4).join(', ')}${notInPack.size() > 4 ? ', ...' : ''})"
    }
    if (shortBy) {
        detail << "its build manifest records ${considered - shortBy} genomes against " +
                  "${considered} references with a path in the pack, which says how many " +
                  "are missing but not which"
    }
    if (newer) {
        // toUnique(), NOT unique(): Groovy's no-arg List.unique() deduplicates IN PLACE and
        // returns the same list, so naming the genomes would silently rewrite the count
        // reported beside them. Twelve newer Filoviridae references share one `name`, and
        // this read "1 reference file" until the harness was pointed at the real manifest.
        def uniq = newer.toUnique()
        detail << "${newer.size()} reference file(s) are newer than the sketch " +
                  "(${uniq.take(6).join(', ')}${uniq.size() > 6 ? ', ...' : ''})"
    }
    // "is" when the recorded genome set proves it, "may be" when only mtime suggests it.
    // Overstating a timestamp as proof is how a guard loses its reader.
    def verdict = (notSketched || notInPack) ? 'is STALE' : 'may be STALE'
    def blind = recorded == null
        ? " This run could compare only TIMESTAMPS: there is no ${sidecar.name} beside the" +
          " sketch, so a reference added with its mtime settled would be invisible here." +
          " Write that record from the database with --stamp."
        : ""
    log.warn "The sylph sketch at ${sylphDb} ${verdict} for the ${packName} pack: " +
             "${detail.join('; ')}.${blind} A genome that was never sketched produces no " +
             "sylph row, which is indistinguishable from finding nothing on it, and " +
             "screen_summary.tsv will still record this arm as ok. Check and fix with:\n" +
             "  scripts/build_sylph_db.sh --pack ${packName} --check   " +
             "(read-only; exits 3 when stale and names the genomes)\n" +
             "  scripts/build_sylph_db.sh --pack ${packName} --force   (rebuild)\n" +
             "Or pass --skip_sylph to run without this arm rather than with a partial one."
}

// The `marker_id` column of a marker TSV, or null if the file has no usable header.
//
// A TOP-LEVEL FUNCTION, NOT A CLOSURE IN A LOCAL `def`. Nextflow 26.04's strict syntax
// resolves `name(...)` against declared functions only, so a closure called that way fails
// to COMPILE - which would have taken down every run with markers enabled, the default.
// Found by extracting the guard and running it, not by reading it.
//
// Both files are heavily commented and the comments are interleaved with data rows as
// section banners, so readLines()[0] is not the header. Same rule as
// read_marker_catalogue(): drop '#' lines, then the first line remaining is the header.
def markerIdColumn(tsv) {
    def rows = tsv.readLines().findAll { it.trim() && !it.trim().startsWith('#') }
    if (rows.size() < 2) { return null }
    def header = rows[0].split('\t')*.trim() as List
    def i = header.indexOf('marker_id')
    if (i < 0) { return null }
    def out = []
    rows.drop(1).each { line ->
        def f = line.split('\t') as List
        if (f.size() > i) {
            def v = f[i]?.trim()
            if (v) { out << v }
        }
    }
    return out
}

// A MARKER DATABASE THAT EXISTS IS NOT A DATABASE OF THIS CATALOGUE.
//
// Third time for D45's predicate - a database is current when it covers its INPUTS, not
// when it opens - and this arm's input is packs/<pack>/markers.tsv. Nothing compared them
// until now, so a marker added to the catalogue and never built was invisible at launch.
//
// THE COST IS A LOST CONFIRMATION, AND IT HIDES IN A KNOWN GAP. marker_evidence.py treats
// markers.manifest.tsv as authoritative for tested-vs-untested, so a catalogue marker with
// no manifest row evaluates UNKNOWN. UNKNOWN never promotes a finding, so nothing breaks
// and no call is obviously wrong - the finding simply settles a tier lower. And this pack
// already carries 16 markers that can only ever be UNKNOWN (3 failed, 4 curation_required,
// 9 motif_todo), so an operator reading UNKNOWN attributes it to the curation gap that is
// genuinely there. The one explanation requiring no action is the one that fits.
//
// Two comparisons here, both over the two plain TSVs:
//
//   the marker_id sets   EXACT. A marker added to the catalogue, or a manifest row for one
//                        no longer in it.
//   mtime                catches an EDITED marker, which the sets CANNOT see: change an
//                        accession or a gene under a marker_id the manifest already knows
//                        and every set still agrees while the sequence underneath is the
//                        old one. That is D41's wrong-gene bug, and it is the inverse of
//                        the sylph case, where the sets caught what mtime missed.
//
// NO SIDECAR RECORD, which is the difference from D46. That record exists because a run
// cannot open a .syldb - sylph lives in its own env and is used inside the process, so the
// head node could only compare timestamps. These artifacts are plain text, so there is
// nothing to keep stamped.
//
// The two comparisons --check also makes and this does not are deliberately left there:
// whether the manifest's claimed tiers are really in the combined FASTAs and HMM file, and
// whether an index predates its source. Those are build-integrity questions rather than
// did-you-rebuild questions, they mean streaming ~17 MB of artifact, and the suite already
// runs the full check against the shipped database. Naming --check is not a deferral of
// something this guard should have done.
//
// WARN rather than exit, like the masking and sketch guards: the arm is still worth having
// for the markers it does cover, and the other arms are unaffected.
def warnIfMarkerDbStale(catalogue, markersDir, packName) {
    if (!catalogue || !catalogue.exists() || !markersDir || !markersDir.exists()) {
        return
    }
    def manifest = file("${markersDir}/markers.manifest.tsv")
    if (!manifest.exists()) {
        log.warn "The marker database at ${markersDir} has no markers.manifest.tsv. " +
                 "marker_evidence.py treats that file as authoritative for tested vs " +
                 "untested, so EVERY marker would be reported UNKNOWN and no finding " +
                 "could be confirmed on marker evidence. Build it with:\n" +
                 "  python3 scripts/build_marker_db.py --pack ${packName}"
        return
    }

    def catIds = markerIdColumn(catalogue)
    def manIds = markerIdColumn(manifest)
    if (catIds == null || manIds == null || !catIds) { return }

    def manSet = manIds as Set
    def catSet = catIds as Set
    def notBuilt = catIds.findAll { !manSet.contains(it) }
    def notInCatalogue = manIds.findAll { !catSet.contains(it) }

    def catalogueNewer = false
    try {
        catalogueNewer = catalogue.lastModified() > manifest.lastModified()
    } catch (Exception ignored) {
        catalogueNewer = false
    }

    if (!notBuilt && !notInCatalogue && !catalogueNewer) { return }

    def detail = []
    if (notBuilt) {
        def uniq = notBuilt.toUnique()
        detail << "${notBuilt.size()} catalogue marker(s) have NO manifest row " +
                  "(${uniq.take(6).join(', ')}${uniq.size() > 6 ? ', ...' : ''})"
    }
    if (notInCatalogue) {
        def uniq = notInCatalogue.toUnique()
        detail << "${notInCatalogue.size()} manifest row(s) are no longer in the " +
                  "catalogue (${uniq.take(6).join(', ')}${uniq.size() > 6 ? ', ...' : ''})"
    }
    if (catalogueNewer) {
        detail << "markers.tsv is NEWER than the manifest, so a marker may have been " +
                  "edited under an id the database already knows - a changed accession " +
                  "or gene is invisible to the set comparison above"
    }
    // "is" when a set proves it, "may be" when only a timestamp suggests it. Overstating a
    // timestamp as proof is how a guard loses its reader.
    def verdict = (notBuilt || notInCatalogue) ? 'is STALE' : 'may be STALE'
    // THE CONSEQUENCE FOLLOWS THE SIGNAL THAT FIRED, because the two differ in kind and
    // one is worse. A missing manifest row costs a confirmation: UNKNOWN never promotes a
    // finding, so the call settles a tier lower. An EDITED marker is still searched, and
    // searched against the sequence the previous build retrieved - so it can return a
    // confident PRESENT or ABSENT for the wrong gene, which is D41's bug and is not a
    // downgrade at all. Printing the mild consequence for the severe case is how a guard
    // teaches its reader to misjudge it.
    def consequence = notBuilt
        ? "A catalogue marker with no manifest row evaluates UNKNOWN, which never " +
          "promotes a finding - the call simply settles a tier lower, and reads exactly " +
          "like the ${packName} pack's existing unevaluable markers."
        : (catalogueNewer && !notInCatalogue
            ? "An edited marker is still SEARCHED, against the sequence the previous " +
              "build retrieved, so it can report a confident PRESENT or ABSENT for the " +
              "wrong gene rather than declining to answer."
            : "A manifest row with no catalogue entry is evidence the database and the " +
              "catalogue were not built together.")
    // --markers-dir, because this guard KNOWS where the database is and the builder's
    // --db-root defaults to the SHARED root while a pack's marker database is built on the
    // project volume. Advice that omits it reports a healthy pack as having no database at
    // all - a false alarm produced by the guard's own instruction, which is the fastest way
    // to teach an operator to ignore both.
    log.warn "The marker database at ${markersDir} ${verdict} for the ${packName} " +
             "pack: ${detail.join('; ')}. ${consequence} Check and fix with:\n" +
             "  python3 scripts/build_marker_db.py --pack ${packName} " +
             "--markers-dir ${markersDir} --check   " +
             "(read-only; exits 3 when stale and names the markers)\n" +
             "  python3 scripts/build_marker_db.py --pack ${packName} " +
             "--db-root <the volume holding bioThreatFinder/>   (rebuild)\n" +
             "Or pass --skip_markers to run without this arm rather than with a partial one."
}

def assertNcbiTaxonomy(dbDir, label) {
    def manifest = file("${dbDir}/BTF_BUILD_MANIFEST.txt")
    if (!manifest.exists()) {
        log.warn "${label}: no BTF_BUILD_MANIFEST.txt under ${dbDir}, so the taxonomy " +
                 "namespace cannot be verified. If this index is GTDB-based, its taxids " +
                 "will not join the NCBI-keyed watchlist and this arm will silently " +
                 "contribute nothing. Build it with scripts/build_metabuli_db.sh to get " +
                 "a manifest, or pass --skip_metabuli to be explicit."
        return
    }
    def line = manifest.readLines().find { it.startsWith('taxonomy_namespace=') }
    if (!line) {
        log.warn "${label}: BTF_BUILD_MANIFEST.txt under ${dbDir} has no " +
                 "taxonomy_namespace entry; namespace could not be verified."
        return
    }
    def ns = line.split('=', 2)[1].trim()
    def nsUpper = ns.toUpperCase()

    // "NCBI" ALONE IS NOT ENOUGH OF AN ANSWER.
    //
    // build_metabuli_db.sh used to derive this line from --mode, so a GTDB prebuilt
    // unpacked under --mode refseq was stamped NCBI and sailed through here. The
    // builder now checks the written taxonomy instead, and says when it COULD NOT
    // verify ("NCBI (UNVERIFIED ...)") or when the artifact and the request disagree
    // ("NCBI (CONTRADICTION ...)"). Both still start with NCBI, so accepting any
    // NCBI-prefixed string would re-open the hole one layer up.
    if (nsUpper.contains('UNVERIFIED') || nsUpper.contains('CONTRADICTION')) {
        exit 1, """
        ERROR: the ${label} index at ${dbDir} records an UNCONFIRMED taxonomy namespace:

          ${ns}

        A GTDB-namespaced index joined to an NCBI-keyed watchlist yields ZERO matches -
        not an error, just an empty result indistinguishable from a clean negative for
        every organism on the watchlist. An unverified namespace is not evidence that
        the index is usable, so it is refused rather than assumed.

        Resolve by one of:
          - rebuild so the taxonomy can be checked: scripts/build_metabuli_db.sh --mode refseq
          - point --metabuli_db at an index whose manifest verifies as NCBI
          - drop the arm deliberately:              --skip_metabuli
        """.stripIndent()
    }
    if (!nsUpper.startsWith('NCBI')) {
        exit 1, """
        ERROR: the ${label} index at ${dbDir} is built in the '${ns}' taxonomy namespace.

        This pipeline joins classifier output to the watchlist by NCBI taxid. A
        non-NCBI index yields ZERO watchlist matches - not an error, just an empty
        result, which is indistinguishable from a clean negative for every organism
        on the watchlist. That is refused rather than warned about.

        Resolve by one of:
          - build an NCBI-namespaced index:  scripts/build_metabuli_db.sh --mode refseq
          - point --metabuli_db at an NCBI-namespaced index
          - drop the arm deliberately:       --skip_metabuli
        """.stripIndent()
    }
}

def parseSamplesheet(path) {
    return Channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def required = ['sample', 'platform', 'fastq_1']
            required.each { c ->
                if (!row.containsKey(c)) exit 1, "Samplesheet is missing column '${c}'"
            }
            def id = row.sample?.trim()
            if (!id) exit 1, "Samplesheet has a row with an empty 'sample' value"

            def platform = (row.platform ?: '').trim().toLowerCase()
            if (!(platform in ['illumina', 'ont'])) {
                exit 1, "Sample '${id}': platform must be 'illumina' or 'ont', got '${row.platform}'"
            }

            def role = ((row.role ?: 'sample').trim().toLowerCase()) ?: 'sample'
            if (!(role in ['sample', 'control'])) {
                exit 1, "Sample '${id}': role must be 'sample' or 'control', got '${row.role}'"
            }

            def r1 = file(row.fastq_1.trim(), checkIfExists: true)
            def r2 = (row.fastq_2?.trim()) ? file(row.fastq_2.trim(), checkIfExists: true) : null

            if (platform == 'ont' && r2) {
                exit 1, "Sample '${id}': platform 'ont' cannot have fastq_2"
            }

            def meta = [id: id, platform: platform, role: role, single_end: (r2 == null)]
            return tuple(meta, r2 ? [r1, r2] : [r1])
        }
        // A HEADER-ONLY SAMPLESHEET MUST NOT LOOK LIKE A SUCCESSFUL RUN.
        //
        // tests/data/samplesheet_test.csv was one line - the header - so `-profile
        // test` yielded no samples: every per-sample channel was empty, SCORE_EVIDENCE
        // ran on .ifEmpty([]) inputs, MAKE_REPORT emitted an empty report and the run
        // exited 0. The profile documented as existing "so the channel wiring can be
        // exercised in minutes" exercised none of it, and a green result was
        // indistinguishable from a working pipeline.
        .ifEmpty { exit 1, """
        ERROR: samplesheet ${params.input} contains no data rows.

        Nothing would be analysed and the run would still exit 0, which is
        indistinguishable from a clean negative for every organism on the watchlist.
        Add at least one sample row: sample,platform,fastq_1,fastq_2,role
        """.stripIndent() }
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    if (params.help) { helpMessage(); exit 0 }
    if (!params.input) { helpMessage(); exit 1, "ERROR: --input is required" }

    // `--skip_metabuli false` IS A SILENT NO-OP ON THE NEXTFLOW CLI, AND THAT HAS
    // ALREADY COST THIS PROJECT A FALSIFICATION TEST. See lib/BtfBool.groovy for the
    // mechanism and the history. params is read-only at runtime - assigning params[k]
    // is silently ignored - so the coerced values live here and every read below uses
    // SKIP, never params.skip_*.
    def SKIP = [:]
    ['skip_krakenuniq', 'skip_metabuli', 'skip_kraken2', 'skip_sylph',
     'skip_host_depletion', 'skip_markers', 'skip_marker_hmm',
     'skip_pavian_export', 'pavian_html'].each { k ->
        try { SKIP[k] = BtfBool.of(k, params[k]) }
        catch (IllegalArgumentException e) { exit 1, "ERROR: ${e.message}" }
    }

    def pack = loadPack(params.pack)

    def arms = []
    if (!SKIP.skip_krakenuniq) arms << 'krakenuniq'
    if (!SKIP.skip_metabuli)   arms << 'metabuli'
    if (!SKIP.skip_kraken2)    arms << 'kraken2/bracken'
    if (!SKIP.skip_sylph)      arms << 'sylph'
    if (arms.isEmpty()) {
        exit 1, "ERROR: every stage-1 classifier is disabled. Nothing would be screened."
    }
    if (SKIP.skip_kraken2) {
        log.warn "Kraken2 is disabled. It is the PRIMARY nucleotide arm and the source " +
                 "of distinct-minimizer evidence (--report-minimizer-data), so without " +
                 "it reads piled on a few conserved k-mers cannot be told apart from " +
                 "reads spread across a genome at stage 1."
    }

    log.info """
    bioThreatFinder v${workflow.manifest.version}
      pack            : ${pack.name}
      stage-1 arms    : ${arms.join(', ')}  (any-classifier OR; consensus recorded, not required)
      marker gating   : ${SKIP.skip_markers ? 'DISABLED' : 'enabled'}
      outdir          : ${params.outdir}
      report          : ${params.outdir}/report/biothreatfinder_report.html
      findings TSV    : ${params.outdir}/score_evidence/findings.tsv

    Findings are in-silico screening evidence only. Every one requires orthogonal
    laboratory confirmation and interpretation by a qualified laboratorian before
    any clinical or public-health action.
    """.stripIndent()

    ch_reads = parseSamplesheet(params.input)

    // -- Watchlist resolution --------------------------------------------
    // Runs once per pipeline invocation, not per sample. Fails the whole run if
    // any watchlist name is unresolvable, before any compute is spent.
    RESOLVE_WATCHLIST(
        file("${pack.dir}/${pack.reporting.watchlist}", checkIfExists: true),
        file(params.taxdump, checkIfExists: true)
    )
    ch_watchlist   = RESOLVE_WATCHLIST.out.resolved
    ch_descendants = RESOLVE_WATCHLIST.out.descendants

    // -- Stage 0: QC ------------------------------------------------------
    ch_reads.branch { meta, _reads ->
        illumina: meta.platform == 'illumina'
        ont:      meta.platform == 'ont'
    }.set { ch_by_platform }

    FASTP(ch_by_platform.illumina)
    CHOPPER(ch_by_platform.ont)

    ch_trimmed = FASTP.out.reads.mix(CHOPPER.out.reads)

    // -- Stage 0b: host depletion -----------------------------------------
    // Host reads are the majority of most clinical/veterinary specimens and a
    // major source of spurious classification via repeats.
    // A run with --skip_host_depletion must not inherit the previous run's QC file.
    // storeDir writes outside the work dir and is not cleaned between runs, so without a
    // sentinel, run B (skipped) would present run A's healthy depletion figures as its own.
    if (SKIP.skip_host_depletion) {
        ch_clean = ch_trimmed
        ch_host_stats = Channel.empty()
        Channel.of("sample\tplatform\talignment_records\thost_matching_records\t" +
                   "host_fraction_of_records\tsurviving_records\n" +
                   "ALL_SAMPLES\tNA\tNA\tNA\tSKIPPED\tNA\n")
            .collectFile(name: 'host_depletion.tsv', storeDir: "${params.outdir}/qc")
            .subscribe { f -> log.warn "[QC] host depletion SKIPPED; eukaryotic references " +
                               "recruit host reads in proportion to what survives it " +
                               "(docs/decisions.md D14) -> ${f}" }
    } else {
        ch_host_idx = ch_trimmed.map { meta, reads ->
            def idx = meta.platform == 'illumina'
                    ? pack.host.minimap2_index_sr
                    : pack.host.minimap2_index_lr
            tuple(meta, reads, file(idx, checkIfExists: true))
        }
        MINIMAP2_HOST_DEPLETE(ch_host_idx)
        ch_clean      = MINIMAP2_HOST_DEPLETE.out.reads
        ch_host_stats = MINIMAP2_HOST_DEPLETE.out.stats
        // Run-level QC artifact. Previously ch_host_stats was produced and consumed by
        // nothing; D14 showed the number is load-bearing, so it is now written out.
        MINIMAP2_HOST_DEPLETE.out.efficiency
            .collectFile(name: 'host_depletion.tsv',
                         storeDir: "${params.outdir}/qc",
                         keepHeader: true, skip: 1, sort: true)
            .subscribe { f -> log.info "[QC] host depletion -> ${f}" }
    }

    SEQKIT_STATS(ch_clean)

    // -- Stage 1: broad screen (multi-classifier) --------------------------
    // Independent arms, because they fail in different directions. Kraken2 is
    // PRIMARY (decisions.md D2); KrakenUniq is opt-in and OFF by default:
    //   kraken2     minimizer counts + duplicity -> specificity at stage 1, primary
    //   metabuli    DNA + amino acid             -> divergent / novel agents
    //   bracken     abundance re-estimation      -> abundance, continuity
    //   sylph       ANI / containment            -> orthogonal, low FP
    //   krakenuniq  HyperLogLog cardinality      -> OPT-IN second opinion
    // Any arm raises a candidate; consensus is recorded, not required. Requiring
    // consensus would discard the sensitivity that running four arms buys.
    // AN ABSENT OPTIONAL PATH IS AN EMPTY LIST, NOT A PLACEHOLDER FILE.
    //
    // This used to be `def no_file = file("assets/NO_FILE")`, and every padded arm got
    // that same file. Nextflow stages inputs by basename, so a process receiving two or
    // more of them died with "input file name collision" - which is what PARSE_SCREEN did
    // the moment two stage-1 arms were skipped. stageAs papered over it by giving each a
    // private directory; passing `[]` removes the cause, because nothing is staged at all
    // and there is no basename to collide. It is also the nf-core idiom, and the script
    // side becomes `x ? ... : ''` instead of comparing a filename to a magic string.

    // pack.yaml supplies the database paths; a CLI --*_db overrides it.
    def ku_db_path = params.krakenuniq_db ?: pack.screen.krakenuniq_db
    def mb_db_path = params.metabuli_db   ?: pack.screen.metabuli_db
    def k2_db_path = params.kraken2_db    ?: pack.screen.kraken2_db
    def sy_db_path = params.sylph_db      ?: pack.screen.sylph_db

    // Every arm is padded out to all samples by fillMissing(), so a skipped or
    // failed classifier never silently drops a sample from the join below.
    ch_arm_base = ch_clean.map { meta, _r -> tuple(meta, []) }

    if (SKIP.skip_krakenuniq) {
        ch_ku = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        KRAKENUNIQ(ch_clean, file(ku_db_path, checkIfExists: true))
        ch_ku = fillMissing(ch_arm_base, KRAKENUNIQ.out.report)
    }

    if (SKIP.skip_metabuli) {
        ch_mb = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        def mb_db = file(mb_db_path, checkIfExists: true)
        // Refuses a GTDB-namespaced index outright - see assertNcbiTaxonomy().
        assertNcbiTaxonomy(mb_db, 'Metabuli')
        METABULI(ch_clean, mb_db)
        ch_mb = fillMissing(ch_arm_base, METABULI.out.report)
    }

    if (SKIP.skip_kraken2) {
        ch_k2      = ch_clean.map { meta, _r -> tuple(meta, []) }
        ch_bracken = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        KRAKEN2(ch_clean, file(k2_db_path, checkIfExists: true))
        ch_k2 = fillMissing(ch_arm_base, KRAKEN2.out.report)
        // Bracken is built on a fixed read length, so it applies to short reads
        // only; ONT read lengths violate its model.
        BRACKEN(KRAKEN2.out.report.filter { meta, _r -> meta.platform == 'illumina' },
                file(k2_db_path))
        ch_bracken = fillMissing(ch_arm_base, BRACKEN.out.report)
    }

    if (SKIP.skip_sylph) {
        ch_sylph = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        def sy_db = file(sy_db_path, checkIfExists: true)
        // checkIfExists only asks whether the sketch is THERE - see warnIfSketchStale().
        warnIfSketchStale(
            file(pack.confirm.refs_manifest ?: "${pack.confirm.refs_dir}/manifest.tsv"),
            sy_db, pack.name)
        SYLPH_PROFILE(ch_clean, sy_db)
        ch_sylph = fillMissing(ch_arm_base, SYLPH_PROFILE.out.profile)
    }

    ch_screen = ch_ku.join(ch_k2).join(ch_bracken).join(ch_mb).join(ch_sylph)

    // WHICH ARMS WERE ENABLED, CARRIED ALONGSIDE THE REPORTS THEMSELVES.
    //
    // fillMissing() pads every absent arm to `[]`, which is the right fix for the
    // basename collision it replaced and which also erases the only distinction that
    // matters to a reader: a disabled arm, an arm that crashed under its
    // `errorStrategy 'ignore'`, and an arm that emitted an empty report all arrive at
    // PARSE_SCREEN as the same missing argument. Only this scope knows which arms were
    // asked for, so it says so and the screen summary can report an ignored failure as
    // a gap in coverage rather than as silence.
    //
    // bracken is handled for non-Illumina HERE rather than in the module, because this
    // is where that rule is made - BRACKEN is fed a channel filtered on
    // `meta.platform == 'illumina'` thirty lines up, since its model takes one fixed
    // read length and ONT violates it. Expecting it on an ONT sample would report the
    // pipeline's own deliberate choice as a failed arm.
    //
    // It is passed as INAPPLICABLE rather than simply dropped, and the distinction is
    // not cosmetic. Dropping it made it fall through to `not_run`, "disabled for this
    // run" - which the first real ONT run duly printed on a sample where kraken2 read
    // `ok` and nobody had disabled anything, inviting the reader to conclude an
    // operator had passed --skip_kraken2. "Nobody asked for it" and "it cannot apply
    // here" are different answers to why something was not searched.
    def armsEnabled = []
    if (!SKIP.skip_krakenuniq) armsEnabled << 'krakenuniq'
    if (!SKIP.skip_kraken2)    armsEnabled += ['kraken2', 'bracken']
    if (!SKIP.skip_metabuli)   armsEnabled << 'metabuli'
    if (!SKIP.skip_sylph)      armsEnabled << 'sylph'

    ch_screen_in = ch_screen.map { meta, ku, k2, br, mb, sy ->
        def illumina = meta.platform == 'illumina'
        def expected = armsEnabled.findAll { it != 'bracken' || illumina }
        def na = (!illumina && armsEnabled.contains('bracken')) ? ['bracken'] : []
        tuple(meta, ku, k2, br, mb, sy, expected.join(','), na.join(','))
    }

    PARSE_SCREEN(ch_screen_in, ch_watchlist.first(), ch_descendants.first())

    // Browsable copies of the stage-1 reports, for Pavian. Takes the same
    // joined channel PARSE_SCREEN reads, so what the user browses is exactly
    // what the screen was scored from, and reads nothing else. Purely
    // presentational: it has no output edge into scoring and cannot move a
    // call. Off with --skip_pavian_export.
    if (!SKIP.skip_pavian_export) {
        PAVIAN_EXPORT(ch_screen, file(params.taxdump, checkIfExists: true))

        // Pavian's static HTML, once per run. Opt-in because the pipeline
        // cannot install Pavian - see the process header.
        if (SKIP.pavian_html) {
            // ifEmpty, because collect() on an empty channel emits NOTHING and
            // the process would simply never run: --pavian_html would produce no
            // HTML and no error, which is exactly the silent skip the process
            // header promises not to do.
            PAVIAN_REPORT(
                PAVIAN_EXPORT.out.reports.collect().ifEmpty {
                    error "--pavian_html was requested but no classifier report " +
                          "was convertible for any sample, so there is nothing " +
                          "to render. Check the pavian/ directory and the " +
                          "PAVIAN_EXPORT logs."
                },
                file(params.outdir).name)
        }
    } else if (SKIP.pavian_html) {
        // Otherwise --pavian_html is accepted and quietly does nothing.
        exit 1, "ERROR: --pavian_html needs the Pavian exports, but " +
                "--skip_pavian_export is set. Drop one of the two."
    }

    // -- Stage 2: targeted confirmation -----------------------------------
    // Fan out one alignment job per (sample, candidate taxon). Candidates are
    // few by construction, so this stays cheap even on large runs.
    // NO SECOND READ-COUNT GATE HERE. THIS ONE SILENTLY DISCARDED EVERY SYLPH HIT.
    //
    // There used to be a
    //     .filter { (row.assigned_reads as Integer) >= params.min_screen_reads }
    // which looks like a harmless echo of the stage-1 threshold and is neither
    // harmless nor an echo. `assigned_reads` is a READ COUNT, and only the read-based
    // arms (kraken2, krakenuniq, metabuli) produce one; sylph reports ANI containment
    // and abundance, so a sylph-raised candidate carries assigned_reads=0 BY
    // CONSTRUCTION. With min_screen_reads=10, `0 >= 10` is false and the candidate was
    // dropped before stage 2 ever saw it.
    //
    // Measured on the first end-to-end run: Y. pestis at sylph abundance 6.46 and
    // E. coli at 22.17 were both raised correctly by PARSE_SCREEN, both filtered out
    // here, MINIMAP2_COMPETITIVE received nothing, and the run exited 0 with a
    // header-only findings.tsv - on a positive control whose scenario file says "any
    // outcome other than confirmed is a configuration fault". A clean negative for the
    // entire watchlist, reported as success.
    //
    // bin/parse_screen.py is the single place candidacy is decided: it applies
    // --min-reads to each read-based classifier, requires non-empty support, and
    // documents "Permissive OR: any classifier raises the candidate. Stage 2 decides."
    // Re-filtering on a read count here overrode that and contradicted it.
    ch_candidates = PARSE_SCREEN.out.candidates
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { meta, row -> tuple(meta, row) }

    ch_confirm_in = ch_candidates
        .combine(ch_clean, by: 0)
        .map { meta, row, reads ->
            // Competitive reference bundle: target + near neighbours, built by
            // scripts/build_pack_refs.py.
            //
            // Prefer the prebuilt, preset-matched minimap2 index. k and w are
            // baked in at index time, so an -x sr index performs badly on ONT
            // reads and vice versa - hence one index per preset rather than one
            // shared. Fall back to the FASTA (minimap2 indexes on the fly) so a
            // refs dir built with --skip-index still works.
            def preset = meta.platform == 'illumina' ? 'sr' : 'map-ont'
            def mmi = file("${pack.confirm.refs_dir}/competitive/${row.taxid}.${preset}.mmi")
            def ref = mmi.exists()
                    ? mmi
                    : file("${pack.confirm.refs_dir}/competitive/${row.taxid}.fna.gz")
            tuple(meta + [taxid: row.taxid, taxon: row.name, group: row.group,
                          assigned_reads: row.assigned_reads], reads, ref, row)
        }
        .branch { _meta, _reads, ref, _row ->
            present: ref.exists()
            missing: true
        }

    // A watchlist taxon with no reference genome cannot be confirmed. That is a
    // gap in the database, not a negative result, and it has to reach the report.
    ch_missing_ref = ch_confirm_in.missing.map { meta, _reads, ref, row ->
        tuple(meta, row + [confirm_status: 'no_reference', reference_path: ref.toString()])
    }

    MINIMAP2_COMPETITIVE(ch_confirm_in.present.map { meta, reads, ref, _row -> tuple(meta, reads, ref) })
    // Soft-mask intervals, when the pack has been masked. Absent is tolerated and
    // reported rather than assumed clean - an unmasked reference set is a known
    // false-positive source, not a neutral default.
    def masked_dir = pack.confirm.masked_dir ? file(pack.confirm.masked_dir) : []
    if (!masked_dir || !masked_dir.exists()) {
        log.warn "No masked reference intervals at ${pack.confirm.masked_dir}. " +
                 "References are UNMASKED, and vector, host, low-complexity or " +
                 "cross-contaminated sequence inside a reference will collect reads " +
                 "and inflate breadth. Run scripts/mask_pack_refs.py --pack ${pack.name}."
        masked_dir = []
    } else {
        warnIfPackPartlyMasked(
            file(pack.confirm.refs_manifest ?: "${pack.confirm.refs_dir}/manifest.tsv"),
            masked_dir, pack.name)
    }
    COVERAGE_STATS(MINIMAP2_COMPETITIVE.out.bam, masked_dir)

    // -- Stage 2b: marker tier --------------------------------------------
    // Run per sample, not per candidate: the marker DB is small, so one pass
    // serves every taxon's rule more cheaply than a fan-out.
    if (SKIP.skip_markers) {
        ch_marker_states = Channel.empty()
    } else {
        ch_markers_dir = file(pack.confirm.markers_dir ?: "${pack.confirm.refs_dir}/../markers")
        ch_marker_cat  = file("${pack.dir}/markers.tsv", checkIfExists: true)

        // checkIfExists only asks whether the catalogue is THERE, and the database beside
        // it can be older than the catalogue it claims to cover - see warnIfMarkerDbStale().
        warnIfMarkerDbStale(ch_marker_cat, ch_markers_dir, pack.name)

        MARKER_ALIGN_NT(ch_clean, ch_markers_dir)
        MARKER_SEARCH_AA(ch_clean, ch_markers_dir)

        if (SKIP.skip_marker_hmm) {
            ch_hmm = ch_clean.map { meta, _r -> tuple(meta, []) }
        } else {
            MARKER_SEARCH_HMM(MARKER_SEARCH_AA.out.hit_reads, ch_markers_dir)
            ch_hmm = MARKER_SEARCH_HMM.out.hits
        }

        ch_marker_in = MARKER_ALIGN_NT.out.bam
            .join(MARKER_SEARCH_AA.out.hits, remainder: true)
            .join(ch_hmm, remainder: true)
            .map { meta, bam, bai, aa, hmm ->
                tuple(meta, bam ?: [], bai ?: [], aa ?: [], hmm ?: [])
            }

        MARKER_EVIDENCE(ch_marker_in, ch_marker_cat, ch_markers_dir)
        ch_marker_states = MARKER_EVIDENCE.out.states
    }

    // -- Evidence scoring, control subtraction, tiering -------------------
    def manifest_tsv = file(pack.confirm.refs_manifest
                            ?: "${pack.confirm.refs_dir}/manifest.tsv")
    SCORE_EVIDENCE(
        COVERAGE_STATS.out.stats.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        PARSE_SCREEN.out.candidates.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        ch_missing_ref.map { meta, row -> "${meta.id}\t${meta.taxid}\t${meta.taxon}\tno_reference" }
            .collectFile(name: 'missing_refs.tsv', newLine: true).ifEmpty([]),
        ch_watchlist.first(),
        file("${pack.dir}/pack.yaml"),
        ch_marker_states.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        SKIP.skip_markers ? [] : file("${pack.dir}/marker_rules.yaml", checkIfExists: true),
        // The samplesheet is the only authoritative statement of which samples are
        // CONTROLS. Counting them from pipeline output misses the spotless control.
        file(params.input, checkIfExists: true),
        // NAMES for the competitive bundles' near neighbours. Cosmetic to the gates
        // and load-bearing for the reader: without it `attributed_to_neighbour` says
        // "reads favour neighbour taxid 620" and never that 620 is Shigella. Not
        // checkIfExists - a pack built before the manifest existed must still score.
        manifest_tsv.exists() ? manifest_tsv : [],
        // The sample's MEASURED read length, which turns the alignment-length gate
        // from a constant into a fraction. MAKE_REPORT has had these all along;
        // SCORE_EVIDENCE, the process that actually decides, did not - which is the
        // only reason that gate was ever expressed as an absolute number.
        SEQKIT_STATS.out.stats.map { meta, s -> s }.collect().ifEmpty([])
    )

    MAKE_REPORT(
        SCORE_EVIDENCE.out.findings,
        SEQKIT_STATS.out.stats.map { meta, s -> s }.collect().ifEmpty([]),
        file("${pack.dir}/pack.yaml"),
        file(params.input, checkIfExists: true),
        // WHAT STAGE 1 ACTUALLY SEARCHED. Until now this file was published and
        // nothing more - it reached neither SCORE_EVIDENCE nor MAKE_REPORT - so the
        // two things only it records never reached a reader: the unclassified
        // fraction, which is the one signal that survives a divergent agent being
        // confidently mislabelled as a known taxon, and an enabled arm that produced
        // no output, which is a hole in the screen's coverage rather than a negative.
        // It does not reach SCORE_EVIDENCE: it is coverage, and must not move a call.
        PARSE_SCREEN.out.summary.map { meta, tsv -> tsv }.collect().ifEmpty([])
    )
}

// NOTE: the top-level `workflow.onComplete` handler this file used to carry is not
// permitted by the Nextflow 24+ strict script syntax, which allows only
// declarations at script level. The run notice now prints from inside the entry
// workflow instead (the run-summary block in the entry workflow), and the caveats live where
// they actually matter - on every finding in the HTML report and the findings TSV,
// not in terminal scrollback that nobody keeps.
