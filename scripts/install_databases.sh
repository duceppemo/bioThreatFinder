#!/usr/bin/env bash
#
# install_databases.sh - get a clone from nothing to a deployment that can run.
#
#   scripts/install_databases.sh --route minimal          # 62 GB, sylph stage 1
#   scripts/install_databases.sh --route standard         # + Kraken2 PlusPF, 177 GB
#   scripts/install_databases.sh --route full             # + Metabuli, 286 GB
#
# TWO WAYS TO GET A DATABASE, AND THE SPLIT IS NOT ABOUT SIZE
# -----------------------------------------------------------
# DOWNLOADED: artifacts that cannot be re-created. There are two, listed in
# assets/databases.manifest.tsv with a sha256 each. The taxonomy snapshot is the reason
# this route exists at all - NCBI overwrites new_taxdump.tar.gz in place, so the pinned
# 2026-09-29 snapshot is gone upstream, and every taxid in every result keys to it.
#
# BUILT: everything else, by the scripts in this directory, from inputs that pin what they
# fetch. The reference pack is built from packs/<pack>/refs.lock.tsv, which names every
# accession AND the md5 of its decompressed sequence, so a clone can prove it built the
# pack the measurements in docs/ were made on (D52). That is a better guarantee than a
# 43 GB archive, because the lock is auditable and an archive is not.
#
# WHAT "MINIMAL" MEANS, HONESTLY
# ------------------------------
# It is not small. A pipeline that can screen anything needs the reference pack (43 GB)
# and, for human samples, the host genome (17 GB) - stage 2 maps against the pack and
# sylph profiles it, so neither is optional. What minimal LEAVES OUT is the 115 GB Kraken2
# index and the 110 GB Metabuli index. Stage 1 then runs on sylph alone, which is a real
# configuration rather than a crippled one: D48 found 19 of 90 watchlist taxa have no node
# in the Kraken2 index at all, so sylph is already the only arm that can see them.
#
# EVERY STEP IS RESUMABLE AND EVERY STEP IS CHECKED
# -------------------------------------------------
# The build scripts skip completed work, so an interrupted install is restarted by running
# the same command again. At the end this prints - or with --verify runs - the read-only
# checks that answer whether each database covers its inputs: D45 for the bundles, D46 for
# the sketch, D47 for the markers, D48 for the index's coverage. An install that finishes
# is not the same as an install that is correct, and those checks are the difference.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
die() { printf '[%s] ERROR: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; exit 1; }

BTF_DB="${BTF_DB:-}"
PACK="human"
ROUTE="minimal"
MANIFEST="${REPO_DIR}/assets/databases.manifest.tsv"
FROM_DIR=""
THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)"
DRY_RUN=0
VERIFY=0
DOWNLOAD_ONLY=0
SKIP_DOWNLOAD=0

usage() {
    cat <<EOF
Usage: scripts/install_databases.sh [options]

Options:
  --btf-db DIR         where project databases live; must end in bioThreatFinder
                       (default: \$BTF_DB; required if unset)
  --pack NAME          domain pack (default: ${PACK})
  --route ROUTE        minimal | standard | full (default: ${ROUTE})
                         minimal   taxonomy + markers + reference pack + host genome
                                   + masking + sylph sketch           ~62 GB
                         standard  + Kraken2 PlusPF                 ~177 GB
                                   (needs ~700 GiB FREE while building)
                         full      + Metabuli                       ~286 GB
  --from-dir DIR       install the published archives from a local directory instead of
                       downloading them. Checksums are verified either way.
  --download-only      fetch and verify the published archives, build nothing
  --skip-download      build only; assume the published artifacts are in place
  --verify             run the read-only database checks at the end instead of printing
                       them
  --threads N          threads for the build steps (default: ${THREADS})
  --manifest FILE      published-artifact manifest (default: assets/databases.manifest.tsv)
  --dry-run            print what would happen, touch nothing
  -h, --help           this

The tool environments must exist first; each step says which one it needs:
  conda env create -f envs/btf-refs.yml
  conda env create -f envs/btf-kraken2.yml
  conda env create -f envs/btf-sylph.yml
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --btf-db)        BTF_DB="${2:?--btf-db needs a value}"; shift 2 ;;
        --pack)          PACK="${2:?--pack needs a value}"; shift 2 ;;
        --route)         ROUTE="${2:?--route needs a value}"; shift 2 ;;
        --from-dir)      FROM_DIR="${2:?--from-dir needs a value}"; shift 2 ;;
        --manifest)      MANIFEST="${2:?--manifest needs a value}"; shift 2 ;;
        --threads)       THREADS="${2:?--threads needs a value}"; shift 2 ;;
        --download-only) DOWNLOAD_ONLY=1; shift ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        --verify)        VERIFY=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        -h|--help)       usage; exit 0 ;;
        *)               usage; die "unknown argument: $1" ;;
    esac
done

case "$ROUTE" in
    minimal|standard|full) : ;;
    *) die "--route must be minimal, standard or full (got '${ROUTE}')" ;;
esac

[[ -n "$BTF_DB" ]] || die "--btf-db is not set and \$BTF_DB is not in the environment.
       \$BTF_DB is where the pipeline reads its databases; nextflow.config refuses to
       start without it, and no script here defaults to a path on any particular machine.
         export BTF_DB=/path/to/project/bioThreatFinder
       or: scripts/install_databases.sh --btf-db DIR"

# The --db-root build scripts append bioThreatFinder/ themselves, so $BTF_DB has to end in
# that segment or the pipeline and the builders write to two different trees - the shape
# that makes every watchlist taxon report no_reference, which reads as "checked and clear".
[[ "$(basename "$BTF_DB")" == "bioThreatFinder" ]] \
    || die "\$BTF_DB must end in 'bioThreatFinder' (got '${BTF_DB}'). The build scripts
       that take --db-root append that segment themselves, so a different name makes the
       builders and the pipeline disagree about where the databases are."
BTF_VOL="$(dirname "$BTF_DB")"

[[ -f "$MANIFEST" ]] || die "no manifest at ${MANIFEST}"
[[ -d "${REPO_DIR}/packs/${PACK}" ]] || die "no pack '${PACK}' under ${REPO_DIR}/packs"

LOCK="${REPO_DIR}/packs/${PACK}/refs.lock.tsv"
TAXDUMP="${BTF_DB}/taxdump/2026-09-29"

run() {
    if (( DRY_RUN )); then
        log "[dry-run] $*"
    else
        log "+ $*"
        "$@"
    fi
}

# --------------------------------------------------------------------------
# 1. the published artifacts
# --------------------------------------------------------------------------
# sha256 IS CHECKED BEFORE ANYTHING IS UNPACKED, and a mismatch is fatal. A truncated
# taxdump does not fail loudly later: resolve_watchlist.py would resolve fewer taxa and
# report success, which is this project's recurring failure shape - a degraded result
# indistinguishable from a clean one.
install_published() {
    local n_done=0 n_skip=0
    while IFS=$'\t' read -r artifact route filename bytes sha256 unpack_to url needed_for; do
        [[ "$artifact" == "artifact" || "$artifact" == \#* || -z "$artifact" ]] && continue

        local dest="${BTF_DB}/${unpack_to}"
        local stamp="${dest}/.btf_installed_${artifact}"
        if [[ -f "$stamp" ]]; then
            log "${artifact}: already installed (${stamp})"
            n_skip=$((n_skip + 1))
            continue
        fi

        local src=""
        if [[ -n "$FROM_DIR" ]]; then
            src="${FROM_DIR%/}/${filename}"
            [[ -f "$src" ]] || die "${artifact}: no ${filename} in ${FROM_DIR}"
        else
            if [[ "$url" == "PENDING" || -z "$url" ]]; then
                # A dry run is for seeing the plan, so it describes this rather than
                # dying on it; a real install cannot proceed and says what to do.
                if (( DRY_RUN )); then
                    log "[dry-run] ${artifact}: no URL in the manifest yet - a real"
                    log "[dry-run]   install needs the Zenodo deposit, or --from-dir DIR"
                    continue
                fi
                die \
"${artifact}: this manifest has no URL yet - the Zenodo deposit is not published.

       Three ways forward, in order of preference:

       1. If you have the archives (from a colleague, or a mirror):
            scripts/install_databases.sh --from-dir DIR ...
          The sha256 in assets/databases.manifest.tsv is checked either way, so a
          copy from anywhere is as trustworthy as a download.

       2. BUILD BOTH FROM SCRATCH, accepting one real difference. The marker
          database rebuilds from packs/<pack>/markers.tsv over NCBI Entrez - hours,
          and 16 of 166 markers are known-unbuildable for recorded reasons:
            python3 scripts/build_marker_db.py --pack ${PACK} --taxdump DIR
          The taxonomy snapshot CANNOT be rebuilt: NCBI overwrites
          new_taxdump.tar.gz in place, so the pinned 2026-09-29 snapshot is gone
          upstream. Take TODAY'S snapshot instead and re-resolve the watchlist
          against it - see step 1 of the README's quick start.
          WHAT THAT COSTS: your taxids are not the published ones. Clade
          membership moves between snapshots (a copy 6.5 months older measured
          +288,764 taxids), so every identity, breadth and limit-of-detection
          figure in docs/ describes a slightly different deployment than yours.
          Everything downstream must then use YOUR snapshot - the Kraken2 index
          especially, or the index and the watchlist disagree silently.

       3. Ask: the deposit is the project's own next step, and these two archives
          are 166 MB in total."
            fi
            src="${BTF_DB}/.btf_downloads/${filename}"
            run mkdir -p "${BTF_DB}/.btf_downloads"
            # -C resumes a partial transfer: these are hundreds of MB over a link that
            # may drop, and restarting from zero is how people give up and copy the
            # database from a colleague instead.
            run curl -fL --retry 5 --retry-delay 5 -C - -o "$src" "$url"
        fi

        if (( DRY_RUN )); then
            log "[dry-run] would verify sha256 of ${src} and unpack into ${dest}"
            continue
        fi

        if [[ "$sha256" == "PENDING" ]]; then
            log "WARNING ${artifact}: the manifest records no sha256, so the archive is"
            log "        unpacked UNVERIFIED. Regenerate the manifest with"
            log "        scripts/package_databases.sh before publishing."
        else
            log "${artifact}: verifying sha256"
            local got
            got="$(sha256sum "$src" | awk '{print $1}')"
            [[ "$got" == "$sha256" ]] || die \
"${artifact}: checksum mismatch on ${filename}
         manifest: ${sha256}
         file:     ${got}
       The download is corrupt or the archive is not the one this manifest describes.
       Delete it and retry; do not unpack it."
        fi

        log "${artifact}: unpacking into ${dest}"
        mkdir -p "$dest"
        tar xzf "$src" -C "$dest"
        date -u +%Y-%m-%dT%H:%M:%SZ > "$stamp"
        printf 'sha256\t%s\nfilename\t%s\n' "$sha256" "$filename" >> "$stamp"
        n_done=$((n_done + 1))
    done < "$MANIFEST"
    log "published artifacts: ${n_done} installed, ${n_skip} already present"
}

# --------------------------------------------------------------------------
# 2. the built artifacts, in dependency order
# --------------------------------------------------------------------------
build_refs() {
    [[ -f "$LOCK" ]] || die "no reference lock at ${LOCK}. Without it the pack resolves to
       whatever NCBI calls latest today, which is a different pack from the one every
       measurement in docs/ was made on (D52)."
    log "reference pack: building from the lock (needs env btf-refs; 43 GB, hours)"
    run python3 "${SCRIPT_DIR}/build_pack_refs.py" --pack "$PACK" --btf-db "$BTF_DB" \
        --from-lock "$LOCK" --taxdump "$TAXDUMP" --threads "$THREADS"
}

build_host() {
    log "host genome: fetching (needs env btf-refs; 17 GB)"
    run bash "${SCRIPT_DIR}/fetch_host_genome.sh" "$PACK" --db-root "$BTF_VOL" \
        --threads "$THREADS"
}

build_masking() {
    # Not cosmetic and not optional: D5 and D34. An unmasked vector or host stretch in a
    # reference recruits reads that are then attributed to the watchlist taxon, and a
    # reference with NO masking result is scored against a mask map that does not contain
    # it - indistinguishable from examined-and-clean.
    log "reference masking: (needs env btf-refs, and the host genome above)"
    run python3 "${SCRIPT_DIR}/mask_pack_refs.py" --pack "$PACK" --btf-db "$BTF_DB" \
        --threads "$THREADS"
}

build_sketch() {
    # 6 s for 375 genomes, and it must come AFTER the refs: the sketch records the genome
    # set it covers, and a sketch built early is a sketch that silently omits whatever
    # arrived later (D46).
    log "sylph sketch: (needs env btf-sylph; seconds)"
    run bash "${SCRIPT_DIR}/build_sylph_db.sh" --pack "$PACK" --btf-db "$BTF_DB" \
        --threads "$THREADS"
}

build_kraken2() {
    log "Kraken2 PlusPF: (needs env btf-kraken2; 115 GB finished, most of a day)"
    log "  AND ~700 GiB FREE WHILE IT BUILDS - several times the finished size."
    log "  build_kraken2_db.sh preflights that and refuses below it, so a disk sized"
    log "  for the route total declines here rather than earlier."
    log "  NOTE the prebuilt indexes upstream are NOT equivalent: the capped variants"
    log "  (k2_pluspf_16_GB and friends) are built with --max-db-size and lose"
    log "  sensitivity, and any upstream index carries its own taxonomy snapshot rather"
    log "  than the pinned one everything else here keys to."
    run bash "${SCRIPT_DIR}/build_kraken2_db.sh" --db-root "$BTF_VOL" \
        --taxdump "$TAXDUMP" --threads "$THREADS"
}

build_metabuli() {
    # Opt-in for a measured reason: across every scenario it changed zero findings at
    # ~22x Kraken2's CPU (D15). It is kept because its one untested case - an agent
    # absent from the index - cannot be simulated from a database.
    log "Metabuli: (110 GB) - off by default at run time; install only for a"
    log "  novel-agent question, where the amino-acid arm is the point"
    run bash "${SCRIPT_DIR}/build_metabuli_db.sh" --btf-db "$BTF_DB" --threads "$THREADS"
}

# --------------------------------------------------------------------------
# 3. did it work?
# --------------------------------------------------------------------------
CHECKS=(
  "python3 scripts/build_marker_db.py --pack ${PACK} --check"
  "bash scripts/build_sylph_db.sh --pack ${PACK} --check"
  "python3 scripts/mask_pack_refs.py --pack ${PACK} --dry-run"
  "python3 scripts/check_index_coverage.py --pack ${PACK}"
)

verify() {
    local rc=0
    for c in "${CHECKS[@]}"; do
        log "check: ${c}"
        if (cd "$REPO_DIR" && BTF_DB="$BTF_DB" eval "$c" >/dev/null 2>&1); then
            log "  ok"
        else
            local code=$?
            # Exit 3 is a verdict ("no, and here is what is wrong"), exit 1 is "could not
            # tell". They must not be conflated: the index check exits 3 BY DESIGN on this
            # pack, because 19 of 90 watchlist taxa have no node in PlusPF (D48).
            log "  exit ${code} - run it directly for the detail"
            rc=1
        fi
    done
    return $rc
}

log "==================== bioThreatFinder database install ===================="
log "pack=${PACK} route=${ROUTE} btf_db=${BTF_DB} threads=${THREADS}"
[[ -n "$FROM_DIR" ]] && log "published archives from: ${FROM_DIR}"
(( DRY_RUN )) && log "DRY RUN - nothing will be downloaded, built or written"

if (( ! SKIP_DOWNLOAD )); then
    install_published
fi
if (( DOWNLOAD_ONLY )); then
    log "--download-only: stopping before the build steps"
    exit 0
fi

build_refs
build_host
build_masking
build_sketch
case "$ROUTE" in
    standard) build_kraken2 ;;
    full)     build_kraken2; build_metabuli ;;
esac

log "-------------------------------------------------------------------------"
if (( DRY_RUN )); then
    log "dry run complete"
    exit 0
fi
if (( VERIFY )); then
    log "running the read-only database checks"
    verify || log "at least one check did not return 0 - see above, and docs/operations.md"
else
    log "An install that finished is not an install that is correct. Run these:"
    for c in "${CHECKS[@]}"; do log "    ${c}"; done
    log "  (0 = yes, 3 = no and it names what is wrong, 1 = could not tell.)"
    log "  check_index_coverage exits 3 on the human pack BY DESIGN: 19 of 90 watchlist"
    log "  taxa have no node in PlusPF and are reachable through sylph instead - D48."
fi
log "Read DISCLAIMER.md before using any output. No threshold here is a measured"
log "operating point, and this is not a validated assay."
log "========================================================================="
