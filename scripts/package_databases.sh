#!/usr/bin/env bash
#
# package_databases.sh - build the archives that get deposited, and the manifest that
# describes them. Run once per database change, by whoever makes the deposit.
#
#   scripts/package_databases.sh --out-dir /somewhere/with/room
#
# WHAT IT PACKAGES, AND WHY ONLY THIS
# -----------------------------------
# Two artifacts, because only two cannot be re-created:
#   taxdump   NCBI overwrites new_taxdump.tar.gz in place, so the pinned snapshot is gone
#             upstream. Every taxid in every result keys to it.
#   markers   61 MB of curation, including 16 markers whose unbuildability is itself a
#             recorded finding. Rebuilding gets a near miss, not this database.
# The 43 GB reference pack is deliberately NOT packaged: packs/<pack>/refs.lock.tsv pins
# every accession and the md5 of its decompressed sequence, so a clone builds the same
# pack and can prove it (D52).
#
# tar.gz, not tar.zst, ON PURPOSE. zstd compresses the taxdump better and faster, but the
# install side then needs zstd installed before it can get a database - and the whole
# point of this route is that it works on a machine with nothing on it. `tar xzf` does not
# need anything.
#
# WHAT GOES INSIDE EVERY ARCHIVE
# ------------------------------
# DISCLAIMER.md and a PROVENANCE record. These are derived data redistributed from public
# NCBI records, snapshotted on purpose and therefore out of date by design; an archive
# that travels without saying so will be read as a curated reference resource, which it is
# not.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
die() { printf '[%s] ERROR: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; exit 1; }

BTF_DB="${BTF_DB:-}"
PACK="human"
OUT_DIR=""
MANIFEST="${REPO_DIR}/assets/databases.manifest.tsv"
TAXDUMP_DATE="2026-09-29"
DRY_RUN=0

usage() {
    cat <<EOF
Usage: scripts/package_databases.sh --out-dir DIR [options]

  --out-dir DIR     where the archives are written (needs ~1 GB free)
  --btf-db DIR      the built databases to package (default: \$BTF_DB)
  --pack NAME       domain pack (default: ${PACK})
  --manifest FILE   manifest to update in place (default: assets/databases.manifest.tsv)
  --dry-run         print what would be packaged
  -h, --help        this

Then: deposit the archives on Zenodo, and put each file's download URL into the manifest's
url column. Nothing else in the manifest needs editing by hand.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir)  OUT_DIR="${2:?--out-dir needs a value}"; shift 2 ;;
        --btf-db)   BTF_DB="${2:?--btf-db needs a value}"; shift 2 ;;
        --pack)     PACK="${2:?--pack needs a value}"; shift 2 ;;
        --manifest) MANIFEST="${2:?--manifest needs a value}"; shift 2 ;;
        --dry-run)  DRY_RUN=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *)          usage; die "unknown argument: $1" ;;
    esac
done

[[ -n "$OUT_DIR" ]] || { usage; die "--out-dir is required"; }
[[ -n "$BTF_DB" ]] || die "--btf-db is not set and \$BTF_DB is not in the environment"
[[ -f "$MANIFEST" ]] || die "no manifest at ${MANIFEST}"

TAXDUMP_DIR="${BTF_DB}/taxdump/${TAXDUMP_DATE}"
MARKERS_DIR="${BTF_DB}/packs/${PACK}/markers"
[[ -d "$TAXDUMP_DIR" ]] || die "no taxdump at ${TAXDUMP_DIR}"
[[ -d "$MARKERS_DIR" ]] || die "no marker database at ${MARKERS_DIR}"

# A marker database that is not current must not be published: a reader cannot tell a
# stale archive from a fresh one, and the manifest is authoritative for tested-vs-untested
# (D47). Refuse rather than ship it.
if ! (cd "$REPO_DIR" && BTF_DB="$BTF_DB" python3 scripts/build_marker_db.py \
        --pack "$PACK" --check >/dev/null 2>&1); then
    die "the marker database does not pass its own freshness check (D47). Publishing it
       would ship a database whose catalogue it no longer covers, and nothing downstream
       could tell. Run:
         python3 scripts/build_marker_db.py --pack ${PACK} --check"
fi
log "marker database passes --check"

mkdir -p "$OUT_DIR"
MARKERS_STAMP="$(date -u -r "${MARKERS_DIR}/markers.manifest.tsv" +%Y-%m-%d)"
TAX_ARCHIVE="btf-taxdump-${TAXDUMP_DATE}.tar.gz"
MARK_ARCHIVE="btf-markers-${PACK}-${MARKERS_STAMP}.tar.gz"

# PROVENANCE travels WITH the data. An archive that arrives without it is a pile of .dmp
# files whose snapshot date nobody can recover, which is exactly the situation that made
# publishing the taxdump necessary in the first place.
write_provenance() {
    local dest="$1" what="$2" extra="$3"
    cat > "${dest}/BTF_PROVENANCE.txt" <<EOF
bioThreatFinder published database artifact
===========================================
artifact:        ${what}
packaged (UTC):  $(date -u +%Y-%m-%dT%H:%M:%SZ)
pipeline commit: $(cd "$REPO_DIR" && git rev-parse --short HEAD 2>/dev/null || echo unknown)
pack:            ${PACK}

${extra}

This is DERIVED DATA: a dated snapshot of public NCBI records, redistributed so that a
specific analysis is reproducible. It is not an authoritative or curated reference
resource and not a substitute for the primary records, which remain the citable source.
It is out of date by design and progressively more so - the pinned taxonomy is frozen on
purpose so that results stay comparable.

READ DISCLAIMER.md IN THIS ARCHIVE BEFORE USING ANY RESULT DERIVED FROM IT. This software
is not a validated assay, no threshold in it is a measured operating point, and nothing it
produces should be used to decide that a sample is safe.
EOF
}

package() {
    local archive="$1" src_parent="$2" src_name="$3" what="$4" extra="$5"
    local out="${OUT_DIR%/}/${archive}"
    if (( DRY_RUN )); then
        log "[dry-run] would package ${src_parent}/${src_name} -> ${out}"
        return 0
    fi
    local staging
    staging="$(mktemp -d "${OUT_DIR%/}/.btf_pkg_XXXXXX")"
    # Copy the two travelling documents in beside the payload rather than into it, so
    # unpacking into $BTF_DB/<unpack_to> leaves them next to the data they describe.
    cp "${REPO_DIR}/DISCLAIMER.md" "$staging/"
    write_provenance "$staging" "$what" "$extra"
    log "packaging ${archive}"
    tar czf "$out" \
        -C "$src_parent" "$src_name" \
        -C "$staging" DISCLAIMER.md BTF_PROVENANCE.txt
    rm -rf "$staging"
    log "  $(du -h "$out" | cut -f1)  ${out}"
}

package "$TAX_ARCHIVE" "${BTF_DB}/taxdump" "$TAXDUMP_DATE" \
    "NCBI new_taxdump snapshot, pinned ${TAXDUMP_DATE}" \
    "Unpacks to taxdump/${TAXDUMP_DATE}/ under \$BTF_DB.
Why this is published at all: NCBI overwrites new_taxdump.tar.gz in place, so this exact
snapshot cannot be fetched from NCBI any more. The watchlist is resolved against it, the
Kraken2 index is built against it, and a different snapshot silently moves clade
membership - a copy 6.5 months older was measured at +288,764 taxids."

package "$MARK_ARCHIVE" "${BTF_DB}/packs/${PACK}" "markers" \
    "Marker-gene database for the ${PACK} pack, built ${MARKERS_STAMP}" \
    "Unpacks to packs/${PACK}/markers/ under \$BTF_DB.
markers.manifest.tsv is authoritative for which markers were actually built: a marker
absent from it evaluates UNKNOWN, which never promotes a finding. 16 of 166 markers can
only ever be UNKNOWN (3 failed, 4 curation_required, 9 motif_todo) and that is recorded
rather than hidden. The marker tier addresses gaps for which no published method exists,
so it is unvalidated - see DISCLAIMER.md."

(( DRY_RUN )) && { log "dry run complete"; exit 0; }

# --------------------------------------------------------------------------
# the manifest, rewritten from what was actually produced
# --------------------------------------------------------------------------
log "updating ${MANIFEST}"
python3 - "$MANIFEST" "$OUT_DIR" "$TAX_ARCHIVE" "$MARK_ARCHIVE" <<'PY'
import hashlib, pathlib, sys

manifest, out_dir, tax_archive, mark_archive = sys.argv[1:5]
by_artifact = {"taxdump": tax_archive, "markers": mark_archive}

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

path = pathlib.Path(manifest)
lines = path.read_text().splitlines()
out = []
for line in lines:
    if line.startswith("#") or not line.strip() or line.startswith("artifact\t"):
        out.append(line)
        continue
    f = line.split("\t")
    name = f[0]
    if name in by_artifact:
        archive = pathlib.Path(out_dir) / by_artifact[name]
        f[2] = by_artifact[name]
        f[3] = str(archive.stat().st_size)
        f[4] = sha256(archive)
        print(f"  {name}: {f[2]} {int(f[3]):,} bytes sha256 {f[4][:16]}...")
    out.append("\t".join(f))
path.write_text("\n".join(out) + "\n")
PY

cat <<EOF

Next, by hand, because it needs an account:

1. Create ONE record holding both files. Attach DISCLAIMER.md to the record as well as
   inside the archives.
2. Put each file's download URL into the url column of
   ${MANIFEST}
   and commit that - the installer refuses to download while it says PENDING.
3. Add the DOI to CITATION.cff and to the README's install section.

Then any clone installs with:
   scripts/install_databases.sh --route minimal

Record metadata - title, keywords, licence, related identifiers and the full
description - is ready to paste in:

    assets/zenodo_record.md

ONE record holding both files, not two: they install together and version together.
Cite the CONCEPT doi in CITATION.cff and the README; put the VERSION record's file
URLs in the manifest, because the installer verifies the sha256 it recorded and
"latest" would eventually stop matching.
EOF
