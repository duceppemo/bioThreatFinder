#!/usr/bin/env python3
"""
Marker-rule engine: three-valued boolean logic over marker evidence.

WHY THREE-VALUED LOGIC

A marker can be in one of three states, and collapsing them to two is a
correctness bug with real consequences:

  PRESENT   tested, found, met coverage/identity thresholds
  ABSENT    tested, and confidently not found
  UNKNOWN   never tested - not in the marker DB, or the sequence could not be
            retrieved at build time

"Absent" is evidence. "Untested" is not. If they are treated the same, then a
gap in the marker database silently reads as a negative result - which is how a
tool ends up reporting "not BCG" for a sample it never checked, or failing to
confirm because a marker it never had was scored as missing.

So evaluation uses Kleene three-valued logic:

  not:  T->F, F->T, U->U
  and:  min(a,b) with F < U < T
  or:   max(a,b) with F < U < T

The useful consequence: `UNKNOWN and FALSE` is FALSE (one confident failure is
enough), but `UNKNOWN and TRUE` is UNKNOWN (cannot claim success on untested
evidence). A rule evaluating to UNKNOWN never promotes a finding to
confirmed_in_silico - it holds it at `supported` with an explicit
"marker evidence incomplete" reason naming the untested markers.

GRAMMAR (recursive descent; no eval(), so a rules file cannot execute code)

  expr    := or_expr
  or_expr := and_expr ('or' and_expr)*
  and_expr:= unary ('and' unary)*
  unary   := 'not' unary | atom
  atom    := '(' expr ')' | ('present'|'absent') '(' ident ')' | ident
            | ('any_of'|'all_of'|'n_of') '(' [int ','] ident (',' ident)* ')'

A bare identifier means present(identifier).

Run `marker_rules.py --self-test` to exercise the logic.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# Ordered so that min/max implement Kleene and/or directly.
FALSE, UNKNOWN, TRUE = 0, 1, 2
STATE_NAME = {FALSE: "FALSE", UNKNOWN: "UNKNOWN", TRUE: "TRUE"}

TOKEN_RE = re.compile(r"\s*(\(|\)|,|\band\b|\bor\b|\bnot\b|[A-Za-z_][A-Za-z0-9_.\-]*|\d+)")


class RuleError(ValueError):
    pass


def tokenize(text: str) -> list[str]:
    tokens, pos = [], 0
    while pos < len(text):
        m = TOKEN_RE.match(text, pos)
        if not m:
            if text[pos:].strip() == "":
                break
            raise RuleError(f"unparseable text at offset {pos}: {text[pos:pos + 20]!r}")
        tokens.append(m.group(1))
        pos = m.end()
    return tokens


class Parser:
    """Recursive descent. Produces a closure tree, not code."""

    FUNCS = {"present", "absent", "any_of", "all_of", "n_of"}

    def __init__(self, tokens: list[str]):
        self.t = tokens
        self.i = 0
        self.referenced: set[str] = set()

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None

    def take(self, expected=None):
        tok = self.peek()
        if tok is None:
            raise RuleError("unexpected end of expression")
        if expected and tok.lower() != expected.lower():
            raise RuleError(f"expected {expected!r}, found {tok!r}")
        self.i += 1
        return tok

    def parse(self):
        node = self.or_expr()
        if self.peek() is not None:
            raise RuleError(f"trailing tokens from {self.peek()!r}")
        return node

    def or_expr(self):
        node = self.and_expr()
        while self.peek() and self.peek().lower() == "or":
            self.take()
            rhs = self.and_expr()
            node = (lambda a, b: lambda st: max(a(st), b(st)))(node, rhs)
        return node

    def and_expr(self):
        node = self.unary()
        while self.peek() and self.peek().lower() == "and":
            self.take()
            rhs = self.unary()
            node = (lambda a, b: lambda st: min(a(st), b(st)))(node, rhs)
        return node

    def unary(self):
        if self.peek() and self.peek().lower() == "not":
            self.take()
            inner = self.unary()
            return lambda st: TRUE - inner(st) if inner(st) != UNKNOWN else UNKNOWN
        return self.atom()

    def _idents(self) -> tuple[list[str], int | None]:
        """Parse `( [n,] id, id, ... )`."""
        self.take("(")
        n = None
        if self.peek() and self.peek().isdigit():
            n = int(self.take())
            self.take(",")
        names = [self.take()]
        while self.peek() == ",":
            self.take(",")
            names.append(self.take())
        self.take(")")
        for name in names:
            self.referenced.add(name)
        return names, n

    def atom(self):
        tok = self.peek()
        if tok is None:
            raise RuleError("unexpected end of expression")

        if tok == "(":
            self.take("(")
            node = self.or_expr()
            self.take(")")
            return node

        low = tok.lower()
        if low in self.FUNCS and self.i + 1 < len(self.t) and self.t[self.i + 1] == "(":
            self.take()
            names, n = self._idents()

            if low == "present":
                if len(names) != 1:
                    raise RuleError("present() takes exactly one marker")
                name = names[0]
                return lambda st: st.get(name, UNKNOWN)

            if low == "absent":
                if len(names) != 1:
                    raise RuleError("absent() takes exactly one marker")
                name = names[0]
                # Untested stays untested: absent(x) is only TRUE when x was
                # tested and confidently not found.
                return lambda st: {TRUE: FALSE, FALSE: TRUE}.get(
                    st.get(name, UNKNOWN), UNKNOWN)

            if low == "any_of":
                return lambda st: max((st.get(x, UNKNOWN) for x in names), default=FALSE)

            if low == "all_of":
                return lambda st: min((st.get(x, UNKNOWN) for x in names), default=TRUE)

            if low == "n_of":
                if n is None:
                    raise RuleError("n_of() requires a leading count, e.g. n_of(2, a, b, c)")
                need = n

                def n_of(st, names=names, need=need):
                    vals = [st.get(x, UNKNOWN) for x in names]
                    n_true = sum(1 for v in vals if v == TRUE)
                    n_unk = sum(1 for v in vals if v == UNKNOWN)
                    if n_true >= need:
                        return TRUE
                    # Could still reach the threshold if untested markers resolve.
                    return UNKNOWN if n_true + n_unk >= need else FALSE

                return n_of

        if low in ("and", "or", "not", ")", ","):
            raise RuleError(f"unexpected token {tok!r}")

        self.take()
        self.referenced.add(tok)
        return lambda st, name=tok: st.get(name, UNKNOWN)


def compile_rule(expr: str):
    """Return (callable(state_dict) -> tri-state, referenced marker ids)."""
    p = Parser(tokenize(expr))
    fn = p.parse()
    return fn, p.referenced


# ---------------------------------------------------------------------------
# Rules file (minimal YAML subset - stdlib only, same reason as elsewhere)
# ---------------------------------------------------------------------------

RULE_KEYS = {"confirm_requires", "support_requires", "veto_if", "vaccine_if",
             "report_always", "rationale", "no_markers_available"}

# Keys whose value is a flag rather than an expression. They are coerced to a real
# bool here, and ANYTHING ELSE IS AN ERROR, because the string "false" is truthy in
# Python and this particular flag lifts the cap that stops a taxon with no marker
# evidence reaching confirmed_in_silico. `no_markers_available: flase` must not read
# as permission; silence must not grant it and neither must a typo.
BOOL_KEYS = {"no_markers_available"}
_TRUE = {"true", "yes", "on", "1"}
_FALSE = {"false", "no", "off", "0"}


def load_rules(path: Path) -> dict[str, dict]:
    """Parse the restricted YAML shape used by marker_rules.yaml.

    Shape:  "<taxon name>":\n    key: value\n    key: >\n      folded text
    """
    rules: dict[str, dict] = {}
    current: str | None = None
    pending_key: str | None = None
    fold: list[str] = []

    def flush():
        nonlocal pending_key, fold
        if current and pending_key and fold:
            rules[current][pending_key] = " ".join(fold).strip()
        pending_key, fold = None, []

    for raw in path.read_text().splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip())
        line = raw.strip()

        if indent == 0:
            flush()
            if not line.endswith(":"):
                raise RuleError(f"expected a taxon key ending in ':', got {line!r}")
            current = line[:-1].strip().strip("'\"")
            # DUPLICATE TAXON BLOCKS ARE AN ERROR, NOT LAST-ONE-WINS. This file is 452
            # lines of hand-curated rules; appending a second block for a taxon that
            # already has one silently DELETED the first - including its veto_if. The
            # B. thailandensis -> B. pseudomallei veto the README highlights could vanish
            # from an ordinary editing mistake with no error, and --validate cannot see it
            # because it validates only what survived the parse. bin/btf_yaml.py refuses
            # duplicate keys for exactly this reason and says so in its docstring; this
            # loader is separate and did not.
            if current in rules:
                raise RuleError(
                    f"duplicate taxon block {current!r}: a second block silently replaces "
                    f"the first, discarding its rules. Merge them into one block.")
            rules[current] = {}
            continue

        if current is None:
            raise RuleError(f"indented line before any taxon key: {line!r}")

        if pending_key and indent >= 4 and ":" not in line.split(" ")[0]:
            fold.append(line)
            continue

        flush()
        if ":" not in line:
            raise RuleError(f"expected 'key: value' in {current!r}, got {line!r}")
        key, val = line.split(":", 1)
        key, val = key.strip(), val.strip()
        if key not in RULE_KEYS:
            raise RuleError(
                f"unknown rule key {key!r} for {current!r}; allowed: {sorted(RULE_KEYS)}")
        # Same reasoning as the duplicate taxon block: last-one-wins on a repeated key
        # drops a curated rule without a word.
        if key in rules[current] or (pending_key == key):
            raise RuleError(
                f"duplicate key {key!r} in block {current!r}: the later value would "
                f"silently replace the earlier rule.")
        if key in BOOL_KEYS:
            low = val.strip("'\"").strip().lower()
            if low in _TRUE:
                rules[current][key] = True
            elif low in _FALSE:
                rules[current][key] = False
            else:
                raise RuleError(
                    f"{key!r} in block {current!r} must be a boolean "
                    f"(true/false), got {val!r}. This flag lifts the cap that keeps "
                    f"a taxon with no marker evidence out of confirmed_in_silico, so "
                    f"an unrecognised value is refused rather than read as true.")
        elif val in (">", "|", ">-", "|-", ""):
            pending_key = key
            fold = []
        else:
            rules[current][key] = val.strip("'\"")
    flush()
    return rules


def validate_rules(rules: dict[str, dict], known_markers: set[str]) -> list[str]:
    """Compile every expression and check referenced markers exist.

    A rule referring to a marker that is not in the catalogue would evaluate to
    UNKNOWN forever, permanently holding findings at `supported` with no obvious
    cause. Better to fail at load.
    """
    problems = []
    for taxon, spec in rules.items():
        for key in ("confirm_requires", "support_requires", "veto_if", "vaccine_if"):
            expr = spec.get(key)
            if not expr:
                continue
            try:
                _, refs = compile_rule(expr)
            except RuleError as exc:
                problems.append(f"{taxon}.{key}: {exc}")
                continue
            unknown = sorted(refs - known_markers)
            if unknown:
                problems.append(
                    f"{taxon}.{key}: references unknown marker(s) {', '.join(unknown)}")
    return problems


def evaluate(taxon: str, spec: dict, state: dict[str, int]) -> dict:
    """Apply one taxon's rules to its marker state. Returns a verdict dict."""
    out = {
        "taxon": taxon,
        "confirm": None, "support": None, "veto": None, "vaccine": None,
        "untested": [], "reasons": [],
    }

    for key, field in (("confirm_requires", "confirm"), ("support_requires", "support"),
                       ("veto_if", "veto"), ("vaccine_if", "vaccine")):
        expr = spec.get(key)
        if not expr:
            continue
        fn, refs = compile_rule(expr)
        val = fn(state)
        out[field] = val
        untested = sorted(r for r in refs if state.get(r, UNKNOWN) == UNKNOWN)
        if val == UNKNOWN and untested:
            out["untested"].extend(untested)
            out["reasons"].append(
                f"{key} is indeterminate - untested marker(s): {', '.join(untested)}")
        else:
            out["reasons"].append(f"{key} = {STATE_NAME[val]}")

    # report_always: markers surfaced regardless of whether any rule fires.
    #
    # This key was accepted by the parser, documented in marker_rules.yaml's header as
    # "markers always surfaced (resistance, subtyping)", and populated for real taxa -
    # V. cholerae's vchol_wbeO1/vchol_wbfZ/vchol_tcpA serogroup and pilus typing, and
    # A. fumigatus's cyp51A resistance alleles, which the README calls "the finding" for
    # that organism. evaluate() never read it and nothing downstream surfaced it, so
    # those markers were curated, built, scored into marker_states.tsv and then dropped
    # before the report: a documented capability that silently produced nothing.
    always = spec.get("report_always") or []
    if isinstance(always, str):
        always = [x.strip() for x in always.replace(",", " ").split() if x.strip()]
    # PRESENT/ABSENT rather than TRUE/FALSE: these are marker states being surfaced to
    # an analyst, not the truth values of a rule.
    marker_state_name = {TRUE: "PRESENT", FALSE: "ABSENT", UNKNOWN: "UNKNOWN"}
    reported = []
    for mid in always:
        reported.append({"marker_id": mid,
                         "state": marker_state_name[state.get(mid, UNKNOWN)]})
    out["report_always"] = reported
    if reported:
        out["reasons"].append(
            "always-reported: "
            + ", ".join(f"{r['marker_id']}={r['state']}" for r in reported))

    out["untested"] = sorted(set(out["untested"]))
    return out


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def self_test() -> int:
    failures = []

    def check(expr, state, expected, label=""):
        fn, _ = compile_rule(expr)
        got = fn(state)
        if got != expected:
            failures.append(
                f"{label or expr!r} with {state} -> {STATE_NAME[got]}, "
                f"expected {STATE_NAME[expected]}")

    P, A = TRUE, FALSE

    # Kleene basics
    check("a and b", {"a": P, "b": P}, TRUE)
    check("a and b", {"a": P, "b": A}, FALSE)
    check("a and b", {"a": P}, UNKNOWN, "untested cannot confirm")
    check("a and b", {"a": A}, FALSE, "one confident failure decides")
    check("a or b", {"a": P}, TRUE, "one confident success decides")
    check("a or b", {"a": A}, UNKNOWN, "untested alternative leaves it open")
    check("a or b", {"a": A, "b": A}, FALSE)
    check("not a", {"a": P}, FALSE)
    check("not a", {"a": A}, TRUE)
    check("not a", {}, UNKNOWN, "negating untested stays untested")

    # present / absent distinction - the core of the design
    check("absent(a)", {"a": A}, TRUE, "tested and not found")
    check("absent(a)", {"a": P}, FALSE)
    check("absent(a)", {}, UNKNOWN, "NEVER TESTED must not read as absent")
    check("present(a)", {"a": P}, TRUE)
    check("present(a)", {}, UNKNOWN)

    # Real anthrax rule: needs a pXO1 toxin gene AND a pXO2 capsule gene
    anthrax = "(pagA or lef or cya) and (capB or capA)"
    check(anthrax, {"pagA": P, "capB": P}, TRUE, "anthrax both plasmids")
    check(anthrax, {"pagA": P, "capB": A, "capA": A}, FALSE, "pXO1 only -> not confirmed")
    check(anthrax, {"pagA": A, "lef": A, "cya": A, "capB": P}, FALSE, "pXO2 only")
    check(anthrax, {"pagA": P}, UNKNOWN, "pXO2 untested -> must not confirm")

    # BCG: MTBC signal with RD1 deleted is the vaccine strain, not disease
    check("absent(mtb_RD1)", {"mtb_RD1": A}, TRUE, "BCG detected")
    check("absent(mtb_RD1)", {"mtb_RD1": P}, FALSE, "wild-type RD1 intact")
    check("absent(mtb_RD1)", {}, UNKNOWN, "RD1 untested - cannot claim BCG or wild type")

    # Aggregates
    check("n_of(2, a, b, c)", {"a": P, "b": P, "c": A}, TRUE)
    check("n_of(2, a, b, c)", {"a": P, "b": A, "c": A}, FALSE)
    check("n_of(2, a, b, c)", {"a": P}, UNKNOWN, "could still reach 2")
    check("any_of(a, b)", {"a": A, "b": P}, TRUE)
    check("all_of(a, b)", {"a": P, "b": A}, FALSE)
    check("n_of(2, a, b, c)", {"a": A, "b": A, "c": A}, FALSE)

    # Precedence: and binds tighter than or
    check("a or b and c", {"a": A, "b": P, "c": A}, FALSE)
    check("(a or b) and c", {"a": P, "b": A, "c": A}, FALSE)
    check("a and (b or c)", {"a": P, "b": A, "c": P}, TRUE)

    # Parser must reject junk rather than misinterpret it
    for bad in ["a and", "and a", "(a", "a)", "present()", "n_of(a, b)",
                "present(a, b)", "a b", ""]:
        try:
            compile_rule(bad)
            failures.append(f"parser accepted invalid expression {bad!r}")
        except RuleError:
            pass

    # No code execution
    for evil in ["__import__('os')", "a or __import__"]:
        try:
            fn, _ = compile_rule(evil)
            fn({})  # treated as plain identifiers, never executed
        except RuleError:
            pass

    if failures:
        sys.stderr.write(f"SELF-TEST FAILED ({len(failures)}):\n")
        for f in failures:
            sys.stderr.write(f"  {f}\n")
        return 1
    sys.stderr.write("marker_rules self-test: all checks passed\n")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--rules", type=Path, help="marker_rules.yaml")
    ap.add_argument("--markers", type=Path, help="markers.tsv, for validation")
    ap.add_argument("--watchlist", type=Path,
                    help="watchlist.resolved.tsv. With --validate, cross-checks that "
                         "rule-block names and catalogue taxon names are watchlist "
                         "entries. The join between these three files is held together "
                         "by hand-maintained literal strings with nothing checking it.")
    ap.add_argument("--validate", action="store_true",
                    help="compile all rules and check marker references")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    if not args.rules:
        ap.error("--rules is required unless --self-test")

    rules = load_rules(args.rules)
    sys.stderr.write(f"[rules] loaded {len(rules)} taxon rule set(s) from {args.rules}\n")

    if args.validate:
        known: set[str] = set()
        if args.markers:
            with args.markers.open() as fh:
                lines = [ln for ln in fh if not ln.startswith("#")]
            known = {r["marker_id"] for r in csv.DictReader(lines, delimiter="\t")
                     if r.get("marker_id")}
            sys.stderr.write(f"[rules] {len(known)} markers in catalogue\n")
        problems = validate_rules(rules, known) if known else []

        # THE THREE FILES AGREE ON LITERAL STRINGS AND NOTHING CHECKS THEM.
        # markers.tsv.taxon, marker_rules.yaml's keys and watchlist.resolved.tsv.name
        # must name the same organisms. They do today - 58 exact matches, verified
        # byte-for-byte - but a typo in any one of them degrades SILENTLY: a rule block
        # whose key no longer matches a watchlist name simply stops gating, which looks
        # exactly like a taxon with no rules, which is D17a's permissive path. A marker
        # whose taxon string drifts cannot ever raise or gate anything for its organism.
        # Reported, not fatal, for the two asymmetries that exist on purpose.
        if args.watchlist and args.watchlist.exists():
            with args.watchlist.open() as fh:
                wl = {r["name"].strip()
                      for r in csv.DictReader(fh, delimiter="\t") if r.get("name")}
            cat_taxa: set[str] = set()
            if args.markers:
                with args.markers.open() as fh:
                    lines = [ln for ln in fh if not ln.startswith("#")]
                cat_taxa = {(r.get("taxon") or "").strip()
                            for r in csv.DictReader(lines, delimiter="\t")
                            if (r.get("taxon") or "").strip()}
            # Confuser taxa (B. cereus, Y. pseudotuberculosis, ...) carry markers on
            # purpose and are deliberately NOT watchlist entries - their markers exist
            # to veto, never to call. So a catalogue taxon outside the watchlist is
            # normal; a RULE BLOCK outside it is not, because a block keyed on a name
            # no watchlist row carries can never fire.
            orphan_blocks = sorted(set(rules) - wl)
            for t in orphan_blocks:
                problems.append(
                    f"rule block {t!r} is not a watchlist name - it can never fire")
            no_own_marker = sorted(set(rules) - cat_taxa) if cat_taxa else []
            for t in no_own_marker:
                sys.stderr.write(
                    f"[rules] NOTE: {t!r} has a rule block but no marker of its own; "
                    f"its rules reference another taxon's markers, so a PRESENT call "
                    f"there raises this taxon too\n")
            no_block = sorted((cat_taxa & wl) - set(rules)) if cat_taxa else []
            for t in no_block:
                sys.stderr.write(
                    f"[rules] NOTE: {t!r} is on the watchlist and has marker(s) but no "
                    f"rule block, so marker evidence for it gates nothing (D17a)\n")
        if problems:
            sys.stderr.write(f"\nERROR: {len(problems)} rule problem(s):\n")
            for p in problems:
                sys.stderr.write(f"  {p}\n")
            return 1
        sys.stderr.write("[rules] all rules compile and reference known markers\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
