#!/usr/bin/env python3
"""
Restricted YAML subset parser - stdlib only.

WHY THIS EXISTS, AND WHY IT IS DELIBERATELY SMALL

The runtime scripts in bin/ are stdlib-only so they run in whatever environment
the pipeline hands them, and PyYAML is not guaranteed present. The pack, rule and
benchmark-scenario files need comments (they encode reasoning that reviewers must
read), which rules out JSON. Three ad-hoc line scrapers were already accumulating
across this codebase; this replaces them with one parser that has a test suite.

The governing design rule: **fail loudly rather than mis-parse quietly.**

A config parser that silently misreads a threshold, drops a duplicated key, or
ignores a construct it does not understand is far more dangerous here than one
that refuses to run - these files gate high-consequence calls. So every YAML
feature outside the supported subset raises BtfYamlError naming the line, instead
of being skipped or guessed at.

SUPPORTED
  nested block mappings and sequences (indentation-based)
  sequence items: `- scalar`, `- key: value`, `-` + nested block
  flow collections on one line: [a, b, c] and {a: 1, b: 2}
  block scalars: `>` (folded) and `|` (literal), with `-` / `+` chomping
  scalars: int, float, true/false, null/~, single- and double-quoted strings
  comments: whole-line and trailing (quote-aware)

REJECTED, with an explicit error
  tab indentation, anchors (&), aliases (*), tags (!), merge keys (<<),
  multiple documents (---/...), explicit complex keys (?), duplicate mapping keys

Duplicate keys are an error rather than last-one-wins because a duplicated key in
a rules or scenario file means a silently discarded rule.
"""

from __future__ import annotations

import re
from pathlib import Path

__all__ = ["load", "loads", "BtfYamlError"]

KEY_RE = re.compile(r"^([A-Za-z_][\w.\-]*|'[^']*'|\"[^\"]*\")\s*:(\s|$)")
BLOCK_SCALAR_RE = re.compile(r"^([>|])([+-]?)$")


class BtfYamlError(ValueError):
    """Raised for anything outside the supported subset, with line context."""


class _Line:
    __slots__ = ("indent", "text", "no")

    def __init__(self, indent: int, text: str, no: int):
        self.indent, self.text, self.no = indent, text, no

    def __repr__(self):
        return f"_Line({self.indent}, {self.text!r}, line {self.no})"


def _split_comment(s: str) -> str:
    """Strip a trailing comment, ignoring '#' inside quotes.

    AN APOSTROPHE IS NOT AN OPENING QUOTE. The previous version toggled its
    in-string flag on every `\'`, so one apostrophe in a prose value left the
    scanner "inside a string" for the rest of the line and the trailing comment was
    never removed: `notes: Yersinia's plasmid content # see D3` parsed as the value
    "Yersinia's plasmid content # see D3". These pack files are full of prose
    notes:/rationale: fields that flow into the report's caveat column, and on a line
    whose value was numeric the glued result would become a string - a string where a
    threshold was expected. A quote only opens a scalar where a scalar may start, so
    that is what is tested here.
    """
    in_s = in_d = False
    openers = " \t:,[{-"
    for i, c in enumerate(s):
        may_open = (i == 0) or (s[i - 1] in openers)
        if c == "'" and not in_d:
            if in_s or may_open:
                in_s = not in_s
        elif c == '"' and not in_s:
            if in_d or may_open:
                in_d = not in_d
        elif c == "#" and not in_s and not in_d:
            if i == 0 or s[i - 1] in " \t":
                return s[:i].rstrip()
    return s.rstrip()


def _err(msg: str, line: _Line | None = None):
    where = f" (line {line.no})" if line else ""
    raise BtfYamlError(f"{msg}{where}")


def _split_flow(body: str, line: _Line) -> list[str]:
    """Split a flow collection body on top-level commas, quote- and nest-aware."""
    parts, buf, depth = [], [], 0
    in_s = in_d = False
    for c in body:
        if c == "'" and not in_d:
            in_s = not in_s
        elif c == '"' and not in_s:
            in_d = not in_d
        if not in_s and not in_d:
            if c in "[{":
                depth += 1
            elif c in "]}":
                depth -= 1
                if depth < 0:
                    _err("unbalanced brackets in flow collection", line)
            elif c == "," and depth == 0:
                parts.append("".join(buf))
                buf = []
                continue
        buf.append(c)
    if in_s or in_d:
        _err("unterminated quoted string in flow collection", line)
    if depth != 0:
        _err("unbalanced brackets in flow collection", line)
    tail = "".join(buf).strip()
    if tail:
        parts.append(tail)
    return [p.strip() for p in parts if p.strip() != ""]


def _scalar(text: str, line: _Line):
    t = text.strip()
    if t == "":
        return None

    if t[0] in "&*!":
        kind = {"&": "anchors", "*": "aliases", "!": "tags"}[t[0]]
        _err(f"YAML {kind} are not supported", line)
    if t.startswith("<<"):
        _err("YAML merge keys (<<) are not supported", line)

    if t[0] in "'\"":
        if len(t) >= 2 and t[-1] == t[0]:
            return t[1:-1]
        _err("unterminated quoted string", line)

    if t.startswith("["):
        if not t.endswith("]"):
            _err("flow sequence must close on the same line", line)
        return [_scalar(p, line) for p in _split_flow(t[1:-1], line)]

    if t.startswith("{"):
        if not t.endswith("}"):
            _err("flow mapping must close on the same line", line)
        out = {}
        for part in _split_flow(t[1:-1], line):
            if ":" not in part:
                _err(f"flow mapping entry {part!r} has no ':'", line)
            k, v = part.split(":", 1)
            key = _scalar(k, line)
            if key in out:
                _err(f"duplicate key {key!r} in flow mapping", line)
            out[key] = _scalar(v, line)
        return out

    if t in ("null", "Null", "NULL", "~"):
        return None
    if t in ("true", "True", "TRUE"):
        return True
    if t in ("false", "False", "FALSE"):
        return False

    # Ints and floats, but never strip a leading zero (e.g. a barcode "007"
    # must stay a string) and never accept YAML 1.1 sexagesimals.
    if re.fullmatch(r"[+-]?\d+", t):
        if not (len(t.lstrip("+-")) > 1 and t.lstrip("+-").startswith("0")):
            return int(t)
    if re.fullmatch(r"[+-]?(\d+\.\d*|\.\d+)([eE][+-]?\d+)?", t) or \
       re.fullmatch(r"[+-]?\d+[eE][+-]?\d+", t):
        return float(t)

    return t


def _prepare(text: str) -> list[_Line]:
    lines = []
    for no, raw in enumerate(text.splitlines(), start=1):
        if raw.strip() in ("---", "..."):
            _err("multi-document YAML is not supported", _Line(0, raw.strip(), no))
        stripped_full = _split_comment(raw)
        if stripped_full.strip() == "":
            continue
        expanded = stripped_full.expandtabs(8)
        indent = len(expanded) - len(expanded.lstrip(" "))
        if "\t" in raw[:len(raw) - len(raw.lstrip())]:
            _err("tab characters in indentation are not supported; use spaces",
                 _Line(indent, stripped_full.strip(), no))
        lines.append(_Line(indent, expanded.strip(), no))
    return lines


class _Parser:
    def __init__(self, text: str):
        self.raw_lines = text.splitlines()
        self.lines = _prepare(text)
        self.i = 0

    def at_end(self) -> bool:
        return self.i >= len(self.lines)

    def cur(self) -> _Line:
        return self.lines[self.i]

    # -- block scalars ---------------------------------------------------
    def _block_scalar(self, style: str, chomp: str, parent_indent: int) -> str:
        """Consume the raw source lines belonging to a > or | block."""
        start_no = self.lines[self.i - 1].no if self.i > 0 else 1
        body: list[str] = []
        base: int | None = None
        n = start_no  # 1-based index into raw_lines of the key line
        while n < len(self.raw_lines):
            raw = self.raw_lines[n]
            if raw.strip() == "":
                body.append("")
                n += 1
                continue
            expanded = raw.expandtabs(8)
            indent = len(expanded) - len(expanded.rstrip("\n").lstrip(" "))
            if indent <= parent_indent:
                break
            if base is None:
                base = indent
            body.append(expanded[base:].rstrip())
            n += 1

        # Re-sync the logical line cursor past the consumed block.
        while not self.at_end() and self.cur().no < n + 1:
            self.i += 1

        while body and body[-1] == "":
            body.pop()
        if style == "|":
            out = "\n".join(body)
        else:
            # Folded: blank lines become paragraph breaks, others join with a
            # space. This is what makes the `rationale:` and `notes:` fields in
            # pack files read as prose.
            paras, cur = [], []
            for ln in body:
                if ln == "":
                    if cur:
                        paras.append(" ".join(cur))
                        cur = []
                else:
                    cur.append(ln)
            if cur:
                paras.append(" ".join(cur))
            out = "\n\n".join(paras)
        if chomp == "-":
            return out
        if chomp == "+":
            return out + "\n"
        return out

    # -- blocks ----------------------------------------------------------
    def parse_block(self, indent: int):
        if self.at_end():
            return None
        if self.cur().indent < indent:
            return None
        if self.cur().text.startswith("- ") or self.cur().text == "-":
            return self.parse_sequence(indent)
        return self.parse_mapping(indent)

    def parse_sequence(self, indent: int) -> list:
        out = []
        while not self.at_end():
            line = self.cur()
            if line.indent < indent:
                break
            if line.indent > indent:
                _err(f"unexpected indentation inside sequence (expected {indent})", line)
            if not (line.text.startswith("- ") or line.text == "-"):
                break

            content = line.text[1:].strip()
            item_indent = indent + 2
            self.i += 1

            if content == "":
                out.append(self.parse_block(item_indent))
                continue

            if KEY_RE.match(content):
                # `- key: value` opens a mapping whose first key is inline.
                out.append(self.parse_mapping(item_indent, inline=content,
                                              inline_line=line))
                continue

            # A ONE-ELEMENT MAPPING MUST NOT DEGRADE INTO A STRING.
            #
            # KEY_RE does not admit unquoted spaces, so `- Bacillus anthracis: 5`
            # fell through to the plain-scalar branch and silently became the STRING
            # "Bacillus anthracis: 5". A consumer doing item.get(...) then finds
            # nothing, or raises. parse_mapping() already errors on the identical
            # construct, so the parser disagreed with itself depending on context.
            # This file's contract is to fail loudly rather than mis-parse quietly.
            if ": " in content or content.endswith(":"):
                _err(f"sequence item looks like a mapping with an unquoted key "
                     f"containing spaces: {content!r}. A plain YAML scalar cannot "
                     f"contain ': '. Quote the key, e.g. "
                     f'- "Bacillus anthracis": 5', line)

            m = BLOCK_SCALAR_RE.match(content)
            if m:
                out.append(self._block_scalar(m.group(1), m.group(2), indent))
                continue

            out.append(_scalar(content, line))
        return out

    def parse_mapping(self, indent: int, inline: str | None = None,
                      inline_line: _Line | None = None) -> dict:
        out: dict = {}

        def add(key, value, line):
            if key in out:
                _err(f"duplicate key {key!r} in mapping", line)
            out[key] = value

        if inline is not None:
            assert inline_line is not None
            key, rest = self._split_kv(inline, inline_line)
            add(key, self._value_for(rest, indent, inline_line), inline_line)

        while not self.at_end():
            line = self.cur()
            if line.indent < indent:
                break
            if line.indent > indent:
                _err(f"unexpected indentation (expected {indent}, got {line.indent})", line)
            if line.text.startswith("- ") or line.text == "-":
                break
            if line.text.startswith("?"):
                _err("explicit complex keys (?) are not supported", line)
            if not KEY_RE.match(line.text):
                _err(f"expected 'key: value', got {line.text!r}", line)

            key, rest = self._split_kv(line.text, line)
            self.i += 1
            add(key, self._value_for(rest, indent, line), line)
        return out

    def _split_kv(self, text: str, line: _Line) -> tuple:
        m = KEY_RE.match(text)
        if not m:
            _err(f"expected 'key: value', got {text!r}", line)
        raw_key = m.group(1)
        key = raw_key[1:-1] if raw_key[0] in "'\"" else raw_key
        rest = text[len(raw_key):].lstrip()
        if rest.startswith(":"):
            rest = rest[1:].strip()
        return key, rest

    def _value_for(self, rest: str, indent: int, line: _Line):
        m = BLOCK_SCALAR_RE.match(rest)
        if m:
            return self._block_scalar(m.group(1), m.group(2), indent)
        if rest == "":
            # Nested block, or an explicitly empty value.
            if not self.at_end() and self.cur().indent > indent:
                return self.parse_block(self.cur().indent)
            return None
        return _scalar(rest, line)


def loads(text: str):
    p = _Parser(text)
    if p.at_end():
        return {}
    value = p.parse_block(p.cur().indent)
    if not p.at_end():
        _err("could not parse to end of document", p.cur())
    return value


def load(path):
    path = Path(path)
    try:
        text = path.read_text()
    except OSError as exc:
        raise BtfYamlError(f"cannot read {path}: {exc}") from exc
    try:
        return loads(text)
    except BtfYamlError as exc:
        raise BtfYamlError(f"{path}: {exc}") from None


# ---------------------------------------------------------------------------

def _self_test() -> int:
    import sys
    fails = []

    def eq(text, expected, label):
        try:
            got = loads(text)
        except BtfYamlError as exc:
            fails.append(f"{label}: unexpected error {exc}")
            return
        if got != expected:
            fails.append(f"{label}: got {got!r}, expected {expected!r}")

    def rejects(text, label):
        try:
            loads(text)
        except BtfYamlError:
            return
        fails.append(f"{label}: should have been rejected but parsed")

    eq("a: 1\nb: two\n", {"a": 1, "b": "two"}, "flat map")
    eq("a:\n  b:\n    c: 3\n", {"a": {"b": {"c": 3}}}, "nested maps")
    eq("a: [1, 2, 3]\n", {"a": [1, 2, 3]}, "flow seq")
    eq("a: {x: 1, y: b}\n", {"a": {"x": 1, "y": "b"}}, "flow map")
    eq("a:\n  - 1\n  - 2\n", {"a": [1, 2]}, "block seq")
    eq("a:\n  - x: 1\n    y: 2\n  - x: 3\n    y: 4\n",
       {"a": [{"x": 1, "y": 2}, {"x": 3, "y": 4}]}, "seq of maps")
    eq("a:\n  - x: 1\n    z: [p, q]\n", {"a": [{"x": 1, "z": ["p", "q"]}]},
       "seq of maps with flow")
    eq("t: true\nf: false\nn: null\ntil: ~\n",
       {"t": True, "f": False, "n": None, "til": None}, "specials")
    eq("f: 1.5\ng: 1e3\nh: -2\n", {"f": 1.5, "g": 1000.0, "h": -2}, "numbers")
    eq('q: "has: colon"\nr: \'and # hash\'\n',
       {"q": "has: colon", "r": "and # hash"}, "quotes protect colon and hash")
    eq("a: 1  # trailing\n# whole line\nb: 2\n", {"a": 1, "b": 2}, "comments")
    eq("z: '007'\n", {"z": "007"}, "leading zero stays a string")

    # Block scalars - the reason JSON was not an option
    eq("a: >\n  one\n  two\nb: 2\n", {"a": "one two", "b": 2}, "folded scalar")
    eq("a: |\n  one\n  two\nb: 2\n", {"a": "one\ntwo", "b": 2}, "literal scalar")
    eq("a: >\n  p1 line1\n  p1 line2\n\n  p2\n", {"a": "p1 line1 p1 line2\n\np2"},
       "folded paragraphs")

    # The real shapes this parser must handle
    eq('"Bacillus anthracis":\n    confirm_requires: (a or b) and (c or d)\n'
       '    rationale: >\n      why it matters\n      across lines\n',
       {"Bacillus anthracis": {"confirm_requires": "(a or b) and (c or d)",
                               "rationale": "why it matters across lines"}},
       "marker_rules shape")
    eq("scenarios:\n  - id: s1\n    spikes:\n      - taxon: X\n        levels_rpm: [1, 10]\n"
       "    truth:\n      present: [X]\n",
       {"scenarios": [{"id": "s1", "spikes": [{"taxon": "X", "levels_rpm": [1, 10]}],
                       "truth": {"present": ["X"]}}]},
       "scenario shape")

    # Must refuse rather than guess
    rejects("a: 1\na: 2\n", "duplicate key")
    rejects("a: {x: 1, x: 2}\n", "duplicate key in flow map")
    rejects("a:\n\t b: 1\n", "tab indentation")
    rejects("a: &anchor 1\n", "anchor")
    rejects("b: *alias\n", "alias")
    rejects("c: !!str 1\n", "tag")
    rejects("---\na: 1\n", "multi-document")
    rejects("? complex\n", "complex key")
    rejects("a: [1, 2\n", "unclosed flow seq")
    rejects("a: 'unterminated\n", "unterminated quote")
    rejects("a: 1\n  b: 2\n", "bad indentation")
    rejects("just a scalar\n", "top-level scalar is not a mapping")
    rejects("a:\n  - 1\n   - 2\n", "ragged sequence indentation")

    if fails:
        sys.stderr.write(f"btf_yaml self-test FAILED ({len(fails)}):\n")
        for f in fails:
            sys.stderr.write(f"  {f}\n")
        return 1
    sys.stderr.write("btf_yaml self-test: all checks passed\n")
    return 0


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        sys.exit(_self_test())
    if len(sys.argv) > 1:
        import json
        print(json.dumps(load(sys.argv[1]), indent=2, default=str))
    else:
        sys.stderr.write(__doc__ + "\nUsage: btf_yaml.py [--self-test | FILE]\n")
        sys.exit(2)
