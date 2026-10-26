#!/usr/bin/env python3
"""Merge multi-sample eccDNA calls from CoReID, CReSIL, NanoCircle, Flec, or FLED.

The program uses one shared coordinate-normalisation, coordinate-clustering,
gene-annotation, and matrix-writing workflow.  Only the input parser changes
with ``--method``.

Manifest format (tab-delimited, no header):
    SampleID<TAB>file_role<TAB>path

``file_role`` is ``main`` for CoReID, CReSIL, and Flec; it is ``simple`` or
``complex`` for FLED and NanoCircle.  Blank lines and lines beginning with
``#`` are ignored.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial
from typing import DefaultDict, Sequence

import pandas as pd


METHOD_ALIASES = {
    "coreid": "CoReID",
    "cresil": "CReSIL",
    "nanocircle": "NanoCircle",
    "flec": "Flec",
    "fled": "FLED",
}
SINGLE_FILE_METHODS = {"CoReID", "CReSIL", "Flec"}
TWO_FILE_METHODS = {"FLED", "NanoCircle"}


@dataclass(frozen=True)
class Segment:
    """One genomic segment in a circular eccDNA structure."""

    chrom: str
    start: int
    end: int
    strand: str


EccDNA = tuple[Segment, ...]


def canonical_chromosome(chrom: object) -> str:
    """Return chromosome names consistently as ``chr*``."""
    value = str(chrom).strip()
    if not value:
        raise ValueError("empty chromosome name")
    if value.lower().startswith("chr"):
        value = value[3:]
    return f"chr{value}"


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    """Natural sort key for standard and non-standard contigs."""
    clean = canonical_chromosome(chrom)[3:]
    if clean.isdigit():
        return (0, int(clean))
    standard = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if clean.upper() in standard:
        return (0, standard[clean.upper()])
    return (1, clean)


def flip_strand(strand: str) -> str:
    return {"+": "-", "-": "+"}.get(strand, strand)


def canonicalise_eccdna(segments: Sequence[Segment]) -> EccDNA:
    """Normalise circular rotation and reverse-complement orientation.

    The rule matches the previous merge scripts: rotate to the segment with
    the smallest chromosome/start key.  If that first segment is on the minus
    strand, retain it as the start and reverse the remaining circular order.
    Unstranded calls (``.``) are intentionally kept unchanged.
    """
    if not segments:
        raise ValueError("eccDNA has no valid segments")

    cleaned: list[Segment] = []
    for segment in segments:
        if segment.start >= segment.end:
            raise ValueError(f"invalid interval: {segment.chrom}:{segment.start}-{segment.end}")
        strand = segment.strand if segment.strand in {"+", "-"} else "."
        cleaned.append(Segment(canonical_chromosome(segment.chrom), segment.start, segment.end, strand))

    base = min(range(len(cleaned)), key=lambda i: (chrom_sort_key(cleaned[i].chrom), cleaned[i].start))
    rotated = cleaned[base:] + cleaned[:base]
    if rotated[0].strand == "-":
        rotated = [
            Segment(rotated[0].chrom, rotated[0].start, rotated[0].end, "+"),
            *[
                Segment(s.chrom, s.start, s.end, flip_strand(s.strand))
                for s in reversed(rotated[1:])
            ],
        ]
    return tuple(rotated)


def format_eccdna(eccdna: EccDNA) -> str:
    """Human-readable, stable identifier used as the coordinate-matrix row."""
    return ",".join(
        f"{s.chrom}:{s.start}-{s.end}_{s.strand}" if s.strand != "." else f"{s.chrom}:{s.start}-{s.end}_."
        for s in eccdna
    )


def parse_coreid_or_cresil(path: str, has_header: bool) -> list[tuple[EccDNA, int]]:
    """Parse CoReID/CReSIL merged calls (``merge_region`` and ``numreads``)."""
    header = 0 if has_header else None
    frame = pd.read_csv(path, sep="\t", header=header)
    results: list[tuple[EccDNA, int]] = []
    for _, row in frame.iterrows():
        try:
            region = row["merge_region"] if has_header else row.iloc[1]
            reads = int(float(row["numreads"] if has_header else row.iloc[5]))
            segments = []
            for text in str(region).split(","):
                chrom_part, strand = text.strip().rsplit("_", 1)
                chrom, positions = chrom_part.split(":", 1)
                start, end = positions.split("-", 1)
                segments.append(Segment(chrom, int(start), int(end), strand))
            results.append((canonicalise_eccdna(segments), reads))
        except (IndexError, KeyError, TypeError, ValueError):
            continue
    return results


def parse_flec(path: str, has_header: bool) -> list[tuple[EccDNA, int]]:
    """Parse Flec calls, retaining records with ``Nfullpass > 0``.

    Each retained Flec record contributes one read, as in the original script.
    """
    header = 0 if has_header else None
    frame = pd.read_csv(path, sep="\t", header=header)
    results: list[tuple[EccDNA, int]] = []
    for _, row in frame.iterrows():
        try:
            if float(row.iloc[1]) <= 0:
                continue
            segments = []
            for text in str(row.iloc[5]).split(","):
                chrom, position_and_strand = text.strip().split(":", 1)
                start, end_and_strand = position_and_strand.split("-", 1)
                end, strand = end_and_strand.split("(", 1)
                segments.append(Segment(chrom, int(start), int(end), strand.rstrip(")")))
            results.append((canonicalise_eccdna(segments), 1))
        except (IndexError, TypeError, ValueError):
            continue
    return results


def parse_fled_simple(path: str) -> list[tuple[EccDNA, int]]:
    """Parse FLED simple calls (reads are columns 7 + 8, one-based)."""
    frame = pd.read_csv(path, sep="\t", header=None)
    results: list[tuple[EccDNA, int]] = []
    for _, row in frame.iterrows():
        try:
            segment = Segment(row.iloc[0], int(row.iloc[1]), int(row.iloc[2]), str(row.iloc[3]))
            reads = int(float(row.iloc[6])) + int(float(row.iloc[7]))
            results.append((canonicalise_eccdna([segment]), reads))
        except (IndexError, TypeError, ValueError):
            continue
    return results


def parse_fled_complex(path: str) -> list[tuple[EccDNA, int]]:
    """Parse FLED complex calls (reads are columns 4 + 5, one-based).

    FLED complex output does not retain strand.  For backward compatibility
    with the supplied FLED merger, these calls are represented internally with
    a ``+`` strand.
    """
    frame = pd.read_csv(path, sep="\t", header=None)
    results: list[tuple[EccDNA, int]] = []
    for _, row in frame.iterrows():
        try:
            reads = int(float(row.iloc[3])) + int(float(row.iloc[4]))
            segments = []
            for text in str(row.iloc[0]).split("|"):
                chrom, positions = text.strip().split(":", 1)
                start, end = positions.split("-", 1)
                segments.append(Segment(chrom, int(start), int(end), "+"))
            results.append((canonicalise_eccdna(segments), reads))
        except (IndexError, TypeError, ValueError):
            continue
    return results


def parse_nanocircle_simple(path: str) -> list[tuple[EccDNA, int]]:
    """Parse NanoCircle simple calls (read count is column 5, one-based)."""
    frame = pd.read_csv(path, sep="\t", header=None)
    results: list[tuple[EccDNA, int]] = []
    for _, row in frame.iterrows():
        try:
            start, end, length = int(row.iloc[1]), int(row.iloc[2]), int(row.iloc[3])
            if length <= 0 or start >= end:
                continue
            results.append((canonicalise_eccdna([Segment(row.iloc[0], start, end, "+")]), int(float(row.iloc[4]))))
        except (IndexError, TypeError, ValueError):
            continue
    return results


def parse_nanocircle_complex(path: str) -> list[tuple[EccDNA, int]]:
    """Parse NanoCircle complex calls (read count is the third-last column)."""
    frame = pd.read_csv(path, sep="\t", header=None)
    results: list[tuple[EccDNA, int]] = []
    for _, row in frame.iterrows():
        try:
            reads = int(float(row.iloc[-3]))
            segments = []
            # Segment columns repeat as chrom, start, end, length before 3 trailing columns.
            for index in range(0, len(row) - 3, 4):
                chrom = str(row.iloc[index])
                if chrom in {".", "nan", ""}:
                    break
                start, end, length = int(row.iloc[index + 1]), int(row.iloc[index + 2]), int(row.iloc[index + 3])
                if length > 0 and start < end:
                    segments.append(Segment(chrom, start, end, "+"))
            if segments:
                results.append((canonicalise_eccdna(segments), reads))
        except (IndexError, TypeError, ValueError):
            continue
    return results


def parse_file(method: str, role: str, path: str, has_header: bool) -> list[tuple[EccDNA, int]]:
    """Dispatch to the correct method- and role-specific parser."""
    if not os.path.isfile(path):
        print(f"[WARN] Missing file; skipped: {path}", file=sys.stderr)
        return []
    if os.path.getsize(path) == 0:
        print(f"[WARN] Empty file; skipped: {path}", file=sys.stderr)
        return []
    try:
        if method in {"CoReID", "CReSIL"}:
            return parse_coreid_or_cresil(path, has_header)
        if method == "Flec":
            return parse_flec(path, has_header)
        if method == "FLED" and role == "simple":
            return parse_fled_simple(path)
        if method == "FLED" and role == "complex":
            return parse_fled_complex(path)
        if method == "NanoCircle" and role == "simple":
            return parse_nanocircle_simple(path)
        if method == "NanoCircle" and role == "complex":
            return parse_nanocircle_complex(path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
        print(f"[WARN] Cannot parse {path}: {exc}", file=sys.stderr)
        return []
    raise ValueError(f"unsupported method/role pair: {method}/{role}")


def parse_manifest_entry(
    entry: tuple[str, str, str], method: str, has_header: bool
) -> tuple[str, list[tuple[EccDNA, int]]]:
    """Pickle-safe worker entry point for parallel manifest parsing."""
    sample, role, path = entry
    return sample, parse_file(method, role, path, has_header)


def load_manifest(path: str, method: str) -> list[tuple[str, str, str]]:
    """Read and validate the three-column input manifest."""
    manifest_dir = os.path.dirname(os.path.abspath(path))
    entries: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    allowed_roles = {"main"} if method in SINGLE_FILE_METHODS else {"simple", "complex"}
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 3:
                raise ValueError(f"manifest line {line_number}: expected exactly 3 tab-delimited columns")
            sample, role, filename = (x.strip() for x in fields)
            if not sample or role not in allowed_roles or not filename:
                roles = ", ".join(sorted(allowed_roles))
                raise ValueError(f"manifest line {line_number}: role must be one of: {roles}")
            if (sample, role) in seen:
                raise ValueError(f"manifest line {line_number}: duplicate SampleID/role pair: {sample}/{role}")
            seen.add((sample, role))
            if not os.path.isabs(filename):
                filename = os.path.join(manifest_dir, filename)
            entries.append((sample, role, os.path.abspath(filename)))
    if not entries:
        raise ValueError("manifest contains no usable records")

    by_sample: DefaultDict[str, set[str]] = defaultdict(set)
    for sample, role, _ in entries:
        by_sample[sample].add(role)
    if method in TWO_FILE_METHODS:
        incomplete = sorted(sample for sample, roles in by_sample.items() if roles != {"simple", "complex"})
        if incomplete:
            raise ValueError("FLED/NanoCircle require both simple and complex entries for each sample: " + ", ".join(incomplete))
    return entries


def interval_similarity(a: Segment, b: Segment) -> float:
    overlap = min(a.end, b.end) - max(a.start, b.start)
    denominator = max(a.end - a.start, b.end - b.start)
    return overlap / denominator if denominator > 0 else -1.0


def cluster_eccdnas(calls: dict[EccDNA, dict[str, int]], threshold: float, max_start_distance: int) -> dict[EccDNA, EccDNA]:
    """Cluster coordinate-near calls with identical segment topology.

    A call joins the first nearby representative for which every corresponding
    segment has overlap / max(segment length) >= ``threshold``.  This retains
    the clustering criterion used by the five original scripts.
    """
    groups: DefaultDict[tuple[tuple[str, str], ...], list[EccDNA]] = defaultdict(list)
    for eccdna in calls:
        groups[tuple((s.chrom, s.strand) for s in eccdna)].append(eccdna)

    mapping: dict[EccDNA, EccDNA] = {}
    for topology, candidates in groups.items():
        candidates.sort(key=lambda item: item[0].start)
        representatives: list[EccDNA] = []
        for candidate in candidates:
            matched: EccDNA | None = None
            for representative in reversed(representatives):
                if candidate[0].start - representative[0].start > max_start_distance:
                    break
                if all(interval_similarity(x, y) >= threshold for x, y in zip(candidate, representative)):
                    matched = representative
                    break
            if matched is None:
                representatives.append(candidate)
                matched = candidate
            mapping[candidate] = matched
    return mapping


def load_gene_annotation(path: str) -> dict[str, list[tuple[int, int, str]]]:
    """Load the first four columns of a BED/BED12 gene annotation file."""
    genes: DefaultDict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"annotation line {line_number}: expected at least four BED columns")
            try:
                genes[canonical_chromosome(fields[0])].append((int(fields[1]), int(fields[2]), fields[3]))
            except ValueError as exc:
                raise ValueError(f"annotation line {line_number}: invalid BED coordinates") from exc
    for chromosome in genes:
        genes[chromosome].sort(key=lambda item: item[0])
    return dict(genes)


def overlapping_genes(eccdna: EccDNA, genes: dict[str, list[tuple[int, int, str]]]) -> set[str]:
    hits: set[str] = set()
    for segment in eccdna:
        for gene_start, gene_end, gene_name in genes.get(segment.chrom, []):
            if gene_start >= segment.end:
                break
            if max(segment.start, gene_start) < min(segment.end, gene_end):
                hits.add(gene_name)
    return hits


def write_matrix(
    merged_calls: dict[EccDNA, dict[str, int]],
    sample_labels: list[str],
    level: str,
    annotation: str | None,
    out_path: str,
) -> None:
    """Write either an eccDNA-coordinate or gene-by-sample count matrix."""
    output: DefaultDict[str, dict[str, int]] = defaultdict(dict)
    if level == "eccDNA":
        for eccdna, sample_counts in merged_calls.items():
            output[format_eccdna(eccdna)] = dict(sample_counts)
        index_name = "eccDNA_ID"
    else:
        if annotation is None:
            raise ValueError("gene-level output requires --annotation")
        gene_db = load_gene_annotation(annotation)
        for eccdna, sample_counts in merged_calls.items():
            for gene_name in overlapping_genes(eccdna, gene_db):
                for sample, count in sample_counts.items():
                    output[gene_name][sample] = output[gene_name].get(sample, 0) + count
        index_name = "Gene_Name"

    frame = pd.DataFrame.from_dict(output, orient="index")
    frame = frame.reindex(columns=sample_labels, fill_value=0).fillna(0).astype(int)
    frame.index.name = index_name
    frame.sort_index().to_csv(out_path, sep="\t")


def run(args: argparse.Namespace) -> None:
    method = METHOD_ALIASES[args.method.lower()]
    if not 0 < args.overlap <= 1:
        raise ValueError("--overlap must be in (0, 1]")
    if args.level == "gene" and not args.annotation:
        raise ValueError("--annotation is required when --level gene is selected")
    if args.annotation and not os.path.isfile(args.annotation):
        raise ValueError(f"annotation file does not exist: {args.annotation}")

    entries = load_manifest(args.input, method)
    sample_labels = sorted({sample for sample, _, _ in entries})
    print(f"[INFO] Method: {method}; samples: {len(sample_labels)}; input files: {len(entries)}")

    all_calls: DefaultDict[EccDNA, dict[str, int]] = defaultdict(dict)
    # ``executor.map`` preserves manifest order and reports failures directly.
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        worker = partial(parse_manifest_entry, method=method, has_header=args.header == "T")
        for sample, parsed_calls in executor.map(worker, entries):
            for eccdna, reads in parsed_calls:
                all_calls[eccdna][sample] = all_calls[eccdna].get(sample, 0) + reads
    print(f"[INFO] Parsed {len(all_calls)} unique coordinate records before clustering.")

    mapping = cluster_eccdnas(dict(all_calls), args.overlap, args.max_start_distance)
    merged: DefaultDict[EccDNA, dict[str, int]] = defaultdict(dict)
    for eccdna, sample_counts in all_calls.items():
        representative = mapping[eccdna]
        for sample, count in sample_counts.items():
            merged[representative][sample] = merged[representative].get(sample, 0) + count
    print(f"[INFO] Retained {len(merged)} coordinate clusters.")

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    suffix = "gene_matrix.tsv" if args.level == "gene" else "eccDNA_matrix.tsv"
    out_path = os.path.join(out_dir, f"{args.prefix}_{method}_{suffix}")
    write_matrix(dict(merged), sample_labels, args.level, args.annotation, out_path)
    print(f"[INFO] Matrix written: {out_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Merge eccDNA calls from CoReID, CReSIL, NanoCircle, Flec, or FLED."
    )
    parser.add_argument("-m", "--method", required=True, choices=sorted(METHOD_ALIASES.values()),
                        help="Method that generated every file in the manifest.")
    parser.add_argument("-i", "--input", required=True, help="Three-column, tab-delimited manifest.")
    parser.add_argument("-l", "--level", choices=["eccDNA", "gene"], default="eccDNA",
                        help="Output level (default: eccDNA).")
    parser.add_argument("-a", "--annotation", help="BED/BED12 gene annotation; required for --level gene.")
    parser.add_argument("-b", "--header", choices=["T", "F"], default="T",
                        help="Whether CoReID/CReSIL/Flec result files have a header (default: T).")
    parser.add_argument("-th", "--overlap", type=float, default=0.9,
                        help="Minimum per-segment overlap ratio for clustering (default: 0.9).")
    parser.add_argument("--max-start-distance", type=int, default=1_000_000,
                        help="Maximum first-segment start distance considered during clustering (default: 1000000).")
    parser.add_argument("-p", "--prefix", default="merged_eccDNA", help="Output filename prefix.")
    parser.add_argument("-o", "--out-dir", default=".", help="Output directory (default: current directory).")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Parser worker processes (default: 1).")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.threads < 1:
        raise SystemExit("ERROR: --threads must be at least 1")
    try:
        run(args)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc


if __name__ == "__main__":
    main()
