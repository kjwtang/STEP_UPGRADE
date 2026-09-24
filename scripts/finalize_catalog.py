#!/usr/bin/env python3
"""Combine chronological STEP_UPGRADE chunk/month tables into one catalog."""

import argparse
import csv
import json
import os
from pathlib import Path
import shutil


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "parts", nargs="+", type=Path,
        help="Chronological output directories from run_npy_validation.py",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    parts = [part.resolve() for part in args.parts]
    for part in parts:
        for name in ("objects.csv", "edges.csv", "events.csv", "family_map.csv"):
            if not (part / name).is_file():
                raise FileNotFoundError(part / name)

    family_rows = read_rows(parts[-1] / "family_map.csv")
    family_map = {
        int(row["branch_id"]): int(row["family_id"]) for row in family_rows
    }
    objects, edges, events = [], [], []
    for part in parts:
        objects.extend(read_rows(part / "objects.csv"))
        edges.extend(read_rows(part / "edges.csv"))
        events.extend(read_rows(part / "events.csv"))

    node_ids = [int(row["node_id"]) for row in objects]
    if len(node_ids) != len(set(node_ids)):
        raise ValueError("duplicate node_id across input parts")
    known_nodes = set(node_ids)
    for row in edges:
        if int(row["parent_node_id"]) not in known_nodes or int(row["child_node_id"]) not in known_nodes:
            raise ValueError("edge references a node absent from the supplied parts")
    for row in objects:
        branch = int(row["branch_id"])
        if branch not in family_map:
            raise ValueError(f"branch {branch} is absent from final family map")
        row["family_id"] = str(family_map[branch])

    objects.sort(key=lambda row: (int(row["time"]), int(row["node_id"])))
    edges.sort(key=lambda row: (
        int(row["time"]), int(row["parent_node_id"]), int(row["child_node_id"])
    ))
    events.sort(key=lambda row: (int(row["time"]), int(row["event_id"])))

    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    if temporary.exists():
        raise FileExistsError(temporary)
    temporary.mkdir()
    try:
        write_rows(temporary / "objects.csv", objects, list(objects[0]) if objects else ["time", "node_id", "branch_id", "family_id"])
        write_rows(temporary / "edges.csv", edges, list(edges[0]) if edges else ["time", "parent_node_id", "child_node_id", "event"])
        write_rows(temporary / "events.csv", events, list(events[0]) if events else ["event_id", "time", "event", "parent_node_ids", "child_node_ids"])
        write_rows(
            temporary / "family_map.csv",
            [{"branch_id": branch, "family_id": family} for branch, family in sorted(family_map.items())],
            ["branch_id", "family_id"],
        )
        with (temporary / "catalog_manifest.json").open("w") as handle:
            json.dump({
                "schema_version": 1,
                "parts": [str(part) for part in parts],
                "nodes": len(objects),
                "edges": len(edges),
                "events": len(events),
                "families": len(set(family_map.values())),
            }, handle, indent=2, sort_keys=True)
        os.replace(temporary, output)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        raise
    print(f"Final catalog: {output}")


if __name__ == "__main__":
    main()
