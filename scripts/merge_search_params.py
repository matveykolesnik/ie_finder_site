#!/usr/bin/env python3
"""Write the search parameters a run uses: the base file plus an optional override.

Keys in the override win; keys it leaves out keep their base values. Each file
is checked on its own first, so an error names the file it came from. The
output lists every parameter, so the report shows the full set.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from v3ps_filters import (
    FilterThresholds,
    read_params_file,
    thresholds_from_params,
    thresholds_to_params,
)


def merge_search_params(base: Path, override: Path | None) -> FilterThresholds:
    """Thresholds from ``base`` with the keys of ``override`` on top.

    Raises:
        OSError, ValueError, yaml.YAMLError: on a missing file, an unknown key
            or a bad value; the message names the file.
    """
    params = read_params_file(base)
    thresholds_from_params(params, str(base))
    if override is not None:
        extra = read_params_file(override)
        thresholds_from_params(extra, str(override))
        params.update(extra)
    return thresholds_from_params(params)


def main() -> None:
    """Parse arguments, merge, and write the parameters used."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", required=True, type=Path, help="search_params.yaml")
    parser.add_argument("--override", type=Path, default=None, help="File whose keys win")
    parser.add_argument("--out", required=True, type=Path, help="Merged parameters to write")
    args = parser.parse_args()
    try:
        thresholds = merge_search_params(args.base, args.override)
    except (OSError, ValueError, yaml.YAMLError) as exc:
        sys.exit(f"Search parameters: {exc}")
    args.out.write_text(yaml.safe_dump(thresholds_to_params(thresholds), sort_keys=False))


if __name__ == "__main__":
    main()
