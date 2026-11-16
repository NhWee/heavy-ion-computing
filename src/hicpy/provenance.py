"""Provenance helpers shared by the event-generation scripts."""
from __future__ import annotations

import pathlib
import subprocess

# What "the code that produced a file" means: outputs in data/ and results/
# change on every run (they record their own production time), so they must
# not make the tree count as dirty.
CODE_PATHS = ("src", "analysis", "scripts", "config", "generators",
              "environment/environment.yml", "environment/requirements-ci.txt")


def git_describe(repo: pathlib.Path) -> str:
    """Short commit hash, plus "-dirty" if code or config differ from it.

    A sample produced from uncommitted code must say so: "abc1234" alone would
    claim the file can be regenerated from that commit when it cannot.
    """
    try:
        head = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL, text=True).strip()
        status = subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain",
                                          "--", *CODE_PATHS],
                                         stderr=subprocess.DEVNULL, text=True).strip()
        return head + ("-dirty" if status else "")
    except Exception:                                   # noqa: BLE001
        return "unknown"


__all__ = ["CODE_PATHS", "git_describe"]
