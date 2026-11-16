"""Require the release tag to match the version frozen into installers."""
import sys
import tomllib
from pathlib import Path

version = tomllib.loads((Path(__file__).resolve().parents[1]/"pyproject.toml").read_text())["project"]["version"]
if sys.argv[1:] != ["v"+version]:
    raise SystemExit(f"Release tag must be v{version}; update pyproject.toml before tagging")
print("Release version:", version)
