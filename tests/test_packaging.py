import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_debian_archive_contains_valid_control_data_and_launcher(tmp_path):
    spec = importlib.util.spec_from_file_location("build_release", ROOT/"scripts/build_release.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    bundle = tmp_path/"bundle"
    bundle.mkdir()
    executable = bundle/"poreckit-gui"
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(0o755)
    out = tmp_path/"test.deb"
    module.debian_package(bundle, out, "0.1.0", "amd64")
    members = {}
    with out.open("rb") as stream:
        assert stream.read(8) == b"!<arch>\n"
        while header := stream.read(60):
            assert len(header) == 60 and header[-2:] == b"`\n"
            name = header[:16].decode().strip().rstrip("/")
            size = int(header[48:58])
            members[name] = stream.read(size)
            if size % 2:
                assert stream.read(1) == b"\n"
    assert members["debian-binary"] == b"2.0\n"
    with tarfile.open(fileobj=io.BytesIO(members["control.tar.gz"])) as tar:
        control = tar.extractfile("control").read().decode()
        assert "Architecture: amd64" in control
        assert "Depends: libc6" in control
    with tarfile.open(fileobj=io.BytesIO(members["data.tar.gz"])) as tar:
        assert tar.getmember("usr/bin/poreckit-gui").linkname == "/opt/poreckit-gui/poreckit-gui"
        assert tar.extractfile("opt/poreckit-gui/poreckit-gui").read() == executable.read_bytes()
        assert tar.getmember("usr/share/applications/poreckit-gui.desktop").uid == 0


def test_ci_builds_all_native_installers_before_publishing():
    workflow = yaml.safe_load((ROOT/".github/workflows/release.yml").read_text())
    job = workflow["jobs"]["build"]
    targets = {v["target"] for v in job["strategy"]["matrix"]["include"]}
    assert targets == {"windows-x86_64", "linux-deb-x86_64", "macos-arm64"}
    runs = "\n".join(step.get("run", "") for step in job["steps"])
    assert "uv sync --locked --group build" in runs
    assert "scripts/build_release.py" in runs
    assert "pytest -q" in runs
    assert "apt-get install -y ./dist/releases/*.deb" in runs
    publisher = workflow["jobs"]["publish"]
    assert set(publisher["needs"]) == {"build", "build-arch"}
    assert publisher["permissions"] == {"contents": "write"}
    publication = publisher["steps"][-1]["run"]
    assert publication.index("sha256sum -c") < publication.index("gh release upload")
    assert publication.index("gh release upload") < publication.index("--draft=false")
    assert len((ROOT/".github/poreckit-ref").read_text().strip()) == 40
    arch = workflow["jobs"]["build-arch"]
    assert arch["container"] == "archlinux:base-devel"
    commands = "\n".join(s.get("run", "") for s in arch["steps"])
    assert "runuser -u builder" in commands
    assert "scripts/build_release.py --format arch" in commands
    assert "pacman -U --noconfirm" in commands
    assert "AppImage" not in (ROOT/"scripts/build_release.py").read_text()


def test_linux_format_follows_the_build_distribution(monkeypatch):
    spec = importlib.util.spec_from_file_location("build_release", ROOT/"scripts/build_release.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for release, expected in [({"ID": "arch"}, "arch"), ({"ID": "endeavouros", "ID_LIKE": "arch"}, "arch"),
                              ({"ID": "ubuntu", "ID_LIKE": "debian"}, "deb")]:
        monkeypatch.setattr(module.platform, "freedesktop_os_release", lambda release=release: release)
        assert module.native_linux_format() == expected


def test_selftest_command_reports_real_gui_worker_results(tmp_path):
    report = tmp_path/"report.json"
    result = subprocess.run([sys.executable, "-m", "poreckit_gui", "--self-test", str(report)],
                            env={**os.environ, "QT_QPA_PLATFORM": "offscreen"}, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, (report.read_text() if report.exists() else result.stdout+result.stderr)
    data = json.loads(report.read_text())
    assert data["ok"]
    assert {
        "native reader and GUI rendering",
        "persistent worker reuse and signal track",
        "whole-genome sparse aggregation",
        "standalone signal, reference, reads and variant queries",
        "assembly inversion and lifted contact rendering",
        "dynamic plot catalog and workbench subprocess",
    }.issubset(data["checks"])
    assert data["catalog_entries"] > 100


def test_version_tag_validation_rejects_mismatched_release():
    result = subprocess.run([sys.executable, ROOT/"scripts/check_release_tag.py", "v999.0.0"], capture_output=True, text=True)
    assert result.returncode != 0 and "Release tag must be" in result.stderr
