"""Freeze, test and wrap the native application into platform installers."""
import argparse
import gzip
import hashlib
import importlib.metadata
import io
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import shlex
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = tomllib.loads((ROOT/"pyproject.toml").read_text())["project"]["version"]


def run(args, **kwargs):
    print("Running:", " ".join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), check=True, **kwargs)


def run_with_retry(args, attempts=4, delay=5, **kwargs):
    """hdiutil on CI runners intermittently fails with 'Resource busy' while
    DiskImages2 or Spotlight still holds a stale mount; retry until it clears."""
    import time
    for attempt in range(attempts):
        try:
            return run(args, **kwargs)
        except subprocess.CalledProcessError:
            if attempt == attempts - 1:
                raise
            print(f"hdiutil attempt {attempt + 1} failed; retrying in {delay}s", flush=True)
            time.sleep(delay)


def digest(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def ar_member(stream, name, data):
    size = len(data) if isinstance(data, bytes) else data.stat().st_size
    header = f"{name+'/':<16}{0:<12}{0:<6}{0:<6}{'100644':<8}{size:<10}`\n".encode("ascii")
    assert len(header) == 60
    stream.write(header)
    if isinstance(data, bytes):
        stream.write(data)
    else:
        with data.open("rb") as source:
            shutil.copyfileobj(source, stream)
    if size % 2:
        stream.write(b"\n")


def debian_package(bundle, output, version, arch):
    """Create the standard ar/control.tar.gz/data.tar.gz Debian layout."""
    epoch = int(os.environ.get("SOURCE_DATE_EPOCH", "0"))
    def root_owner(info):
        info.uid = info.gid = 0
        info.uname = info.gname = "root"
        info.mtime = epoch
        return info
    glibc = platform.libc_ver()[1] or "2.35"
    installed = sum(p.stat().st_size for p in bundle.rglob("*") if p.is_file() and not p.is_symlink())//1024
    control = (f"Package: poreckit-gui\nVersion: {version}\nSection: science\nPriority: optional\nArchitecture: {arch}\n"
               f"Maintainer: poreckit developers\nInstalled-Size: {installed}\n"
               f"Depends: libc6 (>= {glibc}), libglib2.0-0, libgl1, libegl1, libx11-6, libxcb1, libxkbcommon0, "
               "libxkbcommon-x11-0, libxcb-cursor0, libxcb-icccm4, libxcb-image0, libxcb-keysyms1, "
               "libxcb-render-util0, libxcb-xinerama0, libxcb-xkb1, libdbus-1-3, libfontconfig1, libfreetype6\n"
               "Homepage: https://github.com/imjiaoyuan/poreckit\nDescription: Native Hi-C and Pore-C contact map browser\n"
               " Browse contact maps, genomic tracks, annotations and analysis figures.\n").encode()
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=epoch) as compressed:
        # dpkg rejects PAX extended headers ('x' type), so force the GNU format.
        with tarfile.open(fileobj=compressed, mode="w", format=tarfile.GNU_FORMAT) as tar:
            info = tarfile.TarInfo("control")
            info.size, info.mode = len(control), 0o644
            tar.addfile(root_owner(info), io.BytesIO(control))
    with tempfile.TemporaryDirectory(prefix="poreckit-deb-") as folder:
        data = Path(folder)/"data.tar.gz"
        with data.open("wb") as stream, gzip.GzipFile(fileobj=stream, mode="wb", mtime=epoch) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.GNU_FORMAT) as tar:
                tar.add(bundle, arcname="opt/poreckit-gui", filter=root_owner)
                tar.add(ROOT/"packaging/poreckit-gui.desktop", arcname="usr/share/applications/poreckit-gui.desktop", filter=root_owner)
                tar.add(ROOT/"packaging/poreckit.svg", arcname="usr/share/icons/hicolor/scalable/apps/poreckit-gui.svg", filter=root_owner)
                info = tarfile.TarInfo("usr/bin/poreckit-gui")
                info.type, info.linkname, info.mode = tarfile.SYMTYPE, "/opt/poreckit-gui/poreckit-gui", 0o755
                tar.addfile(root_owner(info))
        with output.open("wb") as stream:
            stream.write(b"!<arch>\n")
            ar_member(stream, "debian-binary", b"2.0\n")
            ar_member(stream, "control.tar.gz", buf.getvalue())
            ar_member(stream, "data.tar.gz", data)


def native_linux_format():
    try:
        release = platform.freedesktop_os_release()
    except OSError:
        return "deb"
    families = (release.get("ID", "")+" "+release.get("ID_LIKE", "")).split()
    return "arch" if "arch" in families else "deb"


def arch_package(bundle, output_dir, version, arch):
    """Build a native pacman package with makepkg, including its ALPM metadata."""
    import re
    if not re.fullmatch(r"[A-Za-z0-9_.+]+", version) or arch not in {"x86_64", "aarch64"}:
        raise ValueError("Invalid Arch package version or architecture")
    for tool in ("makepkg", "fakeroot", "bsdtar", "zstd", "pacman"):
        if not shutil.which(tool):
            raise RuntimeError(f"Arch packaging requires {tool}; build it on Arch Linux with base-devel installed")
    if os.geteuid() == 0:
        raise RuntimeError("makepkg must run as an unprivileged user")
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="poreckit-arch-") as folder:
        stage = Path(folder)
        source = stage/"poreckit-gui.tar"
        with tarfile.open(source, "w") as tar:
            tar.add(bundle, arcname="poreckit-gui")
        for name in ("poreckit-gui.desktop", "poreckit.svg"):
            shutil.copyfile(ROOT/"packaging"/name, stage/name)
        recipe = (ROOT/"packaging/arch/PKGBUILD").read_text()
        replacements = dict(VERSION=version, ARCH=arch, GLIBC=platform.libc_ver()[1],
                            BUNDLE_SHA256=digest(source), DESKTOP_SHA256=digest(stage/"poreckit-gui.desktop"),
                            ICON_SHA256=digest(stage/"poreckit.svg"))
        for key, value in replacements.items():
            recipe = recipe.replace("@"+key+"@", value)
        (stage/"PKGBUILD").write_text(recipe)
        config = stage/"makepkg.conf"
        paths = {"BUILDDIR": stage/"build", "SRCDEST": stage, "SRCPKGDEST": stage/"sources", "LOGDEST": stage/"logs", "PKGDEST": output_dir}
        for path in paths.values():
            path.mkdir(parents=True, exist_ok=True)
        config.write_text("source /etc/makepkg.conf\n" + "".join(f"{k}={shlex.quote(str(v))}\n" for k, v in paths.items()) +
                          "PKGEXT='.pkg.tar.zst'\nBUILDENV=(!distcc !color !ccache !check !sign)\n"
                          "COMPRESSZST=(zstd -c -T2 -3 -)\nPACKAGER='poreckit builders'\n")
        run(["makepkg", "--config", config, "--force", "--noconfirm", "--nosign", "--clean"], cwd=stage)
    out = output_dir/f"poreckit-gui-{version}-1-{arch}.pkg.tar.zst"
    if not out.is_file():
        raise RuntimeError("makepkg did not create the expected Arch package")
    return out


def selftest(executable, report):
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    with tempfile.TemporaryDirectory(prefix="poreckit-smoke-") as folder:
        env["MPLCONFIGDIR"] = str(Path(folder)/"matplotlib")
        run([executable, "--self-test", report], cwd=folder, env=env, timeout=240)
    data = json.loads(report.read_text())
    if not data.get("ok") or not data.get("frozen"):
        raise RuntimeError(f"Frozen self-test failed: {data}")


def verify_deb(package, report):
    with tempfile.TemporaryDirectory(prefix="poreckit-deb-test-") as folder, package.open("rb") as stream:
        if stream.read(8) != b"!<arch>\n":
            raise RuntimeError("Invalid Debian archive")
        found = False
        while header := stream.read(60):
            if len(header) != 60 or header[-2:] != b"`\n":
                raise RuntimeError("Invalid Debian member header")
            name = header[:16].decode().strip().rstrip("/")
            size = int(header[48:58])
            if name == "data.tar.gz":
                with tarfile.open(fileobj=stream, mode="r|gz") as tar:
                    selected = (v for v in tar if v.name.startswith("opt/poreckit-gui/"))
                    tar.extractall(folder, members=selected, filter="data")
                found = True
                break
            stream.seek(size+size%2, 1)
        if not found:
            raise RuntimeError("Missing Debian data payload")
        selftest(Path(folder)/"opt/poreckit-gui/poreckit-gui", report)


def verify_arch(package, report):
    run(["pacman", "-Qip", package])
    with tempfile.TemporaryDirectory(prefix="poreckit-arch-test-") as folder:
        run(["bsdtar", "-xf", package, "-C", folder])
        selftest(Path(folder)/"usr/bin/poreckit-gui", report)


def sign_windows(path):
    thumbprint = os.environ.get("WINDOWS_SIGN_CERT_SHA1")
    if not thumbprint:
        return
    candidates = sorted(Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")).glob("Windows Kits/10/bin/*/x64/signtool.exe"))
    tool = os.environ.get("SIGNTOOL_PATH") or shutil.which("signtool") or (str(candidates[-1]) if candidates else "signtool.exe")
    run([tool, "sign", "/sha1", thumbprint, "/s", "My", "/fd", "SHA256", "/tr", "https://timestamp.digicert.com", "/td", "SHA256", path])
    run([tool, "verify", "/pa", path])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-freeze", action="store_true", help="Package an existing dist tree")
    parser.add_argument("--verify-installed", action="store_true", help="CI: run the Windows installer; changes per-user installation registration")
    parser.add_argument("--format", nargs="+", choices=["deb", "arch", "dmg", "exe", "zip"], help="Defaults to native installers")
    parser.add_argument("--output", type=Path, default=ROOT/"dist/releases")
    args = parser.parse_args()
    system = {"Darwin": "macos", "Windows": "windows", "Linux": "linux"}[platform.system()]
    arch = {"amd64": "x86_64", "x86_64": "x86_64", "arm64": "arm64", "aarch64": "arm64"}[platform.machine().lower()]
    forms = args.format or {"linux": [native_linux_format()], "windows": ["exe", "zip"], "macos": ["dmg"]}[system]
    allowed = {"linux": {"deb", "arch"}, "windows": {"exe", "zip"}, "macos": {"dmg", "zip"}}[system]
    if not set(forms) <= allowed:
        parser.error("Installers must be built on their target operating system")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PORECKIT_VERSION": VERSION, "MPLBACKEND": "Agg"}
    run([sys.executable, ROOT/"packaging/create_icons.py"], env=env)
    if not args.skip_freeze:
        run([sys.executable, "-m", "PyInstaller", "--clean", "--noconfirm", "--workpath", ROOT/"build/pyinstaller",
             "--distpath", ROOT/"dist", ROOT/"packaging/poreckit.spec"], cwd=ROOT, env=env)
    bundle = ROOT/"dist/poreckit-gui"
    executable = bundle/("poreckit-gui.exe" if system == "windows" else "poreckit-gui")
    if system == "macos":
        bundle = ROOT/"dist/poreckit-gui.app"
        executable = bundle/"Contents/MacOS/poreckit-gui"
    target = system+"-"+"-".join(forms) if system == "linux" else system
    stem = f"poreckit-gui-{VERSION}-{target}-{arch}"
    if system == "windows":
        if arch != "x86_64":
            parser.error("The Windows installer currently targets x86_64")
        sign_windows(executable)
    report = args.output/(stem+"-selftest.json")
    selftest(executable.resolve(), report)
    assets = []
    for form in forms:
        if form == "deb":
            out = args.output/(stem+".deb")
            debian_package(bundle, out, VERSION, "amd64" if arch == "x86_64" else "arm64")
            verify_deb(out, args.output/(stem+"-deb-extracted-selftest.json"))
        elif form == "arch":
            out = arch_package(bundle, args.output, VERSION, "x86_64" if arch == "x86_64" else "aarch64")
            verify_arch(out, args.output/(stem+"-arch-extracted-selftest.json"))
        elif form == "exe":
            candidates = [p for key in ("ProgramFiles(x86)", "ProgramFiles") for p in
                          Path(os.environ.get(key, r"C:\Program Files (x86)")).glob("Inno Setup*/ISCC.exe")]
            compiler = shutil.which("ISCC") or (str(sorted(candidates)[-1]) if candidates else "ISCC.exe")
            run([compiler, f"/DAppVersion={VERSION}", f"/DSourceDir={bundle}", f"/DOutputDir={args.output}", ROOT/"packaging/windows.iss"])
            out = args.output/(stem+"-setup.exe")
            sign_windows(out)
            if args.verify_installed:
                installed = ROOT/"build/installed-windows"
                run([out, "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART", "/SP-", f"/DIR={installed}"])
                selftest(installed/"poreckit-gui.exe", args.output/(stem+"-installed-selftest.json"))
        elif form == "dmg":
            out = args.output/(stem+".dmg")
            with tempfile.TemporaryDirectory(prefix="poreckit-dmg-") as folder:
                shutil.copytree(bundle, Path(folder)/bundle.name, symlinks=True)
                (Path(folder)/"Applications").symlink_to("/Applications")
                run_with_retry(["hdiutil", "create", "-volname", "poreckit GUI", "-srcfolder", folder, "-ov", "-format", "UDZO", out])
            run(["hdiutil", "verify", out])
            with tempfile.TemporaryDirectory(prefix="poreckit-dmg-test-") as folder:
                run(["hdiutil", "attach", "-readonly", "-nobrowse", "-mountpoint", folder, out])
                try:
                    selftest(Path(folder)/bundle.name/"Contents/MacOS/poreckit-gui", args.output/(stem+"-mounted-selftest.json"))
                finally:
                    run(["hdiutil", "detach", folder])
            identity = os.environ.get("MACOS_CODESIGN_IDENTITY")
            if identity:
                run(["codesign", "--force", "--sign", identity, "--timestamp", out])
            profile = os.environ.get("MACOS_NOTARY_PROFILE")
            if profile:
                if not identity:
                    raise RuntimeError("Notarization requires MACOS_CODESIGN_IDENTITY")
                command = ["xcrun", "notarytool", "submit", out, "--keychain-profile", profile, "--wait", "--timeout", "30m"]
                if os.environ.get("MACOS_NOTARY_KEYCHAIN"):
                    command += ["--keychain", os.environ["MACOS_NOTARY_KEYCHAIN"]]
                run(command)
                run(["xcrun", "stapler", "staple", out])
        else:
            out = Path(shutil.make_archive(str(args.output/stem), "zip", root_dir=bundle.parent, base_dir=bundle.name))
        if not out.is_file():
            raise RuntimeError(f"Installer was not created: {out}")
        assets.append(out)
        out.with_name(out.name+".sha256").write_text(digest(out)+"  "+out.name+"\n")
    manifest = dict(version=VERSION, os=system, formats=forms, architecture=arch, python=sys.version, platform=platform.platform(),
                    glibc=platform.libc_ver(), gui_commit=os.environ.get("GITHUB_SHA", "local"),
                    core_ref=(ROOT/".github/poreckit-ref").read_text().strip(),
                    pyinstaller=importlib.metadata.version("pyinstaller"),
                    signing={"windows": bool(os.environ.get("WINDOWS_SIGN_CERT_SHA1")),
                             "macos_developer_id": bool(os.environ.get("MACOS_CODESIGN_IDENTITY")),
                             "macos_notarized": bool(os.environ.get("MACOS_NOTARY_PROFILE"))},
                    dependencies={name: importlib.metadata.version(name) for name in ("poreckit", "numpy", "matplotlib", "h5py", "PySide6-Essentials")},
                    assets=[dict(name=p.name, sha256=digest(p), bytes=p.stat().st_size) for p in assets])
    (args.output/(stem+"-build.json")).write_text(json.dumps(manifest, indent=2))
    print("Installers:", *assets, sep="\n")


if __name__ == "__main__":
    main()
