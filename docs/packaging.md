# 三端安装包、CI 与发布

PyInstaller 将 Python、Qt、poreckit Rust 扩展和绘图依赖冻结到应用目录，再由平台工具生成安装包。最终用户不需要安装 Python、uv 或 Rust。

| 平台 | CI runner | 安装包 |
| --- | --- | --- |
| Windows x86_64 | windows-2022 | `*-windows-x86_64-setup.exe`，附便携 ZIP |
| macOS Intel | macos-15-intel | `*-macos-x86_64.dmg` |
| macOS Apple Silicon | macos-15 | `*-macos-arm64.dmg` |
| Debian / Ubuntu x86_64 | ubuntu-22.04 | `*-linux-deb-x86_64.deb` |
| Arch Linux x86_64 | Arch 原生容器 | `poreckit-gui-*-1-x86_64.pkg.tar.zst` |

构建入口为 [build_release.py](../scripts/build_release.py)，冻结配置为 [poreckit.spec](../packaging/poreckit.spec)，CI 为 [release.yml](../.github/workflows/release.yml)。这些配置不进行跨系统交叉编译，而是在目标系统上构建。[PyInstaller 平台说明](https://pyinstaller.org/en/stable/usage.html#supporting-multiple-operating-systems)

## 本地构建

源代码保持以下目录结构：

```text
work/
  poreckit/       # Rust + Python 后端
  poreckit_gui/   # 当前项目
```

安装 uv、Rust stable。Windows 需要 Visual Studio C++ Build Tools；macOS 需要 Xcode Command Line Tools。构建固定使用 Python 3.12。

```sh
uv sync --locked --group build
uv run --no-sync python scripts/build_release.py
```

Windows 还需要 Inno Setup 6 或更新版本，编译器 `ISCC.exe` 放在 PATH 或默认安装目录。CI 会自动安装。安装包按用户安装到 `%LOCALAPPDATA%\Programs\poreckit-gui`，提供开始菜单、可选桌面快捷方式和卸载程序，无需管理员权限。

Linux 只生成 DEB 和 Arch 原生包，分别在对应发行版中冻结和打包。脚本在 Arch 及 Arch 衍生发行版默认选择 `arch`，在 Debian / Ubuntu 默认选择 `deb`。Debian 包由标准 ar/tar 格式生成，本地打包不要求 dpkg，也不会修改系统安装目录。

Arch 使用 [PKGBUILD](../packaging/arch/PKGBUILD) 和系统 `makepkg`，生成 pacman 可读取的 `.PKGINFO`、`.BUILDINFO` 和 `.MTREE`。安装 `base-devel`、zstd 及 PKGBUILD 中列出的 Qt 系统依赖，以普通用户执行构建；不要用 root 运行 makepkg。打包关闭 strip/debug/LTO，防止破坏 PyInstaller 内嵌归档。源码归档和桌面资源均检查 SHA-256，且运行时不依赖系统 Python。

macOS 使用 `hdiutil` 生成 DMG，包含 `.app` 和 Applications 快捷方式。Intel 与 Apple Silicon 分别构建，不混装两种 CPU 的 Python/Qt 二进制。

输出目录为 `dist/releases/`，包括安装包、每个安装包的 `.sha256`、冻结应用自检报告和构建信息 JSON。`dist/poreckit-gui/` 是完整冻结目录；macOS 还会生成 `dist/poreckit-gui.app/`。不要只复制可执行文件而遗漏依赖目录。

开发时可以只重打安装包，或者指定本平台格式：

```sh
uv run --no-sync python scripts/build_release.py --skip-freeze
uv run --no-sync python scripts/build_release.py --format deb
uv run --no-sync python scripts/build_release.py --format arch
uv run --no-sync python scripts/build_release.py --format dmg
uv run --no-sync python scripts/build_release.py --format exe zip
```

`--skip-freeze` 复用上一次冻结结果；修改源码后应运行不带此选项的完整命令。`--format` 仅接受当前系统可生成的格式。DEB 应在 Debian / Ubuntu 构建，Arch 包应在 Arch 构建，避免将另一个发行版冻结的依赖误当作原生构建。

## 构建验证

每次打包都会从临时工作目录运行冻结应用，移除 `PYTHONPATH` 和 `PYTHONHOME` 的影响，确认它独立于源码和开发环境。自检覆盖：

1. 原生数据读取与 GUI 热图。
2. 持久绘图进程复用和信号轨道。
3. 全基因组稀疏汇总。
4. 动态绘图目录及 Plot Workbench 子进程。

CI 对 Windows 使用 `--verify-installed`，静默安装到 `build/installed-windows/` 并验证已安装程序；默认本地构建不会安装程序或修改安装注册表。该选项会改变用户安装注册信息，因此仅在一次性测试环境中使用。DEB 和 Arch 包都会提取到临时目录并运行包内程序；Arch 还使用 `pacman -Qip` 检查包信息和入口。CI 会分别用 apt 和 pacman 实际安装，再从已安装位置自检。DMG 会执行 `hdiutil verify`、只读挂载并从盘内运行 `.app` 自检。自检使用独立设置文件，不覆盖日常使用的窗口配置。所有检查通过后才进入发布阶段。

手动复查冻结应用：

```sh
# Linux
QT_QPA_PLATFORM=offscreen dist/poreckit-gui/poreckit-gui --self-test report.json

# macOS
QT_QPA_PLATFORM=offscreen dist/poreckit-gui.app/Contents/MacOS/poreckit-gui --self-test report.json
```

```powershell
# Windows PowerShell
$env:QT_QPA_PLATFORM = 'offscreen'
& .\dist\poreckit-gui\poreckit-gui.exe --self-test report.json
```

检查 `report.json` 中的 `ok` 与 `frozen` 均为 `true`。冻结入口恢复 Windows GUI 子进程的管道标准流；Matplotlib 明确包含 Agg、SVG 和 PDF 后端，避免源码环境可运行、冻结后无法绘图或导出的情况。

## CI 与版本固定

PR、主分支 push、手动 workflow_dispatch 均运行五个目标的源码测试、冻结自检和原生安装包构建，然后上传 Actions artifacts。Arch 作业在 `archlinux:base-devel` 容器内安装工具链，由非 root 的 builder 用户构建，安装测试才使用容器 root。DEB 与 Arch 的报告和构建信息使用不同文件名，合并发布资产时不会覆盖。固定后端提交在 [.github/poreckit-ref](../.github/poreckit-ref)，避免后端 main 分支变化影响相同 GUI 版本的构建。

更新后端时，先更新 `../poreckit`，将通过 GUI 测试的完整 40 位提交写入 `.github/poreckit-ref`，更新 `uv.lock` 后一并提交。`pyproject.toml`、`uv.lock`、PKGBUILD 和后端提交应一起审查。

如果后端仓库是私有仓库，在 GUI 仓库的 Actions secrets 添加 `PORECKIT_READ_TOKEN`：它只需要对 `imjiaoyuan/poreckit` 的 Contents 只读权限。公开仓库可以使用默认 `GITHUB_TOKEN`。来自外部 fork 的 PR 无法取得私有仓库 secret，应由维护者在有权限的可信分支验证；不要把私有仓库 token 暴露给 fork PR。

DEB 发布 runner 固定 Ubuntu 22.04，Arch 包使用当前 Arch 环境。PyInstaller 不打包 glibc，因此在 Arch 等较新的系统上构建的本地安装包不能据此宣称兼容旧发行版。两种包均按实际构建主机记录 glibc 版本依赖；构建信息 JSON 同时记录主机、架构和依赖版本。[Linux 二进制兼容说明](https://pyinstaller.org/en/stable/usage.html#making-gnu-linux-apps-forward-compatible)

## 发布到 GitHub Releases

先更新 `pyproject.toml` 中版本，再更新并提交锁文件。只有与项目版本精确一致的 `v` 标签可以发布：

```sh
uv lock
git add pyproject.toml uv.lock .github packaging scripts src docs README.md
git commit -m "Prepare desktop release"
git tag v0.1.0
git push origin main
git push origin v0.1.0
```

示例标签必须改成当前项目版本。分支构建和手动构建只上传 artifacts，不自动创建 Release。标签构建要求五个目标全部成功，核对安装包校验和，再创建草稿 Release、上传全部资产，最后发布。带 `-` 的预发行标签会标记 prerelease。已发布的版本不能被工作流覆盖；失败的草稿可以重跑，正式更新需新版本号。

## 可选代码签名

默认 Windows 安装包未做 Authenticode 签名；macOS 使用 PyInstaller 的 ad-hoc 签名，不等同于 Developer ID 签名或 Apple 公证。构建信息会记录签名状态。

配置以下 Actions secrets 后，可信分支和标签构建可自动签名；PR 构建不使用这些证书：

| Windows secret | 内容 |
| --- | --- |
| `WINDOWS_CERTIFICATE_PFX_BASE64` | PFX 证书的 Base64 |
| `WINDOWS_CERTIFICATE_PASSWORD` | PFX 密码 |

CI 将证书导入临时 runner 的用户证书库，签名应用 EXE 与安装 EXE，并用 SHA-256 和时间戳验证。私钥密码不会作为打包脚本参数打印。

| macOS secret | 内容 |
| --- | --- |
| `MACOS_CERTIFICATE_P12_BASE64` | Developer ID Application P12 的 Base64 |
| `MACOS_CERTIFICATE_PASSWORD` | P12 密码 |
| `MACOS_CODESIGN_IDENTITY` | Developer ID Application identity |
| `APPLE_ID` | Apple 开发者账号 |
| `APPLE_APP_SPECIFIC_PASSWORD` | App-specific password |
| `APPLE_TEAM_ID` | Apple Team ID |

CI 创建临时钥匙串，签名 `.app` 与 DMG，提交 Apple 公证并 staple DMG，之后才计算校验和。macOS 证书一旦配置，其他相关 secret 也必须齐备，缺失时构建会失败，不会把未公证包标为已公证。临时证书和钥匙串会在结束时清理。

本地签名可设置 `WINDOWS_SIGN_CERT_SHA1`（已导入用户证书库的证书指纹）、`MACOS_CODESIGN_IDENTITY`；公证使用本机钥匙串中预先保存的 `MACOS_NOTARY_PROFILE`，自定义钥匙串通过 `MACOS_NOTARY_KEYCHAIN` 指定。

## 用户安装

- Windows：运行 `*-setup.exe`，卸载使用系统“已安装的应用”。ZIP 是独立的便携版本。
- macOS：打开与 CPU 匹配的 DMG，将应用拖入 Applications。
- Debian / Ubuntu：`sudo apt install ./poreckit-gui-*.deb`；卸载 `sudo apt remove poreckit-gui`。
- Arch Linux：`sudo pacman -U ./poreckit-gui-*.pkg.tar.zst`；卸载 `sudo pacman -R poreckit-gui`。

开发者应在不同目标系统上检查实际运行结果；配置了 CI 不代表 Windows/macOS 已经通过。冻结自检不替代全部真实数据、显示驱动或签名分发验证。
