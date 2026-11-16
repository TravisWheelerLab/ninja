#!/usr/bin/env python3
"""Write the body of a GitHub release: the CHANGELOG entry for the version,
then a table of the downloads and how to run them.

usage: release_notes.py <owner/repo> <version> <archive>...

Each archive is named ninja-<platform>.tar.gz or .zip. The notes go to
standard output; the table has one row per archive, in the order of
PLATFORMS. The CHANGELOG entry is the section headed "## <version>" in
CHANGELOG.md, up to the next "## " heading; without one the notes carry
only the table.
"""

import os
import re
import sys

# platform in the archive name -> (OS, architecture), in table order
PLATFORMS = {
    "macos-universal": ("macOS", "Universal (Intel and Apple Silicon)"),
    "x86_64-apple-darwin": ("macOS", "Intel 64 bit"),
    "aarch64-apple-darwin": ("macOS", "ARM 64 bit (Apple Silicon)"),
    "x86_64-unknown-linux-gnu": ("Linux", "Intel/AMD 64 bit"),
    "aarch64-unknown-linux-gnu": ("Linux", "ARM 64 bit"),
    "x86_64-pc-windows-msvc": ("Windows", "Intel/AMD 64 bit"),
}


def platform_of(archive: str) -> str:
    """'ninja-x86_64-apple-darwin.tar.gz' -> 'x86_64-apple-darwin'"""
    for ext in (".tar.gz", ".zip"):
        if archive.endswith(ext):
            return archive[: -len(ext)].split("-", 1)[1]
    sys.exit(f"{archive}: not a .tar.gz or .zip archive")


def changelog_entry(version: str) -> str:
    """The body of the CHANGELOG section for `version`, or an empty string."""
    try:
        with open("CHANGELOG.md", encoding="utf-8") as f:
            text = f.read()
    except OSError:
        return ""
    m = re.search(rf"^## {re.escape(version)}\b[^\n]*\n(.*?)(?=^## |\Z)", text, re.M | re.S)
    return m.group(1).strip() if m else ""


def main() -> None:
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    repo, version = sys.argv[1], sys.argv[2]
    tag = f"v{version}"
    archives = {platform_of(a): a for a in map(os.path.basename, sys.argv[3:])}

    unknown = sorted(set(archives) - set(PLATFORMS))
    missing = [p for p in PLATFORMS if p not in archives]
    if unknown or missing:
        sys.exit(f"unexpected platforms: {unknown}; no archive for: {missing}")

    entry = changelog_entry(version)
    if entry:
        print(entry)
        print()

    print("### Downloads")
    print("| OS | Architecture | Link |")
    print("|---------|----------|-------------|")
    for platform, (system, arch) in PLATFORMS.items():
        url = f"https://github.com/{repo}/releases/download/{tag}/{archives[platform]}"
        print(f"| {system} | {arch} | [Download]({url}) |")

    print(
        """
Each download unpacks to a folder that holds the program, `ninja`, and its license. Copy `ninja` to a folder on your `PATH`. The Linux program needs glibc 2.17 or newer.

On macOS, you may see the following message when you try to run `ninja`:

```
"ninja" cannot be opened because the developer cannot be verified.
```

To fix this, run the following command:

```
xattr -dr com.apple.quarantine <path to the unpacked folder>
```"""
    )


if __name__ == "__main__":
    main()
