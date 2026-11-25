#!/usr/bin/env python3
"""
Select a ZIP of Sanger *unclipped.seq files, concatenate them, and extract barcodes.

Workflow:
  1. File dialog: pick the .zip archive.
  2. Parameter form: edit upstream/downstream flanks, length range, max errors
     (pre-filled with the defaults in BarcodeParams).
  3. The archive is extracted to a temporary folder, all *unclipped.seq files are
     concatenated into <zip_folder>/<zip_name>_combined_unclipped.seq, the temporary
     folder is removed, and extract_barcodes() runs on the combined file.

Falls back to terminal prompts if no GUI/display is available.

Usage:
    uv run main.py          # dialogs
    uv run main.py --cli    # terminal prompts
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path

from utils.barcode_extractor import extract_barcodes


# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class BarcodeParams:
    upstream_seq: str = "GGTCTCT"
    downstream_seq: str = "CGTACGC"
    min_len: int = 10
    max_len: int = 23
    max_errors: int = 2


# Label shown in the form / CLI prompt for each parameter (order = display order)
FIELD_LABELS = {
    "upstream_seq": "Upstream sequence",
    "downstream_seq": "Downstream sequence",
    "min_len": "Min barcode length",
    "max_len": "Max barcode length",
    "max_errors": "Max errors (mismatches)",
}

IUPAC_DNA = set("ACGTRYSWKMBDHVN")


def parse_params(raw: dict[str, str]) -> BarcodeParams:
    """Validate user-entered strings and build a BarcodeParams. Raises ValueError."""
    seqs = {}
    for key in ("upstream_seq", "downstream_seq"):
        seq = raw[key].strip().upper()
        if not seq:
            raise ValueError(f"{FIELD_LABELS[key]} cannot be empty.")
        bad = set(seq) - IUPAC_DNA
        if bad:
            raise ValueError(
                f"{FIELD_LABELS[key]} contains non-DNA characters: {''.join(sorted(bad))}"
            )
        seqs[key] = seq

    ints = {}
    for key in ("min_len", "max_len", "max_errors"):
        try:
            ints[key] = int(raw[key].strip())
        except ValueError:
            raise ValueError(f"{FIELD_LABELS[key]} must be a whole number.") from None
        if ints[key] < 0:
            raise ValueError(f"{FIELD_LABELS[key]} cannot be negative.")

    if ints["min_len"] < 1:
        raise ValueError("Min barcode length must be at least 1.")
    if ints["min_len"] > ints["max_len"]:
        raise ValueError("Min barcode length cannot exceed max barcode length.")

    return BarcodeParams(**seqs, **ints)


# --------------------------------------------------------------------------- #
# Core logic (GUI-independent)
# --------------------------------------------------------------------------- #
def safe_extract(
    zip_path: Path,
    dest_dir: Path,
    password: bytes | None = None,
    progress=print,
) -> list[Path]:
    """Extract zip_path into dest_dir, refusing entries that escape dest_dir."""
    dest_dir = dest_dir.resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []

    with zipfile.ZipFile(zip_path) as zf:
        members = zf.infolist()
        n = len(members)
        for i, member in enumerate(members, start=1):
            target = (dest_dir / member.filename).resolve()
            if not target.is_relative_to(dest_dir):  # zip-slip protection
                raise ValueError(f"Unsafe path in archive, aborting: {member.filename!r}")
            zf.extract(member, path=dest_dir, pwd=password)
            extracted.append(target)
            progress(f"[{i}/{n}] {member.filename}")

    return extracted


def needs_password(zip_path: Path) -> bool:
    """True if any entry has the 'encrypted' flag bit set."""
    with zipfile.ZipFile(zip_path) as zf:
        return any(info.flag_bits & 0x1 for info in zf.infolist())


def concat_unclipped(directory: Path, output: Path) -> Path:
    """Concatenate all *unclipped.seq files under `directory` (recursively) into `output`."""
    files = sorted(
        p
        for p in directory.rglob("*unclipped.seq")
        if p.is_file()
        and "__MACOSX" not in p.parts      # macOS Finder metadata folder
        and not p.name.startswith("._")    # macOS resource-fork stubs
        and p.resolve() != output.resolve()
    )
    if not files:
        raise FileNotFoundError("No *unclipped.seq files found in the archive.")

    with output.open("wb") as out:
        for f in files:
            with f.open("rb") as src:
                shutil.copyfileobj(src, out)
                if src.tell() > 0:  # non-empty: ensure it ends with a newline
                    src.seek(-1, 2)
                    if src.read(1) != b"\n":
                        out.write(b"\n")

    print(f"Concatenated {len(files)} files -> {output}")
    return output


def run_pipeline(
    zip_path: Path,
    params: BarcodeParams,
    password: bytes | None = None,
) -> Path:
    """Extract -> concatenate -> clean up -> extract barcodes. Returns the combined file."""
    combined = zip_path.parent / f"{zip_path.stem}_combined_unclipped.seq"

    # Extracted files are temporary: TemporaryDirectory deletes them even on error,
    # and never touches a folder that existed before.
    with tempfile.TemporaryDirectory(prefix=f"{zip_path.stem}_") as tmp:
        safe_extract(zip_path, Path(tmp), password=password)
        concat_unclipped(Path(tmp), output=combined)

    print(f"Running barcode extraction with {asdict(params)}")
    extract_barcodes(
        input_fasta=combined,
        output_prefix=zip_path.stem,
        **asdict(params),
    )
    return combined


# --------------------------------------------------------------------------- #
# GUI front-end
# --------------------------------------------------------------------------- #
def ask_params_gui(root, zip_path: Path, defaults: BarcodeParams) -> BarcodeParams | None:
    """Modal form pre-filled with defaults. Returns None if the user cancels."""
    import tkinter as tk
    from tkinter import messagebox, ttk

    win = tk.Toplevel(root)
    win.title("Barcode extraction parameters")
    win.resizable(False, False)
    win.attributes("-topmost", True)

    frame = ttk.Frame(win, padding=12)
    frame.grid(sticky="nsew")

    ttk.Label(frame, text=f"Archive: {zip_path.name}", font=("TkDefaultFont", 0, "bold")).grid(
        row=0, column=0, columnspan=2, sticky="w", pady=(0, 8)
    )

    vars_: dict[str, tk.StringVar] = {}
    for row, (name, label) in enumerate(FIELD_LABELS.items(), start=1):
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10), pady=3)
        var = tk.StringVar(value=str(getattr(defaults, name)))
        entry = ttk.Entry(frame, textvariable=var, width=24)
        entry.grid(row=row, column=1, sticky="ew", pady=3)
        if row == 1:
            entry.focus_set()
        vars_[name] = var

    result: dict[str, BarcodeParams | None] = {"params": None}

    def on_run(_event=None):
        try:
            result["params"] = parse_params({k: v.get() for k, v in vars_.items()})
        except ValueError as e:
            messagebox.showerror("Invalid parameter", str(e), parent=win)
            return
        win.destroy()

    def on_reset():
        for name, var in vars_.items():
            var.set(str(getattr(defaults, name)))

    def on_cancel(_event=None):
        win.destroy()

    buttons = ttk.Frame(frame)
    buttons.grid(row=len(FIELD_LABELS) + 1, column=0, columnspan=2, sticky="e", pady=(10, 0))
    ttk.Button(buttons, text="Reset defaults", command=on_reset).pack(side="left", padx=(0, 6))
    ttk.Button(buttons, text="Cancel", command=on_cancel).pack(side="left", padx=(0, 6))
    ttk.Button(buttons, text="Run", command=on_run).pack(side="left")

    win.bind("<Return>", on_run)
    win.bind("<Escape>", on_cancel)
    win.protocol("WM_DELETE_WINDOW", on_cancel)

    win.grab_set()
    win.focus_force()
    root.wait_window(win)
    return result["params"]


def run_gui(root) -> int:
    from tkinter import filedialog, messagebox, simpledialog

    root.withdraw()
    root.attributes("-topmost", True)

    zip_name = filedialog.askopenfilename(
        title="Select a ZIP file",
        filetypes=[("ZIP archives", "*.zip"), ("All files", "*.*")],
    )
    if not zip_name:
        print("No file selected. Exiting.")
        return 1
    zip_path = Path(zip_name)

    if not zipfile.is_zipfile(zip_path):
        messagebox.showerror("Not a ZIP", f"{zip_path.name} is not a valid ZIP archive.")
        return 1

    params = ask_params_gui(root, zip_path, BarcodeParams())
    if params is None:
        print("Cancelled.")
        return 1

    pwd: bytes | None = None
    if needs_password(zip_path):
        p = simpledialog.askstring("Password", "This archive is encrypted. Password:", show="*")
        if p is None:
            return 1
        pwd = p.encode()

    try:
        run_pipeline(zip_path, params, password=pwd)
    except (RuntimeError, zipfile.BadZipFile, ValueError, NotImplementedError,
            FileNotFoundError) as e:
        # RuntimeError: bad password; NotImplementedError: unsupported compression (AES)
        messagebox.showerror("Processing failed", str(e))
        return 1
    finally:
        root.destroy()

    return 0


# --------------------------------------------------------------------------- #
# Terminal fallback
# --------------------------------------------------------------------------- #
def ask_params_cli(defaults: BarcodeParams) -> BarcodeParams:
    """Prompt for each parameter; Enter keeps the default. Re-asks on invalid input."""
    while True:
        raw = {}
        for name, label in FIELD_LABELS.items():
            default = str(getattr(defaults, name))
            raw[name] = input(f"{label} [{default}]: ").strip() or default
        try:
            return parse_params(raw)
        except ValueError as e:
            print(f"Invalid parameter: {e}\nPlease try again.\n")


def run_cli() -> int:
    import getpass

    zip_path = Path(input("Path to ZIP file: ").strip().strip('"').strip("'")).expanduser()
    if not zipfile.is_zipfile(zip_path):
        print(f"Error: {zip_path} is not a valid ZIP archive.", file=sys.stderr)
        return 1

    params = ask_params_cli(BarcodeParams())
    pwd = getpass.getpass("Password: ").encode() if needs_password(zip_path) else None

    try:
        run_pipeline(zip_path, params, password=pwd)
    except (RuntimeError, zipfile.BadZipFile, ValueError, NotImplementedError,
            FileNotFoundError) as e:
        print(f"Processing failed: {e}", file=sys.stderr)
        return 1
    return 0


# --------------------------------------------------------------------------- #
def main() -> int:
    parser = argparse.ArgumentParser(description="Extract barcodes from a ZIP of .seq files.")
    parser.add_argument("--cli", action="store_true", help="use terminal prompts instead of dialogs")
    args = parser.parse_args()

    if args.cli:
        return run_cli()

    # Only GUI *start-up* failures trigger the fallback; errors in the pipeline itself
    # are reported normally instead of being mistaken for "no display".
    try:
        import tkinter as tk
        root = tk.Tk()
    except Exception as e:  # ImportError (no tkinter) or TclError (no display)
        print(f"GUI unavailable ({e}); falling back to terminal prompts.\n")
        return run_cli()

    return run_gui(root)


if __name__ == "__main__":
    sys.exit(main())