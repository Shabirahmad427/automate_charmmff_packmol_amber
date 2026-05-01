#!/usr/bin/env python3
#1. Manual box:
#python3 prepare_charmm_packmol_amber.py --pdb input.pdb --pqr input_ph55.pqr --outdir my_system --box-side 110 --neutralize-only
#2. Auto box from protein size:
#python3 prepare_charmm_packmol_amber.py --pdb input.pdb --pqr input_ph55.pqr --outdir my_system--neutralize-only --padding 15

from __future__ import annotations
import argparse
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

CHARMM_TOP = "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf"
CHARMM_PAR = "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmpar1.5/par_all36_prot.prm"
CHARMM_WATER_IONS = "/usr/local/lib/vmd/plugins/noarch/tcl/trunctraj1.5/toppar/stream/toppar_water_ions.str"
DEFAULT_PACKMOL = shutil.which("packmol") or "packmol"
DEFAULT_PSFGEN = "/home/shabir/Downloads/NAMD3.1/psfgen"
DEFAULT_PARMED = "/home/shabir/Downloads/ambertools25/bin/parmed"
@dataclass(frozen=True)
class ResidueKey:
    chain: str
    resid: int
    icode: str
def read_lines(path: Path) -> List[str]:
    return path.read_text().splitlines()
def write_text(path: Path, text: str) -> None:
    path.write_text(text)
def is_atom_record(line: str) -> bool:
    return line.startswith("ATOM") or line.startswith("HETATM")
def pdb_resname(line: str) -> str:
    return line[17:20].strip()
def pdb_chain(line: str) -> str:
    return line[21:22]
def pdb_resid(line: str) -> int:
    return int(line[22:26])
def pdb_icode(line: str) -> str:
    return line[26:27]
def pdb_atom_name(line: str) -> str:
    return line[12:16]
def pdb_coords(line: str) -> Tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])
def replace_resname(line: str, resname: str) -> str:
    return f"{line[:17]}{resname:>3}{line[20:]}"
def replace_atom_name(line: str, atom_name: str) -> str:
    return f"{line[:12]}{atom_name:>4}{line[16:]}"
def residue_key(line: str) -> ResidueKey:
    return ResidueKey(pdb_chain(line), pdb_resid(line), pdb_icode(line))
def map_pqr_resname_to_psfgen(name: str) -> str:
    return {
        "HID": "HSD",
        "HIE": "HSE",
        "HIP": "HSP",
        "ASH": "ASP",
        "GLH": "GLU",
    }.get(name, name)

def load_protonation_map(pqr_path: Path) -> Dict[ResidueKey, str]:
    protonation: Dict[ResidueKey, str] = {}
    for line in read_lines(pqr_path):
        if not is_atom_record(line):
            continue
        key = residue_key(line)
        protonation.setdefault(key, map_pqr_resname_to_psfgen(pdb_resname(line)))
    return protonation

def map_pdb_from_pqr(pdb_path: Path, pqr_path: Path) -> List[str]:
    protonation = load_protonation_map(pqr_path)
    out: List[str] = []
    for line in read_lines(pdb_path):
        if is_atom_record(line):
            key = residue_key(line)
            if key in protonation:
                line = replace_resname(line, protonation[key])
        out.append(line)
    return out

def generate_patch_lines(pqr_path: Path) -> List[str]:
    patches: List[str] = []
    seen: set[ResidueKey] = set()
    for line in read_lines(pqr_path):
        if not is_atom_record(line):
            continue
        key = residue_key(line)
        if key in seen:
            continue
        seen.add(key)
        resname = pdb_resname(line)
        if resname == "GLH":
            patches.append(f"patch GLUP {key.chain}:{key.resid}")
        elif resname == "ASH":
            patches.append(f"patch ASPP {key.chain}:{key.resid}")
    return patches

def protein_only_lines(mapped_lines: Sequence[str]) -> List[str]:
    atoms = [line for line in mapped_lines if line.startswith("ATOM")]
    if not atoms:
        return []
    max_resid_by_chain: Dict[str, int] = {}
    for line in atoms:
        chain = pdb_chain(line)
        max_resid_by_chain[chain] = max(max_resid_by_chain.get(chain, pdb_resid(line)), pdb_resid(line))
    fixed: List[str] = []
    for line in atoms:
        chain = pdb_chain(line)
        if pdb_resid(line) == max_resid_by_chain[chain]:
            atom = pdb_atom_name(line)
            if atom == " O  ":
                line = replace_atom_name(line, "OT1")
            elif atom == " OXT":
                line = replace_atom_name(line, "OT2")
        fixed.append(line)
    return fixed

def bounding_box(lines: Sequence[str]) -> Tuple[float, float, float, float, float, float]:
    coords = [pdb_coords(line) for line in lines if line.startswith("ATOM")]
    xs = [x for x, _, _ in coords]
    ys = [y for _, y, _ in coords]
    zs = [z for _, _, z in coords]
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)

def centered_rmax(lines: Sequence[str]) -> Tuple[float, Tuple[float, float, float]]:
    xmin, xmax, ymin, ymax, zmin, zmax = bounding_box(lines)
    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    cz = (zmin + zmax) / 2.0
    max_r = 0.0
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        x, y, z = pdb_coords(line)
        dx = x - cx
        dy = y - cy
        dz = z - cz
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        if r > max_r:
            max_r = r
    return max_r, (cx, cy, cz)

def net_charge_from_pqr(pqr_path: Path) -> float:
    total = 0.0
    for line in read_lines(pqr_path):
        if not is_atom_record(line):
            continue
        total += float(line[54:62])
    return total

def water_count_for_cube(box_side: float) -> int:
    volume = box_side ** 3
    return int(round(0.0334 * volume))

def ion_counts(net_charge: float, box_side: float, molarity: float | None) -> Tuple[int, int]:
    na = 0
    cl = 0
    if molarity and molarity > 0:
      volume_a3 = box_side ** 3
      pairs = int(round(molarity * 6.02214076e23 * volume_a3 * 1e-27))
      na += pairs
      cl += pairs
    rounded = int(round(net_charge))
    if rounded < 0:
        na += abs(rounded)
    elif rounded > 0:
        cl += rounded
    return na, cl

def packmol_input(box_side: float, waters: int, na: int, cl: int) -> str:
    half = box_side / 2.0
    lines = [
        "tolerance 2.0",
        "output solvated.pdb",
        "add_box_sides 1.0",
        "filetype pdb",
        "seed -1",
        "packall",
        f"pbc {-half:.1f} {-half:.1f} {-half:.1f}, {half:.1f} {half:.1f} {half:.1f}",
        "",
        "structure step1_pdbreader.pdb",
        "    number 1",
        "    center",
        "    fixed 0. 0. 0. 0. 0. 0.",
        "end structure",
        "",
        "structure WATER.pdb",
        f"    number {waters}",
        "end structure",
        "",
        "structure SOD.pdb",
        f"    number {na}",
        "end structure",
    ]
    if cl > 0:
        lines.extend([
            "",
            "structure CLA.pdb",
            f"    number {cl}",
            "end structure",
        ])
    return "\n".join(lines) + "\n"

def write_template_pdbs(outdir: Path) -> None:
    write_text(outdir / "WATER.pdb", "\n".join([
        "HETATM    1  OH2 TIP A   1       0.000   0.000   0.000  1.00  0.00           O",
        "HETATM    2  H1  TIP A   1       0.957   0.000   0.000  1.00  0.00           H",
        "HETATM    3  H2  TIP A   1      -0.239   0.927   0.000  1.00  0.00           H",
        "END",
        "",
    ]))
    write_text(outdir / "SOD.pdb", "\n".join([
        "HETATM    1 SOD  SOD A   1       0.000   0.000   0.000  1.00  0.00          NA",
        "END",
        "",
    ]))
    write_text(outdir / "CLA.pdb", "\n".join([
        "HETATM    1 CLA  CLA A   1       0.000   0.000   0.000  1.00  0.00          CL",
        "END",
        "",
    ]))

def write_splitters(outdir: Path) -> None:
    write_text(outdir / "split_solvated_packmol.awk", "\n".join([
        'BEGIN { water_out = "waters_packmol.pdb"; ion_out = "ions_packmol.pdb"; prot_out = "protein_from_packmol.pdb" }',
        '/^(ATOM|HETATM)/ {',
        '  res = substr($0, 18, 4); gsub(/ /, "", res)',
        '  if (res == "TIP" || res == "TIP3") print >> water_out',
        '  else if (res == "SOD" || res == "CLA") print >> ion_out',
        '  else print >> prot_out',
        '}',
        'END { print "END" >> water_out; print "END" >> ion_out; print "END" >> prot_out }',
        "",
    ]))
    write_text(outdir / "split_packmol_waters_by_chain.awk", "\n".join([
        '/^(ATOM|HETATM)/ { f="waters_packmol_" substr($0,22,1) ".pdb"; print >> f; seen[f]=1 }',
        'END { for (f in seen) print "END" >> f }',
        "",
    ]))
    write_text(outdir / "fix_psf_header.awk", 'NR == 1 { print "PSF CMAP XPLOR"; next } { print }\n')

def run(cmd: Sequence[str] | str, cwd: Path, stdin_text: str | None = None) -> None:
    subprocess.run(
        cmd,
        cwd=cwd,
        input=stdin_text,
        text=True,
        check=True,
        shell=isinstance(cmd, str),
    )

def write_build_scripts(outdir: Path, patches: Sequence[str], water_chain_files: Sequence[str] | None = None) -> None:
    write_text(outdir / "patches_from_pqr.tcl", "\n".join(patches) + ("\n" if patches else ""))
    protein_script = [
        "package require psfgen",
        f"topology {CHARMM_TOP}",
        f"topology {CHARMM_WATER_IONS}",
        "pdbalias atom ILE CD1 CD",
        "pdbalias atom PHE OXT OT2",
        "segment A {",
        "  first NTER",
        "  last CTER",
        "  pdb 9upt_psfgen_ph55_protein.pdb",
        "}",
        "coordpdb 9upt_psfgen_ph55_protein.pdb A",
        "source patches_from_pqr.tcl",
        "guesscoord",
        "writepsf x-plor cmap atbgl1a_protein.psf",
        "writepdb atbgl1a_protein.pdb",
        "exec awk -f fix_psf_header.awk atbgl1a_protein.psf > atbgl1a_protein.psf.tmp",
        "exec mv atbgl1a_protein.psf.tmp atbgl1a_protein.psf",
        "",
    ]
    write_text(outdir / "build_protein_psf.tcl", "\n".join(protein_script))

    solvated = [
        "package require psfgen",
        f"topology {CHARMM_TOP}",
        f"topology {CHARMM_WATER_IONS}",
        "pdbalias atom ILE CD1 CD",
        "pdbalias atom PHE OXT OT2",
        "pdbalias residue TIP TIP3",
        "segment A {",
        "  first NTER",
        "  last CTER",
        "  pdb 9upt_psfgen_ph55_protein.pdb",
        "}",
        "coordpdb 9upt_psfgen_ph55_protein.pdb A",
        "source patches_from_pqr.tcl",
    ]
    if water_chain_files:
        for idx, filename in enumerate(water_chain_files):
            seg = f"WT{idx:02d}"
            solvated.extend([
                f"segment {seg} {{",
                "  auto none",
                f"  pdb {filename}",
                "}",
                f"coordpdb {filename} {seg}",
            ])
    else:
        solvated.extend([
            "segment WT1 {",
            "  auto none",
            "  pdb waters_packmol.pdb",
            "}",
            "coordpdb waters_packmol.pdb WT1",
        ])
    solvated.extend([
        "segment ION {",
        "  auto none",
        "  pdb ions_packmol.pdb",
        "}",
        "coordpdb ions_packmol.pdb ION",
        "guesscoord",
        "writepsf x-plor cmap atbgl1a_solvated.psf",
        "writepdb atbgl1a_solvated.pdb",
        "exec awk -f fix_psf_header.awk atbgl1a_solvated.psf > atbgl1a_solvated.psf.tmp",
        "exec mv atbgl1a_solvated.psf.tmp atbgl1a_solvated.psf",
        "",
    ])
    write_text(outdir / "build_solvated_psf.tcl", "\n".join(solvated))

    parmed_script = "\n".join([
        "chamber \\",
        f"  -top {CHARMM_TOP} \\",
        f"  -param {CHARMM_PAR} \\",
        f"  -str {CHARMM_WATER_IONS} \\",
        "  -psf atbgl1a_solvated.psf \\",
        "  -crd atbgl1a_solvated.pdb \\",
        "  -box bounding",
        "",
        "parmout atbgl1a_charmm.parm7 atbgl1a_charmm.rst7",
        "go",
        "",
    ])
    write_text(outdir / "to_amber_from_charmm.parmed.in", parmed_script)

def write_amber_inputs(outdir: Path, temperatures: Iterable[int]) -> None:
    write_text(outdir / "min1.in", "\n".join([
        "Initial minimization with heavy restraints",
        "&cntrl",
        "  imin=1, maxcyc=10000, ncyc=5000, cut=10.0, ntb=1, ntx=1, irest=0,",
        "  ntpr=100, ntr=1, restraint_wt=10.0, restraintmask='!:WAT & !@H=',",
        "/",
        "",
    ]))
    write_text(outdir / "min2.in", "\n".join([
        "Unrestrained minimization",
        "&cntrl",
        "  imin=1, maxcyc=15000, ncyc=7500, cut=10.0, ntb=1, ntx=1, irest=0,",
        "  ntpr=100, ntr=0,",
        "/",
        "",
    ]))
    for temp in temperatures:
        heat = "\n".join([
            f"Heating to {temp} K under NVT with protein heavy-atom restraints",
            "&cntrl",
            "  imin=0, ntx=1, irest=0, nstlim=250000, dt=0.002, ntc=2, ntf=2,",
            "  cut=10.0, ntb=1, ntp=0, tempi=50.0,",
            f"  temp0={temp}.0, ntt=3, gamma_ln=2.0, ntpr=5000, ntwx=5000, ntwr=5000,",
            "  ig=-1, ioutfm=1, ntr=1, restraint_wt=5.0, restraintmask='!:WAT & !@H=',",
            "/",
            f"&wt type='TEMP0', istep1=1, istep2=250000, value1=50.0, value2={temp}.0 /",
            "&wt type='END' /",
            "",
        ])
        write_text(outdir / f"heat_{temp}K.in", heat)
    write_text(outdir / "eq_npt.in", "\n".join([
        "NPT equilibration at target temperature with lighter restraints",
        "&cntrl",
        "  imin=0, ntx=5, irest=1, nstlim=500000, dt=0.002, ntc=2, ntf=2, cut=10.0,",
        "  ntb=2, ntp=1, pres0=1.0, taup=2.0, temp0=343.0, ntt=3, gamma_ln=2.0,",
        "  ntpr=5000, ntwx=5000, ntwr=5000, ig=-1, ioutfm=1, ntr=1,",
        "  restraint_wt=1.0, restraintmask='!:WAT & !@H=',",
        "/",
        "",
    ]))
    write_text(outdir / "eq_relax.in", "\n".join([
        "Unrestrained equilibration under NPT",
        "&cntrl",
        "  imin=0, ntx=5, irest=1, nstlim=2500000, dt=0.002, ntc=2, ntf=2, cut=10.0,",
        "  ntb=2, ntp=1, pres0=1.0, taup=2.0, temp0=343.0, ntt=3, gamma_ln=2.0,",
        "  ntpr=5000, ntwx=5000, ntwr=5000, ig=-1, ioutfm=1, ntr=0,",
        "/",
        "",
    ]))
    write_text(outdir / "prod.in", "\n".join([
        "Production MD under NPT",
        "&cntrl",
        "  imin=0, ntx=5, irest=1, nstlim=50000000, dt=0.002, ntc=2, ntf=2, cut=10.0,",
        "  ntb=2, ntp=1, pres0=1.0, taup=2.0, temp0=343.0, ntt=3, gamma_ln=2.0,",
        "  ntpr=5000, ntwx=5000, ntwr=5000, ig=-1, ioutfm=1, ntr=0,",
        "/",
        "",
    ]))

def split_waters_by_chain(outdir: Path) -> List[str]:
    waters = outdir / "waters_packmol.pdb"
    chain_files: Dict[str, List[str]] = {}
    for line in read_lines(waters):
        if not is_atom_record(line):
            continue
        chain = pdb_chain(line)
        chain_files.setdefault(chain, []).append(line)
    filenames: List[str] = []
    for chain, lines in sorted(chain_files.items()):
        filename = f"waters_packmol_{chain}.pdb"
        write_text(outdir / filename, "\n".join(lines + ["END", ""]))
        filenames.append(filename)
    return filenames


def prepare_only(args: argparse.Namespace) -> None:
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    mapped = map_pdb_from_pqr(args.pdb, args.pqr)
    protein = protein_only_lines(mapped)
    write_text(outdir / "mapped_from_pqr.pdb", "\n".join(mapped) + "\n")
    write_text(outdir / "9upt_psfgen_ph55_protein.pdb", "\n".join(protein) + "\n")
    shutil.copyfile(outdir / "9upt_psfgen_ph55_protein.pdb", outdir / "step1_pdbreader.pdb")

    patches = generate_patch_lines(args.pqr)
    write_splitters(outdir)
    write_template_pdbs(outdir)
    write_build_scripts(outdir, patches)
    write_amber_inputs(outdir, args.temperatures)

    net = net_charge_from_pqr(args.pqr)
    xmin, xmax, ymin, ymax, zmin, zmax = bounding_box(protein)
    rmax, center = centered_rmax(protein)
    box_side = args.box_side
    if box_side is None:
        box_side = 2.0 * (rmax + args.padding)
        if args.round_box and args.round_box > 0:
            box_side = float(math.ceil(box_side / args.round_box) * args.round_box)
    waters = water_count_for_cube(box_side)
    na, cl = ion_counts(net, box_side, None if args.neutralize_only else args.salt_molarity)

    prep_summary = "\n".join([
        f"net_charge {net:.3f}",
        f"protein_dimensions {xmax - xmin:.3f} {ymax - ymin:.3f} {zmax - zmin:.3f}",
        f"bounding_box_center {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}",
        f"rmax {rmax:.3f}",
        f"padding {args.padding:.3f}",
        f"box_side {box_side:.3f}",
        f"water_count {waters}",
        f"sodium_count {na}",
        f"chloride_count {cl}",
        "",
    ])
    write_text(outdir / "prep_summary.txt", prep_summary)
    write_text(outdir / "packmol_system.inp", packmol_input(box_side, waters, na, cl))

def full_run(args: argparse.Namespace) -> None:
    outdir = args.outdir.resolve()
    prepare_only(args)

    run(f"{args.packmol} < packmol_system.inp", cwd=outdir)
    run(["awk", "-f", "split_solvated_packmol.awk", "solvated.pdb"], cwd=outdir)
    water_chain_files = split_waters_by_chain(outdir)
    write_build_scripts(outdir, generate_patch_lines(args.pqr), water_chain_files)
    run([args.psfgen, "build_protein_psf.tcl"], cwd=outdir)
    run([args.psfgen, "build_solvated_psf.tcl"], cwd=outdir)
    run([args.parmed, "-n", "-O", "-i", "to_amber_from_charmm.parmed.in"], cwd=outdir)

def parse_temperatures(text: str) -> List[int]:
    return [int(part.strip()) for part in text.split(",") if part.strip()]

def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare a pdb2pqr -> Packmol -> psfgen -> ParmEd -> AMBER workflow."
    )
    parser.add_argument("--pdb", type=Path, required=True, help="Original PDB file.")
    parser.add_argument("--pqr", type=Path, required=True, help="pdb2pqr output with desired protonation.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--box-side", type=float, default=None, help="Cubic box side length in A. If omitted, compute from rmax and padding.")
    parser.add_argument("--padding", type=float, default=15.0, help="Target solvent padding in A when auto-computing box size.")
    parser.add_argument("--round-box", type=float, default=10.0, help="Round auto-computed box side up to this increment in A. Use 0 to disable rounding.")
    parser.add_argument("--neutralize-only", action="store_true", help="Neutralize only; do not add bulk salt.")
    parser.add_argument("--salt-molarity", type=float, default=0.15, help="NaCl molarity if not neutralize-only.")
    parser.add_argument("--temperatures", type=parse_temperatures, default=[300, 338, 343, 353], help="Comma-separated AMBER heat temperatures.")
    parser.add_argument("--packmol", default=DEFAULT_PACKMOL, help="Packmol executable.")
    parser.add_argument("--psfgen", default=DEFAULT_PSFGEN, help="psfgen executable.")
    parser.add_argument("--parmed", default=DEFAULT_PARMED, help="ParmEd executable.")
    parser.add_argument("--prepare-only", action="store_true", help="Write workflow files only; do not run Packmol/psfgen/ParmEd.")
    return parser

def main() -> None:
    parser = make_parser()
    args = parser.parse_args()
    if args.prepare_only:
        prepare_only(args)
    else:
        full_run(args)


if __name__ == "__main__":
    main()
