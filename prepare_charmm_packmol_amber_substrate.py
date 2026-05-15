#!/usr/bin/env python3
#1. Manual box:
#python3 prepare_charmm_packmol_amber.py --pdb input.pdb --pqr input_ph5.5.pqr --outdir my_system --box-side 110 --neutralize-only
#2. Auto box from protein size:
#python3 prepare_charmm_packmol_amber.py --pdb input.pdb --pqr input_ph5.5.pqr --outdir my_system--neutralize-only --padding 15

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
CHARMM_CGENFF = "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_cgenff.rtf"
CHARMM_CGENFF_PAR = "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmpar1.5/par_all36_cgenff.prm"
CHARMM_CARB = "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_carb.rtf"
CHARMM_CARB_PAR = "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmpar1.5/par_all36_carb.prm"

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
    keep = []

    for line in mapped_lines:
        if not is_atom_record(line):
            continue

        resname = pdb_resname(line)

        # Exclude solvent/ions
        if resname in {"TIP3", "TIP3", "WAT", "HOH", "SOD", "CLA"}:
            continue

        keep.append(line)

    if not keep:
        return []

    # Fix terminal oxygens ONLY for protein atoms
    max_resid_by_chain: Dict[str, int] = {}
    for line in keep:
        if line.startswith("ATOM"):  # only protein
            chain = pdb_chain(line)
            max_resid_by_chain[chain] = max(
                max_resid_by_chain.get(chain, pdb_resid(line)),
                pdb_resid(line),
            )

    fixed: List[str] = []
    for line in keep:
        if line.startswith("ATOM"):
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

        # ONLY protein atoms
        if line.startswith("ATOM"):
            parts = line.split()
            try:
                total += float(parts[-2])
            except (IndexError, ValueError):
                raise ValueError(f"Bad PQR line: {line}")

    return round(total, 3)


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
        "HETATM    1  OH2 TIP3A   1       0.000   0.000   0.000  1.00  0.00           O",
        "HETATM    2  H1  TIP3A   1       0.957   0.000   0.000  1.00  0.00           H",
        "HETATM    3  H2  TIP3A   1      -0.239   0.927   0.000  1.00  0.00           H",
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
        'BEGIN {',
        '  water_out = "waters_packmol.pdb";',
        '  ion_out = "ions_packmol.pdb";',
        '  prot_out = "protein_from_packmol.pdb";',
        '  sub_out = "substrate_from_packmol.pdb";',
        '}',
        '',
        '/^(ATOM|HETATM)/ {',
        '  res = substr($0,18,4);',
        '  gsub(/ /,"",res);',
        '',
        '  # ---------------- WATER ----------------',
        '  if (res == "TIP" || res == "TIP3") {',
        '    print >> water_out;',
        '    next;',
        '  }',
        '',
        '  # ---------------- IONS ----------------',
        '  else if (res == "SOD" || res == "CLA") {',
        '    print >> ion_out;',
        '    next;',
        '  }',
        '',
        '  # ---------------- SUBSTRATE ----------------',
        '  else if (res == "GLC" || res == "STR") {',
        '    print >> sub_out;',
        '    next;',
        '  }',
        '',
        '  # ---------------- PROTEIN ----------------',
        '  else {',
        '    print >> prot_out;',
        '  }',
        '}',
        '',
        'END {',
        '  print "END" >> water_out;',
        '  print "END" >> ion_out;',
        '  print "END" >> prot_out;',
        '  print "END" >> sub_out;',
        '}',
        ""
    ]))

    write_text(outdir / "fix_psf_header.awk",
        'NR == 1 { print "PSF CMAP x-plor"; next } { print }\n'
    )

    write_text(outdir / "split_waters.py", "\n".join([
        '"""Split waters_packmol.pdb into 3 equal PDB files with unique residue numbers.',
        '',
        'Run after:  awk -f split_solvated_packmol.awk solvated.pdb',
        'Produces:   waters_A.pdb  waters_B.pdb  waters_C.pdb',
        '"""',
        'import math',
        'from pathlib import Path',
        '',
        'atoms = [l for l in Path("waters_packmol.pdb").read_text().splitlines()',
        '         if l.startswith(("ATOM", "HETATM"))]',
        'assert len(atoms) % 3 == 0, "Water atom count not divisible by 3"',
        'n_waters = len(atoms) // 3',
        'n = math.ceil(n_waters / 9999)',
        'sizes = [9999] * (n - 1) + [n_waters - 9999 * (n - 1)]',
        'wi = 0',
        'for i, count in enumerate(sizes):',
        '    fname = f"waters_{chr(65 + i)}.pdb"',
        '    out = []',
        '    for resid in range(1, count + 1):',
        '        for k in range(3):',
        '            line = atoms[wi * 3 + k]',
        '            out.append(line[:22] + f"{resid:4d}" + line[26:])',
        '        wi += 1',
        '    Path(fname).write_text("\\n".join(out + ["END", ""]))',
        '    print(f"Wrote {fname}: {count} waters")',
        '',
    ]))


def run(cmd: Sequence[str] | str, cwd: Path, stdin_text: str | None = None) -> None:
    subprocess.run(
        cmd,
        cwd=cwd,
        input=stdin_text,
        text=True,
        check=True,
        shell=isinstance(cmd, str),
    )


def write_build_scripts(
    outdir: Path,
    patches: Sequence[str],
    water_chain_files=None,
    n_water_segments: int = 4,
    ligand_str: str | None = None,
    substrate_pdb: str = "substrate_from_packmol_h.pdb",
    use_carb: bool = False,
) -> None:
    water_chain_files = water_chain_files or []

    write_text(
        outdir / "patches_from_pqr.tcl",
        "\n".join(patches) + ("\n" if patches else "")
    )

    # -------------------------
    # PROTEIN PSF
    # -------------------------
    protein_script = [
        "package require psfgen",
        f"topology {CHARMM_TOP}",
        f"topology {CHARMM_WATER_IONS}",
        "pdbalias atom ILE CD1 CD",
        "pdbalias atom PHE OXT OT2",

        "segment A {",
        "  first NTER",
        "  last CTER",
        "  pdb protein_from_packmol.pdb",
        "}",
        "coordpdb protein_from_packmol.pdb A",

        "source patches_from_pqr.tcl",
        "guesscoord",

        "writepsf x-plor cmap atbgl1a_protein.psf",
        "writepdb atbgl1a_protein.pdb",

        "exec awk -f fix_psf_header.awk atbgl1a_protein.psf > atbgl1a_protein.psf.tmp",
        "exec mv atbgl1a_protein.psf.tmp atbgl1a_protein.psf",
        "",
    ]
    write_text(outdir / "build_protein_psf.tcl", "\n".join(protein_script))

    # -------------------------
    # SOLVATED + SUBSTRATE SYSTEM
    # -------------------------
    solvated_top = [
        "package require psfgen",
        f"topology {CHARMM_TOP}",
        f"topology {CHARMM_WATER_IONS}",
        f"topology {CHARMM_CGENFF}",
        f"topology {CHARMM_CARB}",
    ]
    if ligand_str:
        solvated_top += [
            f"topology {ligand_str}",
        ]
    solvated_top += [
        "pdbalias atom ILE CD1 CD",
        "pdbalias atom PHE OXT OT2",
        "pdbalias residue TIP3 TIP3",
    ]

    solvated = solvated_top + [
        "# --------------------",
        "# PROTEIN",
        "# --------------------",
        "segment A {",
        "  first NTER",
        "  last CTER",
        "  pdb protein_from_packmol.pdb",
        "}",
        "coordpdb protein_from_packmol.pdb A",

        "# --------------------",
        "# SUBSTRATE (GLC / ligands)",
        "# --------------------",
        "segment SUB {",
        "  auto none",
        f"  pdb {substrate_pdb}",
        "}",
        f"coordpdb {substrate_pdb} SUB",
    ]

    # -------------------------
    # WATER HANDLING
    # Always use pre-split files; full_run() produces them via
    # split_waters_n_segments(), prepare-only users run split_waters.py first.
    # -------------------------
    if not water_chain_files:
        water_chain_files = [f"waters_{chr(65 + i)}.pdb" for i in range(n_water_segments)]
    for idx, filename in enumerate(water_chain_files):
        seg = f"WT{idx:02d}"
        solvated.extend([
            f"segment {seg} {{",
            "  auto none",
            f"  pdb {filename}",
            "}",
            f"coordpdb {filename} {seg}",
        ])

    # -------------------------
    # IONS + FINALIZE
    # -------------------------
    solvated.extend([
        "segment ION {",
        "  auto none",
        "  pdb ions_packmol.pdb",
        "}",
        "coordpdb ions_packmol.pdb ION",

        "source patches_from_pqr.tcl",
        "guesscoord",

        "writepsf x-plor cmap atbgl1a_solvated.psf",
        "writepdb atbgl1a_solvated.pdb",

        "exec awk -f fix_psf_header.awk atbgl1a_solvated.psf > atbgl1a_solvated.psf.tmp",
        "exec mv atbgl1a_solvated.psf.tmp atbgl1a_solvated.psf",
        "",
    ])

    write_text(outdir / "build_solvated_psf.tcl", "\n".join(solvated))

    parmed_lines = [
        "chamber \\",
        f"  -top {CHARMM_TOP} \\",
        f"  -param {CHARMM_PAR} \\",
        f"  -param {CHARMM_CGENFF_PAR} \\",
        f"  -param {CHARMM_CARB_PAR} \\",
    ]
    if ligand_str:
        parmed_lines += [
            f"  -str {ligand_str} \\",
        ]
    parmed_lines += [
        f"  -str {CHARMM_WATER_IONS} \\",
        "  -psf atbgl1a_solvated.psf \\",
        "  -crd atbgl1a_solvated.pdb \\",
        "  -box bounding",
        "",
        "parmout atbgl1a_charmm.parm7 atbgl1a_charmm.rst7",
        "go",
        "",
    ]
    write_text(outdir / "to_amber_from_charmm.parmed.in", "\n".join(parmed_lines))


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


def split_waters_n_segments(outdir: Path) -> List[str]:
    """Split waters_packmol.pdb into PDB files of at most 9999 residues each.

    waters_packmol.pdb is the direct output of split_solvated_packmol.awk,
    which passes water records through unchanged from solvated.pdb (standard
    PDB column layout, 4-char resnum at cols 23-26).  This function assigns
    local sequential residue numbers (1..count) per output file, keeping
    coordinates at their standard column positions.
    Number of segments = ceil(n_waters / 9999).
    """
    atoms = [l for l in read_lines(outdir / "waters_packmol.pdb") if is_atom_record(l)]
    assert len(atoms) % 3 == 0, "Water atom count not divisible by 3 — check waters_packmol.pdb"
    n_waters = len(atoms) // 3

    n = math.ceil(n_waters / 9999)
    sizes = [9999] * (n - 1) + [n_waters - 9999 * (n - 1)]

    filenames: List[str] = []
    water_idx = 0
    for i, count in enumerate(sizes):
        label = chr(ord("A") + i)
        fname = f"waters_{label}.pdb"
        out_lines: List[str] = []
        for local_resid in range(1, count + 1):
            for k in range(3):  # OH2, H1, H2
                line = atoms[water_idx * 3 + k]
                line = line[:22] + f"{local_resid:4d}" + line[26:]
                out_lines.append(line)
            water_idx += 1
        write_text(outdir / fname, "\n".join(out_lines + ["END", ""]))
        filenames.append(fname)
    return filenames


def prepare_only(args: argparse.Namespace) -> None:
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    mapped = map_pdb_from_pqr(args.pdb, args.pqr)
    protein = protein_only_lines(mapped)
    write_text(outdir / "mapped_from_pqr.pdb", "\n".join(mapped) + "\n")
    write_text(outdir / "9upt_psfgen_ph5.5_protein.pdb", "\n".join(protein) + "\n")
    shutil.copyfile(outdir / "9upt_psfgen_ph5.5_protein.pdb", outdir / "step1_pdbreader.pdb")

    net = net_charge_from_pqr(args.pqr)
    print(f"[DEBUG] Raw net charge = {net}")
    xmin, xmax, ymin, ymax, zmin, zmax = bounding_box(protein)
    rmax, center = centered_rmax(protein)
    box_side = args.box_side
    if box_side is None:
        box_side = 2.0 * (rmax + args.padding)
        if args.round_box and args.round_box > 0:
            box_side = float(math.ceil(box_side / args.round_box) * args.round_box)
    waters = water_count_for_cube(box_side)
    print(f"[DEBUG] Box side = {box_side}")
    print(f"[DEBUG] Waters = {waters}")
    na, cl = ion_counts(net, box_side, None if args.neutralize_only else args.salt_molarity)

    patches = generate_patch_lines(args.pqr)
    write_splitters(outdir)
    write_template_pdbs(outdir)
    ligand_str = str(args.ligand_str.resolve()) if args.ligand_str else None
    n_water_segments = math.ceil(waters / 9999)
    write_build_scripts(
        outdir, patches,
        n_water_segments=n_water_segments,
        ligand_str=ligand_str,
        substrate_pdb=args.substrate_pdb,
        use_carb=args.carb,
    )
    write_amber_inputs(outdir, args.temperatures)
    print(f"[DEBUG] Sodium (Na) = {na}")
    print(f"[DEBUG] Chloride (Cl) = {cl}")
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
    water_chain_files = split_waters_n_segments(outdir)
    patches = generate_patch_lines(args.pqr)
    ligand_str = str(args.ligand_str.resolve()) if args.ligand_str else None
    write_build_scripts(
        outdir, patches,
        water_chain_files=water_chain_files,
        ligand_str=ligand_str,
        substrate_pdb=args.substrate_pdb,
        use_carb=args.carb,
    )
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
    parser.add_argument(
        "--ligand-str", type=Path, default=None,
        help="CGenFF .str file for the non-standard ligand/substrate. "
             "Loads top_all36_cgenff.rtf + this .str in psfgen and "
             "par_all36_cgenff.prm + this .str in ParmEd.",
    )
    parser.add_argument(
        "--substrate-pdb", default="substrate_from_packmol_h.pdb",
        help="Substrate PDB filename used in the psfgen build script "
             "(default: substrate_from_packmol_h.pdb, the CGenFF-named hydrogenated file). "
             "Must have atom names matching the CGenFF topology STR.",
    )
    parser.add_argument(
        "--carb", action="store_true",
        help="Load CHARMM36 carbohydrate force field (top_all36_carb.rtf / "
             "par_all36_carb.prm) in psfgen and ParmEd. Use for standard "
             "sugars (AGLC, BGLC, BGAL, etc.) from the CHARMM carb topology.",
    )
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
