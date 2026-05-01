
# automate_charmm_packmol_amber

Automates a `pdb2pqr -> Packmol -> psfgen -> ParmEd -> AMBER` preparation workflow for protein MD systems. The repository also keeps a concrete AtBgl1A example at `pH 5.5`, including CHARMM-based system building, AMBER inputs.

## What This Repository Contains

- `prepare_charmm_packmol_amber.py`
  Main helper script. It maps protonation states from a `.pqr`, writes Packmol and `psfgen` inputs, generates AMBER MD input files, and can optionally run the whole build.
  
## Workflow Summary

The Python script uses the protonation pattern from `pdb2pqr` as the source of truth and prepares a CHARMM-to-AMBER pipeline:

1. Read the original `.pdb` and the protonated `.pqr`.
2. Map `pdb2pqr` residue names into `psfgen`/CHARMM-compatible names.
3. Generate protonation patch commands for protonated acids.
4. Write a protein-only PDB and Packmol input files.
5. Write `psfgen` scripts for unsolvated and solvated systems.
6. Write AMBER minimization, heating, equilibration, and production inputs.
7. Optionally run `Packmol`, `psfgen`, and `ParmEd` to produce AMBER-ready topology and coordinates.

## Protonation Mapping

The script converts common `pdb2pqr` names into CHARMM naming conventions:

- `HID -> HSD`
- `HIE -> HSE`
- `HIP -> HSP`
- `ASH -> ASP` plus `ASPP` patch
- `GLH -> GLU` plus `GLUP` patch

This matters because protonated acids should not remain as `ASH` or `GLH` in the `psfgen` input PDB. They are converted to standard residue names and then protonated through patches.

## Requirements

The current defaults in `prepare_charmm_packmol_amber.py` assume:

- `packmol` is available in `PATH`
- `psfgen` is available at `/home/shabir/Downloads/NAMD3.1/psfgen`
- `parmed` is available at `/home/shabir/Downloads/ambertools25/bin/parmed`
- CHARMM36 topology and parameter files are available from the local VMD plugin installation:
  - `/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf`
  - `/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmpar1.5/par_all36_prot.prm`
  - `/usr/local/lib/vmd/plugins/noarch/tcl/trunctraj1.5/toppar/stream/toppar_water_ions.str`

If your tools live elsewhere, pass `--packmol`, `--psfgen`, and `--parmed` explicitly.

## Usage

### 1. Prepare files only

This writes the workflow files without running external tools.

```bash
python3 prepare_charmm_packmol_amber.py \
  --pdb receptor.pdb \
  --pqr receptor_ph5.5.pqr \
  --outdir generic_receptor_autobox \
  --padding 15 \
  --neutralize-only \
  --prepare-only
```

### 2. Prepare and run the full build

```bash
python3 prepare_charmm_packmol_amber.py \
  --pdb receptor.pdb \
  --pqr receptor_ph5.5.pqr \
  --outdir my_system \
  --box-side 110 \
  --neutralize-only
```

### 3. Common options

- `--box-side 110`
  Use a fixed cubic box size in angstrom.
- `--padding 15`
  When `--box-side` is omitted, set solvent padding around the protein.
- `--round-box 10`
  Round the auto-computed box side up to this increment.
- `--neutralize-only`
  Add only counterions needed to neutralize the system.
- `--salt-molarity 0.15`
  Add bulk NaCl if `--neutralize-only` is not used.
- `--temperatures 300,338,343,353`
  Write heating inputs for these target temperatures.

## Generated Outputs

Each output directory can contain:

- `mapped_from_pqr.pdb`
  Input PDB with protonation states translated from the `.pqr`.
- `receptor_psfgen_ph55_protein.pdb`
  Protein-only PDB for `psfgen`.
- `patches_from_pqr.tcl`
  `ASPP` and `GLUP` protonation patches.
- `packmol_system.inp`
  Packmol system definition.
- `build_protein_psf.tcl`, `build_solvated_psf.tcl`
  `psfgen` build scripts.
- `to_amber_from_charmm.parmed.in`
  ParmEd conversion script.
- `min1.in`, `min2.in`, `heat_*K.in`, `eq_npt.in`, `eq_relax.in`, `prod.in`
  AMBER MD templates.
- `prep_summary.txt`
  Computed box, charge, and solvent-count summary.

If the full run is executed, the workflow also produces:

- `solvated.pdb`
- `protein.psf`, `atbgl1a_protein.pdb`
- `solvated.psf`, `atbgl1a_solvated.pdb`
- `receptor_charmm.parm7`, `receptor_charmm.rst7`

## Example

The current example setup is based on:

- input structure: `receptor.pdb`
- protonation source: `receptor_ph5.5.pqr`
- target study pH: `5.5`
- temperature set: `300 K`, `338 K`, `343 K`, `353 K`
## Notes

- The final `parm7/rst7` files are AMBER-format outputs derived from a CHARMM36-based build. The force field does not become native AMBER protein parameters just because the file format changes.
- The generated `eq_npt.in`, `eq_relax.in`, and `prod.in` currently use `temp0=343.0` by default. If you run other temperatures, adjust those files or generate separate stage files as needed.
- The current build scripts assume the protein segment is chain `A`. If your system uses multiple protein chains, update the generated `psfgen` scripts accordingly.
