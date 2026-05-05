
# automate_CHARMMff_packmol_amber

This workflow uses your `pdb2pqr` protonation assignment at `pH 5.5` as the source of truth, then maps it into CHARMMff/psfgen naming.

Pipeline:

1. `pdb2pqr` protonation at `pH 5.5`
2. map residue names into `psfgen`-compatible PDB naming
3. generate protonation patch commands for `ASPP` and `GLUP`
4. build `PSF/PDB` with `psfgen`
5. solvate with `Packmol`
6. build final solvated `PSF/PDB` with `psfgen`
7. convert CHARMMff `PSF/PDB` to AMBER package `parm7/rst7` with ParmEd `chamber`

## Important naming translation

`pdb2pqr` names are not the same as CHARMMff names:

- `HID` -> `HSD`
- `HIE` -> `HSE`
- `HIP` -> `HSP`
- `ASH` -> `ASP` plus `ASPP` patch
- `GLH` -> `GLU` plus `GLUP` patch

For this reason, protonated acids must **not** remain as `ASH/GLH` in the psfgen PDB. They must be converted to `ASP/GLU` and then patched inside `psfgen`.

## Files

- `receptor_psfgen_ph5.5.pdb`
  Full mapped PDB with `pdb2pqr` protonation mapped into CHARMMff residue names.
- `receptor_psfgen_ph5.5_protein.pdb`
  Protein-only version used by `psfgen` for the protein segment.
- `step1_pdbreader.pdb`
  Working protein file for your Packmol-style input.
- `WATER.pdb`, `SOD.pdb`, `CLA.pdb`
  Minimal template molecules for Packmol using CHARMMff-compatible names.
- `patches_from_pqr.tcl`
  Generated patch commands for protonated acidic residues.
- `build_protein_psf.tcl`
  Builds the unsolvated protein PSF/PDB.
- `build_solvated_psf.tcl`
  Example solvated-system psfgen build script for protein + waters + ions.
- `packmol_water_ions_template.inp`
  Template for Packmol solvation.
- `to_amber_from_CHARMMff.parmed.in`
  ParmEd script that converts the CHARMMff-built system to `parm7/rst7`.
- `split_solvated_packmol.awk`
  Splits `solvated.pdb` into `protein_from_packmol.pdb`, `waters_packmol.pdb`, and `ions_packmol.pdb`.

## Local topology/parameter files used

- `top_all36_prot.rtf`
  `/usr/local/lib/vmd/plugins/noarch/tcl/readCHARMMfftop1.2/top_all36_prot.rtf`
- `par_all36_prot.prm`
  `/usr/local/lib/vmd/plugins/noarch/tcl/readCHARMMffpar1.5/par_all36_prot.prm`
- `toppar_water_ions.str`
  `/usr/local/lib/vmd/plugins/noarch/tcl/trunctraj1.5/toppar/stream/toppar_water_ions.str`

These files define the actual force field used in the workflow:

- `CHARMMff36` protein topology/parameters
- CHARMMff water/ion definitions from `toppar_water_ions.str`

The PSF format written by the scripts is requested explicitly as:

- `X-PLOR`
- with `CMAP`

using:

```tcl
writepsf x-plor cmap receptor_protein.psf
```

and similarly for the solvated system.

Because this `psfgen` build still writes `PSF CMAP` on the first line while indicating `x-plor psf file` in the remarks, the workflow applies a final header normalization step so the first line is strictly:

```text
PSF CMAP XPLOR
```

## Generated protonation state examples

From your `pdb2pqr` file at `pH 5.5`:

- catalytic `E166` is protonated in `pdb2pqr`, so this becomes:
  - `GLU 166` in the PDB
  - `patch GLUP A:166` in psfgen
- catalytic `E355` remains standard `GLU`
- histidines are mapped to `HSD/HSE/HSP`

## Typical usage

### 1. Solvate with Packmol

For receptor in this workspace, a concrete starting file is provided:

- `packmol_receptor_110A_cube.inp`

It uses:

- a `110 A` cubic box
- protein centering by Packmol
- `2` sodium ions for neutralization only
- about `44,500` waters

This corresponds to:

- protein net charge `-2` at `pH 5.5`
- no added bulk salt
- neutralization only

Run Packmol from this directory:

```bash
packmol < packmol_receptor_110A_cube.inp
```
Then split the resulting `solvated.pdb`:

```bash
awk -f split_solvated_packmol.awk solvated.pdb
```
then split the resulting `waters_packmol.pb`:
```bash
awk -f split_packmol_waters_by_chain.awk waters_packmol.pdb
```
### 2. Build protein-only PSF/PDB

```bash
cd /media/shabir/Coaraci/GH/psfgen_receptor
/home/shabir/Downloads/NAMD3.1/psfgen build_protein_psf.tcl
```

### 3. Build solvated PSF/PDB

```bash
/home/shabir/Downloads/NAMD3.1/psfgen build_solvated_psf.tcl
```

### 4. Convert to AMBER package `parm7/rst7`

```bash
/home/shabir/Downloads/ambertools25/bin/parmed -n -O -i to_amber_from_CHARMMff.parmed.in
```

The ParmEd conversion script uses:

- `-box bounding`

This is appropriate when the solvated PDB from your `Packmol + psfgen` workflow does not already carry reliable periodic box metadata and you want ParmEd to infer a rectangular bounding box from the coordinates.

## Important force-field note

After conversion, `parm7/rst7` are AMBER-format files, but the force field remains **CHARMMff-derived**. This is a format conversion, not a change to AMBER package protein parameters.

## Practical note on waters and ions

If you use Packmol, the water and ion residue names and atom names must match the CHARMMff topology exactly. For the supplied example:

- water residue: `TIP3`
- water atoms: `OH2`, `H1`, `H2`
- sodium residue/atom: `SOD`
- chloride residue/atom: `CLA`

If your Packmol output uses different names, rename them before `psfgen`.

## Usage

### 1. Prepare files only

This writes the workflow files without running external tools.

```bash
python3 prepare_CHARMMff_packmol_amber.py \
  --pdb receptor.pdb \
  --pqr receptor_ph5.5.pqr \
  --outdir generic_receptor_autobox \
  --padding 15 \
  --neutralize-only \
  --prepare-only
```

### 2. Prepare and run the full build

```bash
python3 prepare_CHARMMff_packmol_amber.py \
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
- `receptor_psfgen_ph5.5_protein.pdb`
  Protein-only PDB for `psfgen`.
- `patches_from_pqr.tcl`
  `ASPP` and `GLUP` protonation patches.
- `packmol_system.inp`
  Packmol system definition.
- `build_protein_psf.tcl`, `build_solvated_psf.tcl`
  `psfgen` build scripts.
- `to_amber_from_CHARMMff.parmed.in`
  ParmEd conversion script.
- `min1.in`, `min2.in`, `heat_*K.in`, `eq_npt.in`, `eq_relax.in`, `prod.in`
  AMBER MD templates.
- `prep_summary.txt`
  Computed box, charge, and solvent-count summary.

If the full run is executed, the workflow also produces:

- `solvated.pdb`
- `protein.psf`, `receptor_protein.pdb`
- `solvated.psf`, `receptor_solvated.pdb`
- `receptor_CHARMMff.parm7`, `receptor_CHARMMff.rst7`

## Notes

- The final `parm7/rst7` files are AMBER-format outputs derived from a CHARMMff36-based build. The force field does not become native AMBER protein parameters just because the file format changes.
- The generated `eq_npt.in`, `eq_relax.in`, and `prod.in` currently use `temp0=343.0` by default. If you run other temperatures, adjust those files or generate separate stage files as needed.
- The current build scripts assume the protein segment is chain `A`. If your system uses multiple protein chains, update the generated `psfgen` scripts accordingly.
