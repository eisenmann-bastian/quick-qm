# Molecular aligner

Select molecular interaction sites from SDF files, derive candidate axes, and
generate aligned molecular pairs for geometric scanning and optional torsion
optimization. It does not run ORCA calculations.

## Features

- Identify named functional groups and interaction sites using RDKit.
- Derive configurable axes and export selection CSV, alignment CSV, or annotated SDF.
- Align rigid and flexible molecule pairs using Cartesian or spherical scans.
- Optionally optimize rotation around the flexible axis using distance or MMFF94 scores.

## Getting Started

PubChem SDF downloading is no longer included in this repository and will be
provided by a separate GitHub project. Supply SDF files to the selector as
described below.

The selector and shared logging are provided by the separate
[`molecular-selector`](https://github.com/eisenmann-bastian/molecular-selector)
and [`molecular-utilities`](https://github.com/eisenmann-bastian/molecular-utilities)
packages. Create and activate a Python environment, then install this project:

```bash
python -m pip install -e .
```

Both dependencies are fetched from GitHub as specified in
`pyproject.toml`; installation requires Git and network access. Shared logging
lives in `molecular_utilities.logging`; CLI log messages use the
`molecular_aligner` namespace and standard error. The old `core.logging`
compatibility module has been removed.

The selector dependency requires a published selector revision whose
distribution name is `molecular-selector`.

The installed commands are `molecule-axes`, `molecular-aligner`, and
`molecular-aligner-select`. The previous `molecule-aligner` and
`molecule-aligner-select` commands remain available as compatibility aliases.

### Select functional groups and interaction sites

Install the selector dependency and process a list of SDF paths:

```bash
molecule-selector --csv examples/reference_molecules.csv -o matches.csv
molecule-axes matches.csv -o axes.csv
# Compatibility combined workflow:
python3 aligner/selector.py --csv examples/reference_molecules.csv --format alignment -o axes.csv
python3 aligner/selector.py --csv examples/reference_molecules.csv --format sdf -o axes.sdf
```

The standalone selector produces selection CSV. The compatibility command
`aligner/selector.py` additionally supports alignment CSV or annotated SDF
without an intermediate file. To list selectors, use
`molecule-selector --list-selectors`.

Selection definitions and matching live in the separate selector repository.
Axis calculation, configuration and axis output live under `aligner/`;
`aligner/selector.py` provides the combined entry point.

The batch CSV has one column with a header and one SDF path per row:

```csv
sdf
molecules/sulfide.sdf
molecules/aromatic.sdf
```

Relative SDF paths are resolved relative to the CSV file. Use `--delimiter ';'`
to read older multi-column input with a different delimiter. A batch error
returns a nonzero exit without partial output. The `--csv` option cannot be
combined with a positional SDF, `--name`, or `--smarts`.

For compatibility, multi-column batch CSVs with selector headers and marked
cells still select only the marked selectors. Single-SDF `--name` and
`--smarts` queries remain available and can also use any output format.

Results are CSV on stdout by default. `--output`/`-o` writes a new file
instead and refuses to overwrite an existing file. Output always uses commas,
regardless of the batch input `--delimiter`. Each match is one row; a
molecule with no matches has a row with empty match fields. The columns are
`csv_line` (batch input line; empty for single-SDF input), `sdf_file`,
`record_index`, `molecule_name`, `selector`, `smarts`, `atom_index_base`,
`match_index`, `atom_indices`, and `site_atom_indices`. Output `sdf_file`
paths are relative to the command's working directory, even when `--output`
writes to another directory. Both record and match
indices are zero-based. Atom-index cells contain JSON arrays of zero-based
indices, for example `"[0,1,2]"`, from the original SDF record; CSV quoting
preserves arrays containing commas. Explicit hydrogens are retained.
`atom_indices` identifies the whole SMARTS match and `site_atom_indices`
identifies the highlighted interaction site.
For example, `thioether` matches C–S–C but highlights only S. Custom SMARTS
highlights all matched atoms. Named matches sharing the same highlighted sites
are reported once. Names are case-insensitive; spaces, hyphens, and underscores
are interchangeable.

### Derive alignment axes

With `--format alignment`, the program reads the original SDF record and
highlighted atom indices for each match, retaining explicit hydrogen atoms.
It writes CSV to stdout or, with `-o`, to a new file without overwriting.
Each output row retains the selector
columns and adds `axis_type`, `axis_name`, `center_atom_index`, `center`, and `axis`.
The latter two are JSON arrays of three coordinates in the SDF coordinate
frame; `axis` is a unit vector pointing outward from `center`. Atom indices
are zero-based. Ring axes use the ring centroid and leave
`center_atom_index` empty. A match can produce multiple rows, including
duplicate geometries from different selectors.

#### Configure axes

The bundled `aligner/axes.config` is a YAML file loaded automatically for
alignment CSV and SDF output. Every named definition has an `enabled`
boolean and may have a `directions` mapping to enable or disable individual
numbered directions:

```yaml
axes:
  aromatic_ring_normal:
    enabled: true
    directions: {1: true, 2: false}
```

This is an excerpt: a custom config must contain **all** axis definitions
from the bundled file. Copy and edit it, then use
`--axes-config path/to/axes.config` with either the unified selector or the
legacy alignment command. Configurations with unknown or missing names,
duplicate keys, non-boolean switches, or invalid direction numbers are
errors. Setting `enabled: false` disables the entire definition; unspecified
direction numbers inherit `enabled`. Selection CSV is unaffected.

The output `axis_name` is the definition name followed by a 1-based direction
number, for example `aromatic_ring_normal_1` and `aromatic_ring_normal_2`.
Numbers restart for each site and are global across molecules and selectors:
disabling direction 2 applies to every instance of that definition. Direction
order follows the SDF atom/bond ordering and RDKit ring ordering, not a
universal spatial orientation. Renumbering atoms or changing input ordering
can swap which physical direction is numbered 1 or 2. Axis calculation and
geometry validation occur before filtering; disabling a direction does not
repair invalid or degenerate input geometry.

Alignment CSV rows also include `group_type`, `group_atom_indices` and
`group_charge`. For a carboxylate oxygen these identify the delocalized
carboxylate group, its zero-based carbon and oxygen indices, and its total
formal charge of −1; they are empty for other sites. The charge is a
group-level annotation, not a calculated partial charge.

Supported structural geometries are:

- Carboxylate: both oxygens use the same geometric rules regardless of the
  input resonance form. Each has an outward C–O extension and two in-plane
  candidate directions at 120 degrees to the O–C bond. The plane is defined
  by the carbon and both oxygens.
- Singly bonded O, S, Se or Te: an outward bond extension centered on the
  chalcogen for each bond to a non-hydrogen atom. For example, a thioether
  has one axis for each of its two C–S bonds. Chalcogen–H bonds only produce
  such a chalcogen-centered axis for thiols as described below.
- Thioethers (C–S–C): in addition to the two bond extensions, two opposite
  unit normals to the C–S–C plane and the outward in-plane C–S–C angle
  bisector, all centered on sulfur. These are geometric candidate directions
  for lone-pair interactions, not calculated electron-pair positions.
- Thiols (C–S–H): the same two plane normals and outward angle bisector,
  centered on sulfur. Both C–S and H–S have sulfur-centered extensions
  opposite their respective bonded atom. The S–H donor axis remains centered
  on hydrogen. The H–S extension and plane axes require an explicitly
  positioned hydrogen and are skipped with a warning otherwise.
- Double-bonded O, S or Se: the outward bond extension; when the other atom is
  carbon with two explicitly positioned substituents, also two in-plane
  lone-pair directions at 120 degrees to the bond from the chalcogen to carbon.
- Singly bonded halogens and hydrogen: the outward bond extension, centered on
  the halogen or hydrogen. For selected N/O/S/Se sites, explicit bonded
  hydrogens also produce hydrogen axes.
- Neutral nitrile nitrogen: the outward C–N extension. Neutral two-coordinate
  aromatic or unsaturated ring nitrogen without hydrogen: the outward
  in-plane angle bisector. Neutral C=N–H imine nitrogen: the outward
  in-plane C–N–H angle bisector.
- Named aromatic rings: both directions of the ring-plane normal at the
  centroid, including separate axes for fused rings.

These are geometric approximations, not calculated electron densities or
evidence of an interaction. A site is classified from its actual structure,
not inferred from the selector name, so axes do not by themselves establish
whether the selected site is an interaction donor or acceptor. Unsupported
sites are skipped with a warning. A hydrogen direction requires an explicitly
positioned hydrogen; missing hydrogens are not generated or guessed.
SDF records used for matches must have 3D coordinates; invalid CSV, missing
SDFs, malformed records, invalid indices, and degenerate geometry are errors.
Input errors produce no output file.

To visualize the axes as pseudo-atoms, use SDF output instead:

```bash
python3 aligner/selector.py --csv examples/reference_molecules.csv --format sdf -o axes.sdf
```

The input is the one-column SDF-path CSV. The output contains one SDF molecule per
referenced source record, including records without supported axes (which are
written unchanged with a warning). Each distinct axis adds an unbonded RDKit
dummy atom (`*`, atomic number 0), positioned at `center + axis`, one
ångström from the center. For aromatic rings, the center is the ring centroid.
Identical center and direction coordinates within 1e-8 are deduplicated.
Original atoms, bonds, coordinates and properties are retained. The SDF
properties `quick_qm_source` and `quick_qm_alignment_axes` contain JSON source
and axis metadata; each axis lists the zero-based pseudo-atom index and all
selector rows that produced it, their axis names, and carboxylate group metadata where
applicable. Original bond orders and atom formal charges are not altered.
Existing source properties with these names
are rejected rather than overwritten. SDF output requires `-o` and does not
overwrite an existing file.

Use `molecule-axes matches.csv -o axes.csv` to process a selection-result CSV.
The `aligner/alignment.py` command is also available. Input SDF paths
are resolved relative to the command's working directory, as before.

Named patterns include common functional groups and hydrogen-bond donor and
acceptor sites. The latter use RDKit's Lipinski SMARTS definitions. These are
structural heuristics, not a prediction that an interaction actually occurs.
Named functional-group SMARTS are intentionally limited definitions, not an
exhaustive chemical classification. Invalid input or an unreadable SDF record
causes a nonzero exit; no matches is a successful result.

The `carboxylate` selector identifies deprotonated carboxylic acids
(`R–COO⁻`, including formate) and highlights both oxygen atoms as one group.
It recognizes either assignment of the single and double C–O bonds. It does
not match neutral carboxylic acids, esters, carbonates or other resonance
groups. The `carbonyl` and `negative charge` selectors retain their literal
bond-order and formal-charge semantics; use `carboxylate` to select both
oxygens in one match.

Additional interaction-site selectors include:

- `halogen bond donor`: carbon-bound, neutral Cl, Br or I; the halogen is
  highlighted. Fluorine and free halide ions are excluded.
- `positive charge` and `negative charge`: atoms with formal charges from +1
  to +4 and -1 to -4 respectively, including charged metals. Protonation and
  charge must already be represented in the input.
- `aromatic ring` (also `pi-pi` or `cation-pi ring`): each RDKit-perceived
  aromatic ring is returned as a set of its atoms, including heteroaromatic
  rings. Fused rings can share atoms. The CSV `smarts` field is empty
  because this selector uses RDKit ring perception rather than a SMARTS.
- `chalcogen bond donor`: neutral divalent S, Se or Te bound to two carbon
  atoms, highlighting the chalcogen.
- `metal center`: positively charged Li, Na, K, Mg, Ca, Mn, Fe, Co, Ni, Cu,
  Zn, Cd, Hg, Pt or Pd; `metal coordination ligand` selects potential
  coordinating atoms.

`halogen bond acceptor`, `chalcogen bond acceptor`, and
`metal coordination ligand` reuse RDKit's Lipinski hydrogen-bond acceptor
pattern as an approximate candidate filter. These definitions are deliberately
incomplete; for example, they do not capture every metal ligand or chalcogen
donor. They identify potential sites only and do not test geometry, ring
orientation, protonation equilibria, interaction strength or actual bonding.

### Align pairs of molecular axes

Install the packages as described above, then use an alignment CSV generated by
`molecule-axes` or `aligner/selector.py --format alignment`:

```bash
python3 aligner/aligner.py axes.csv \
  --rigid '^241$' '^aromatic ring$' '^aromatic_normal$' \
  --flexible '^878$' '^thiol$' '^hydrogen_sigma$'
```

The CLI triplets contain Python regular expressions for `molecule_name`,
`selector` and `axis_type`. Matching uses regex search; use `^` and `$` for
full-field matching. Each matching rigid row is paired with each matching
flexible row, including duplicate or self pairs. Repeat either option to
specify alternative triplets; a row matching several alternatives is selected
once. For either side, include and exclude lists are independent:

```yaml
rigid:
  include:
    - {molecule_name: "^241$"}
    - {molecule_name: "^175$"}
  exclude:
    - {molecule_name: "^175$", selector: "carbonyl"}
flexible:
  exclude:
    - {molecule_name: "^241$"}
```

Each YAML filter must specify at least one of `molecule_name`, `selector` and
`axis_type`; omitted fields impose no restriction. All fields specified in a
filter must match. Each list can contain any number of filters; alternatives
within each list are ORed. A row is retained if it matches at least one include
filter (or the include list is absent or empty) and matches no exclude filter.
With neither list, all rows are selected. After exclusions, an empty side is an
error. Legacy direct mappings or lists remain accepted as include filters.

On the CLI, repeat `--rigid` or `--flexible` for includes, and repeat
`--rigid-exclude` or `--flexible-exclude` for exclusions, each followed by
three regexes (use `.*` for an unrestricted field). CLI include options replace
the YAML include list for that side; CLI exclude options independently replace
the YAML exclude list.

The rigid molecule stays fixed. The flexible molecule is rotated so its axis
points opposite the rigid axis, then translated so its axis center lies along
the rigid axis. By default, the separation is the sum of the two center atoms'
RDKit van der Waals radii. An aromatic-ring centroid contributes zero, so a
ring–atom pair uses the atom's radius and two ring centroids coincide.
`--distance 3.5` overrides this separation in ångströms, including for rings.
The rotation is deterministic but rotation around the aligned axis is not
optimized. No clash detection, relaxation or interaction-energy calculation
is performed; these are starting geometries, not predicted complexes.

For multiple distances per axis pair, select `distance-scan`:

```bash
python3 aligner/aligner.py axes.csv \
  --rigid '^241$' '^aromatic ring$' '^aromatic_normal$' \
  --flexible '^878$' '^thiol$' '^hydrogen_sigma$' \
  --alignment-mode distance-scan --scan-start 2 --scan-end 4 \
  --scan-step 0.1 --no-include-start --include-end
```

The range is in ångströms, with non-negative start and end greater than start.
Both boundaries are included by default; `--no-include-start` and
`--no-include-end` exclude them independently. Specify exactly one of
`--scan-step` or `--scan-count`. A step must divide the full boundary span
exactly, even when boundaries are excluded, and must produce at least one
pose. With a count, positions form an evenly spaced grid: excluded boundaries
are grid points, not poses. The count specifies the positive number of points
strictly between the boundaries; each included boundary adds one more pose.
For example, count 1 gives the midpoint, plus either or both boundaries if
included. `--distance` is only valid in `default` mode. Invalid scan
settings abort before output is created.

#### Lateral X and Z offsets

The original rigid axis defines local Y. Optionally add independent X and Z
offsets; by default the flexible molecule's orientation remains unchanged:

```bash
python3 aligner/aligner.py axes.csv \
  --rigid '^241$' '^aromatic ring$' '^aromatic_normal$' \
  --flexible '^878$' '^thiol$' '^hydrogen_sigma$' \
  --distance 3 --x-alignment-mode distance-scan \
  --x-scan-start -1 --x-scan-end 1 --x-scan-step 0.5 --z-distance -0.5
```

Each axis supports a fixed offset or an independent scan with the same
endpoint and step/count rules as Y. Use `--x-` or `--z-` prefixes for
`alignment-mode`, `distance`, `scan-start`, `scan-end`, `include-start`,
`include-end`, `scan-step`, and `scan-count`, including `--no-x-include-start`
and similar exclusions. X and Z accept negative values and default to zero
when omitted. Y remains non-negative. The Cartesian product of the Y, X, and Z
positions is generated for each axis pair (Z varies fastest).

Set `flexible_orientation: toward-center` in YAML or pass
`--flexible-orientation toward-center` to rotate the flexible molecule around
its placed center so its axis points toward the rigid center. This direction
is the line between the centers, not necessarily parallel to the rigid axis
when X or Z is nonzero. If the centers coincide, the existing opposite-to-Y
orientation is kept. The default, `stable`, always retains that orientation.
CLI overrides YAML. The total rotation and translation in the CSV reflect the
selected orientation, and the mode is recorded in both CSV and SDF.

In YAML, group the axis settings under `grid_based_scan`. A range automatically
selects scan mode for that axis; use `distance` instead for a fixed offset:

```yaml
grid_based_scan:
  x:
    scan_start: -1
    scan_end: 1
    scan_count: 1
  y:
    scan_start: 2
    scan_end: 4
    scan_count: 3
    include_start: false
    include_end: false
  z:
    scan_start: -0.5
    scan_end: 0.5
    scan_step: 0.5
reference_axis:
  priority: [A, B, C]
  enabled: {A: true, B: true, C: true}
```

For example, replace the `z` scan with `distance: 0.5` for a fixed Z offset.
An omitted Y axis uses the van der Waals distance; omitted X and Z use zero.
Y must remain non-negative. CLI values override individual YAML values;
a CLI step/count replaces the alternative sampling choice. A CLI fixed
distance replaces that axis's YAML scan, and CLI scan settings replace a
YAML fixed distance. The former `distance_scan` key, top-level scan keys and
top-level `x`/`z` mappings are no longer accepted; move them into
`grid_based_scan`.

An X or Z configuration requires a reference even for a zero offset. The
reference point determines X by projecting the center-to-reference vector
onto the plane perpendicular to Y and normalizing it; Z = X × Y. The first
enabled, usable method in `priority` wins:

* **A:** lowest RDKit atom index among all neighbors of the center atom's
  neighbors, excluding the center itself.
* **B:** highest such atom index.
* **C:** the global coordinate origin `(0, 0, 0)`, even for atom centers.

All methods are enabled by default; omit methods from `priority` or set their
`enabled` value to `false` to disable them. A/B consider all atoms, including
explicit hydrogens, and the IDs are zero-based SDF/RDKit indices. A method
is unusable if its vector is zero, parallel to Y, or non-finite. If no enabled
method is usable, the run fails before writing output. Method C depends on
the input coordinate origin and is not translation-invariant. With no X/Z
configuration no reference is needed, preserving Y-only behavior.

#### Spherical coordinate scan

Use `polar_coordinate_based_scan` instead of `grid_based_scan` to specify a
radius in ångströms, phi and theta. Each component accepts `distance` for a
fixed value or `scan_start`, `scan_end`, and either `scan_step` or
`scan_count`, plus optional endpoint inclusion flags. Radius defaults to the
van der Waals separation; phi and theta default to zero. Values are combined
in radius, phi, theta order, with theta varying fastest.

```yaml
polar_coordinate_based_scan:
  radius: {distance: 3.0}
  phi: {scan_start: 0, scan_end: 360, scan_step: 90}
  phi_unit: degrees
  theta: {scan_start: 0, scan_end: 180, scan_step: 90}
  theta_unit: degrees
  merge_tolerance: 0.05
```

`phi_unit` and `theta_unit` independently accept `degrees` (default) or
`radians`. Theta ranges from 0 to 180 degrees (π radians) measured from +Y;
phi ranges from 0 to 360 degrees (2π radians), measured from +X toward +Z.
For local unit vectors X, Y, Z, displacement is
`r * (sin(theta)*cos(phi)*X + cos(theta)*Y + sin(theta)*sin(phi)*Z)`.
The polar scan requires a usable `reference_axis` to define X and Z.

Positions within `merge_tolerance` ångströms (default 0.05, strictly
positive) of an already retained position are merged per axis pair, with a
warning. The first point in scan order is retained, including its original
radius and angles. This removes repeated pole, zero-radius and 0/360-degree
positions; nearby distinct positions can also be merged if they fall within
the tolerance. No merging occurs within a Cartesian grid.

For CLI use, select `--coordinate-system polar` or supply any `--polar-*`
option. Fixed values use `--polar-radius`, `--polar-phi`, `--polar-theta`;
each supports `--polar-<component>-scan-start`, `-scan-end`, `-scan-step` or
`-scan-count`, and `--[no-]polar-<component>-include-start/end`.
`--polar-phi-unit`, `--polar-theta-unit`, and `--polar-merge-tolerance`
override their YAML values. CLI-only grid and polar options cannot be combined.
The CSV and SDF include the retained spherical coordinates and units;
`distance_angstrom`, X and Z offsets remain local Cartesian components.
`stable` and `toward-center` orientations are both supported.

#### Multiple scans in one run

Use an ordered `scans` list to run any number of grid and polar scans for
each axis pair. Each entry contains exactly one scan section:

```yaml
scan_merge_tolerance: 0.05
scans:
  - grid_based_scan:
      y: {scan_start: 2, scan_end: 3, scan_step: 0.5}
  - polar_coordinate_based_scan:
      radius: {distance: 3}
      phi: {scan_start: 0, scan_end: 360, scan_step: 90}
      theta: {distance: 90}
  - grid_based_scan:
      y: {distance: 4}
```

Alternatively, `grid_based_scan` and `polar_coordinate_based_scan` can each
be a mapping or a nonempty list of mappings, and both can appear at the top
level. In this form, sections run in YAML key order and entries run in list
order. Do not mix these top-level sections with `scans`.

All scans use the same filters, reference axis, orientation, and torsion
settings. Poses are appended to the existing per-pair or single SDF and CSV
outputs. Their 1-based `scan_index` is recorded in the CSV and as
`quick_qm_scan_index` in the SDF; pose numbering remains continuous per axis
pair. For each axis pair, a position within `scan_merge_tolerance` of a
retained position from an earlier scan is skipped; the earlier scan wins.
The tolerance is in ångströms, must be positive, and defaults to 0.05.
Set it in YAML or override it with `--scan-merge-tolerance`. Within an
individual scan the existing rules still apply: grid positions are not merged,
and polar positions use that scan's own `merge_tolerance`.

Grid CLI scan options override every grid scan, and `--polar-*` options
override every polar scan. An override for a coordinate system absent from
the configuration is an error. When both systems are configured, omit
`--coordinate-system`; specifying a single system conflicts with the mixed
configuration. With no configured scans, CLI-only behavior is unchanged.

#### Optimize rotation around the flexible axis

Optional `torsion_optimization` runs after placement and orientation, once
per retained pose in either scan system:

```yaml
torsion_optimization:
  enabled: true
  score: minimum_distance     # average_distance or mmff94_interaction
  objective: maximize         # or minimize
  include_vdw: true
  exclude_flexible_site: true
  exclude_rigid_site: false
  exclude_flexible_atoms: []  # additional 0-based SDF atom indices
  exclude_rigid_atoms: []
  steps_per_iteration: 12
  iterations: 3
```

By default optimization is disabled; when enabled it maximizes minimum
distance. Average and minimum scores use every remaining rigid–flexible atom
pair, including hydrogens, in ångströms. With `include_vdw`, each pair's
score is its center distance minus the two RDKit van der Waals radii, so
negative values indicate overlap. The flexible selector's `site_atom_indices`
are excluded by default; rigid site exclusion is optional. If site indices are
missing, the center atom is used; a ring center without site indices is an
error. Additional explicit exclusions are combined with site exclusions.
Invalid indices or no remaining atoms cause an error.

`mmff94_interaction` uses RDKit MMFF94 complex energy with interfragment
interactions enabled, minus isolated rigid and flexible molecule energies,
in kcal/mol. It scores **all atoms**, ignoring distance-score exclusions;
`include_vdw` is not applicable. Missing MMFF parameters are an error.
Minimize this energy to seek a low-energy orientation. Neither this model nor
distance scores include solvent or relaxation.

The deterministic search samples the full circle (including the original
orientation), then refines around the best angle. The defaults are 12 samples
per pass and three passes; this is a discrete approximation, not a guarantee
of a global optimum. Only the flexible molecule rotates about its placed axis;
the axis and both centers stay fixed. CSV and SDF record the selected angle,
score, score type and objective; the output transform includes this rotation.
CLI overrides include `--optimize-torsion`, `--torsion-score`,
`--torsion-objective`, `--torsion-steps`, `--torsion-iterations`,
`--[no-]torsion-include-vdw`, `--[no-]torsion-exclude-flexible-site`,
`--[no-]torsion-exclude-rigid-site`, and the `--torsion-exclude-*-atoms`
lists. Use `--no-optimize-torsion` to disable a YAML setting.

By default, `aligned_pairs/` contains one SDF per rigid–flexible axis pair,
named `axis-pair-000001.sdf` and so on. Each file contains one record per
pose for that pair; `aligned_pairs.csv` describes the results. Each SDF record holds both
molecules as disconnected components, rigid first. `--output NAME` changes
the new directory name. Use `--output-format single --output pairs.sdf` for
one multi-record SDF and a corresponding `pairs.csv`. Existing destinations
are rejected by default. `--overwrite-output` permits replacing the CSV and
matching SDF files; other files in an existing output directory, including
SDFs left over from earlier runs, are not removed. Only the current CSV lists
the current results. Output is staged before replacement, and input files
cannot be overwritten. Use `--no-overwrite-output` to override a YAML setting.
The CSV records each source row, transform, center,
direction, distance and output record index; `axis_pair_id` groups poses for
the same axis pair and `pose_index` is one-based within that pair. SDF
properties also include these identifiers, the Y distance, lateral offsets,
and, when used, the reference method and X/Z unit vectors. CSV records the
reference point and atom index as well. `distance_angstrom` remains the Y
component, not the Euclidean center-to-center distance. Input row fields are JSON
objects. All paths in it are relative to the command's working directory.

Parameters may instead be supplied via `--config examples/configs/aligner.config`;
the input CSV remains a required positional argument. The example YAML holds
`rigid` and `flexible` mappings with include and exclude filter lists, plus
optional `output`, `output_format`, `overwrite_output`, `grid_based_scan`,
`polar_coordinate_based_scan`, `scans`, `scan_merge_tolerance`, `torsion_optimization`, and
`reference_axis`, plus `flexible_orientation`. Each axis accepts `distance`, or `scan_start`, `scan_end`,
`include_start`, `include_end` and either `scan_step` or `scan_count`.
Do not mix a fixed `distance` with scan settings on the same axis.
Relative source SDF paths in the input CSV are
resolved against the current working directory, matching selector output.
Run the aligner from the same directory used to generate that CSV.

The `examples/configs/` directory contains runnable examples:
`y_step.yaml` and `y_count.yaml` scan Y by step or interior-point count;
`x_step.yaml` and `z_count.yaml` scan lateral offsets; `xyz_grid.yaml`
combines scans into a Cartesian grid; `polar_grid.yaml` scans spherical
coordinates; `multiple_scans.yaml` combines grid and polar scans;
`torsion_optimization.yaml` optimizes torsions. They use the supplied
`examples/test_out_aligned.csv` and its relative SDF paths. Run an example
from the project root (or use the script's absolute path from another directory):

```bash
./examples/run-example.sh xyz_grid.yaml
./examples/run-example.sh polar_grid.yaml
```

The configuration name defaults to `xyz_grid.yaml`. The script writes results
to a new temporary directory and prints its path. To choose a directory, pass
it as the second argument; existing output is never overwritten:

```bash
./examples/run-example.sh torsion_optimization.yaml /tmp/torsion-example
```

Set `PYTHON` to the desired interpreter (use an absolute path when invoking
the example script). The example runs from `examples/` so the CSV's relative
SDF paths resolve correctly.

## Development and tests

Install the three local packages as described above, then run:

```sh
./test.sh
# Alternatively:
./examples/run-tests.sh
```

Set `PYTHON` to the desired interpreter, for example
`PYTHON=.venv/bin/python ./test.sh`. The test suite uses temporary output
directories; no network access is required. `aligner/test-selection.sh` and
`aligner/test.sh` run the corresponding subsets.

## Repository structure

- `aligner/`: selector and aligner CLIs, axis input, geometry, scan planning,
  pose generation, torsion optimization, transactional output, and tests.
- `examples/configs/`: runnable scan configurations.
- `examples/run-example.sh` and `examples/run-tests.sh`: example and test runners.
- `examples/reference_molecules/`, `examples/reference_molecules.csv` and
  `examples/test_out_aligned.csv`: sample input structures and axes.
- `test.sh`: unified offline test runner.

The installed CLI command names and `quick_qm` SDF property names are retained
for compatibility. Source paths and Python imports previously under `selection`
or `alignment` now live under `aligner`; import selector classes from
`molecular_selector` and axis classes from `aligner`. Generated output and
Python caches are ignored by Git.

## License

[BSD 3-Clause License](LICENSE).
