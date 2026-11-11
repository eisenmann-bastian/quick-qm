<h1 align="center">
  <br>
  [Project's Logo]
  <br>
  QUICK-QM
  <br>
</h1>

<h4 align="center">This Python based project serves as a wrapper around the ORCA tool suite to quickly retrieve molecular properties and interaction energies of two molecular systems.</h4>

<!--<p align="center">
  [Project's badges]
</p>-->

<p align="center">
  <a href="#about">About</a> •
  <a href="#key-features">Key Features</a> •
  <a href="#getting-started">Getting Started</a> •
  <a href="#faq">FAQ</a> •
  <a href="#roadmap">Roadmap</a> •
  <a href="#support">Support</a> •
  <a href="#license">License</a>
</p>

<!--![screenshot](screenshots/1.jpg)-->

## About

The ORCA tool suite offers various functions to either retrieve properties of a single molecular (A) or retrieving the interaction energy forming an adduct of n molecules (B), which both are commonly used approaches in various QM based studies.

While (A) requires a simple run and adequate calling of sub procedures, (B) requires first the calculation of the individual system energies and then retrieval of the interaction energy by forming an adduct of n molecules and running a calculation afterwards. The final interaction energy is given by the difference of the adduct system energy and the sum of all single system energies.

As (A) is a key component of any single molecular property retrieval and (B) could be used to model any molecule-molecule interaction, this approach can overall be widely used.

This package provides an easy-to-use wrapper, especially targeting two typical conducting simulation approaches: Single-Molecule-QM, to retrieve properties of a single molecule and Two-Molecule-QM, allowing to simulate the interaction of two molecules.

The tool-suite reduces the complexity and speeds up the process by defining the ORCA input files and allows a selector based alignment, quickly creating the alignments needed based on molecular functional groups. 

## Key Features

- Feature 1
- Feature 2
  - Feature 2.1
  - Feature 2.2
- Feature 3

## Getting Started

### Download structures from PubChem

Install the downloader dependencies and run it from this repository:

```bash
python3 -m pip install -r pubchem/requirements.txt
python3 pubchem/downloader.py compounds.csv --output-dir sdf
```

The CSV must have at least one of the case-insensitive headers `cid`, `name`,
`inchi`, `smiles`, or `formula`. Other columns are ignored. For each row the
first populated identifier is used in that order; lower-priority values are not
used as fallbacks if the lookup fails.

```csv
name,inchi,smiles,formula
benzene,,,
,,CCO,
```

`pubchem/sample.csv` contains all five identifiers for each example compound;
the downloader selects the CID in every row. `pubchem/sample-mixed.csv` uses a
different identifier in each row. Its formula row may be rejected as
ambiguous, since formulas can match multiple CIDs.

Each successful match produces `row-<CSV line>-cid-<CID>.sdf`. Files are not
overwritten. CID, name and InChI rows are downloaded directly with PubChemPy;
their responses must contain exactly one SDF record.

Formula and SMILES searches use `--result-selection`:

- `exact` (default): require exactly one CID.
- `first`: download the first CID returned by PubChem.
- `all`: download every distinct CID returned by PubChem.

For SMILES, `--filter-mode` chooses the PubChem search type: `identity`
(default), `substructure`, `superstructure`, or `similarity`. Search modes
other than identity can return different compounds and many results. For
example:

```bash
python3 pubchem/downloader.py compounds.csv --filter-mode similarity --result-selection first
```

`--result-selection` does not apply to direct CID, name or InChI downloads.
Rows with no results, ambiguous results in `exact` mode, or download failures
are reported as errors. Processing continues, including remaining CIDs in
`all` mode, and the command exits nonzero if any row fails.

The default is a 3D SDF. Pass `--record-type 2d` or `--record-type 3d` to
choose explicitly, and `--delimiter ';'` for a semicolon-separated file. A
missing 3D structure is reported as an error; it is not automatically
replaced with a 2D structure. PubChem searches may be slow and `all` may
download many files.

Python tools can share the logging configuration in `core.logging`: call
`configure_logging()` at the application entry point and use
`get_logger("tool_name")` in modules. The downloader requires Python 3.10 or
newer.

### Select functional groups and interaction sites

Install the selector dependency and process a list of SDF paths:

```bash
python3 -m pip install -r selection/requirements.txt
python3 selection/selector.py --csv selection/reference_molecules.csv -o matches.csv
python3 selection/selector.py --csv selection/reference_molecules.csv --format alignment -o axes.csv
python3 selection/selector.py --csv selection/reference_molecules.csv --format sdf -o axes.sdf
```

The single program detects all named selectors for each input molecule and
produces selection CSV (the default), alignment CSV, or annotated SDF. No
intermediate selection file is needed for alignment or SDF output. To list
the selectors, use `python3 selection/selector.py --list-selectors`.

The selector code separates named definitions (`selection/definitions.py`),
SDF matching (`selection/matching.py`), batch input (`selection/batch.py`),
alignment (`selection/alignment.py`), and pipeline orchestration
(`selection/pipeline.py`); `selection/selector.py` is the CLI entry point.

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

The bundled `selection/axes.config` is a YAML file loaded automatically for
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
python3 selection/selector.py --csv selection/reference_molecules.csv --format sdf -o axes.sdf
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

The former `selection/alignment.py` command remains available for processing
an existing selection-result CSV. Its input SDF paths are resolved relative
to the command's working directory, as before.

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

This is an example of how you may give instructions on setting up your project locally. To get a local copy up and running follow these simple example steps.

### Prerequisites

This is an example of how to list things you need to use the software and how to install them.

- npm

### Installing and Running

Below is an example of how you can instruct your audience on installing and setting up your app. This template doesn't rely on any external dependencies or services.

```bash
# Clone this repository
$ git clone https://github.com/[your-github-user]/[your-repository-name]

# Go into the repository
$ cd [your-repository-name]
```

## FAQ

### Is it any good?

[yes.](https://news.ycombinator.com/item?id=3067434)

## Definitions

- This are constraints you give to your project
- Here's another constraint: use this README!

## Roadmap

- [x] Add this README.
- [ ] Finish the unicorn project.

## Emailware

<!--[Project's name] is an [emailware](https://en.wiktionary.org/wiki/emailware). Meaning, if you liked using this app or it has helped you in any way, I'd like you send me an email at <your-email@email.com> about anything you'd want to say about this software. I'd really appreciate it!-->

## Contributing

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the project
2. Clone it open the repository in command line
3. Create your feature branch (`git checkout -b feature/amazing-feature`)
4. Commit your changes (`git commit -m 'Add some amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a pull request from your feature branch from your repository into this repository main branch, and provide a description of your changes

## Support

You can also support us by:

<p align="left">
  <a href="https://buymeacoffee.com/bastian.eisenmann" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/purple_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a> &nbsp &nbsp
  <!--<a href="https://www.patreon.com">
    <img src="https://c5.patreon.com/external/logo/become_a_patron_button@2x.png" width="160">
  </a>-->
</p>

## License (BSD 3-Clause License)

BSD 3-Clause License

Copyright (c) 2026, Eisenmann, Bastian

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

<!--
## Acknowledgments

Inspiration, code snippets, etc.

- [Markdownify's README](https://github.com/amitmerchant1990/electron-markdownify#readme)

## You may also like...

List of apps or libs that do similar stuff as your project.

- [Best-README-Template](https://github.com/othneildrew/Best-README-Template)
- [Simple README.md template](https://gist.github.com/DomPizzie/7a5ff55ffa9081f2de27c315f5018afc)

---
-->

> [hidek.xyz](https://hidek.xyz) &nbsp;&middot;&nbsp;
> GitHub [@abehidek](https://github.com/abehidek) &nbsp;&middot;&nbsp;
> Twitter [@guilhermehabe](https://twitter.com/guilhermehabe)
