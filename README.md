# ReShiki Importer

[![MoleditPy](https://img.shields.io/badge/MoleditPy->=4.0.0-3577F7)](https://github.com/HiroYokoyama/python_molecular_editor)

A [MoleditPy](https://github.com/HiroYokoyama/python_molecular_editor) plugin that opens
[ReShiki](https://github.com/Ameyanagi/ReShiki) drawings (`.rsk`) in the 2D editor. Draw a
structure in ReShiki, open it in MoleditPy, and carry on to 3D conversion and input
generation.

## Use

- **File → Import → Import .rsk (ReShiki Importer)...**
- Drag a `.rsk` file onto the MoleditPy window.
- Or open one from the command line: `moleditpy drawing.rsk`.

The drawing replaces the current 2D structure as one undo step. The status bar says
what was imported and what was left out.

### Paste from ReShiki

Copy in ReShiki (Ctrl+C), then in MoleditPy choose **Edit → Paste from ReShiki** or
press **Ctrl+Alt+V**. The copied structure is added to the current drawing at the centre
of the view, as one undo step; the same rules as opening a file apply to what is kept.

This reads ReShiki's own clipboard data: its native format on Windows and macOS, and its
text copy on Linux. On macOS the shortcut is **Cmd+Option+V**.

## What is imported

| ReShiki | MoleditPy |
|---|---|
| Atoms, charges, radicals | Same |
| Single, double, triple bonds | Same |
| Aromatic bonds (including circle rings) | Kekulé structure, assigned with RDKit |
| Wedge, hollow wedge, bold single bond | Wedge |
| Hash, hashed single bond | Dash |
| Abbreviations (Ph, OMe, …) | Their full structures |
| Dative, partial bonds | Single bond (reported) |
| Quadruple bonds | Triple bond (reported) |

Coordinates are scaled to MoleditPy's bond length and centred in the view.

**Left out, and reported:** reaction arrows, text, shapes, reaction schemes, isotope
labels, centroid atoms of π-complexes and their bonds, and hydrogen-interaction bonds.
The editor holds molecules, not figures.

## Install

Download `reshiki_importer.py` from the latest release and install it through
**Plugins → Plugin Manager**, or copy it into your plugins directory
(`%USERPROFILE%\.moleditpy\plugins` on Windows, `~/.moleditpy/plugins` on macOS/Linux)
and restart MoleditPy.

## Format support

Written against ReShiki document format version 15 (ReShiki 0.8). A newer file still
imports, with a note that its version is newer than the importer knows.

## Development

```bash
python -m pytest tests/ -v
```

With `python_molecular_editor` and `ReShiki` checked out beside this repository, the
suite also checks the plugin against the host API and imports every `.rsk` file in the
ReShiki repository. Without them those tests skip.

## Licence

GPL-3.0 (see [LICENSE](LICENSE)). ReShiki is a separate project by Ameyanagi,
licensed MIT OR Apache-2.0; this plugin reads its file format and contains none of its code.
