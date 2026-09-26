"""Open ReShiki drawings (.rsk) in the MoleditPy 2D editor, or paste them.

A .rsk file is ReShiki's document as JSON. The molecule is imported: atoms,
bonds, charges, radicals and wedge/hash stereo. Reaction arrows, text, shapes
and other figure content have no place in the editor and are skipped.

Edit > Paste from ReShiki (Ctrl+Alt+V) adds what ReShiki last copied to the
current drawing.
"""

import json
import logging
import math
import os
import sys
from collections import Counter

PLUGIN_NAME = "ReShiki Importer"
PLUGIN_VERSION = "0.2.1"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Open ReShiki drawings (.rsk) in the 2D editor, or paste them from ReShiki "
    "with Ctrl+Alt+V, with their atoms, bonds, charges, radicals and wedge/hash stereo."
)
PLUGIN_CATEGORY = "Import"
PLUGIN_TAGS = ["Import", "ReShiki"]
PLUGIN_DEPENDENCIES = []
PLUGIN_OPTIONAL_DEPENDENCIES = []
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"

EXTENSION = ".rsk"
# The newest ReShiki document version this importer was written against.
KNOWN_VERSION = 15
# MoleditPy's standard 2D bond length, in scene units.
BOND_LENGTH = 75.0
# ReShiki's default bond length, used when a drawing has no bonds to measure.
RESHIKI_BOND_LENGTH = 42.0

PASTE_SHORTCUT = "Ctrl+Alt+V"
# ReShiki's own clipboard format: the registered Windows format name and the
# macOS pasteboard type.
NATIVE_TYPE = "dev.reshiki.drawing"
# The same cap ReShiki puts on clipboard data.
CLIPBOARD_LIMIT = 64 * 1024 * 1024
# The native format as Qt names it.
_CLIPBOARD_FORMATS = (
    f'application/x-qt-windows-mime;value="{NATIVE_TYPE}"',
    NATIVE_TYPE,
)
# Where ReShiki has no native clipboard (Linux) it copies the document as text.
_TEXT_PREFIXES = ("RESHIKI_DRAWING_V1\n", "MORUNO_DRAWING_V1\n")

WEDGE = 1
DASH = 2
_WEDGE_DISPLAYS = {"wedge", "hollow_wedge", "bold"}
_DASH_DISPLAYS = {"hash", "hashed"}
AROMATIC = "aromatic"
# ReShiki bond orders MoleditPy cannot draw, and what they become.
_APPROXIMATED = {
    5: (1, "dative bonds imported as single"),
    6: (3, "quadruple bonds imported as triple"),
    7: (1, "partial bonds imported as single"),
}

_context = None


class RskError(ValueError):
    """The file is not a ReShiki drawing this importer can read."""


def parse_rsk(document):
    """Turn a decoded .rsk document into atoms and bonds for the editor.

    Returns ``(atoms, bonds, notes)``. Atoms carry ``id, symbol, x, y, charge,
    radical, explicit_h``; bonds carry ``a, b, order, stereo`` where ``order``
    may still be ``AROMATIC``. ``notes`` counts what was skipped or changed.
    """
    if not isinstance(document, dict) or not isinstance(document.get("atoms"), list):
        raise RskError("This is not a ReShiki drawing: it has no atom list.")
    version = document.get("version")
    if not isinstance(version, int):
        raise RskError("This is not a ReShiki drawing: it has no document version.")

    notes = Counter()
    if version > KNOWN_VERSION:
        notes[f"file is ReShiki format {version}; importer knows up to {KNOWN_VERSION}"] += 1

    atoms = []
    kept = set()
    for raw in document["atoms"]:
        try:
            atom_id = int(raw["id"])
            symbol = str(raw["element"])
            x = float(raw["position"]["x"])
            y = float(raw["position"]["y"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RskError(f"An atom in this drawing is incomplete: {exc}") from exc
        if symbol in ("", "*") or raw.get("centroid"):
            notes["centroid and placeholder atoms skipped"] += 1
            continue
        if raw.get("isotope"):
            notes["isotope labels dropped"] += 1
        atoms.append(
            {
                "id": atom_id,
                "symbol": symbol,
                "x": x,
                "y": y,
                "charge": int(raw.get("charge", 0) or 0),
                "radical": int(raw.get("radical_electrons", 0) or 0),
                "explicit_h": int(raw.get("explicit_h", 0) or 0),
            }
        )
        kept.add(atom_id)

    bonds = []
    for raw in document.get("bonds", []):
        try:
            a, b, order = int(raw["a"]), int(raw["b"]), int(raw["order"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RskError(f"A bond in this drawing is incomplete: {exc}") from exc
        if a not in kept or b not in kept or a == b:
            notes["bonds to skipped atoms dropped"] += 1
            continue
        if order == 0:
            notes["hydrogen-interaction bonds dropped"] += 1
            continue
        if order == 4:
            order = AROMATIC
        elif order in _APPROXIMATED:
            order, note = _APPROXIMATED[order]
            notes[note] += 1
        elif order not in (1, 2, 3):
            notes["unknown bond orders imported as single"] += 1
            order = 1
        stereo = 0
        if order == 1 and not raw.get("projection"):
            display = raw.get("display", "plain")
            if display in _WEDGE_DISPLAYS:
                stereo = WEDGE
            elif display in _DASH_DISPLAYS:
                stereo = DASH
        bonds.append({"a": a, "b": b, "order": order, "stereo": stereo})

    for key, label in (
        ("arrows", "arrows"),
        ("annotations", "text annotations"),
        ("graphics", "shapes"),
        ("reactions", "reaction schemes"),
    ):
        count = len(document.get(key) or [])
        if count:
            notes[f"{label} skipped"] += count
    return atoms, bonds, notes


def kekulize(atoms, bonds, notes):
    """Give aromatic bonds single/double orders, in place.

    Uses RDKit, which MoleditPy ships. When the ring cannot be kekulized the
    aromatic bonds become single so the drawing still loads.
    """
    aromatic = [bond for bond in bonds if bond["order"] == AROMATIC]
    if not aromatic:
        return
    try:
        from rdkit import Chem

        mol = Chem.RWMol()
        index = {}
        for atom in atoms:
            rd_atom = Chem.Atom(atom["symbol"])
            rd_atom.SetFormalCharge(atom["charge"])
            rd_atom.SetNumRadicalElectrons(atom["radical"])
            rd_atom.SetNumExplicitHs(atom["explicit_h"])
            index[atom["id"]] = mol.AddAtom(rd_atom)
        for bond in bonds:
            if bond["order"] == AROMATIC:
                bond_type = Chem.BondType.AROMATIC
                mol.GetAtomWithIdx(index[bond["a"]]).SetIsAromatic(True)
                mol.GetAtomWithIdx(index[bond["b"]]).SetIsAromatic(True)
            else:
                bond_type = {
                    1: Chem.BondType.SINGLE,
                    2: Chem.BondType.DOUBLE,
                    3: Chem.BondType.TRIPLE,
                }[bond["order"]]
            mol.AddBond(index[bond["a"]], index[bond["b"]], bond_type)
            if bond_type == Chem.BondType.AROMATIC:
                mol.GetBondBetweenAtoms(index[bond["a"]], index[bond["b"]]).SetIsAromatic(True)
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        for bond in aromatic:
            rd_bond = mol.GetBondBetweenAtoms(index[bond["a"]], index[bond["b"]])
            bond["order"] = 2 if rd_bond.GetBondType() == Chem.BondType.DOUBLE else 1
    except Exception as exc:  # noqa: BLE001 - any RDKit failure falls back to single bonds
        logging.warning("%s: could not kekulize aromatic bonds: %s", PLUGIN_NAME, exc)
        for bond in aromatic:
            bond["order"] = 1
        notes["aromatic bonds that could not be kekulized imported as single"] += len(aromatic)


def layout(atoms, bonds):
    """Scene positions for the atoms: MoleditPy bond length, centred on the origin."""
    if not atoms:
        return {}
    position = {atom["id"]: (atom["x"], atom["y"]) for atom in atoms}
    lengths = sorted(
        math.dist(position[bond["a"]], position[bond["b"]])
        for bond in bonds
        if math.dist(position[bond["a"]], position[bond["b"]]) > 1e-6
    )
    typical = lengths[len(lengths) // 2] if lengths else RESHIKI_BOND_LENGTH
    scale = BOND_LENGTH / typical
    cx = sum(x for x, _ in position.values()) / len(position)
    cy = sum(y for _, y in position.values()) / len(position)
    # Both programs draw with y pointing down, so no flip is needed.
    return {key: ((x - cx) * scale, (y - cy) * scale) for key, (x, y) in position.items()}


def read_rsk(path):
    """Read and convert a .rsk file. Returns ``(atoms, bonds, positions, notes)``."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RskError(f"Could not read {os.path.basename(path)}: {exc}") from exc
    return convert(document)


def convert(document):
    """Convert a decoded document. Returns ``(atoms, bonds, positions, notes)``."""
    atoms, bonds, notes = parse_rsk(document)
    if not atoms:
        raise RskError("This ReShiki drawing contains no atoms to import.")
    kekulize(atoms, bonds, notes)
    return atoms, bonds, layout(atoms, bonds), notes


def summary(atoms, bonds, notes, verb="Imported"):
    text = f"{verb} {len(atoms)} atoms and {len(bonds)} bonds from ReShiki."
    if notes:
        text += " " + "; ".join(f"{count} {label}" for label, count in sorted(notes.items())) + "."
    return text


def macos_pasteboard_data(kind=NATIVE_TYPE):
    """Bytes of ``kind`` on the macOS general pasteboard, or None.

    Qt only lists pasteboard types it has a converter for, so ReShiki's own
    type never reaches ``QMimeData`` on macOS; read it from AppKit instead.
    Every call has its exact prototype: ``objc_msgSend`` with a wrong one
    crashes the process rather than raising.
    """
    import ctypes

    try:
        objc = ctypes.CDLL("/usr/lib/libobjc.A.dylib")
        ctypes.CDLL("/System/Library/Frameworks/AppKit.framework/AppKit")
    except OSError:
        return None
    ptr = ctypes.c_void_p
    objc.objc_getClass.restype = ptr
    objc.objc_getClass.argtypes = [ctypes.c_char_p]
    objc.sel_registerName.restype = ptr
    objc.sel_registerName.argtypes = [ctypes.c_char_p]
    objc.objc_autoreleasePoolPush.restype = ptr
    objc.objc_autoreleasePoolPush.argtypes = []
    objc.objc_autoreleasePoolPop.restype = None
    objc.objc_autoreleasePoolPop.argtypes = [ptr]
    address = ctypes.cast(objc.objc_msgSend, ptr).value
    send = ctypes.CFUNCTYPE(ptr, ptr, ptr)(address)
    send_object = ctypes.CFUNCTYPE(ptr, ptr, ptr, ptr)(address)
    send_string = ctypes.CFUNCTYPE(ptr, ptr, ptr, ctypes.c_char_p)(address)
    send_length = ctypes.CFUNCTYPE(ctypes.c_ulong, ptr, ptr)(address)

    def sel(name):
        return objc.sel_registerName(name.encode())

    pool = objc.objc_autoreleasePoolPush()
    try:
        pasteboard_class = objc.objc_getClass(b"NSPasteboard")
        string_class = objc.objc_getClass(b"NSString")
        if not pasteboard_class or not string_class:
            return None
        pasteboard = send(pasteboard_class, sel("generalPasteboard"))
        name = send_string(string_class, sel("stringWithUTF8String:"), kind.encode())
        if not pasteboard or not name:
            return None
        data = send_object(pasteboard, sel("dataForType:"), name)
        if not data:
            return None
        length = send_length(data, sel("length"))
        if length > CLIPBOARD_LIMIT:
            raise RskError("The ReShiki drawing on the clipboard is too large to paste.")
        start = send(data, sel("bytes"))
        # Copy out before the pool releases the NSData.
        return ctypes.string_at(start, length) if start and length else None
    finally:
        objc.objc_autoreleasePoolPop(pool)


def clipboard_document(mime, native=None):
    """The ReShiki document on the clipboard, decoded, or None when there is none.

    ``native`` is ReShiki's own format read outside Qt (macOS), used when Qt
    does not offer it.
    """
    if mime is None and native is None:
        return None
    raw = native
    for fmt in _CLIPBOARD_FORMATS if mime is not None else ():
        if mime.hasFormat(fmt):
            raw = bytes(mime.data(fmt))
            break
    if raw is not None:
        # Windows rounds clipboard memory up, so the JSON can trail NULs.
        text = raw.decode("utf-8", errors="replace").rstrip("\x00")
    else:
        text = mime.text() if mime.hasText() else ""
        for prefix in _TEXT_PREFIXES:
            if text.startswith(prefix):
                text = text[len(prefix):]
                break
        else:
            return None
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise RskError(f"The ReShiki drawing on the clipboard could not be read: {exc}") from exc


def load_into_editor(context, atoms, bonds, positions, replace=True):
    from PyQt6.QtCore import QPointF

    scene = context.scene
    if replace:
        context.clear_canvas(push_to_undo=True)
    cx, cy = view_center(context)
    items = {}
    for atom in atoms:
        x, y = positions[atom["id"]]
        atom_id = scene.create_atom(
            atom["symbol"], QPointF(x + cx, y + cy), charge=atom["charge"], radical=atom["radical"]
        )
        items[atom["id"]] = scene.atom_items[atom_id]
    for bond in bonds:
        scene.create_bond(
            items[bond["a"]], items[bond["b"]], bond_order=bond["order"], bond_stereo=bond["stereo"]
        )
    context.refresh_2d_scene()
    context.push_undo_checkpoint()


def view_center(context):
    """Scene point at the middle of the 2D view, where the host places opened files."""
    try:
        view = context.get_main_window().init_manager.view_2d
        center = view.mapToScene(view.viewport().rect().center())
        return center.x(), center.y()
    except AttributeError:
        return 0.0, 0.0


def _warn(message):
    try:
        from PyQt6.QtWidgets import QMessageBox

        QMessageBox.warning(_context.get_main_window(), PLUGIN_NAME, message)
    except Exception:  # noqa: BLE001 - fall back to the status bar when no dialog can open
        _context.show_status_message(message, 8000)


def open_rsk(path):
    """File opener: replace the 2D drawing with the ReShiki drawing at ``path``."""
    try:
        atoms, bonds, positions, notes = read_rsk(path)
    except RskError as exc:
        _warn(str(exc))
        raise
    load_into_editor(_context, atoms, bonds, positions)
    _context.show_status_message(summary(atoms, bonds, notes), 8000)


def paste_from_clipboard():
    """Add the drawing ReShiki last copied to the 2D editor, at the view centre."""
    from PyQt6.QtWidgets import QApplication

    try:
        native = macos_pasteboard_data() if sys.platform == "darwin" else None
        document = clipboard_document(QApplication.clipboard().mimeData(), native)
        if document is None:
            _context.show_status_message(
                "The clipboard holds no ReShiki drawing. Copy one in ReShiki first.", 5000
            )
            return
        atoms, bonds, positions, notes = convert(document)
    except RskError as exc:
        _warn(str(exc))
        return
    load_into_editor(_context, atoms, bonds, positions, replace=False)
    _context.show_status_message(summary(atoms, bonds, notes, "Pasted"), 8000)


def _handle_drop(path):
    if not str(path).lower().endswith(EXTENSION):
        return False
    try:
        open_rsk(path)
    except RskError:
        pass
    return True


def initialize(context):
    global _context
    _context = context
    context.register_file_opener(EXTENSION, open_rsk)
    context.register_drop_handler(_handle_drop)
    context.add_menu_action(
        "Edit/Paste from ReShiki", paste_from_clipboard, shortcut=PASTE_SHORTCUT
    )
