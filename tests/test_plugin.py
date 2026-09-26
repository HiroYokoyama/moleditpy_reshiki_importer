import importlib.util
import json
import math
import os
import re
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = ROOT / "reshiki_importer.py"
RESHIKI = ROOT.parent / "ReShiki"


@pytest.fixture
def plugin():
    spec = importlib.util.spec_from_file_location("reshiki_importer", PLUGIN)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FakeContext:
    def __init__(self):
        self.openers = {}
        self.drop_handlers = []
        self.menu_actions = {}
        self.messages = []
        self.calls = []
        self.scene = FakeScene()

    def register_file_opener(self, extension, callback, priority=0):
        self.openers[extension] = callback

    def register_drop_handler(self, callback, priority=0):
        self.drop_handlers.append(callback)

    def add_menu_action(self, path, callback, text=None, icon=None, shortcut=None):
        self.menu_actions[path] = (callback, shortcut)

    def show_status_message(self, message, timeout=3000):
        self.messages.append(message)

    def get_main_window(self):
        return None

    def clear_canvas(self, push_to_undo=True):
        self.calls.append("clear")
        self.scene.atoms.clear()
        self.scene.bonds.clear()

    def refresh_2d_scene(self):
        self.calls.append("refresh")

    def fit_2d_view(self):
        self.calls.append("fit")

    def push_undo_checkpoint(self):
        self.calls.append("undo")


class FakeScene:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.atom_items = {}

    def create_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = len(self.atoms)
        self.atoms.append((symbol, pos.x(), pos.y(), charge, radical))
        self.atom_items[atom_id] = atom_id
        return atom_id

    def create_bond(self, start, end, bond_order=None, bond_stereo=None):
        self.bonds.append((start, end, bond_order, bond_stereo))


class FakeMime:
    def __init__(self, formats=None, text=None):
        self.formats = formats or {}
        self._text = text

    def hasFormat(self, fmt):
        return fmt in self.formats

    def data(self, fmt):
        return self.formats[fmt]

    def hasText(self):
        return self._text is not None

    def text(self):
        return self._text


class FakePoint:
    def __init__(self, x, y):
        self._x, self._y = x, y

    def x(self):
        return self._x

    def y(self):
        return self._y


@pytest.fixture
def qt_stub(monkeypatch):
    """QPointF without a GUI toolkit, so placement runs headless."""
    core = types.ModuleType("PyQt6.QtCore")
    core.QPointF = FakePoint
    widgets = types.ModuleType("PyQt6.QtWidgets")
    widgets.QMessageBox = None
    package = types.ModuleType("PyQt6")
    package.QtCore, package.QtWidgets = core, widgets
    clipboard = types.SimpleNamespace(mime=None)
    clipboard.mimeData = lambda: clipboard.mime
    widgets.QApplication = types.SimpleNamespace(clipboard=lambda: clipboard)
    monkeypatch.setitem(sys.modules, "PyQt6", package)
    monkeypatch.setitem(sys.modules, "PyQt6.QtCore", core)
    monkeypatch.setitem(sys.modules, "PyQt6.QtWidgets", widgets)
    return clipboard


def atom(atom_id, element="C", x=0.0, y=0.0, **extra):
    return {"id": atom_id, "element": element, "position": {"x": x, "y": y}, **extra}


def bond(a, b, order=1, display="plain", **extra):
    return {"a": a, "b": b, "order": order, "display": display, **extra}


def drawing(atoms, bonds=(), **extra):
    return {"version": 15, "atoms": list(atoms), "bonds": list(bonds), **extra}


def benzene(order=4):
    atoms = [
        atom(i + 1, x=42 * math.cos(math.radians(60 * i)), y=42 * math.sin(math.radians(60 * i)))
        for i in range(6)
    ]
    return drawing(atoms, [bond(i + 1, (i + 1) % 6 + 1, order) for i in range(6)])


def test_metadata_matches_the_registry_contract(plugin):
    for name in ("PLUGIN_NAME", "PLUGIN_VERSION", "PLUGIN_AUTHOR", "PLUGIN_DESCRIPTION"):
        assert getattr(plugin, name).strip()
    assert re.fullmatch(r"\d+\.\d+\.\d+", plugin.PLUGIN_VERSION)
    assert not hasattr(plugin, "run") and not hasattr(plugin, "autorun")


def test_importing_needs_no_qt_or_rdkit():
    source = PLUGIN.read_text(encoding="utf-8")
    top_level = [line for line in source.splitlines() if line.startswith(("import ", "from "))]
    assert not any("PyQt6" in line or "rdkit" in line for line in top_level)


def test_initialize_registers_the_opener_drop_handler_and_paste(plugin):
    context = FakeContext()
    plugin.initialize(context)
    assert set(context.openers) == {".rsk"}
    assert len(context.drop_handlers) == 1
    assert context.menu_actions == {
        "Edit/Paste from ReShiki": (plugin.paste_from_clipboard, "Ctrl+Alt+V")
    }


def test_atoms_charges_and_radicals(plugin):
    atoms, bonds, notes = plugin.parse_rsk(
        drawing(
            [atom(1, "N", charge=1), atom(2, "C", 42, radical_electrons=1), atom(3, "O", 84)],
            [bond(1, 2), bond(2, 3, 2)],
        )
    )
    assert [(a["symbol"], a["charge"], a["radical"]) for a in atoms] == [
        ("N", 1, 0),
        ("C", 0, 1),
        ("O", 0, 0),
    ]
    assert [(b["a"], b["b"], b["order"]) for b in bonds] == [(1, 2, 1), (2, 3, 2)]
    assert not notes


@pytest.mark.parametrize(
    "display, stereo",
    [("wedge", 1), ("hollow_wedge", 1), ("bold", 1), ("hash", 2), ("hashed", 2), ("wavy", 0), ("plain", 0)],
)
def test_single_bond_displays_become_stereo(plugin, display, stereo):
    _, bonds, _ = plugin.parse_rsk(drawing([atom(1), atom(2, x=42)], [bond(1, 2, 1, display)]))
    assert bonds[0]["stereo"] == stereo


def test_bold_aromatic_edges_and_projections_carry_no_stereo(plugin):
    _, bonds, _ = plugin.parse_rsk(
        drawing(
            [atom(1), atom(2, x=42), atom(3, x=84)],
            [bond(1, 2, 4, "bold"), bond(2, 3, 1, "wedge", projection=True)],
        )
    )
    assert [b["stereo"] for b in bonds] == [0, 0]


def test_unusual_orders_are_approximated_and_reported(plugin):
    atoms = [atom(i, x=42 * i) for i in range(1, 7)]
    bonds = [bond(1, 2, 5), bond(2, 3, 6), bond(3, 4, 7), bond(4, 5, 0), bond(5, 6, 9)]
    _, parsed, notes = plugin.parse_rsk(drawing(atoms, bonds))
    assert [b["order"] for b in parsed] == [1, 3, 1, 1]
    assert notes["dative bonds imported as single"] == 1
    assert notes["quadruple bonds imported as triple"] == 1
    assert notes["hydrogen-interaction bonds dropped"] == 1
    assert notes["unknown bond orders imported as single"] == 1


def test_centroids_and_their_bonds_are_skipped(plugin):
    atoms, bonds, notes = plugin.parse_rsk(
        drawing(
            [atom(1, "Fe"), atom(2, "*", centroid=[3]), atom(3, x=42)],
            [bond(1, 2), bond(2, 3)],
        )
    )
    assert [a["id"] for a in atoms] == [1, 3]
    assert bonds == []
    assert notes["centroid and placeholder atoms skipped"] == 1
    assert notes["bonds to skipped atoms dropped"] == 2


def test_figure_content_is_reported(plugin):
    _, _, notes = plugin.parse_rsk(
        drawing([atom(1)], arrows=[{}, {}], annotations=[{}], graphics=[], reactions=[{}])
    )
    assert notes == {"arrows skipped": 2, "text annotations skipped": 1, "reaction schemes skipped": 1}


def test_abbreviation_members_are_imported_expanded(plugin):
    atoms, bonds, _ = plugin.parse_rsk(
        drawing(
            [atom(1), atom(2, "O", 42), atom(3, x=84)],
            [bond(1, 2), bond(2, 3)],
            abbreviations=[{"label": "OMe", "anchor": 2, "members": [2, 3]}],
        )
    )
    assert len(atoms) == 3 and len(bonds) == 2


def test_newer_versions_load_with_a_note(plugin):
    atoms, _, notes = plugin.parse_rsk({"version": 99, "atoms": [atom(1)]})
    assert len(atoms) == 1
    assert any("format 99" in label for label in notes)


@pytest.mark.parametrize(
    "document",
    [[], {"atoms": []}, {"version": 15}, {"version": "15", "atoms": []}, {"version": 15, "atoms": [{"id": 1}]}],
)
def test_malformed_documents_are_rejected(plugin, document):
    with pytest.raises(plugin.RskError):
        plugin.parse_rsk(document)


def test_layout_uses_the_editor_bond_length_and_centres(plugin):
    atoms, bonds, _ = plugin.parse_rsk(drawing([atom(1, x=100, y=50), atom(2, x=142, y=50)], [bond(1, 2)]))
    positions = plugin.layout(atoms, bonds)
    assert positions[1] == pytest.approx((-37.5, 0.0))
    assert positions[2] == pytest.approx((37.5, 0.0))


def test_layout_without_bonds_uses_reshiki_default_length(plugin):
    atoms, bonds, _ = plugin.parse_rsk(drawing([atom(1), atom(2, x=42)]))
    positions = plugin.layout(atoms, bonds)
    assert math.dist(positions[1], positions[2]) == pytest.approx(75.0)


def test_kekulize_benzene(plugin):
    pytest.importorskip("rdkit")
    atoms, bonds, notes = plugin.parse_rsk(benzene())
    plugin.kekulize(atoms, bonds, notes)
    assert sorted(b["order"] for b in bonds) == [1, 1, 1, 2, 2, 2]
    assert not notes


def test_kekulize_pyrrole_nh(plugin):
    pytest.importorskip("rdkit")
    atoms = [atom(1, "N", explicit_h=1)] + [atom(i, x=42 * i) for i in range(2, 6)]
    bonds = [bond(i, i % 5 + 1, 4) for i in range(1, 6)]
    atoms, bonds, notes = plugin.parse_rsk(drawing(atoms, bonds))
    plugin.kekulize(atoms, bonds, notes)
    assert sorted(b["order"] for b in bonds) == [1, 1, 1, 2, 2]
    assert not notes


def test_unkekulizable_rings_fall_back_to_single(plugin, monkeypatch):
    atoms, bonds, notes = plugin.parse_rsk(benzene())
    monkeypatch.setitem(sys.modules, "rdkit", None)
    plugin.kekulize(atoms, bonds, notes)
    assert all(b["order"] == 1 for b in bonds)
    assert notes["aromatic bonds that could not be kekulized imported as single"] == 6


def test_open_replaces_the_drawing(plugin, qt_stub, tmp_path):
    path = tmp_path / "ethanol.rsk"
    path.write_text(
        json.dumps(
            drawing(
                [atom(1), atom(2, x=42), atom(3, "O", 84)],
                [bond(1, 2, 1, "wedge"), bond(2, 3)],
                arrows=[{}],
            )
        ),
        encoding="utf-8",
    )
    context = FakeContext()
    plugin.initialize(context)
    context.openers[".rsk"](str(path))
    assert [a[0] for a in context.scene.atoms] == ["C", "C", "O"]
    assert context.scene.bonds == [(0, 1, 1, 1), (1, 2, 1, 0)]
    assert context.calls == ["clear", "refresh", "undo"]
    assert "Imported 3 atoms and 2 bonds" in context.messages[-1]
    assert "1 arrows skipped" in context.messages[-1]


def test_open_centres_the_molecule_in_the_current_view(plugin, qt_stub, tmp_path):
    class View:
        def viewport(self):
            return self

        def rect(self):
            return self

        def center(self):
            return "viewport centre"

        def mapToScene(self, point):
            assert point == "viewport centre"
            return FakePoint(500.0, -200.0)

    window = types.SimpleNamespace(init_manager=types.SimpleNamespace(view_2d=View()))
    path = tmp_path / "ethane.rsk"
    path.write_text(json.dumps(drawing([atom(1), atom(2, x=42)], [bond(1, 2)])), encoding="utf-8")
    context = FakeContext()
    context.get_main_window = lambda: window
    plugin.initialize(context)
    context.openers[".rsk"](str(path))
    placed = [(x, y) for _, x, y, _, _ in context.scene.atoms]
    assert placed == [pytest.approx((462.5, -200.0)), pytest.approx((537.5, -200.0))]


def test_open_rejects_a_broken_file_without_touching_the_drawing(plugin, qt_stub, tmp_path):
    path = tmp_path / "broken.rsk"
    path.write_text("{not json", encoding="utf-8")
    context = FakeContext()
    plugin.initialize(context)
    with pytest.raises(plugin.RskError):
        context.openers[".rsk"](str(path))
    assert context.calls == []
    assert "broken.rsk" in context.messages[-1]


def test_drop_handler_claims_only_rsk_files(plugin, qt_stub, tmp_path):
    context = FakeContext()
    plugin.initialize(context)
    drop = context.drop_handlers[0]
    assert drop(str(tmp_path / "molecule.mol")) is False
    broken = tmp_path / "Broken.RSK"
    broken.write_text("[]", encoding="utf-8")
    assert drop(str(broken)) is True
    assert context.calls == []


@pytest.mark.skipif(not RESHIKI.is_dir(), reason="ReShiki checkout not found beside this repository")
def test_every_drawing_in_the_reshiki_repository_imports():
    """Parity check against real files, when ReShiki is checked out next to this repo."""
    spec = importlib.util.spec_from_file_location("reshiki_importer", PLUGIN)
    plugin = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plugin)
    files = [
        p
        for p in RESHIKI.rglob("*.rsk")
        if not {"target", "node_modules", ".git"} & set(p.parts)
    ]
    assert files
    imported = 0
    for path in files:
        document = json.loads(path.read_text(encoding="utf-8"))
        if not plugin.parse_rsk(document)[0]:
            with pytest.raises(plugin.RskError):
                plugin.read_rsk(str(path))
            continue
        atoms, bonds, positions, notes = plugin.read_rsk(str(path))
        imported += 1
        assert len(positions) == len(atoms), path
        assert all(b["order"] in (1, 2, 3) for b in bonds), path
        if "rdkit" in sys.modules:
            assert not any("kekulized" in label for label in notes), (path, notes)
    assert imported


ETHANOL = drawing([atom(1), atom(2, x=42), atom(3, "O", 84)], [bond(1, 2), bond(2, 3)])


@pytest.mark.parametrize(
    "fmt", ['application/x-qt-windows-mime;value="dev.reshiki.drawing"', "dev.reshiki.drawing"]
)
def test_clipboard_reads_reshikis_native_format(plugin, fmt):
    # Windows hands back the whole rounded-up allocation, NULs included.
    raw = json.dumps(ETHANOL).encode("utf-8") + b"\x00\x00\x00"
    mime = FakeMime({fmt: raw}, text="CCO")
    assert plugin.clipboard_document(mime) == ETHANOL


@pytest.mark.parametrize("prefix", ["RESHIKI_DRAWING_V1\n", "MORUNO_DRAWING_V1\n"])
def test_clipboard_reads_reshikis_text_copy(plugin, prefix):
    mime = FakeMime(text=prefix + json.dumps(ETHANOL))
    assert plugin.clipboard_document(mime) == ETHANOL


@pytest.mark.parametrize("mime", [None, FakeMime(), FakeMime(text="CCO"), FakeMime(text='{"version": 15}')])
def test_clipboard_without_a_reshiki_drawing(plugin, mime):
    assert plugin.clipboard_document(mime) is None


def test_clipboard_with_a_broken_drawing_is_an_error(plugin):
    with pytest.raises(plugin.RskError):
        plugin.clipboard_document(FakeMime(text="RESHIKI_DRAWING_V1\n{not json"))


def test_paste_adds_to_the_drawing(plugin, qt_stub):
    qt_stub.mime = FakeMime({"dev.reshiki.drawing": json.dumps(ETHANOL).encode()})
    context = FakeContext()
    context.scene.atoms.append(("N", 0.0, 0.0, 0, 0))
    plugin.initialize(context)
    context.menu_actions["Edit/Paste from ReShiki"][0]()
    assert [a[0] for a in context.scene.atoms] == ["N", "C", "C", "O"]
    assert context.calls == ["refresh", "undo"]
    assert context.messages[-1].startswith("Pasted 3 atoms and 2 bonds from ReShiki.")


def test_paste_with_nothing_to_paste_leaves_the_drawing(plugin, qt_stub):
    qt_stub.mime = FakeMime(text="CCO")
    context = FakeContext()
    plugin.initialize(context)
    plugin.paste_from_clipboard()
    assert context.calls == [] and context.scene.atoms == []
    assert "no ReShiki drawing" in context.messages[-1]


@pytest.mark.parametrize("text", ["{not json", json.dumps(drawing([]))])
def test_paste_of_an_unusable_drawing_warns(plugin, qt_stub, text):
    qt_stub.mime = FakeMime(text="RESHIKI_DRAWING_V1\n" + text)
    context = FakeContext()
    plugin.initialize(context)
    plugin.paste_from_clipboard()
    assert context.calls == []
    assert "ReShiki" in context.messages[-1]
