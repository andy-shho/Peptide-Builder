import os
import pytest
import numpy as np
from numpy.testing import assert_array_equal
from residue import ShapeError

# Fix your moodule name
from residue import Residue, create_residue

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


# =====================================================================
# STAGE 1: Basic Constructor
# Students should start here - just getting the class to accept arguments
# =====================================================================

@pytest.mark.stage1
def test_constructor_basic_interface():
    """Residue must accept required arguments and expose basic attributes."""

    atoms = {
        "atom_ids": ["N", "CA", "C"],
        "symbols": ["N", "C", "C"],
        "coordinates": [
            [0.0, 1.0, 2.0],
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
        ],
    }

    bonds = [
        ("N", "CA", "SING"),
        ("CA", "C",  "SING"),
    ]

    res = Residue(
        name="MockResidue",
        atoms=atoms,
        bonds=bonds,
        one_letter_code="M",
        three_letter_code="MCK",
        residue_type="amino_acid",
    )

    assert res.name == "MockResidue"
    assert res.one_letter_code == "M"
    assert res.three_letter_code == "MCK"
    assert res.residue_type == "amino_acid"


@pytest.mark.stage1
def test_atom_ids_symbols_coordinates_alignment():
    """atom_ids, symbols, and coordinates must align in order and length."""

    atoms = {
        "atom_ids": ["A1", "A2", "A3"],
        "symbols": ["C", "O", "N"],
        "coordinates": [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.0, 1.5, 0.0],
        ],
    }

    res = Residue(
        name="Aligned",
        atoms=atoms,
        bonds=[],
        one_letter_code=None,
        three_letter_code="ALG",
        residue_type="amino_acid",
    )

    assert len(res.atom_ids) == 3
    assert len(res.symbols) == 3
    assert len(res.coordinates) == 3

    for i, atom_id in enumerate(res.atom_ids):
        assert isinstance(atom_id, str)
        assert isinstance(res.symbols[i], str)
        coord = res.coordinates[i]
        assert len(coord) == 3
        assert all(isinstance(x, (int, float)) for x in coord)


# =====================================================================
# STAGE 2: Properties
# Adding @property for computed attributes
# =====================================================================

@pytest.mark.stage2
def test_n_atoms_property():
    """n_atoms must be a @property and equal to len(atom_ids)."""

    atoms = {
        "atom_ids": ["X", "Y"],
        "symbols": ["C", "C"],
        "coordinates": [[0, 0, 0], [1, 1, 1]],
    }

    res = Residue(
        name="CountTest",
        atoms=atoms,
        bonds=[],
        one_letter_code=None,
        three_letter_code="CNT",
        residue_type="amino_acid",
    )

    # Ensure n_atoms is a property on the class
    assert isinstance(type(res).n_atoms, property)
    assert res.n_atoms == 2
    assert res.n_atoms == len(res.atom_ids)


# =====================================================================
# STAGE 3: Bond Validation
# Ensuring bonds reference valid atoms
# =====================================================================

@pytest.mark.stage3
def test_bonds_reference_valid_atom_ids():
    """Each bond must reference valid atom_ids and have a string order."""

    atoms = {
        "atom_ids": ["A", "B", "C"],
        "symbols": ["C", "N", "O"],
        "coordinates": [[0,0,0], [1,0,0], [0,1,0]],
    }

    bonds = [
        ("A", "B", "SING"),
        ("B", "C", "DOUB"),
    ]

    res = Residue(
        name="BondCheck",
        atoms=atoms,
        bonds=bonds,
        one_letter_code=None,
        three_letter_code="BDC",
        residue_type="amino_acid",
    )

    valid = set(res.atom_ids)

    for a1, a2, order in res.bonds:
        assert a1 in valid
        assert a2 in valid
        assert isinstance(order, str)
        assert len(order) > 0


# =====================================================================
# STAGE 4: CIF Parsing with from_cif
# Implementing the classmethod to read CIF files
# =====================================================================

@pytest.mark.stage4
def test_from_cif_exists_and_returns_residue():
    """Residue.from_cif must exist and return a Residue instance."""

    assert hasattr(Residue, "from_cif")
    assert callable(Residue.from_cif)

    path = os.path.join(DATA_DIR, "ALA.cif")
    res = Residue.from_cif(path)

    assert isinstance(res, Residue)


@pytest.mark.stage4
def test_from_cif_invalid_path_raises():
    """Invalid paths must raise FileNotFoundError."""

    bogus = os.path.join(DATA_DIR, "DOES_NOT_EXIST.cif")

    with pytest.raises(FileNotFoundError):
        Residue.from_cif(bogus)


@pytest.mark.stage4
@pytest.mark.parametrize(
    "filename, expected_three_letter, expected_type, expected_n_atoms",
    [
        ("ALA.cif", "ALA", "amino_acid", 13),
        ("LYS.cif", "LYS", "amino_acid", 25),
        ("DA.cif",  "DA",  "nucleic_acid", 36),
    ],
)
def test_from_cif_basic_fields(
    filename, expected_three_letter, expected_type, expected_n_atoms
):
    """Check essential fields extracted correctly from CIF files."""

    path = os.path.join(DATA_DIR, filename)
    res = Residue.from_cif(path)

    assert res.three_letter_code == expected_three_letter
    assert res.residue_type == expected_type
    assert res.n_atoms == expected_n_atoms
    assert len(res.atom_ids) == expected_n_atoms
    assert len(res.symbols) == expected_n_atoms
    assert len(res.coordinates) == expected_n_atoms

    for coord in res.coordinates:
        assert len(coord) == 3
        assert all(isinstance(x, (int, float)) for x in coord)


# =====================================================================
# STAGE 5: Factory Method
# Implementing create_residue factory function
# =====================================================================

@pytest.mark.stage5
def test_create_residue_exists():
    """create_residue function must exist and be callable."""

    assert callable(create_residue)


@pytest.mark.stage5
def test_create_residue_returns_residue():
    """create_residue must return a Residue instance."""

    res = create_residue("ALA")
    assert isinstance(res, Residue)


@pytest.mark.stage5
@pytest.mark.parametrize(
    "code, expected_type",
    [
        ("ALA", "amino_acid"),
        ("LYS", "amino_acid"),
        ("DA", "nucleic_acid"),
    ],
)
def test_create_residue_loads_correct_residue(code, expected_type):
    """create_residue must load the correct CIF file based on code."""

    res = create_residue(code)
    assert res.three_letter_code == code
    assert res.residue_type == expected_type

# =====================================================================
# STAGE 6: Setting Coordinates Mechanism
# Implementing getter and setter for coordinates attribute
# =====================================================================

@pytest.mark.stage6
def test_coordinates():
    atoms = {
        "atom_ids": ["X", "Y"],
        "symbols": ["C", "C"],
        "coordinates": [[0, 0, 0], [1, 1, 1]],
    }

    res = Residue(
        name="CountTest",
        atoms=atoms,
        bonds=[],
        one_letter_code=None,
        three_letter_code="CNT",
        residue_type="amino_acid",
    )
    assert_array_equal(res.coordinates, atoms["coordinates"])

    new_coordinates = [[2, 2, 2], [4, 4, 4]]

    res.coordinates = new_coordinates

    assert_array_equal(res.coordinates, new_coordinates)

    with pytest.raises(ShapeError):
        res.coordinates = [[2, 2, 2]]

# =====================================================================
# STAGE 7: Removing atoms
# Implementing the remove_atom method
# =====================================================================

@pytest.mark.stage7
def test_remove_atom_atom_ids():
    atoms = {
        "atom_ids": ["X", "Y", "Z"],
        "symbols": ["C", "C", "H"],
        "coordinates": [[0, 0, 0], [1, 1, 1], [2, 2, 2]],
    }

    res = Residue(
        name="CountTest",
        atoms=atoms,
        bonds=[
            ("X", "Y", "SING"),
            ("Y", "Z", "SING")
        ],
        one_letter_code=None,
        three_letter_code="CNT",
        residue_type="amino_acid",
    )

    res.remove_atom("X")

    assert res.n_atoms == 2
    assert "X" not in res.atom_ids

@pytest.mark.stage7
def test_remove_atom_coordinates():
    atoms = {
        "atom_ids": ["X", "Y", "Z"],
        "symbols": ["C", "C", "H"],
        "coordinates": [[0, 0, 0], [1, 1, 1], [2, 2, 2]],
    }

    res = Residue(
        name="CountTest",
        atoms=atoms,
        bonds=[
            ("X", "Y", "SING"),
            ("Y", "Z", "SING")
        ],
        one_letter_code=None,
        three_letter_code="CNT",
        residue_type="amino_acid",
    )

    res.remove_atom("X")

    assert len(res.coordinates) == 2
    assert_array_equal(res.coordinates, atoms["coordinates"][1:])

@pytest.mark.stage7
def test_remove_atom_bonds():
    atoms = {
        "atom_ids": ["X", "Y", "Z"],
        "symbols": ["C", "C", "H"],
        "coordinates": [[0, 0, 0], [1, 1, 1], [2, 2, 2]],
    }

    res = Residue(
        name="CountTest",
        atoms=atoms,
        bonds=[
            ("X", "Y", "SING"),
            ("Y", "Z", "SING")
        ],
        one_letter_code=None,
        three_letter_code="CNT",
        residue_type="amino_acid",
    )

    res.remove_atom("X")

    assert len(res.bonds) == 1
    assert res.bonds == [("Y", "Z", "SING")]

# =====================================================================
# STAGE 8: Computing the backbone N->C vector
# Implementing the mechanism to compute backbone N->C
# =====================================================================

@pytest.mark.stage8
def test_computing_backbone_vector():
    res = create_residue("ALA")

    backbone_vector = res.backbone_vector

    assert len(backbone_vector) == 3
    assert isinstance(backbone_vector, np.ndarray)
    assert backbone_vector[2] == 0