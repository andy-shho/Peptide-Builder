import os
import pytest
import numpy as np
from numpy.testing import assert_array_equal

# Fix your moodule name
from peptide import Peptide, NoInputError

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# =====================================================================
# STAGE 9: Converting one-letter-sequence to three-letter-sequence
# Implementing the mechanism for conversion
# =====================================================================

@pytest.mark.stage9
def test_sequence_conversion():
    pep = Peptide("AGDV")

    correct_three_letter_seq = "Ala-Gly-Asp-Val"
    assert pep.three_letter_sequence == correct_three_letter_seq

# =====================================================================
# STAGE 10: Creating peptide chain
# Implementing the mechanism to create peptide chain
# =====================================================================

@pytest.mark.stage10
def test_peptide_chain_creation():
    pep = Peptide("AGDV")

    assert len(pep.bonds) == 3
    assert len(pep.peptide_chain) == 4
    assert len(pep.unique_residue_index) == 4

# =====================================================================
# STAGE 11: Getting residue
# Implementing the get_residue method to access a residue
# =====================================================================

@pytest.mark.stage11
def test_get_residue():
    pep = Peptide("AGDV")

    assert pep.get_residue(index=0).one_letter_code == "A"
    assert pep.get_residue(unique_name="ALA_1").one_letter_code == "A"
    with pytest.raises(NoInputError):
        pep.get_residue()