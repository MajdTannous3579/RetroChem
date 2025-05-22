from rdkit import Chem
import os
import pytest

from retrochem.preliminary_functions import (
    canonicalize_smiles,
    name_to_smiles,
    structure_to_smiles,
)

"""
Tests for the `canonicalize_smiles` function from `preliminary_functions` module.

The `canonicalize_smiles` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the input SMILES is already in canonical form
Test in the case a non-canonical SMILES is inputed
Test in the case the input is an aromatic lower-case input
Test in the case the input is an invalid smiles string

"""


def test_canonicalize_smiles_already_canonical():
    """
    If the input SMILES is already in canonical form, it should return the same string.
    """
    smi = "CCO"  # ethanol
    assert canonicalize_smiles(smi) == "CCO"

def test_canonicalize_smiles_non_canonical_input():
    """
    A non-canonical SMILES (atom order shuffled) should be converted to the canonical form.
    """
    # 'OC(C)' is chemically the same as 'CCO' but with a different atom order
    assert canonicalize_smiles("OC(C)") == "CCO"

def test_canonicalize_smiles_aromatic_vs_aliphatic():
    """
    Aromatic lower-case input should be canonicalized to uppercase Kekulé form (if RDKit chooses that).
    """
    # benzene aromatic notation
    result = canonicalize_smiles("c1ccccc1")
    # RDKit canonical SMILES for benzene is "c1ccccc1"
    assert result.lower() == "c1ccccc1"

def test_canonicalize_smiles_invalid_raises():
    """
    An invalid SMILES string should raise ValueError.
    """
    with pytest.raises(ValueError):
        canonicalize_smiles("this-is-not-smiles")




"""
Tests for the `name_to_smiles` function from `preliminary_functions` module.

The `name_to_smiles` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test when the API returns a valid SMILES string.
Test when the API returns an empty response (should raise ValueError).
Test when `urlopen` raises an exception (network or other I/O error).
Test that names with spaces/special characters get URL-encoded.
Test that leading/trailing whitespace in the API response is stripped.

"""

class DummyResponse:
    """
    Simulate a urllib response object with a .read() method.
    """
    def __init__(self, data: bytes):
        self._data = data

    def read(self) -> bytes:
        return self._data


def test_name_to_smiles_success(monkeypatch):
    """
    When the API returns a valid SMILES string, name_to_smiles should return it.
    """
    # Simulate a 200-OK response with "C(C)O\n"
    monkeypatch.setattr(
        "retrochem.preliminary_functions.urlopen",
        lambda url: DummyResponse(b"C(C)O\n")
    )
    assert name_to_smiles("ethanol") == "C(C)O"


def test_name_to_smiles_empty_response(monkeypatch):
    """
    An empty response body should trigger a ValueError.
    """
    monkeypatch.setattr(
        "retrochem.preliminary_functions.urlopen",
        lambda url: DummyResponse(b"")
    )
    with pytest.raises(ValueError) as excinfo:
        name_to_smiles("nothing")
    assert "Empty response for name: nothing" in str(excinfo.value)


def test_name_to_smiles_network_error(monkeypatch):
    """
    If urlopen itself raises, name_to_smiles should wrap it in a ValueError.
    """
    def bad_connect(url):
        raise RuntimeError("no network")

    monkeypatch.setattr(
        "retrochem.preliminary_functions.urlopen",
        bad_connect
    )

    with pytest.raises(ValueError) as excinfo:
        name_to_smiles("ethanol")
    assert "Could not convert name 'ethanol' to SMILES" in str(excinfo.value)


"""
Tests for the `structure_to_smiles` function from `preliminary_functions` module.

The `structure_to_smiles` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test that passing an RDKit Mol object returns its SMILES string.
Test that passing a MolBlock string (containing 'M  END') returns the correct SMILES.
Test that reading from a .mol file on disk returns the correct SMILES.
Test that giving a file with an unsupported extension raises ValueError.
Test that passing a string which is neither a file nor a MolBlock raises ValueError.
Test that passing a non-str, non-Mol (e.g. an integer) raises TypeError.
Test that a MolBlock which fails to parse (RDKit returns None) raises ValueError.

"""


MOLBLOCK = """
  Mrv2007 07082006022D          

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.2094    0.0000    0.0000 C   0  0
   -0.6070    1.0475    0.0000 O   0  0
  1  2  1  0
  1  3  1  0
M  END
"""

def test_structure_to_smiles_from_mol_object():
    """Passing an RDKit Mol should return its SMILES."""
    mol = Chem.MolFromSmiles("CCO")
    assert structure_to_smiles(mol) == "CCO"

def test_structure_to_smiles_from_molblock_string():
    """A string containing a MolBlock (with 'M  END') should parse to non-canonical SMILES."""
    assert structure_to_smiles(MOLBLOCK) == "C(C)O"

def test_structure_to_smiles_from_file(tmp_path):
    """A .mol file on disk should be read and converted to non-canonical SMILES."""
    molfile = tmp_path / "ethanol.mol"
    molfile.write_text(MOLBLOCK)
    result = structure_to_smiles(str(molfile), input_format="mol")
    assert result == "C(C)O"

def test_structure_to_smiles_unsupported_file_extension(tmp_path):
    """Files with unsupported extensions should raise ValueError."""
    txt = tmp_path / "foo.txt"
    txt.write_text(MOLBLOCK)
    with pytest.raises(ValueError) as exc:
        structure_to_smiles(str(txt))
    assert "Unsupported file format" in str(exc.value)

def test_structure_to_smiles_bad_string_input():
    """A string that is neither a file nor a MolBlock should raise ValueError."""
    with pytest.raises(ValueError) as exc:
        structure_to_smiles("not a file or block")
    assert "neither a valid file path nor a MolBlock" in str(exc.value)

def test_structure_to_smiles_wrong_type():
    """Passing a non-str, non‐Mol (e.g. an integer) should raise TypeError."""
    with pytest.raises(TypeError):
        structure_to_smiles(12345)

def test_structure_to_smiles_parsing_failure(monkeypatch):
    """
    If RDKit fails to parse a MolBlock (returns None), we should get ValueError.
    """
    # MolBlock branch: string contains 'M  END' but is invalid
    bad_block = "M  END"
    with pytest.raises(ValueError) as exc:
        structure_to_smiles(bad_block)
    assert "Failed to parse molecular structure" in str(exc.value)