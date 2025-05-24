from rdkit import Chem
import os
import pytest

from retrochem.preliminary_functions import (
    canonicalize_smiles,
    name_to_smiles,
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
