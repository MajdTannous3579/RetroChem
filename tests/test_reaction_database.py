import json
import os
import pytest

from retrochem.reaction_database import (
    reverse_reaction_generator,
    load_database,
    register_database,
    list_reactants,
    add_new_smart,
    REACTION_DATABASES,
    REACTION_REVERSERS,
)


# ──────────────────────────────────────────────────────────────────────────────
# Global setup/teardown
# ──────────────────────────────────────────────────────────────────────────────
def setup_function(function):
    """Clear global state before each test."""
    REACTION_DATABASES.clear()
    REACTION_REVERSERS.clear()

"""
Tests for the `reverse_reaction_generator` function from the `reaction_database` module.

The `reverse_reaction_generator` function is tested under different scenarios to ensure its correctness.

Test Cases:
- Test in the case of a valid reaction SMARTS, the reverser returns the expected reactant SMILES and conditions.
- Test in the case the input SMILES does not match the reactant SMARTS, the reverser returns None.
- Test in the case the input SMILES is invalid, the reverser raises ValueError.
"""


def test_reverse_reaction_generator_success():
    """
    Given a reaction SMARTS that splits CC → C.C, calling the returned reverser
    on "CC" should give back ("C.C", conditions).
    """
    reaction_smart = ("CC>>C.C", {"temperature": "25 °C"})
    reverser = reverse_reaction_generator(reaction_smart)

    result = reverser("CC")
    # sort of “C.C” is canonical: both fragments are “C”
    assert result == ("C.C", {"temperature": "25 °C"})


def test_reverse_reaction_generator_no_products():
    """
    If the input SMILES doesn't match the reactant pattern,
    the reverser should return None.
    """
    reaction_smart = ("CC>>C.C", {"foo": "bar"})
    reverser = reverse_reaction_generator(reaction_smart)

    # “C” does not match the left-hand CC pattern → no products
    assert reverser("C") is None


def test_reverse_reaction_generator_invalid_smiles():
    """
    Passing a truly invalid SMILES string should raise ValueError("Invalid SMILES").
    """
    reaction_smart = ("CC>>C.C", {})
    reverser = reverse_reaction_generator(reaction_smart)

    with pytest.raises(ValueError) as excinfo:
        reverser("not-a-smiles")
    assert "Invalid SMILES" in str(excinfo.value)



"""
Tests for the `load_database` function from the `reaction_database` module.

The `load_database` function is tested under different scenarios to ensure its correctness.

Test Cases:
- Test in the case the JSON file exists and contains a valid list of SMARTS–condition pairs.
- Test in the case the file does not exist (should return None).
- Test in the case the file contains malformed JSON (should return None).
"""


def test_load_database_valid(tmp_path):
    """
    Given a JSON file containing a list of SMARTS‐condition pairs,
    load_database should return the same list.
    """
    data = [
        ["CC>>C.C", {"foo": "bar"}],
        ["CCC>>C.C.C", {"temp": "100K"}]
    ]
    db_file = tmp_path / "test.db"
    db_file.write_text(json.dumps(data))

    loaded = load_database(str(db_file))
    assert loaded == data


def test_load_database_nonexistent_file():
    """
    If the file does not exist, load_database should return None.
    """
    assert load_database("this_file_does_not_exist.db") is None


def test_load_database_malformed_json(tmp_path):
    """
    If the file exists but contains invalid JSON, load_database should return None.
    """
    bad_file = tmp_path / "bad.db"
    bad_file.write_text("{ not valid JSON }")

    assert load_database(str(bad_file)) is None



"""
Tests for the `register_database` function from the `reaction_database` module.

The `register_database` function is tested under different scenarios to ensure its correctness.

Test Cases:
- Test in the case of registering a new database: REACTION_DATABASES and REACTION_REVERSERS populate correctly.
- Test in the case of registering multiple entries: dictionaries reflect all entries.
- Test in the case of overwriting an existing database: old entries are replaced by new values.
"""

def test_register_database_single_entry():
    """Registering one SMARTS–conditions pair creates entries in both dicts."""
    values = [("CC>>C.C", {"cond": "A"})]
    register_database(values, "db1")

    assert "db1" in REACTION_DATABASES
    assert REACTION_DATABASES["db1"] == values
    assert "db1" in REACTION_REVERSERS
    assert len(REACTION_REVERSERS["db1"]) == 1
    assert callable(REACTION_REVERSERS["db1"][0])


def test_register_database_multiple_entries():
    """Registering multiple pairs stores them all in the dicts."""
    values = [
        ("CC>>C.C", {"x": 1}),
        ("CCC>>C.C.C", {"y": 2})
    ]
    register_database(values, "db2")

    assert REACTION_DATABASES["db2"] == values
    assert len(REACTION_REVERSERS["db2"]) == 2
    for fn in REACTION_REVERSERS["db2"]:
        assert callable(fn)


def test_register_database_overwrite():
    """Registering the same database name twice replaces the old entries."""
    initial = [("C>>C", {})]
    register_database(initial, "db3")

    new_vals = [("CC>>C.C", {"z": 3})]
    register_database(new_vals, "db3")

    assert REACTION_DATABASES["db3"] == new_vals
    assert len(REACTION_REVERSERS["db3"]) == 1
    # Ensure it’s the new function, not leftover from initial
    assert callable(REACTION_REVERSERS["db3"][0])



"""
Tests for the `list_reactants` function from the `reaction_database` module.

The `list_reactants` function is tested under different scenarios to ensure its correctness.

Test Cases:
- Test in the case the database name is not registered (should return None).
- Test in the case the database is registered but no reactions match the input SMILES (should return an empty list).
- Test in the case the database is registered and some reactions match (should return all matching SMILES–conditions pairs).
"""

def test_list_reactants_unknown_database():
    """If the database does not exist, list_reactants should return None."""
    assert list_reactants("CC", "no_such_db") is None

def test_list_reactants_no_matches():
    """If no reactions match the input SMILES, list_reactants should return an empty list."""
    register_database([("CC>>C.C", {"foo": "bar"})], "dbA")
    assert list_reactants("C", "dbA") == []

def test_list_reactants_with_matches():
    """
    If reactions match, list_reactants should return all SMILES–conditions pairs.
    Note: the CC>>C.C rule also applies when you feed in "CCC", so you get two entries.
    """
    values = [
        ("CC>>C.C", {"foo": "bar"}),
        ("CCC>>C.C.C", {"baz": "qux"})
    ]
    register_database(values, "dbB")

    # Only the first rule matches "CC"
    assert list_reactants("CC", "dbB") == [("C.C", {"foo": "bar"})]

    # Both rules match "CCC" (CC is a substructure of CCC, and the second rule is exact)
    assert list_reactants("CCC", "dbB") == [
        ("C.C", {"foo": "bar"}),
        ("C.C.C", {"baz": "qux"})
    ]


"""
Tests for the `add_new_smart` function from the `reaction_database` module.

The `add_new_smart` function is tested under different scenarios to ensure its correctness.

Test Cases:
- Test in the case of creating a new database: file is created, globals are updated.
- Test in the case of appending a second SMARTS: file contains both entries, globals reflect both.
- Test in the case of no explicit conditions: default empty dict is used.
"""

def test_add_new_smart_creates_file_and_registers(tmp_path, monkeypatch):
    """Creating a new SMARTS entry writes `<db>.db` and updates globals."""
    monkeypatch.chdir(tmp_path)

    # Add first SMARTS entry
    add_new_smart("dbX", product="CC", reactants=["C", "C"], conditions={"a": "b"})

    # The database file must exist
    db_file = tmp_path / "dbX.db"
    assert db_file.exists()

    # Pull out what was actually registered
    smarts0, cond0 = REACTION_DATABASES["dbX"][0]

    # Check globals
    assert cond0 == {"a": "b"}
    assert callable(REACTION_REVERSERS["dbX"][0])

    # Check file contents include the exact SMARTS and JSON-dumped conditions
    content = db_file.read_text()
    assert smarts0 in content
    assert json.dumps(cond0) in content


def test_add_new_smart_appends_second_entry(tmp_path, monkeypatch):
    """Adding a second SMARTS to the same database appends to file and globals."""
    monkeypatch.chdir(tmp_path)

    # First entry
    add_new_smart("dbY", product="CC",  reactants=["C", "C"],   conditions={"x": "1"})
    # Second entry
    add_new_smart("dbY", product="CCC", reactants=["C", "C", "C"], conditions={"y": "2"})

    # There should be two entries in globals
    assert len(REACTION_DATABASES["dbY"]) == 2
    assert len(REACTION_REVERSERS["dbY"]) == 2

    # And the file must mention both SMARTS + both condition-dicts
    db_file = tmp_path / "dbY.db"
    content = db_file.read_text()
    for smarts, cond in REACTION_DATABASES["dbY"]:
        assert smarts in content
        assert json.dumps(cond) in content


def test_add_new_smart_defaults_empty_conditions(tmp_path, monkeypatch):
    """If no conditions dict is passed, an empty dict is recorded."""
    monkeypatch.chdir(tmp_path)

    # No conditions argument here
    add_new_smart("dbZ", product="CO", reactants=["C", "O"])

    smarts0, cond0 = REACTION_DATABASES["dbZ"][0]
    assert cond0 == {}  # default

    # File should include the empty-dict literal
    content = (tmp_path / "dbZ.db").read_text()
    assert json.dumps(cond0) in content  # "{}"