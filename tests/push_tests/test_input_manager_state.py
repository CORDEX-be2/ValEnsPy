"""Regression tests: InputManager / CatalogBuilder must not share mutable state across instances.

Every "dataset" here is a temp dir of empty .nc files - the catalog builder only lists and
regex-matches file names, it never opens them - so no machine-specific paths are needed.
"""

import pytest

from valenspy.input import esm_catalog_builder
from valenspy.input.converter import INPUT_CONVERTORS
from valenspy.input.esm_catalog_builder import CatalogBuilder
from valenspy.input.manager import InputManager


class _StubConvertor:
    """Just the attributes CatalogBuilder reads from an InputConverter."""

    raw_variables = {"T2M"}
    raw_variables_long_names = set()
    cordex_variables = {"tas"}

    def get_CORDEX_variable(self, raw_variable):
        return "tas"


@pytest.fixture(autouse=True)
def _restore_global_convertors():
    """If a test below fails because it leaked into INPUT_CONVERTORS, don't poison the rest of the session."""
    snapshot = dict(INPUT_CONVERTORS)
    yield
    INPUT_CONVERTORS.clear()
    INPUT_CONVERTORS.update(snapshot)


def _make_dataset(root, variable_id="tas", source_id="SRC"):
    root.mkdir(parents=True, exist_ok=True)
    (root / f"{variable_id}_2000.nc").touch()
    return {
        "root": str(root),
        "pattern": "<variable_id>_<time_period>.nc",
        "meta_data": {
            "source_id": source_id,
            "source_type": "RCM",
            "domain_id": "BE",
            "experiment_id": "historical",
            "version": "v1",
            "resolution": "1km",
            "frequency": "mon",
        },
    }


def test_builder_does_not_mutate_builtin_dataset_info(tmp_path, monkeypatch):
    monkeypatch.setitem(
        esm_catalog_builder.DATASET_PATHS,
        "testmachine",
        {"dsA": _make_dataset(tmp_path / "A", source_id="SRC_A")},
    )

    first = CatalogBuilder("testmachine")
    first.add_dataset("dsB", _make_dataset(tmp_path / "B", source_id="SRC_B"))

    assert list(esm_catalog_builder.DATASET_PATHS["testmachine"]) == ["dsA"]
    second = CatalogBuilder("testmachine")
    assert list(second.datasets_info) == ["dsA"]
    assert set(second.df["source_id"]) == {"SRC_A"}


def test_builder_does_not_mutate_callers_datasets_info(tmp_path):
    mine = {"dsA": _make_dataset(tmp_path / "A", source_id="SRC_A")}

    builder = CatalogBuilder("anything", datasets_info=mine)
    builder.add_dataset("dsB", _make_dataset(tmp_path / "B", source_id="SRC_B"))

    assert list(mine) == ["dsA"]


def test_managers_do_not_share_input_convertors(tmp_path):
    m1 = InputManager("m1", datasets_info={"dsA": _make_dataset(tmp_path / "A1")})
    m2 = InputManager("m2", datasets_info={"dsA": _make_dataset(tmp_path / "A2")})

    # A copy that still carries every built-in convertor, not an empty or shared dict.
    assert m1.input_convertors is not INPUT_CONVERTORS
    assert m1.input_convertors is not m2.input_convertors
    assert set(m1.input_convertors) == set(INPUT_CONVERTORS)

    m1.add_input_convertor("ONLY_FOR_M1", _StubConvertor())

    assert "ONLY_FOR_M1" in m1.input_convertors
    assert "ONLY_FOR_M1" not in m2.input_convertors
    assert "ONLY_FOR_M1" not in INPUT_CONVERTORS


def test_managers_do_not_share_xarray_kwargs(tmp_path):
    m1 = InputManager("m1", datasets_info={"dsA": _make_dataset(tmp_path / "A1")})
    m1.xarray_open_kwargs["chunks"] = {"time": 1}
    m1.xarray_combine_by_coords_kwargs["combine_attrs"] = "override"

    # Created after m1's mutation - the shared-default-dict case.
    m2 = InputManager("m2", datasets_info={"dsA": _make_dataset(tmp_path / "A2")})
    assert m2.xarray_open_kwargs == {}
    assert m2.xarray_combine_by_coords_kwargs == {}

    # A dict the caller passes in is copied too, not aliased.
    opts = {"chunks": {}}
    m3 = InputManager(
        "m3",
        datasets_info={"dsA": _make_dataset(tmp_path / "A3")},
        xarray_open_kwargs=opts,
    )
    m3.xarray_open_kwargs["engine"] = "h5netcdf"
    assert opts == {"chunks": {}}


def test_update_catalog_uses_managers_own_convertors(tmp_path):
    """Per-instance convertors must still reach the catalog builder.

    The builder used to read the module-global INPUT_CONVERTORS, and add_input_convertor
    followed by update_catalog_from_dataset_info only worked because the manager's dict
    *was* that global - copying it in the manager alone would silently drop the new
    dataset's variable-name translation.
    """
    manager = InputManager(
        "m", datasets_info={"dsA": _make_dataset(tmp_path / "A", source_id="SRC_A")}
    )
    dsB = _make_dataset(tmp_path / "B", variable_id="T2M", source_id="SRC_B")

    manager.add_input_convertor("dsB", _StubConvertor())
    manager.update_catalog_from_dataset_info(
        "dsB", dsB["root"], dsB["pattern"], dsB["meta_data"]
    )

    row = manager.catalog_builder.df.query("source_id == 'SRC_B'").iloc[0]
    assert row["raw_variable_id"] == ["T2M"]
    assert row["variable_id"] == ["tas"]


def test_custom_input_convertors_used_for_catalog_build(tmp_path):
    """input_convertors passed to InputManager used to affect only `preprocess`, not the catalog."""
    manager = InputManager(
        "m",
        datasets_info={
            "dsT": _make_dataset(tmp_path / "T", variable_id="T2M", source_id="SRC_T")
        },
        input_convertors={"dsT": _StubConvertor()},
    )

    row = manager.catalog_builder.df.iloc[0]
    assert row["raw_variable_id"] == ["T2M"]
    assert row["variable_id"] == ["tas"]
