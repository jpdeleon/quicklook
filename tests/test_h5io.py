"""Tests for quicklook.h5io — round-tripping dicts of arrays + scalars."""

import json
import numpy as np
import pytest
import h5py

from quicklook import h5io


def test_round_trip_arrays_and_scalars(tmp_path):
    data = {
        "flux": np.array([1.0, 2.0, 3.0]),
        "time": np.linspace(0, 10, 5),
        "period": 1.234,
        "sde": 12,
        "pipeline": "spoc",
        "use_star_priors": True,
    }
    path = tmp_path / "out.h5"
    h5io.save(path, data)
    out = h5io.load(path)

    np.testing.assert_array_equal(out["flux"], data["flux"])
    np.testing.assert_array_equal(out["time"], data["time"])
    assert out["period"] == pytest.approx(1.234)
    assert out["sde"] == 12
    assert out["pipeline"] == "spoc"
    assert out["use_star_priors"] is True


def test_numpy_scalars_are_jsonable(tmp_path):
    """np.int64/np.float64/np.bool_ must serialize via _NumpyEncoder."""
    data = {
        "n": np.int64(7),
        "x": np.float64(3.14),
        "flag": np.bool_(False),
    }
    path = tmp_path / "scalars.h5"
    h5io.save(path, data)
    out = h5io.load(path)
    assert out["n"] == 7 and isinstance(out["n"], int)
    assert out["x"] == pytest.approx(3.14)
    assert out["flag"] is False


def test_nested_dict_scalar_only(tmp_path):
    data = {"meta": {"pipeline": "tglc", "size": 50}}
    path = tmp_path / "nested.h5"
    h5io.save(path, data)
    out = h5io.load(path)
    assert out["meta"] == {"pipeline": "tglc", "size": 50}


def test_unserializable_falls_back_to_str(tmp_path):
    """Objects that aren't JSON-serializable are stored via str()."""

    class _Opaque:
        def __repr__(self):
            return "<opaque-thing>"

    data = {"weird": _Opaque(), "n": 1}
    path = tmp_path / "opaque.h5"
    h5io.save(path, data)
    out = h5io.load(path)
    assert out["weird"] == "<opaque-thing>"
    assert out["n"] == 1


def test_load_legacy_deepdish_attrs(tmp_path):
    """Files written before this rewrite stored scalars in file-level attrs
    (no ``_scalar_data`` key). load() must still parse them."""
    path = tmp_path / "legacy.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("flux", data=np.array([1.0, 2.0, 3.0]))
        f.attrs["period"] = 4.5
        f.attrs["pipeline"] = "spoc"
        # Skipped pytables/deepdish bookkeeping must not appear in result.
        f.attrs["TITLE"] = "ignore me"
        f.attrs["DEEPDISH_IO_VERSION"] = 7

    out = h5io.load(path)
    np.testing.assert_array_equal(out["flux"], np.array([1.0, 2.0, 3.0]))
    assert out["period"] == pytest.approx(4.5)
    assert out["pipeline"] == "spoc"
    assert "TITLE" not in out
    assert "DEEPDISH_IO_VERSION" not in out


def test_load_legacy_tuple_group(tmp_path):
    """Legacy deepdish tuples stored as Groups with i0/i1 attrs must
    deserialize back into a Python tuple."""
    path = tmp_path / "legacy_tuple.h5"
    with h5py.File(path, "w") as f:
        g = f.create_group("ephem")
        g.attrs["TITLE"] = "tuple:3"
        g.attrs["i0"] = 1.23
        g.attrs["i1"] = 4.56
        g.attrs["i2"] = 7.89

    out = h5io.load(path)
    assert out["ephem"] == (pytest.approx(1.23), pytest.approx(4.56), pytest.approx(7.89))


def test_save_overwrites_existing(tmp_path):
    """h5io.save uses mode='w' — second save must replace, not append."""
    path = tmp_path / "out.h5"
    h5io.save(path, {"a": np.array([1, 2, 3])})
    h5io.save(path, {"b": np.array([10, 20])})
    out = h5io.load(path)
    assert "a" not in out
    np.testing.assert_array_equal(out["b"], np.array([10, 20]))
