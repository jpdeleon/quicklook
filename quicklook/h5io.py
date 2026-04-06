"""Simple HDF5 dict save/load to replace flammkuchen."""

import json
import numpy as np
import h5py


class _NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy scalar types."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save(path, data):
    """Save a dict to an HDF5 file.

    Supports numpy arrays, scalars, strings, tuples, lists, and nested dicts.
    Non-array values are stored as JSON in an attribute.
    """
    arrays = {}
    attrs = {}

    for key, value in data.items():
        if isinstance(value, np.ndarray):
            arrays[key] = value
        else:
            try:
                json.dumps(value, cls=_NumpyEncoder)
                attrs[key] = value
            except (TypeError, ValueError):
                attrs[key] = str(value)

    with h5py.File(path, "w") as f:
        for key, arr in arrays.items():
            f.create_dataset(key, data=arr)
        f.attrs["_scalar_data"] = json.dumps(attrs, cls=_NumpyEncoder)


def load(path):
    """Load a dict from an HDF5 file.

    Returns a dict with numpy arrays and scalar values.
    """
    result = {}
    with h5py.File(path, "r") as f:
        for key in f.keys():
            result[key] = f[key][()]
        if "_scalar_data" in f.attrs:
            scalars = json.loads(f.attrs["_scalar_data"])
            result.update(scalars)
    return result
