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


def _decode_attr(val):
    """Convert an HDF5 attribute value to a plain Python type."""
    if isinstance(val, np.integer):
        return int(val)
    if isinstance(val, np.floating):
        return float(val)
    if isinstance(val, np.bool_):
        return bool(val)
    if isinstance(val, (bytes, np.bytes_)):
        return val.decode("utf-8", errors="replace")
    if isinstance(val, np.ndarray):
        return val.tolist()
    return val


def _load_legacy_group(group):
    """Load a deepdish/flammkuchen-style Group (tuple or dict stored as attrs)."""
    title = group.attrs.get("TITLE", b"")
    if isinstance(title, (bytes, np.bytes_)):
        title = title.decode("utf-8", errors="replace")

    # Tuple stored as i0, i1, ... with TITLE like "tuple:N"
    if title.startswith("tuple:"):
        n = int(title.split(":")[1])
        items = []
        for i in range(n):
            key = f"i{i}"
            if key in group.attrs:
                items.append(_decode_attr(group.attrs[key]))
            elif key in group:
                items.append(group[key][()])
            else:
                items.append(None)
        return tuple(items)

    return None


def load(path):
    """Load a dict from an HDF5 file.

    Supports both the current h5io format (_scalar_data JSON attribute) and
    the legacy deepdish/flammkuchen format (scalars in file-level attrs,
    tuples as Groups with i0/i1/... attrs).
    """
    result = {}
    with h5py.File(path, "r") as f:
        for key in f.keys():
            obj = f[key]
            if isinstance(obj, h5py.Dataset):
                result[key] = obj[()]
            elif isinstance(obj, h5py.Group):
                val = _load_legacy_group(obj)
                if val is not None:
                    result[key] = val

        if "_scalar_data" in f.attrs:
            scalars = json.loads(f.attrs["_scalar_data"])
            result.update(scalars)
        else:
            # Legacy deepdish format: scalars in file-level attributes
            skip = {
                "CLASS",
                "DEEPDISH_IO_VERSION",
                "PYTABLES_FORMAT_VERSION",
                "TITLE",
                "VERSION",
                "_scalar_data",
            }
            for attr_name, attr_val in f.attrs.items():
                if attr_name in skip:
                    continue
                result[attr_name] = _decode_attr(attr_val)
    return result
