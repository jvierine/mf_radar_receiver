"""Small atomic HDF5 helpers for numerical radar products."""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np


def save_array(
    path: str | os.PathLike,
    values: np.ndarray,
    dataset: str = "data",
) -> None:
    """Atomically replace one HDF5 file containing a numerical dataset."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.create_dataset(
            dataset,
            data=np.asarray(values),
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
    os.replace(temporary, destination)


def load_array(
    path: str | os.PathLike,
    dataset: str = "data",
) -> np.ndarray:
    """Load one numerical dataset from an HDF5 file."""
    with h5py.File(path, "r") as handle:
        return np.asarray(handle[dataset])
