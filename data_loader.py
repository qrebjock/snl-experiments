import numpy as np
import os
import scipy.io


def gather_mat_files(folder):
    """
    Loads and merges .mat files from a folder, filtering out rows with NaNs.

    Returns a dictionary with concatenated data for each variable.
    """

    keys = [
        "dimgt",
        "dimopt",
        "n",
        "ps",
        "noiselevel",
        "ground_truths",
        "num_ps",
        "chunk_size",
        "connectedness",
        "costs_ground_truths",
        "costs_ground_truths_reopt",
        "costs",
        "gradnorms",
        "mineigvals",
        "niters",
        "normprocrustdiff",
        "costs_proj",
        "gradnorms_proj",
        "mineigvals_proj",
        "normprocrustdiff_proj",
        "costs_reopt",
        "gradnorms_reopt",
        "mineigvals_reopt",
        "niters_reopt",
        "normprocrustdiff_reopt",
        "normprocrustdiffysreopt_reopt",
    ]

    files = [file for file in os.listdir(folder) if file.endswith(".mat")]
    mats = [
        scipy.io.loadmat(os.path.join(folder, file), variable_names=keys)
        for file in files
    ]

    # Sanity checks
    for field in ["dimgt", "dimopt", "n", "num_ps"]:
        unique_vals = np.unique([mat[field] for mat in mats])
        if unique_vals.size != 1:
            raise ValueError(f"Inconsistent '{field}' across files: {unique_vals}")

    finalmat = dict()

    for key in keys:
        try:
            concatenated = np.concatenate([mat[key] for mat in mats])
            keeprows = ~np.any(
                np.isnan(concatenated), axis=tuple(range(1, concatenated.ndim))
            )
            finalmat[key] = concatenated[keeprows, :]
        except (KeyError, ValueError) as e:
            print(f"Skipping key '{key}': {e}")

    # Copy static metadata from first file
    for key in ["dimgt", "dimopt", "n", "ps", "num_ps"]:
        val = mats[0][key]
        finalmat[key] = val.item() if val.size == 1 else np.ravel(val)

    return finalmat
