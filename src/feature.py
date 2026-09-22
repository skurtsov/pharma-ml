import numpy as np
from rdkit import DataStructs

from molecule import generate_fingerprint


def fingerprint_to_array(smiles):
    fingerprint = generate_fingerprint(smiles)

    array = np.zeros((2048,), dtype=np.int8)

    DataStructs.ConvertToNumpyArray(
        fingerprint,
        array
    )

    return array