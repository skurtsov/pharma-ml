from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)


def generate_fingerprint(smiles):
    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    return generator.GetFingerprint(molecule)