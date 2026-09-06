from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
#Aspirin
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
#Ibuprofen
smiles_2 = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
#Paracetamol
smiles_3 = "CC(=O)NC1=CC=C(C=C1)O"
##

    
generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)
def generate_fingerprint(smiles):
    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    return generator.GetFingerprint(molecule)

fingerprint = generate_fingerprint(smiles)
fingerprint_2 = generate_fingerprint(smiles_2)
fingerprint_3 = generate_fingerprint(smiles_3)

aspirin_ibuprofen = DataStructs.TanimotoSimilarity(
    fingerprint,
    fingerprint_2
)

aspirin_paracetamol = DataStructs.TanimotoSimilarity(
    fingerprint,
    fingerprint_3
)

ibuprofen_paracetamol = DataStructs.TanimotoSimilarity(
    fingerprint_2,
    fingerprint_3
)
print("Aspirin vs Ibuprofen:", round(aspirin_ibuprofen, 3))
print("Aspirin vs Paracetamol:", round(aspirin_paracetamol, 3))
print("Ibuprofen vs Paracetamol:", round(ibuprofen_paracetamol, 3))