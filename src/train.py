import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import defaultdict
from feature import fingerprint_to_array


df = pd.read_csv("data/egfr_activities.csv")

X = np.array([
    fingerprint_to_array(smiles)
    for smiles in df["smiles"]
])

y = df["label"].to_numpy()

print("X shape:", X.shape)
print("y shape:", y.shape)

def get_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)

    return Chem.MolToSmiles(scaffold)


def scaffold_split(df, test_size=0.2):
    scaffold_groups = defaultdict(list)

    for index, smiles in enumerate(df["smiles"]):
        scaffold = get_scaffold(smiles)
        scaffold_groups[scaffold].append(index)

    groups = sorted(
        scaffold_groups.values(),
        key=len,
        reverse=True
    )

    train_indices = []
    test_indices = []

    target_train_size = int(len(df) * (1 - test_size))

    for group in groups:
        if len(train_indices) + len(group) <= target_train_size:
            train_indices.extend(group)
        else:
            test_indices.extend(group)

    return train_indices, test_indices

train_indices, test_indices = scaffold_split(df)

X_train = X[train_indices]
X_test = X[test_indices]

y_train = y[train_indices]
y_test = y[test_indices]

print("X_train:", X_train.shape)
print("X_test:", X_test.shape)
print("y_train:", y_train.shape)
print("y_test:", y_test.shape)