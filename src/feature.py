import numpy as np
from rdkit import DataStructs
from molecule import generate_fingerprint
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report, roc_auc_score
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import defaultdict

def fingerprint_to_array(smiles):
    fingerprint = generate_fingerprint(smiles)

    array = np.zeros((2048,), dtype=np.int8)

    DataStructs.ConvertToNumpyArray(
        fingerprint,
        array
    )

    return array

test_smiles = "CCO"

features = fingerprint_to_array(test_smiles)

print("Shape:", features.shape)
print("Active bits:", features.sum())
print(features[:20])
active_indices = np.where(features == 1)[0]

print("Active indices:", active_indices)


df = pd.read_csv("data/egfr_activities.csv")

X = np.array([
    fingerprint_to_array(smiles)
    for smiles in df["smiles"]
])

y = df["label"].to_numpy()

print("X shape:", X.shape)
print("y shape:", y.shape)

 #Data separation

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y
)

print("X_train:", X_train.shape)
print("X_test:", X_test.shape)
print("y_train:", y_train.shape)
print("y_test:", y_test.shape)

# Model training - Random Forest 

model = RandomForestClassifier(
    n_estimators=100,
    random_state=42
)

model.fit(X_train, y_train)

y_pred = model.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)

print("Accuracy:", accuracy)

##

y_proba = model.predict_proba(X_test)[:, 1]

print(classification_report(y_test, y_pred))
print("ROC-AUC:", roc_auc_score(y_test, y_proba))
#scaffold 
def get_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)

    return Chem.MolToSmiles(scaffold)
print(get_scaffold(df["smiles"].iloc[0]))
print(get_scaffold(df["smiles"].iloc[1]))
print(get_scaffold(df["smiles"].iloc[2]))
scaffold_groups = defaultdict(list)

for index, smiles in enumerate(df["smiles"]):
    scaffold = get_scaffold(smiles)
    scaffold_groups[scaffold].append(index)

print("Unique scaffolds:", len(scaffold_groups))

largest_groups = sorted(
    scaffold_groups.items(),
    key=lambda item: len(item[1]),
    reverse=True
)


##

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

X_train_scaffold = X[train_indices]
X_test_scaffold = X[test_indices]

y_train_scaffold = y[train_indices]
y_test_scaffold = y[test_indices]

print("Scaffold X_train:", X_train_scaffold.shape)
print("Scaffold X_test:", X_test_scaffold.shape)
print("Scaffold y_train:", y_train_scaffold.shape)
print("Scaffold y_test:", y_test_scaffold.shape)

scaffold_model = RandomForestClassifier(
    n_estimators=100,
    random_state=42
)

scaffold_model.fit(
    X_train_scaffold,
    y_train_scaffold
)

scaffold_pred = scaffold_model.predict(
    X_test_scaffold
)

scaffold_proba = scaffold_model.predict_proba(
    X_test_scaffold
)[:, 1]

print(
    "Scaffold Accuracy:",
    accuracy_score(
        y_test_scaffold,
        scaffold_pred
    )
)

print(
    classification_report(
        y_test_scaffold,
        scaffold_pred
    )
)

print(
    "Scaffold ROC-AUC:",
    roc_auc_score(
        y_test_scaffold,
        scaffold_proba
    )
)
train_scaffolds = {
    get_scaffold(df["smiles"].iloc[index])
    for index in train_indices
}

test_scaffolds = {
    get_scaffold(df["smiles"].iloc[index])
    for index in test_indices
}

overlap = train_scaffolds & test_scaffolds

print("Train scaffolds:", len(train_scaffolds))
print("Test scaffolds:", len(test_scaffolds))
print("Scaffold overlap:", len(overlap))

print(
    "Train ACTIVE ratio:",
    y_train_scaffold.mean()
)

print(
    "Test ACTIVE ratio:",
    y_test_scaffold.mean()
)