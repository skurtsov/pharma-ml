import numpy as np
from rdkit import DataStructs
from molecule import generate_fingerprint
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report, roc_auc_score


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