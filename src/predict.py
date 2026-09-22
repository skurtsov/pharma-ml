import joblib
from molecule import generate_fingerprint
from rdkit import DataStructs
import numpy as np
#Smiles to features function
def smiles_to_features(smiles):
    fingerprint = generate_fingerprint(smiles)

    array = np.zeros((2048,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fingerprint, array)

    return array.reshape(1, -1)

#Prediction function
def predict_activity(smiles):
    features = smiles_to_features(smiles)

    prediction = model.predict(features)[0]
    probability = model.predict_proba(features)[0][1]

    return {
        "smiles": smiles,
        "prediction": "ACTIVE" if prediction == 1 else "INACTIVE",
        "active_probability": float(probability),
    }

model = joblib.load(
    "models/egfr_random_forest.joblib"
)

print(type(model))
print("Trees:", model.n_estimators)
print("Features:", model.n_features_in_)
print("Classes:", model.classes_)

result = predict_activity("CCO")
print(result)