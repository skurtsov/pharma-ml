from chembl_webresource_client.new_client import new_client
from molecule import generate_fingerprint
TARGET_ID = "CHEMBL203"

activity = new_client.activity

# Activity to record
def activity_to_record(activity):
    molecule_id = activity["molecule_chembl_id"]
    molecule_data = new_client.molecule.get(molecule_id)

    smiles = molecule_data["molecule_structures"]["canonical_smiles"]

    ic50 = float(activity["standard_value"])
    label = 1 if ic50 <= 1000 else 0

    return {
        "molecule_chembl_id": molecule_id,
        "smiles": smiles,
        "ic50": ic50,
        "label": label
    }
#Activities filter
activities = new_client.activity.filter(
    target_chembl_id=TARGET_ID,
    standard_type="IC50",
    standard_relation="=",
    standard_units="nM",
    assay_type="B"
)
print("Activities:", len(activities))

molecule = new_client.molecule

first_activity = activities[0]

molecule_id = first_activity["molecule_chembl_id"]
molecule_data = molecule.get(molecule_id)

print("ID:", molecule_id)
#SMILES+fingerprint
smiles = molecule_data["molecule_structures"]["canonical_smiles"]
print(
    "SMILES:",
    smiles
)
fingerprint = generate_fingerprint(smiles)

print("Fingerprint size:", len(fingerprint))
print("Active bits:", fingerprint.GetNumOnBits())
#IC50 label
ic50 = float(first_activity["standard_value"])
label = 1 if ic50 <= 1000 else 0

print("IC50:", ic50)
print("Label:", label)