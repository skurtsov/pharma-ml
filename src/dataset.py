from chembl import get_activities, get_molecule_smiles
import csv
from statistics import median

def activity_to_record(activity):
    molecule_id = activity["molecule_chembl_id"]
    smiles = get_molecule_smiles(molecule_id)

    ic50 = float(activity["standard_value"])
    label = 1 if ic50 <= 1000 else 0

    return {
        "molecule_chembl_id": molecule_id,
        "smiles": smiles,
        "ic50": ic50,
        "label": label
    }


activities = get_activities()
first_activity = activities[0]

print(first_activity.keys())
records = []

for activity in activities[:500]:
    try:
        record = activity_to_record(activity)
        records.append(record)
    except Exception as error:
        print("Skipped:", activity["molecule_chembl_id"], error)

print("Records:", len(records))


def aggregate_records(records):
    grouped_records = {}

    # 1. Group records by molecule_id
    for record in records:
        molecule_id = record["molecule_chembl_id"]

        if molecule_id not in grouped_records:
            grouped_records[molecule_id] = []

        grouped_records[molecule_id].append(record)

    # 2. Create one record per molecule
    aggregated_records = []

    for molecule_id, molecule_records in grouped_records.items():
        ic50_values = [
            record["ic50"]
            for record in molecule_records
        ]

        median_ic50 = median(ic50_values)
        label = 1 if median_ic50 <= 1000 else 0
        measurements_count = len(molecule_records)
        aggregated_record = {
            "molecule_chembl_id": molecule_id,
            "smiles": molecule_records[0]["smiles"],
            "ic50": median_ic50,
            "label": label,
            "measurements_count": measurements_count
        }

        aggregated_records.append(aggregated_record)

    return aggregated_records
    
records = aggregate_records(records)
    
active_count = sum(record["label"] == 1 for record in records)
inactive_count = sum(record["label"] == 0 for record in records)

print("Active:", active_count)
print("Inactive:", inactive_count)
# dedup
molecule_ids = [record["molecule_chembl_id"] for record in records]

unique_molecules = set(molecule_ids)

print("Unique molecules:", len(unique_molecules))
print("Duplicate records:", len(records) - len(unique_molecules))
#Results to CSV
with open("data/egfr_activities.csv", "w", newline="") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "molecule_chembl_id",
            "smiles",
            "ic50",
            "label",
            "measurements_count"
        ]
    )

    writer.writeheader()
    writer.writerows(records)