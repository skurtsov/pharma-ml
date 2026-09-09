from chembl_webresource_client.new_client import new_client

TARGET_ID = "CHEMBL203"


def get_activities():
    return new_client.activity.filter(
        target_chembl_id=TARGET_ID,
        standard_type="IC50",
        standard_relation="=",
        standard_units="nM",
        assay_type="B"
    )


def get_molecule_smiles(molecule_id):
    molecule_data = new_client.molecule.get(molecule_id)

    return molecule_data["molecule_structures"]["canonical_smiles"]