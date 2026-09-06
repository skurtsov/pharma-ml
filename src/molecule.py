from rdkit import Chem
from rdkit.Chem import Descriptors

smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

molecule = Chem.MolFromSmiles(smiles)
molecular_weight = Descriptors.MolWt(molecule)
log_p = Descriptors.MolLogP(molecule)
h_bond_donors = Descriptors.NumHDonors(molecule)
h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
print("Molecular Weight:", round(molecular_weight, 3))
print("LogP:", round(log_p, 4))
print("H-Bond Donors:", h_bond_donors)      
print("H-Bond Acceptors:", h_bond_acceptors)    

num_atoms = molecule.GetNumAtoms()
print(num_atoms)

for atom in molecule.GetAtoms():
    print(
        atom.GetSymbol(),
        atom.GetIdx(),
        "H:",
        atom.GetTotalNumHs()
    )

for bond in molecule.GetBonds():
    print(
        bond.GetBeginAtomIdx(),
        "->",
        bond.GetEndAtomIdx(),
        bond.GetBondType()
    )