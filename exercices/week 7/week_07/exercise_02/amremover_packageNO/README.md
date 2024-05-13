# amremover_package

amremover_package is a Python package for removing atom mapping and canonicalizing reaction SMILES.

## Installation

You can install amremover_package using pip: pip install amremover_package


## Usage

```python
from amremover_package.chemutils.utils import remove_atom_mapping_and_canonicalize_rxn_smiles

rxn_smiles_with_atom_mapping = '[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1.C1CCCCC1>[OH-].[OH-].[Pd+2].CCO>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1'
rxn_smiles_without_atom_mapping = remove_atom_mapping_and_canonicalize_rxn_smiles(rxn_smiles_with_atom_mapping)

print(f"RXN SMILES without atom mapping: {rxn_smiles_without_atom_mapping}")

Contributing

Contributions are welcome! Please submit a pull request if you have any improvements or bug fixes.

License

This package is released under the MIT License. 
