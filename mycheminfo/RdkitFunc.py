from rdkit import Chem
from rdkit.Chem import Descriptors

def get_mw_form_smiles(smile=None):
    """
    return molecular weight of a species
    smile(str): SMILE
    :return: mw(float)
    """
    mw = Descriptors.ExactMolWt(Chem.MolFromSmiles(smile))
    return mw

