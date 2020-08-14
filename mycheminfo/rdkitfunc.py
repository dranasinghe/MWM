from rdkit import Chem
from rdkit.Chem import Descriptors
from pubchempy import get_compounds
import numpy as np
import urllib

def get_mw_form_smiles(smile:str = None)->float:
    """
    return molecular weight of a species
    smile(str): SMILE
    :return: mw(float)
    """
    mw = Descriptors.ExactMolWt(Chem.MolFromSmiles(smile))
    return mw

#alternative to pubchempy we can use cactus API to convert chemical names to smiles.
def get_smiles_from_cactus(name:str = None):
    """ for given compound name try to retrive SMILES from cactus
    Args:
        name(str): chemical name or CAS number
    Returns:
        str: SMILES if sucessful else NAN
    """
    url = "https://cactus.nci.nih.gov/chemical/structure/{0}/smiles".format(urllib.parse.quote(name))
    try:
        f = urllib.request.urlopen(url, timeout=10)
    except:
        print('Invalid identifier for NCI resolver.')
        return np.nan
    smiles = f.read().decode('utf-8')
    return smiles

def get_simles(name:str = None)->str:
    """ convert compound name or CAS number to SMILES. Uses pubchempy or cactus.
    Args:
        name(str): chemical name or CAS number
    Returns:
        str: SMILES if sucessful else NAN
    """
    try:
        smiles = [compound.isomeric_smiles for compound in get_compounds(name, 'name')]
        if len(smiles)>0:
            return smiles[0]
        # if pubchempy fail
        else:
            return get_smiles_from_cactus(name)
    except:
        # incase, try again
        return get_smiles_from_cactus(name)
