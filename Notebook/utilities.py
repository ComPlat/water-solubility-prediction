from rdkit import Chem
from rdkit.Chem import Descriptors,rdMolDescriptors,EState  # For molecular descriptors in RDKit
import pandas as pd
import numpy as np 
from rdkit import DataStructs

def calculate_aromatic_proportion(smiles):
    # Parse SMILES string and generate molecular representation
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None  # Invalid SMILES
    
    # Identify aromatic atoms
    aromatic_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    
    # Calculate aromatic proportion
    total_atoms = mol.GetNumAtoms()
    aromatic_proportion = len(aromatic_atoms) / total_atoms
    
    return aromatic_proportion
####  function for handling single smiles 
def create_feature_multiple(SMILES):
        df125_new=calculate_rdkit_features1(SMILES)
        df125 = df125_new.iloc[:, :125]
        df38=generate_features38(SMILES)
        df128=fingerprint1(SMILES,2,128)
        df7=get_functional_groups1(SMILES)
        return pd.concat([df125, df128, df7, df38], axis=1)
### Function to create descriptors for single smiles
def create_feature_single(SMILES):
    df125_new = calculate_rdkit_features(SMILES)
    df125 = df125_new.iloc[:, :125]
    df38 = generate_features_single38(SMILES)
    df128 = fingerprint1(SMILES, 2, 128)
    df7 = get_functional_groups1(SMILES)

    return pd.concat([df125, df128, df7, df38], axis=1)
    

    return pd.concat([df125, df128, df7, df38], axis=1)
            
def calculate_rdkit_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string.")
    
    # List of RDKit descriptors
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    descriptor_funcs = [desc[1] for desc in Descriptors._descList]
    
    # Compute descriptors
    features = [func(mol) for func in descriptor_funcs]
    
    # Create DataFrame
    df = pd.DataFrame([features], columns=descriptor_names)
    return df
#### Function for list of the smiles 
def calculate_rdkit_features1(smiles_list):
    """
    Calculate RDKit descriptors for a list of SMILES strings.

    Parameters:
        smiles_list (list of str): List of SMILES strings.

    Returns:
        pd.DataFrame: A DataFrame where each row corresponds to a SMILES string
                      and each column corresponds to a descriptor.
    """
    # List of RDKit descriptors
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    descriptor_funcs = [desc[1] for desc in Descriptors._descList]

    data = []
    invalid_smiles = []

    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            invalid_smiles.append((i, smiles))
            continue
        
        # Compute descriptors
        features = [func(mol) for func in descriptor_funcs]
        data.append(features)

    # Create DataFrame
    df = pd.DataFrame(data, columns=descriptor_names, index=[i for i, smiles in enumerate(smiles_list) if smiles not in dict(invalid_smiles)])

    # Log invalid SMILES strings
    if invalid_smiles:
        print(f"Invalid SMILES found at indices: {[idx for idx, sm in invalid_smiles]}\nInvalid SMILES: {[sm for idx, sm in invalid_smiles]}")

    return df
### function for single smiles and list of the smiles 
def fingerprint1(smiles, r, n): 
    # Check if input is a single SMILES string or a list
    if isinstance(smiles, str):
        smiles = [smiles]  # Convert single SMILES string to a list

    # Convert SMILES strings to RDKit molecules
    mols = [Chem.MolFromSmiles(SMILES_string) for SMILES_string in smiles]

    # Compute fingerprints
    fingerprints_array = []
    for m in mols:
        bi = {}
        fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius=r, nBits=n, bitInfo=bi)
        array = np.zeros((1,), dtype=int)
        DataStructs.ConvertToNumpyArray(fingerprint, array)
        fingerprints_array.append(array)

    # Convert to a DataFrame with bits as columns
    fingerprints_df = pd.DataFrame(fingerprints_array)

    # Return the DataFrame for a single SMILES
    return fingerprints_df
def get_charges(smiles):
    if '+' in smiles:
        return 1
    elif '-' in smiles:
        return -1
    else:
        return 0

def get_many_double_bonds(smiles):
    mol = Chem.MolFromSmiles(smiles)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    return 1 if double_bond_count > 4 else 0

def get_atom_degrees(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bonds = mol.GetBonds()
    sum_degree_vector = np.zeros(7)
    for bond in bonds:
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atom_degree_vector = np.array([1 if atom1.GetDegree() == d else 0 for d in range(7)])
        sum_degree_vector += atom_degree_vector
        atom_degree_vector = np.array([1 if atom2.GetDegree() == d else 0 for d in range(7)])
        sum_degree_vector += atom_degree_vector
    return sum_degree_vector.astype(int)

def get_atom_valences(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bonds = mol.GetBonds()
    sum_valence_vector = np.zeros(7)
    for bond in bonds:
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atom_valence_vector = np.array([1 if atom1.GetTotalValence() == v else 0 for v in range(7)])
        sum_valence_vector += atom_valence_vector
        atom_valence_vector = np.array([1 if atom2.GetTotalValence() == v else 0 for v in range(7)])
        sum_valence_vector += atom_valence_vector
    return sum_valence_vector.astype(int)

def get_atom_hybridization(smiles):
    hybridizations = [
        Chem.rdchem.HybridizationType.S,
        Chem.rdchem.HybridizationType.SP, 
        Chem.rdchem.HybridizationType.SP2, 
        Chem.rdchem.HybridizationType.SP3, 
        Chem.rdchem.HybridizationType.SP3D, 
        Chem.rdchem.HybridizationType.SP3D2, 
        Chem.rdchem.HybridizationType.UNSPECIFIED]
    
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bonds = mol.GetBonds()
    sum_hybrid_vector = np.zeros(7)
    for bond in bonds:
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atom_hybrid_vector = np.array([1 if atom1.GetHybridization() == h else 0 for h in hybridizations])
        sum_hybrid_vector += atom_hybrid_vector
        atom_hybrid_vector = np.array([1 if atom2.GetHybridization() == h else 0 for h in hybridizations])
        sum_hybrid_vector += atom_hybrid_vector
    return sum_hybrid_vector.astype(int)

def get_aromatic_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bonds = mol.GetBonds()
    sum_aromatic_vector = np.zeros(1)
    for bond in bonds:
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atom_aromatic_vector = np.array([1 if atom1.GetIsAromatic() else 0])
        sum_aromatic_vector += atom_aromatic_vector
        atom_aromatic_vector = np.array([1 if atom2.GetIsAromatic() else 0])
        sum_aromatic_vector += atom_aromatic_vector
    return sum_aromatic_vector.astype(int)[0]

def get_bond_types(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bonds = mol.GetBonds()
    sum_bond_type_vector = np.zeros(5)
    for bond in bonds:
        bond_type = bond.GetBondType().name
        bond_type_vector = np.array([1 if t == bond_type else 0 for t in ['SINGLE', 'DOUBLE', 'TRIPLE', 'AROMATIC', 'ZERO']])
        sum_bond_type_vector += bond_type_vector
    return sum_bond_type_vector.astype(int)

def is_conjugated(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bonds = mol.GetBonds()
    sum_conjugated = np.zeros(1)
    for bond in bonds:
        is_conjugated = 1 if bond.GetIsConjugated() else 0
        conjugation_vector = np.array([is_conjugated])
        sum_conjugated += conjugation_vector
    return sum_conjugated.astype(int)[0]

def get_bonds_in_ring(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return len(Chem.GetSymmSSSR(mol))

def get_bond_chirality(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bonds = mol.GetBonds()
    sum_chirality = np.zeros(4)
    for bond in bonds:
        chirality = bond.GetStereo()
        chirality_vector = np.array([1 if chirality == c else 0 for c in [Chem.rdchem.BondStereo.STEREONONE,
                                                                          Chem.rdchem.BondStereo.STEREOANY,
                                                                          Chem.rdchem.BondStereo.STEREOZ,
                                                                          Chem.rdchem.BondStereo.STEREOE]])
        sum_chirality += chirality_vector
    return sum_chirality.astype(int)

def get_n_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return mol.GetNumAtoms()

def get_n_bonds(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return mol.GetNumBonds()

def get_n_rings(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return len(Chem.GetSymmSSSR(mol))
### Funntion for creating 38 features for single smiles 
def generate_features_single38(smiles):
    columns = [
        'charge', 'many_double_bonds', 'atoms_degree_0', 'atoms_degree_1',
        'atoms_degree_2', 'atoms_degree_3', 'atoms_degree_4', 'atoms_degree_5',
        'atoms_degree_6', 'atoms_valence_0', 'atoms_valence_1', 'atoms_valence_2',
        'atoms_valence_3', 'atoms_valence_4', 'atoms_valence_5', 'atoms_valence_6',
        'atom_hybridization_S', 'atom_hybridization_SP', 'atom_hybridization_SP2',
        'atom_hybridization_SP3', 'atom_hybridization_SP3D', 'atom_hybridization_SP3D2',
        'atom_hybridization_UNSPECIFIED', 'aromatic_atoms', 'single_bonds', 'double_bonds',
        'triple_bonds', 'aromatic_bonds', 'zero_bonds', 'conjugated_bonds', 'bonds_in_ring',
        'chirality_none', 'chirality_any', 'chirality_z', 'chirality_e', 'n_atoms',
        'n_bonds', 'n_rings'
    ]

    # Calculate features for the single SMILES string
    charge = get_charges(smiles)
    many_double_bonds = get_many_double_bonds(smiles)
    atom_degrees = get_atom_degrees(smiles).tolist()
    atom_valences = get_atom_valences(smiles).tolist()
    atom_hybridizations = get_atom_hybridization(smiles).tolist()
    aromatic_atoms = get_aromatic_atoms(smiles)
    bond_types = get_bond_types(smiles).tolist()
    conjugated_bonds = is_conjugated(smiles)
    bonds_in_ring = get_bonds_in_ring(smiles)
    bond_chirality = get_bond_chirality(smiles).tolist()
    n_atoms = get_n_atoms(smiles)
    n_bonds = get_n_bonds(smiles)
    n_rings = get_n_rings(smiles)

    # Combine all features into a single row
    feature_row = [
        charge, many_double_bonds, atom_degrees[0], atom_degrees[1],
        atom_degrees[2], atom_degrees[3], atom_degrees[4], atom_degrees[5],
        atom_degrees[6], atom_valences[0], atom_valences[1], atom_valences[2],
        atom_valences[3], atom_valences[4], atom_valences[5], atom_valences[6],
        atom_hybridizations[0], atom_hybridizations[1], atom_hybridizations[2],
        atom_hybridizations[3], atom_hybridizations[4], atom_hybridizations[5],
        atom_hybridizations[6], aromatic_atoms, bond_types[0], bond_types[1],
        bond_types[2], bond_types[3], bond_types[4], conjugated_bonds,
        bonds_in_ring, bond_chirality[0], bond_chirality[1], bond_chirality[2],
        bond_chirality[3], n_atoms, n_bonds, n_rings
    ]

    # Convert the feature row into a DataFrame
    return pd.DataFrame([feature_row], columns=columns)
#### function for creating the 38 features for list of the smiles 
def generate_features38(smiles_list):
    columns = [
        'charge', 'many_double_bonds', 'atoms_degree_0', 'atoms_degree_1',
        'atoms_degree_2', 'atoms_degree_3', 'atoms_degree_4', 'atoms_degree_5',
        'atoms_degree_6', 'atoms_valence_0', 'atoms_valence_1', 'atoms_valence_2',
        'atoms_valence_3', 'atoms_valence_4', 'atoms_valence_5', 'atoms_valence_6',
        'atom_hybridization_S', 'atom_hybridization_SP', 'atom_hybridization_SP2',
        'atom_hybridization_SP3', 'atom_hybridization_SP3D', 'atom_hybridization_SP3D2',
        'atom_hybridization_UNSPECIFIED', 'aromatic_atoms', 'single_bonds', 'double_bonds',
        'triple_bonds', 'aromatic_bonds', 'zero_bonds', 'conjugated_bonds', 'bonds_in_ring',
        'chirality_none', 'chirality_any', 'chirality_z', 'chirality_e', 'n_atoms',
        'n_bonds', 'n_rings'
    ]
    
    features = []
    for smiles in smiles_list:
        charge = get_charges(smiles)
        many_double_bonds = get_many_double_bonds(smiles)
        atom_degrees = get_atom_degrees(smiles).tolist()
        atom_valences = get_atom_valences(smiles).tolist()
        atom_hybridizations = get_atom_hybridization(smiles).tolist()
        aromatic_atoms = get_aromatic_atoms(smiles)
        bond_types = get_bond_types(smiles).tolist()
        conjugated_bonds = is_conjugated(smiles)
        bonds_in_ring = get_bonds_in_ring(smiles)
        bond_chirality = get_bond_chirality(smiles).tolist()
        n_atoms = get_n_atoms(smiles)
        n_bonds = get_n_bonds(smiles)
        n_rings = get_n_rings(smiles)

        feature_row = [
            charge, many_double_bonds, atom_degrees[0], atom_degrees[1],
            atom_degrees[2], atom_degrees[3], atom_degrees[4], atom_degrees[5],
            atom_degrees[6], atom_valences[0], atom_valences[1], atom_valences[2],
            atom_valences[3], atom_valences[4], atom_valences[5], atom_valences[6],
            atom_hybridizations[0], atom_hybridizations[1], atom_hybridizations[2],
            atom_hybridizations[3], atom_hybridizations[4], atom_hybridizations[5],
            atom_hybridizations[6], aromatic_atoms, bond_types[0], bond_types[1],
            bond_types[2], bond_types[3], bond_types[4], conjugated_bonds,
            bonds_in_ring, bond_chirality[0], bond_chirality[1], bond_chirality[2],
            bond_chirality[3], n_atoms, n_bonds, n_rings
        ]
        
        features.append(feature_row)
    
    return pd.DataFrame(features, columns=columns)

### Function for 7 groups for single smiles and list of the smiles 
def get_functional_groups1(smiles):
    from rdkit import Chem
    import pandas as pd

    functional_groups = {
        # Polar functional groups
        'Hydroxyl Group': '[OH]',
        'Carbonyl Group': 'C=O',
        'Amide Group': 'C(=O)N',
        'Carboxyl Group': 'C(=O)[OH]',
        # Non-polar functional groups
        'Alkyl': '[R]', 
        'Aromatic Rings': 'c',
        'Alkene': 'C=C'
    }
    
    # Check if input is a single SMILES string
    if isinstance(smiles, str):
        smiles = [smiles]  # Wrap it in a list for consistency
    
    results = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        if mol:  # Ensure valid SMILES
            fg_presence = {fg: 1 if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) else 0 for fg, smarts in functional_groups.items()}
        else:
            fg_presence = {fg: 0 for fg in functional_groups}  # Default to 0 if SMILES is invalid
        #fg_presence['SMILES'] = s
        results.append(fg_presence)
    
    # Convert results to a DataFrame
    data = pd.DataFrame(results)
    
    return data
### for single smiles 
def calc_mol_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolWt(mol)
    return None
#### for smiles list 
def calc_mol_weight1(smiles_list):
    """
    Calculate molecular weights for a list of SMILES strings.

    Parameters:
        smiles_list (list of str): List of SMILES strings.

    Returns:
        list: Molecular weights for valid SMILES strings, or None for invalid ones.
    """
    return [
        Descriptors.MolWt(Chem.MolFromSmiles(smiles)) if Chem.MolFromSmiles(smiles) else None
        for smiles in smiles_list
    ]
def calculate_descriptors_sorkun(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None  # Skip invalid molecules
    saturated_atoms = sum(
        1 for atom in mol.GetAtoms() if atom.GetTotalDegree() == 1 and not atom.GetIsAromatic()
    )
    spiro_atoms = set()
    ring_info = mol.GetRingInfo()
    spiro_atoms = set()
    ring_info = mol.GetRingInfo()
    bridgehead_atoms = set()
    ring_info = mol.GetRingInfo()
    zagreb_index_1 = sum([atom.GetDegree()**2 for atom in mol.GetAtoms()])
    # Iterate over all the rings to detect bridgehead atoms
    for ring in ring_info.AtomRings():
        # If the molecule is bicyclic, consider the atoms at the junction
        if len(ring) == 2:  # Typically for bicyclic systems
            bridgehead_atoms.add(ring[0])
            bridgehead_atoms.add(ring[1])
    for ring in ring_info.AtomRings():
        for atom_idx in ring:
            # Check if the atom is part of exactly two rings
            if sum(1 for other_ring in ring_info.AtomRings() if atom_idx in other_ring) == 2:
                spiro_atoms.add(atom_idx)

    for ring in ring_info.AtomRings():
        for atom_idx in ring:
            # Check if the atom is part of exactly two rings
            if sum(1 for other_ring in ring_info.AtomRings() if atom_idx in other_ring) == 2:
                spiro_atoms.add(atom_idx)
    lipinski_hbd = rdMolDescriptors.CalcNumHBD(mol)

    lipinski_hba = rdMolDescriptors.CalcNumHBA(mol)
    bridgehead_atoms = set()
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        for atom_idx in ring:
            # Check if the atom is part of more than one ring
            if sum(1 for other_ring in ring_info.AtomRings() if atom_idx in other_ring) > 1:
                bridgehead_atoms.add(atom_idx)
    estate_indices = EState.EStateIndices(mol)


    heterocycles = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() not in [6] for idx in ring):  # 6 is carbon (C)
            heterocycles += 1
    num_rings = mol.GetRingInfo().NumRings() 
    descriptors = {
        "nHeavyAtom": mol.GetNumHeavyAtoms(),
        "nHBAcc": rdMolDescriptors.CalcNumHBA(mol),
        "nHBDon": rdMolDescriptors.CalcNumHBD(mol),
        "nRot": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "nBonds": mol.GetNumBonds(),
        "nBondsO": sum(1 for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 8),
        "nBondsS": sum(1 for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 16),
        "nBondsD": sum(1 for bond in mol.GetBonds() if bond.GetBondType().name == "DOUBLE"),
        "TopoPSA(NO)": rdMolDescriptors.CalcTPSA(mol, includeSandP=False),
        "TopoPSA": rdMolDescriptors.CalcTPSA(mol),
        "LabuteASA": rdMolDescriptors.CalcLabuteASA(mol),
        "bpol": rdMolDescriptors.CalcCrippenDescriptors(mol)[1],
        "nAcid": rdMolDescriptors.CalcNumLipinskiHBA(mol),
        "nBase": rdMolDescriptors.CalcNumLipinskiHBD(mol),
        "ECIndex": estate_indices[0],
        "GGI1": Descriptors.FpDensityMorgan1(mol),
        "SLogP": Descriptors.MolLogP(mol),
        "SMR": Descriptors.MolMR(mol),
        "BertzCT": Descriptors.BertzCT(mol),
        "BalabanJ": Descriptors.BalabanJ(mol),
        "WPol": Descriptors.MolWt(mol),
        "Zagreb1": zagreb_index_1,
       # "ABCGG": abcgg_value,
        "nHRing": Descriptors.NumAromaticRings(mol),
        "naHRing": Descriptors.NumAliphaticRings(mol),
        "NsCH3": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "NssCH2": Descriptors.NumRotatableBonds(mol),
        "NaaCH": Descriptors.NumAliphaticCarbocycles(mol),
        "NaaaC": Descriptors.NumAliphaticHeterocycles(mol),
        "NssssC": Descriptors.NumAromaticHeterocycles(mol),
        "SsCH3": Descriptors.NumAromaticCarbocycles(mol),
        "SdCH2": Descriptors.NumSaturatedCarbocycles(mol),
        "SssCH2": Descriptors.NumSaturatedHeterocycles(mol),
        "StCH": Descriptors.NumRadicalElectrons(mol),
        "SdsCH": Descriptors.NumValenceElectrons(mol),
        "SaaCH": Descriptors.NumHeteroatoms(mol),
        "SsssCH": Descriptors.NumHAcceptors(mol),
        "StsC": Descriptors.NumHDonors(mol),
        "SdssC": mol.GetNumHeavyAtoms(),
        "SaasC": sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()),
        "SaaaC": saturated_atoms,
        "SssssC": num_rings,
        "SsNH2": heterocycles,
        "SssNH": len(bridgehead_atoms),
        "SaaN": len(spiro_atoms),
        "SsssN": lipinski_hba,
        "SaasN": lipinski_hbd,
        "SsOH": Descriptors.NumRotatableBonds(mol),
        "SdO": Descriptors.MolWt(mol),
        "SssO": Descriptors.MolLogP(mol),
        "SaaO": Descriptors.MolMR(mol),
        "SsF": Descriptors.FractionCSP3(mol),
        "SdsssP": Descriptors.HallKierAlpha(mol),
        "SdS": len(spiro_atoms),
        "SddssS": len(bridgehead_atoms),
        "SsCl": Descriptors.NumHDonors(mol),
        "SsI": Descriptors.NumHAcceptors(mol),
        "C": mol.GetNumHeavyAtoms(),
    }
    return descriptors
def clean_data(data):
    """Fixes NaN, Infinity, and extreme values in the dataset."""
    data = np.array(data, dtype=np.float64)

    # Replace NaN with the mean (or 0 if all values are NaN)
    nan_mask = np.isnan(data)
    if np.any(nan_mask):
        mean_value = np.nanmean(data) if not np.all(nan_mask) else 0.0
        data[nan_mask] = mean_value  # Replace NaNs

    # Replace Infinity with a large finite number
    inf_mask = np.isinf(data)
    data[inf_mask] = np.sign(data[inf_mask]) * 1e6  # Replace +inf with 1e6, -inf with -1e6

    # Clip excessively large values
    data = np.clip(data, -1e6, 1e6)

    return data