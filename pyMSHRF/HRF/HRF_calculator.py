from .. import pyisopach
import numpy as np
from re import findall

ELECTRON_WEIGHT=0.0005485799


def generate_structure_dict(molecule:str) -> dict:
    parsed = findall(r"([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)",  molecule)
    structure_dict = {}
    for element_details in parsed:
        element = element_details[0]
        if element not in structure_dict:
            structure_dict[element] = 0
            
        element_count = sum(
            [int(x) for x in element_details[1:] if x != ""])
        if element_count > 0:
            structure_dict[element] += element_count
        else:
            structure_dict[element] += 1
    return structure_dict

def mol_to_formula(molstr:str) -> str:
    structure_dict = generate_structure_dict(molstr)
    formula = []
    for key in structure_dict:
        if structure_dict[key] ==1:
            formula.append(key)
        else: 
            formula.append(key + str(structure_dict[key]))
    formula = "".join(formula)
    return formula
    
    
def generate_sub_formulas(formula:str):
    elements_dict = generate_structure_dict(formula)
    def backtrack(idx, current_formula):
        if idx == len(keys):
            sub_formulas.append(current_formula)
            return

        element = keys[idx]
        count = elements_dict[element]

        for i in range(count + 1):
            new_formula = current_formula + element * i
            backtrack(idx + 1, new_formula)

    sub_formulas = []
    keys = list(elements_dict.keys())
    backtrack(0, "")
    sub_formulas.remove("")
    
    sub_formulas_new=[]
    for sub_formula in sub_formulas:
        sub_formulas_new.append(mol_to_formula(sub_formula))
    
    return sub_formulas_new


# calculate HRF
def HRF(formula:str,query_spectrum, delta_ppm=None,delta_da=None): 
    """
     If both delta_ppm and delta_da is defined, delta_da will be used.

    """
    query_spectrum= np.array(query_spectrum,dtype=np.float64)
    if np.sum(query_spectrum[:,1]) == 0:
        return 0
    
    def _get_isotopic_weight(sub_formula):
        try:
            mol = pyisopach.Molecule(sub_formula)
            isotopic_weight = mol.isotopic_distribution()[0] - ELECTRON_WEIGHT
        except:
            mol = pyisopach.Element(sub_formula,1)
            isotopic_weight = np.array(mol.isotopic_weight, dtype=np.float64)  - ELECTRON_WEIGHT
        return isotopic_weight
    
    
    sub_formulas=generate_sub_formulas(formula)
    sub_formulas_dict = {}
    for sub_formula in sub_formulas:
        isotopic_weight = _get_isotopic_weight(sub_formula)
        sub_formulas_dict[sub_formula] = isotopic_weight[0]  # dic of subformulas and its monoisotopic mass  
    
    # generate annotated spec merged
    a = 0
    spec_merged = []     
    while a < len(query_spectrum):
        exp_mass = query_spectrum[a][0]
        closest_formula = None
        closest_mass_diff = float('inf')
        if delta_da is None:
            delta_da = delta_ppm * exp_mass * 1e-6
        
        for formula, mass in sub_formulas_dict.items():
            mass_diff = abs(mass - exp_mass)
            if mass_diff < closest_mass_diff and mass_diff <= delta_da:
                closest_mass_diff = mass_diff
                closest_formula = formula
        
        if closest_formula is not None:
            # peak a annotated
            spec_merged.append(query_spectrum[a])
            if "_iso" not in closest_formula:
                match_isotopic_weight_dis = _get_isotopic_weight(closest_formula)
                iso_No = 1
            else:
                match_isotopic_weight_dis = _get_isotopic_weight(closest_formula.partition("_iso")[0])
                iso_No = int(closest_formula.partition("_iso")[2]) +1
            if len(match_isotopic_weight_dis) > iso_No:
                sub_formulas_dict[closest_formula.partition("_iso")[0]+'_iso'+str(iso_No)] = match_isotopic_weight_dis[iso_No]
        a+=1
    
    
    # calculate score
    if len(spec_merged) >0:
        spec_merged = np.array(spec_merged)
        HRF = np.sum(spec_merged[:,0] * spec_merged[:,1]) /  np.sum(query_spectrum[:,0] * query_spectrum[:,1])
    else:
        HRF=0
        
    return HRF


# calculate RHRF
def RHRF(formula,query_spectrum, ref_spectrum, delta_ppm=None,delta_da=None,match_delta_ppm=None,match_delta_da=0.6): 
    """
     If both delta_ppm and delta_da is defined, delta_da will be used.
    """
    query_spectrum= np.array(query_spectrum,dtype=np.float64)
    ref_spectrum= np.array(ref_spectrum,dtype=np.float64)
    query_spectrum = np.array(sorted(query_spectrum, key=lambda x: x[0]))
    ref_spectrum = np.array(sorted(ref_spectrum, key=lambda x: x[0]))
    
    query_spec_match = [] 
    a=0
    b=0
    
    while a < query_spectrum.shape[0] and b < ref_spectrum.shape[0]:
        if match_delta_da is None:
            match_delta_da = match_delta_ppm * query_spectrum[a, 0] * 1e-6
        mass_delta = query_spectrum[a, 0] - ref_spectrum[b, 0]
        
        if mass_delta < -match_delta_da:
            # Peak only existed in query_spectrum .
            a += 1
        elif mass_delta > match_delta_da:
            # Peak only existed in ref_spectrum.
            b += 1
            
        else:
            # Peak existed in both spec.
            query_spec_match.append(query_spectrum[a])
            a += 1
            b += 1
    
    if query_spec_match:
        query_spec_match = np.array(query_spec_match, dtype=np.float64)
    else:
        query_spec_match = np.array([[0., 0.]], dtype=np.float64)
    
    RHRF = HRF(formula,query_spec_match, delta_ppm, delta_da)
    return RHRF
    