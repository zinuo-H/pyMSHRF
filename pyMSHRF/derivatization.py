from .HRF_calculator import parse_formula, format_formula
    

def derivatization(input_formula, num_tms=0, num_meox=0):
    """
    Calculates the derivatized chemical formula based on the input formula and the number of TMS and MeOX groups added.
    Args:
        input_formula (str): The original chemical formula, e.g., "C8H8O3".
        num_tms (int): Number of TMS groups added.
        num_meox (int): Number of MeOX groups added.
    Returns:
        str: The derivatized chemical formula, e.g., "C12H18NO4Si".
    """
    # Parse the input chemical formula
    elements = parse_formula(input_formula)
    
    # Ensure all necessary elements are initialized
    elements["Si"] = elements.get("Si", 0)  # Default to 0 if not present
    elements["C"] = elements.get("C", 0)
    elements["H"] = elements.get("H", 0)
    elements["N"] = elements.get("N", 0)
    elements["O"] = elements.get("O", 0)

    # Add TMS group elements
    elements["Si"] += num_tms           # Each TMS group adds one Si
    elements["C"] += num_tms * 3        # Each TMS group adds three C
    elements["H"] += num_tms * 9 - num_tms  # Each TMS group adds 9 H and removes 1 H
    
    # Add MeOX group elements
    elements["C"] += num_meox           # Each MeOX group adds one C
    elements["H"] += num_meox * 3 - num_meox  # Each MeOX group adds 3 H and removes 1 H
    elements["N"] += num_meox           # Each MeOX group adds one N
    elements["O"] += num_meox           # Each MeOX group adds one O
    
    # Format and return the new chemical formula
    return format_formula(elements)


