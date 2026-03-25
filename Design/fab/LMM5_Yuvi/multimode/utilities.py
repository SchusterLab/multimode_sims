import numpy as np
def find_E_J(E_c, omega, is_squid=False):
    """
    Find E_J given E_c and omega. 
    If is_squid is True, divideby 2 because parallel addition of two junctions halves Lj  and doubles Ej
    this function is outputting Ej needed per junction
    """
    E_J = ((omega-E_c)**2) / (8 * E_c)
    if is_squid:
        E_J /= 2
    return E_J

def find_transmon_freq(E_J, E_c):
    """
    Calculate transmon qubit frequency omega (in GHz) from Josephson energy E_J and charging energy E_c (both in GHz).
    omega = sqrt(8 * E_J * E_c) - E_c
    """
    omega = np.sqrt(8 * E_J * E_c) - E_c
    return omega

def capacitance_from_Ec(Ec):
    """
    Convert charging energy E_c (in GHz) to capacitance C (in Farads).
    E_c = e^2 / (2C * hbar)
    """
    e = 1.602176634e-19  # C
    hbar = 1.0545718e-34  # J*s
    Ec_J = Ec * 1e9 * hbar  # Convert GHz to Joules
    C = e**2 / (2 * Ec_J)
    return C

def find_L_J(E_J_freq):
    """
    Calculate Josephson inductance L_J from Josephson energy frequency E_J_freq (in GHz).

    L_J = Φ₀² / (4 * π² * h * E_J_freq)
    where:
        Φ₀ = flux quantum = 2.068e-15 Wb
        h = Planck's constant = 2π * ħ
    """
    phi0 = 2.068e-15  # Wb
    h = 6.62607015e-34  # J*s
    E_J_Hz = E_J_freq * 1e9  # Convert GHz to Hz
    L_J = phi0**2 / (4 * np.pi**2 * h * E_J_Hz)
    return L_J

# Conversions from SEM 
def bridge_width_to_junction_length(bridge_width):
    junction_length = 0 
    if bridge_width == 0.2: 
        junction_length = 0.28
    return junction_length

def Ej_from_pin_length(pin_length, critical_current_per_area, gap, sem_correction=True):
    """
    Calculate the Josephson energy E_j from pin length and critical current density.
    pin_length: in um
    gap: in um
    critical_current_per_area: in microA/um^2
    Returns: E_j in GHz
    """
    # Constants
    hbar = 1.0545718e-34  # J*s
    e = 1.602176634e-19   # C
    GHz_to_J = 6.62607015e-34  # Planck's constant (J*s) * 1 Hz

    # Area in um^2
    if sem_correction: 
        print('Applying SEM correction to both bridge and pin when calculating E_j')
        pin_length += 0.07  # um
    area = pin_length * bridge_width_to_junction_length(bridge_width=gap)  # um^2

    # Convert critical_current_per_area to A/um^2
    jc = critical_current_per_area * 1e-6

    # Calculate I_c
    I_c = jc * area  # in Amps

    # Calculate E_j in Joules
    E_j_J = (hbar / (2 * e)) * I_c  # in Joules

    # Convert E_j to GHz
    E_j_GHz = E_j_J / (GHz_to_J * 1e9)

    return E_j_GHz

def pin_size_from_Ej(E_j, critical_current_per_area, bridge_width, sem_correction=True):
    """
    Calculate the pin size (area) needed for a given E_j.
    E_j: Josephson energy in GHz
    critical_current_per_area: in microA/um^2
    gap: in um (not used in calculation, just for reference)
    Returns: area in um^2
    """
    # Constants
    hbar = 1.0545718e-34  # J*s
    e = 1.602176634e-19   # C
    GHz_to_J = 6.62607015e-34  # Planck's constant (J*s) * 1 Hz

    # Convert E_j from GHz to Joules
    E_j_J = E_j * GHz_to_J * 1e9

    # Calculate I_c from E_j
    I_c = (2 * e * E_j_J) / hbar  # in Amps

    # Convert critical_current_per_area to A/um^2
    jc = critical_current_per_area * 1e-6

    # Area in um^2
    area = I_c / jc

    # pin length (SEM correction)
    pin_length = area / bridge_width_to_junction_length(bridge_width)
    if sem_correction: 
        print('Applying SEM correction to both bridge and pin ')
        pin_length -= 0.07  # um
    return area, pin_length 

def calculate_R_n(pin_length, bridge_width, critical_current_per_area, sem_correction=True, is_squid=False):
    """
    Calculate the normal resistance R_n of the junction.

    Formulas used:
    - Area = pin_length × gap
    - Critical current: I_c = critical_current_per_area × Area
    - Josephson relation: E_J = (ħ / 2e) * I_c
    - I_C = (π Δ(0)) / (2e R_n), where Δ(0) = 170 μeV for Aluminum

    Args:
        pin_length: length of the pin in um
        gap: gap in um
        critical_current_per_area: in microA/um^2

    Returns:
        R_n in Ohms
    """
    # Constants
    e = 1.602176634e-19  # C
    Delta_0 = 170e-6 * e  # 170 μeV in Joules

    # Calculate area
    if sem_correction: 
        print('Applying SEM correction to both bridge and pin when calculating R_n')
        pin_length += 0.07  # um
    area = pin_length * bridge_width_to_junction_length(bridge_width)  # um^2

    # Calculate critical current (I_c)
    critical_current_per_area_A = critical_current_per_area * 1e-6  # A/um^2
    I_c = critical_current_per_area_A * area  # A

    # Calculate R_n using I_C = (π Δ(0)) / (2e R_n)
    R_n = (np.pi * Delta_0) / (2 * e * I_c)  # Ohms
    if is_squid:
        R_n /= 2  # For SQUID, R_n is halved

    return R_n