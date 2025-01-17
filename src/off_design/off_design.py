from off_design.interpolation_off_disign import get_sensitivity_solidity_1_5, get_momentum_thickness_off
import numpy as np
import get_class as cl

def get_pitching_OFF_design(compressor, atm, stage_design, blade_cascade) :
    T_tot_1 = atm.Ta
    p_tot_1 = atm.Pa 
    rho_1   = atm.Pa/(atm.R * atm.Ta)
    vm      =  compressor.m_dot / (compressor.area_entrance * atm.rho)

    T_1 = T_tot_1 - (vm**2)/(2*atm.Cp)
    p_1 = p_tot_1*(T_1/T_tot_1)**(atm.gamma/(atm.gamma-1))
    matrix_stage_OFF_design = np.zeros((compressor.number_stage) * 2 + 2, dtype= object) # take into acount the IGV and the state before it
    stage_0 = cl.station_compressor(0, T_1, p_1, T_tot_1, p_tot_1, rho_1,vm, compressor)
    matrix_stage_OFF_design[0] = stage_0
    stage_IGV = get_IGV_off_design(compressor, atm, stage_0 ,stage_design, blade_cascade)

    matrix_stage_OFF_design[1] = stage_IGV
    for i in range(compressor.number_stage-1) :
        rotor_stage = get_rotor_off(compressor, atm, matrix_stage_OFF_design[i * 2 + 1], stage_design[2*i + 2].area, blade_cascade)
        matrix_stage_OFF_design[i*2+2]     = rotor_stage
        stator_stage = get_stator_off(compressor, atm, matrix_stage_OFF_design[2*i + 2], stage_design[2*i + 3].area, blade_cascade)
        matrix_stage_OFF_design[(i*2) + 3] = stator_stage
    rotor_before_OGV = get_rotor_off(compressor, atm, matrix_stage_OFF_design[2 * (compressor.number_stage) - 1], stage_design[2 * (compressor.number_stage)].area, blade_cascade)

    matrix_stage_OFF_design[(compressor.number_stage) * 2 ] = rotor_before_OGV
    OGV_stage = get_OGV_off(compressor, atm, matrix_stage_OFF_design[(compressor.number_stage )*2 ], stage_design[(compressor.number_stage)*2 + 1].area, blade_cascade)
    matrix_stage_OFF_design[(compressor.number_stage)*2 + 1] = OGV_stage

    return matrix_stage_OFF_design

def get_f_mac(mac, atm) :
    f_m = 1 + (atm.gamma - 1)/2 * mac**2
    return f_m

def get_mac(area_design, T_tot, p_tot, angle, atm, m_dot, tol=1e-6, max_iter=100):
    """
    Calculate the Mach number (MAC) using a convergence method.

    Parameters:
    - area_design: section area
    - T_tot: total temperature
    - p_tot: total pressure
    - angle: attack angle in radians
    - atm: object containing atmospheric properties (R, gamma)
    - m_dot: mass flow rate
    - tol: tolerance for convergence (default: 1e-6)
    - max_iter: maximum number of iterations (default: 100)

    Returns:
    - mac: converged Mach number
    """
    
    def residual(mac):
        c_rho = m_dot / (area_design * np.cos(angle))
        left_side = c_rho * np.sqrt(atm.R * T_tot) / p_tot
        f_m = get_f_mac(mac, atm)
        right_side = f_m**(- (atm.gamma + 1) / (2 * (atm.gamma - 1))) * np.sqrt(atm.gamma) * mac
        return left_side - right_side

    def derivative(mac):
        # Numerical approximation of the residual function derivative
        delta = 1e-8
        return (residual(mac + delta) - residual(mac - delta)) / (2 * delta)

    # Initialization
    mac = 0.5  # Initial guess
    for i in range(max_iter):
        res = residual(mac)
        deriv = derivative(mac)
        
        if abs(res) < tol:  # Check for convergence
            return mac
        
        if deriv == 0:
            raise ValueError("Zero derivative encountered, Newton's method fails.")
        
        mac -= res / deriv  # Update MAC
        
        if mac < 0.1 or mac > 0.8:
            raise ValueError(f"Mach number out of bounds {mac}. Check input parameters.")
    
    raise RuntimeError("Convergence not achieved after the maximum number of iterations.")

def get_IGV_off_design(compressor, atm, stage_off_design, stage_design,blade_cascade) :
    p_tot2 = stage_off_design.p_tot
    T_tot2 = stage_off_design.T_tot
    mac = get_mac(stage_design[1].area, T_tot2, p_tot2, blade_cascade.alpha1_design, atm, compressor.m_dot)

    T_stat = T_tot2 / get_f_mac(mac, atm)
    p_stat = p_tot2 *(T_stat/T_tot2)**(atm.gamma/(atm.gamma-1))
    speed_sound = np.sqrt(atm.gamma * atm.R * T_stat)
    v2 = mac * speed_sound
    rho2 = p_stat / (atm.R * T_stat)
    vm = v2 * np.cos(blade_cascade.alpha1_design)
    vu = v2 * np.sin(blade_cascade.alpha1_design)
    wu = vu - blade_cascade.u
    beta_2  = np.arctan(wu/vm)
    alpha_2 = np.arctan(vu/vm)
    stage_IGV = cl.station_compressor(1, T_stat, p_stat, T_tot2, p_tot2, rho2, vm, compressor,wu, vu, beta_2, alpha_2, mac)
    return stage_IGV


def get_rotor_off(compressor, atm, stage_off_design, area, blade_cascade) :    
    
    sensitivity = get_sensitivity_solidity_1_5(stage_off_design.beta)
    aoa         = -stage_off_design.beta - blade_cascade.stragger_rotor
    delta_aoa   = (aoa - blade_cascade.aoa_design_rotor)
    beta_2      = stage_off_design.beta + blade_cascade.deflection_design_rotor + sensitivity * delta_aoa

    if abs(stage_off_design.beta) >= blade_cascade.stall_rotor_pos or abs(stage_off_design.beta) <= blade_cascade.stall_rotor_neg:
        raise ValueError(f'Beta rotor is in the stall region {np.rad2deg(abs(stage_off_design.beta))}')

    delta_aoa = np.abs(aoa - blade_cascade.aoa_design_rotor)

    alph_naca = 0.0117

    Deq = np.cos(-beta_2) / np.cos(-stage_off_design.beta) * (1.12 + alph_naca * (delta_aoa)**1.43 
                                                              + 0.61 * np.cos(-stage_off_design.beta)**2/compressor.solidity 
                                                              * (np.tan(-stage_off_design.beta) - np.tan(-beta_2)))
    
    thinkness = get_momentum_thickness_off(Deq)


    T_tot_rel_1 = stage_off_design.T_stat + (stage_off_design.wu**2 + stage_off_design.vm**2)/(2 * atm.Cp) 
    p_tot_rel_1 = stage_off_design.p_tot  * (T_tot_rel_1/stage_off_design.T_tot)**(atm.gamma/(atm.gamma - 1))
    
    coef_loss   = 2 * thinkness * (compressor.solidity / np.cos(-beta_2)) *(np.cos(-stage_off_design.beta) / np.cos(-beta_2))**2
    w1          = - np.sqrt(stage_off_design.wu**2 + stage_off_design.vm**2)
    delta_p_tot = 1/2 * coef_loss * stage_off_design.rho * w1**2

    T_tot_rel_2  = T_tot_rel_1 # u same 
    p_tot_rel_2  = p_tot_rel_1 - delta_p_tot


    mac = get_mac(area,T_tot_rel_2, p_tot_rel_2,- beta_2, atm, compressor.m_dot)

    T_stat_2 = T_tot_rel_2 / (get_f_mac(mac, atm))
    p_stat_2 = p_tot_rel_2 * (T_stat_2/T_tot_rel_2)**(atm.gamma/(atm.gamma-1))

    w2  = mac * np.sqrt(atm.gamma * atm.R * T_stat_2)

    wu2 = w2 * np.sin(beta_2)

    vm2 = w2 * np.cos(beta_2)

    vu2     = wu2 + stage_off_design.u
    alpha_2 = np.arctan(vu2/vm2)


    rho2 = p_stat_2 / (atm.R * T_stat_2)

    T_tot = T_stat_2 + (vm2**2 + vu2**2)/(2 * atm.Cp)
    p_tot = p_stat_2 * (T_tot/T_stat_2)**(atm.gamma/(atm.gamma - 1))

    stage_rotor = cl.station_compressor(stage_off_design.stage_number + 1, T_stat_2, p_stat_2, T_tot, p_tot, rho2, vm2, compressor,  wu2, vu2, beta_2, alpha_2,mac)

    return stage_rotor

def get_stator_off(compressor, atm, stage_off_design, area, blade_cascade) :

    sensitivity = get_sensitivity_solidity_1_5(stage_off_design.alpha)
    aoa         = stage_off_design.alpha - blade_cascade.stragger_stator
    delta_aoa   = (aoa - blade_cascade.aoa_design_stator)
    alpha_2     = stage_off_design.alpha - blade_cascade.deflection_design_stator - sensitivity * delta_aoa

    if np.abs(stage_off_design.alpha) >= blade_cascade.stall_stator_pos or np.abs(stage_off_design.alpha) <= blade_cascade.stall_stator_neg:
        raise ValueError(f'Alpha stator is in the stall region {np.rad2deg(alpha_2)}')
    

    delta_aoa = np.abs(aoa - blade_cascade.aoa_design_stator)
    alpha_naca = 0.0117
    D_eq = np.cos(alpha_2) / np.cos(stage_off_design.alpha) * (1.12 + alpha_naca * (delta_aoa)**1.43 
                                                               + 0.61 * np.cos(stage_off_design.alpha)**2/compressor.solidity 
                                                               * (np.tan(stage_off_design.alpha) - np.tan(alpha_2)))
    thinkness = get_momentum_thickness_off(D_eq)
    coefficient_loss = 2 * thinkness * (compressor.solidity / np.cos(alpha_2)) * (np.cos(stage_off_design.alpha) / np.cos(alpha_2))**2
    v1 = np.sqrt(stage_off_design.vu**2 + stage_off_design.vm**2)
    delta_p_tot = 1/2 * coefficient_loss * stage_off_design.rho * v1**2

    T_tot_2 = stage_off_design.T_tot
    p_tot_2 = stage_off_design.p_tot - delta_p_tot

    mac = get_mac(area, T_tot_2, p_tot_2, alpha_2, atm, compressor.m_dot)

    T_stat_2 = T_tot_2 / get_f_mac(mac, atm)
    p_stat_2 = p_tot_2 * (T_stat_2/T_tot_2)**(atm.gamma/(atm.gamma-1))

    rho2 = p_stat_2 / (atm.R * T_stat_2)

    speed_sound = np.sqrt(atm.gamma * atm.R * T_stat_2)
    v2 = mac * speed_sound


    vu2    = v2 * np.sin(alpha_2)
    vm2    = v2 * np.cos(alpha_2)
    wu2    = vu2 - stage_off_design.u
    beta_2 = np.arctan(wu2/vm2)
    T_tot = T_stat_2 + (vm2**2 + vu2**2)/(2*atm.Cp)
    p_tot = p_stat_2 * (T_tot/T_stat_2)**(atm.gamma/(atm.gamma - 1))

    stage_stator = cl.station_compressor(stage_off_design.stage_number + 1, T_stat_2, p_stat_2, T_tot, p_tot, rho2, vm2, compressor, wu2, vu2, beta_2, alpha_2, mac)

    return stage_stator



def get_OGV_off(compressor, atm, stage_off_design, area, blade_cascade) :


    sensitivity = get_sensitivity_solidity_1_5(stage_off_design.alpha)
    aoa         = stage_off_design.alpha - blade_cascade.stragger_OGV
    delta_aoa   = (aoa - blade_cascade.aoa_design_OGV)
    alpha_2     = stage_off_design.alpha - blade_cascade.deflection_design_OGV - sensitivity * delta_aoa

    if abs(stage_off_design.alpha) >= blade_cascade.stall_OGV_pos or abs(stage_off_design.alpha) <= blade_cascade.stall_OGV_neg:
        raise ValueError(f'Alpha OGV is in the stall region {np.rad2deg(alpha_2)}')

    alpha_naca = 0.0117
    delta_aoa = np.abs(delta_aoa)
    
    D_eq = np.cos(alpha_2) / np.cos(stage_off_design.alpha) * (1.12 + alpha_naca * (delta_aoa)**1.43 
                                                               + 0.61 * np.cos(stage_off_design.alpha)**2/compressor.solidity 
                                                               * (np.tan(stage_off_design.alpha) - np.tan(alpha_2)))
    
    thinkness        = get_momentum_thickness_off(D_eq)
    coefficient_loss = 2 * thinkness * (compressor.solidity / np.cos(alpha_2)) * (np.cos(stage_off_design.alpha) / np.cos(alpha_2))**2
    
    v1 = np.sqrt(stage_off_design.vu**2 + stage_off_design.vm**2)
    delta_p_tot = 1/2 * coefficient_loss * stage_off_design.rho * v1**2

    T_tot_2 = stage_off_design.T_tot
    p_tot_2 = stage_off_design.p_tot - delta_p_tot

    mac = get_mac(area, T_tot_2, p_tot_2, alpha_2, atm, compressor.m_dot)

    T_stat_2 = T_tot_2 / get_f_mac(mac, atm)
    p_stat_2 = p_tot_2 * (T_stat_2/T_tot_2)**(atm.gamma/(atm.gamma-1))

    speed_sound = np.sqrt(atm.gamma * atm.R * T_stat_2)
    v2 = mac * speed_sound

    rho2 = p_stat_2 / (atm.R * T_stat_2)

    vu2 = v2 * np.sin(alpha_2)
    vm2 = v2 * np.cos(alpha_2)

    wu2 = vu2 - stage_off_design.u

    beta_2 = np.arctan(wu2/vm2)

    T_tot = T_stat_2 + (vm2**2 + vu2**2)/(2*atm.Cp)
    p_tot = p_stat_2 * (T_tot/T_stat_2)**(atm.gamma/(atm.gamma-1))

    stage_OGV = cl.station_compressor(stage_off_design.stage_number + 1, T_stat_2, p_stat_2, T_tot, p_tot, rho2,  vm2, compressor, wu2, vu2, beta_2, alpha_2,mac)
    return stage_OGV