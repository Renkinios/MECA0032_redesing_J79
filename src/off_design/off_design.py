from off_design.interpolation_off_disign import get_sensitivity_solidity_1_5, get_momentum_thickness_off
import numpy as np
import get_class as cl

def get_pitching_OFF_design(compressor, atm, stage_design, blade_cascade) :
    T_tot_1 = atm.Ta
    p_tot_1 = atm.Pa 
    rho_1 = atm.Pa/(atm.R * atm.Ta)
    vm =  compressor.m_dot / (compressor.A * atm.rho)
    print("vm", vm)
    T_1 = T_tot_1 - (vm**2)/(2*atm.Cp)
    p_1 = p_tot_1*(T_1/T_tot_1)**(atm.gamma/(atm.gamma-1))
    matrix_stage_OFF_design = np.zeros((compressor.number_stage) * 2 + 2, dtype= object) # take into acount the IGV and the state before it
    stage_0 = cl.Stage(0, T_1, p_1, T_tot_1, p_tot_1, rho_1, 0, 0, 0, 0, vm, compressor)
    matrix_stage_OFF_design[0] = stage_0
    stage_IGV = get_IGV_off_design(compressor, atm, stage_0 ,stage_design, blade_cascade)
    matrix_stage_OFF_design[1] = stage_IGV
    for i in range(compressor.number_stage-1) :
        rotor_stage = get_rotor_off(compressor, atm, matrix_stage_OFF_design[i * 2 + 1], stage_design[2*i + 2].area, blade_cascade)
        matrix_stage_OFF_design[i*2+2]     = rotor_stage
        # print(rotor_stage)
        stator_stage = get_stator_off(compressor, atm, matrix_stage_OFF_design[2*i + 2], stage_design[2*i + 3].area, blade_cascade)
        matrix_stage_OFF_design[(i*2) + 3] = stator_stage
        # print(stator_stage)
        # if i ==2 :
        
    rotor_before_OGV = get_rotor_off(compressor, atm, matrix_stage_OFF_design[2 * (compressor.number_stage) - 1], stage_design[2 * (compressor.number_stage)].area, blade_cascade)
    print(rotor_before_OGV)
    matrix_stage_OFF_design[(compressor.number_stage) * 2 ] = rotor_before_OGV
    OGV_stage = get_OGV_off(compressor, atm, matrix_stage_OFF_design[(compressor.number_stage )*2 ], stage_design[(compressor.number_stage)*2 + 1].area, blade_cascade)
    matrix_stage_OFF_design[(compressor.number_stage)*2 + 1] = OGV_stage

    return matrix_stage_OFF_design


def get_IGV_off_design(compressor, atm, stage_off_design, stage_design,blade_cascade) :
    p_tot2 = stage_off_design.p_tot
    T_tot2 = stage_off_design.T_tot
    # print("arrea taking", stage_design[1].area)
    # print("##### OFF DESIGN IGV #####")

    # print(f"Area: {stage_design[1].area}\n"
    #     f"T_tot2: {T_tot2}\n"
    #     f"p_tot2: {p_tot2}\n"
    #     f"Alpha1 Design (degrees): {np.degrees(blade_cascade.alpha1_design)}\n"
    #     f"Compressor Mass Flow (m_dot): {compressor.m_dot}")
    mac, c_rho = get_mac(stage_design[1].area, T_tot2, p_tot2, blade_cascade.alpha1_design, atm, compressor.m_dot)
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
    # print("vm", vm)
    # print("Beta 2", np.rad2deg(beta_2))
    # print("Alpha 2", np.rad2deg(alpha_2))
    stage_IGV = cl.Stage(1, T_stat, p_stat, T_tot2, p_tot2, rho2, wu, vu, beta_2, alpha_2, vm, compressor)
    return stage_IGV

def get_mac(area_design, T_tot, p_tot, angle, atm, m_dot) :

    mac = np.arange(0.1,1.2,0.001)

    c_rho      = m_dot / (area_design * np.cos(angle))
    left_side  = c_rho * np.sqrt(atm.R * T_tot) / p_tot
    f_m        = get_f_mac(mac, atm)
    right_side = f_m**(- (atm.gamma+1)/(2*(atm.gamma-1))) * np.sqrt(atm.gamma) * mac
    mac = mac[np.argmin(abs(left_side - right_side))]
    print("Mac", mac)
    print('error', np.min(abs(left_side - right_side)))
    return mac, c_rho



def get_f_mac(mac, atm) :
    f_m = 1 + (atm.gamma-1)/2*mac**2
    return f_m

def get_rotor_off(compressor, atm, stage_off_design, area, blade_cascade) :
    # print("##### OFF DESIGN ROTOR #####")
    # print(f"Beta: {stage_off_design.beta} rad/s {np.rad2deg(stage_off_design.beta)}\n",
    #     f"T_stat: {stage_off_design.T_stat}\n"
    #     f"p_stat: {stage_off_design.p_stat}\n"
    #     f"T_tot: {stage_off_design.T_tot}\n"
    #     f"p_tot: {stage_off_design.p_tot}\n"
    #     f"vm: {stage_off_design.vm}\n"
    #     f"vu: {stage_off_design.vu}\n"
    #     f"Blade Cascade u: {blade_cascade.u}\n"
    #     f"Compressor Mass Flow: {compressor.m_dot}\n"
    #     f"Area: {area}")    
    
    sensitivity = get_sensitivity_solidity_1_5(stage_off_design.beta)
    aoa         = stage_off_design.beta - blade_cascade.deflection_design_rotor
    delta_aoa   =(aoa - blade_cascade.aoa_design_rotor)

    # print("aoa", np.rad2deg(aoa))
    # print("aoa_design_rotor", np.rad2deg(blade_cascade.aoa_design_rotor), blade_cascade.aoa_design_rotor)
    # print("deflection_design_rotor", np.rad2deg(blade_cascade.deflection_design_rotor), blade_cascade.deflection_design_rotor)
    beta_2      = stage_off_design.beta + blade_cascade.deflection_design_rotor + sensitivity * delta_aoa
    beta_2 = beta_2
    delta_aoa = np.abs(aoa - blade_cascade.aoa_design_rotor)

    T_tot_rel_1 = stage_off_design.T_stat + (stage_off_design.wu**2 + stage_off_design.vm**2)/(2*atm.Cp) 
    p_tot_rel_1 = stage_off_design.p_tot  * (T_tot_rel_1/stage_off_design.T_tot)**(atm.gamma/(atm.gamma-1))

    alph_naca = 0.0117

    Deq = np.cos(beta_2) / np.cos(stage_off_design.beta) * (1.12 + alph_naca * (delta_aoa)**1.43 
                                                              + 0.61 * np.cos(stage_off_design.beta)**2/compressor.solidity 
                                                              * (np.tan(stage_off_design.beta) - np.tan(beta_2)))
    thinkness = get_momentum_thickness_off(Deq)

    coef_loss = 2 * thinkness * (compressor.solidity / np.cos(stage_off_design.beta)) *(np.cos(stage_off_design.beta) / np.cos(beta_2))**2
    w_1  = - np.sqrt(stage_off_design.wu**2 + stage_off_design.vm**2)
    delta_p_tot = 1/2 * coef_loss * stage_off_design.rho * w_1**2

    T_tot_rel_2  = T_tot_rel_1 # u same 
    p_tot_rel_2  = p_tot_rel_1 - delta_p_tot


    mac, c_rho = get_mac(area,T_tot_rel_2, p_tot_rel_2, beta_2, atm, compressor.m_dot)

    T_stat_2 = T_tot_rel_2 / (get_f_mac(mac, atm))
    p_stat_2 = p_tot_rel_2 * (T_stat_2/T_tot_rel_2)**(atm.gamma/(atm.gamma-1))

    w2  = mac * np.sqrt(atm.gamma * atm.R * T_stat_2)

    wu2 = w2 * np.sin(beta_2)
    vm2 = w2 * np.cos(beta_2)

    vu2     = wu2 + stage_off_design.u
    alpha_2 = np.arctan(vm2/vu2)

    rho2 = p_stat_2 / (atm.R * T_stat_2)

    T_tot = T_stat_2 + (vm2**2 + vu2**2)/(2*atm.Cp)
    p_tot = p_stat_2 * (T_tot/T_stat_2)**(atm.gamma/(atm.gamma-1))

    stage_rotor = cl.Stage(stage_off_design.stage_number + 1, T_stat_2, p_stat_2, T_tot, p_tot, rho2, wu2, vu2,beta_2, alpha_2, vm2, compressor)

    return stage_rotor

def get_stator_off(compressor, atm, stage_off_design, area, blade_cascade) :

    sensitivity = get_sensitivity_solidity_1_5(stage_off_design.alpha)
    aoa         = stage_off_design.alpha - blade_cascade.deflection_design_stator
    delta_aoa   =(aoa - blade_cascade.aoa_design_stator)
    alpha_2     = stage_off_design.alpha - blade_cascade.deflection_design_stator - sensitivity * delta_aoa
    # print("##### OFF DESIGN STATOR #####")
    # print("alpha 2 : ", np.rad2deg(alpha_2))
    delta_aoa = np.abs(aoa - blade_cascade.aoa_design_stator)
    alpha_naca = 0.0117
    D_eq = np.cos(alpha_2) / np.cos(stage_off_design.alpha) * (1.12 + alpha_naca * (delta_aoa)**1.43 
                                                               + 0.61 * np.cos(stage_off_design.alpha)**2/compressor.solidity 
                                                               * (np.tan(stage_off_design.alpha) - np.tan(alpha_2)))
    thinkness = get_momentum_thickness_off(D_eq)
    coefficient_loss = 2 * thinkness * (compressor.solidity / np.cos(stage_off_design.alpha)) * (np.cos(stage_off_design.alpha) / np.cos(alpha_2))**2
    v1 = np.sqrt(stage_off_design.vu**2 + stage_off_design.vm**2)
    delta_p_tot = 1/2 * coefficient_loss * stage_off_design.rho * v1**2

    T_tot_2 = stage_off_design.T_tot
    p_tot_2 = stage_off_design.p_tot - delta_p_tot

    mac, c_rho = get_mac(area, T_tot_2, p_tot_2, alpha_2, atm, compressor.m_dot)

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

    stage_stator = cl.Stage(stage_off_design.stage_number + 1, T_stat_2, p_stat_2, T_tot, p_tot, rho2, wu2, vu2, beta_2, alpha_2, vm2, compressor)

    return stage_stator



def get_OGV_off(compressor, atm, stage_off_design, area, blade_cascade) :

    print("##### OFF DESIGN OGV #####")
    print("stage_off_design.alpha", np.rad2deg(stage_off_design.alpha))

    sensitivity = get_sensitivity_solidity_1_5(stage_off_design.alpha)
    aoa         = stage_off_design.alpha - blade_cascade.deflection_design_OGV
    delta_aoa   = np.abs(aoa - blade_cascade.aoa_design_OGV)
    alpha_2     = stage_off_design.alpha - blade_cascade.deflection_design_OGV - sensitivity * delta_aoa

    print("alpha_2", np.rad2deg(alpha_2))
    alpha_naca = 0.0117

    D_eq = np.cos(alpha_2) / np.cos(stage_off_design.alpha) * (1.12 + 0.007 * (alpha_naca)**1.43 
                                                               + 0.61 * np.cos(stage_off_design.alpha)**2/compressor.solidity 
                                                               * (np.tan(stage_off_design.alpha) - np.tan(alpha_2)))
    thinkness = get_momentum_thickness_off(D_eq)
    coefficient_loss = 2 * thinkness * (compressor.solidity / np.cos(stage_off_design.alpha)) * (np.cos(stage_off_design.alpha) / np.cos(alpha_2))**2
    v1 = np.sqrt(stage_off_design.vu**2 + stage_off_design.vm**2)
    delta_p_tot = 1/2 * coefficient_loss * stage_off_design.rho * v1**2

    T_tot_2 = stage_off_design.T_tot
    p_tot_2 = stage_off_design.p_tot - delta_p_tot

    mac, c_rho = get_mac(area, T_tot_2, p_tot_2, alpha_2, atm, compressor.m_dot)

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

    stage_OGV = cl.Stage(stage_off_design.stage_number + 1, T_stat_2, p_stat_2, T_tot, p_tot, rho2, wu2, vu2, beta_2, alpha_2, vm2, compressor)
    return stage_OGV