import numpy as np
import VizData as viz
import get_class as cl
import design.difusion_interpolate as di


haller_criterion = 0.68
def get_pitching_design(compresor, atm) :

    triangle_compressor = cl.CompressorTriangle(compresor, atm)
    
    T_tot_1 = atm.Ta
    p_tot_1 = atm.Pa 
    rho_1 = atm.Pa/(atm.R * atm.Ta)
    
    T_1 = T_tot_1 - (triangle_compressor.vm**2)/(2*atm.Cp)
    p_1 = p_tot_1*(T_1/T_tot_1)**(atm.gamma/(atm.gamma-1))
    
    matrix_stage = np.zeros((compresor.number_stage) * 2 + 2, dtype= object) # take into acount the IGV and the state before it
    stage_0 = cl.Stage(0, T_1, p_1, T_tot_1, p_tot_1, rho_1, 0, 0, 0, 0, triangle_compressor.vm, compresor)
    matrix_stage[0] = stage_0
    stage_IGV = get_IGV_outlet(stage_0, triangle_compressor, compresor, atm)
    matrix_stage[1] = stage_IGV
    for i in range(compresor.number_stage-1) :
        rotor_stage = get_rotor_outlet(matrix_stage[i * 2 + 1], triangle_compressor, compresor, atm)
        matrix_stage[i*2+2]     = rotor_stage
        stator_stage = get_stator_outlet(matrix_stage[2*i + 2], triangle_compressor, compresor, atm)
        matrix_stage[(i*2) + 3] = stator_stage

    rotor_before_OGV = get_rotor_outlet(matrix_stage[2 * (compresor.number_stage) - 1], triangle_compressor, compresor, atm)
    matrix_stage[(compresor.number_stage) * 2 ] = rotor_before_OGV
    OGV_stage = get_OGV_outlet(matrix_stage[(compresor.number_stage )*2 ], triangle_compressor, compresor, atm)
    matrix_stage[(compresor.number_stage)*2 + 1] = OGV_stage    
    return matrix_stage, triangle_compressor

def get_IGV_outlet(stage_inlet, triangle_compressor, compressor, atm) :
    # Goal is to fix the velocity of the IGV on the rotor so on the angle beta_1

    wu_2 = triangle_compressor.vm * np.tan(triangle_compressor.beta_1)
    vu_2 = wu_2 + triangle_compressor.u
    T_tot_2 = stage_inlet.T_tot
    p_tot_2 = stage_inlet.p_tot     
    T_stat = T_tot_2 - (triangle_compressor.vm**2 + vu_2**2)/(2*atm.Cp)
    p_stat = p_tot_2 * (T_stat/T_tot_2)**(atm.gamma/(atm.gamma-1))
    rho_2 = p_stat / (atm.R*T_stat)
    a = np.sqrt(atm.gamma*atm.R*T_stat)   
    v = np.sqrt(vu_2**2 + triangle_compressor.vm**2)
    Mw = v/a
    # print("Mac IGV", Mw)
    if Mw > 0.8:
        raise ValueError(f"Mach number is too high in the IGV {Mw} > 0.8")    
    outlet_IGV = cl.Stage(1, T_stat, p_stat, T_tot_2, p_tot_2, rho_2, wu_2, vu_2,triangle_compressor.beta_2, triangle_compressor.alpha_2, triangle_compressor.vm, compressor)
    return outlet_IGV


def get_rotor_outlet(stage_inlet, triangle_compressor, compressor, atm) :
    
    wu_1 = stage_inlet.wu 
    vu_1 = stage_inlet.vu
    wu_2 = triangle_compressor.vm * np.tan(triangle_compressor.beta_2)
    vu_2 = wu_2 + triangle_compressor.u

    w1   =  - np.sqrt(wu_1**2 + triangle_compressor.vm**2)
    w2   =  - np.sqrt(wu_2**2 + triangle_compressor.vm**2)
    if w2/w1 <= haller_criterion :
        raise ValueError(f'The Haller criterion is not respected in the rotor  {w2/w1} < 0.7')
    
    # like rotor use the relative frae of reference
    T_tot_rel_1 = stage_inlet.T_stat + (wu_1**2 + triangle_compressor.vm**2)/(2*atm.Cp) # erreur ici je pense bien 
    p_tot_rel_1 = stage_inlet.p_tot  * (T_tot_rel_1/stage_inlet.T_tot)**(atm.gamma/(atm.gamma-1))
    # Danger ajouter le polytrpoique effisency normalment
    DF = 1 - w2/w1 + (wu_1 - wu_2)/(2 * compressor.solidity * w1)
    if DF > 0.6 :   
        raise ValueError(f'The diffusoin factor is not respected in the rotor {DF} > 0.6')

    moment_thickness = di.get_interpolation_diffusion(DF)

    loss_coeff       = 2 * moment_thickness * (compressor.solidity/(np.cos(triangle_compressor.beta_2)))  * (np.cos(triangle_compressor.beta_1) /np.cos(triangle_compressor.beta_2))**2
    delta_p_tot      = 1/2 * loss_coeff * stage_inlet.rho * w1**2

    T_tot_rel_2  = T_tot_rel_1 # u same 
    p_tot_rel_2  = p_tot_rel_1 - delta_p_tot

    T_stat_2 = T_tot_rel_2 - (triangle_compressor.vm**2 + wu_2**2)/(2*atm.Cp)
    p_stat_2 = p_tot_rel_2 * (T_stat_2/T_tot_rel_2)**(atm.gamma/(atm.gamma-1))
    rho_2 = p_stat_2 / (atm.R*T_stat_2)

    T_tot_2 = T_stat_2 + (triangle_compressor.vm**2 + vu_2**2)/(2*atm.Cp)
    p_tot_2 = p_stat_2 * (T_tot_2/T_stat_2)**(atm.gamma/(atm.gamma-1))

    a = np.sqrt(atm.gamma*atm.R*T_stat_2)
    w = np.sqrt(wu_2**2 + triangle_compressor.vm**2)
    mac = w/a
    # print("MAC rotor", mac)
    if mac > 0.8:
        raise ValueError(f'Mach number is too high in the rotor in the rotor {mac} > 0.8')
    outlet_rotor = cl.Stage(stage_inlet.stage_number + 1, T_stat_2, p_stat_2, T_tot_2, p_tot_2, rho_2, wu_2, vu_2,triangle_compressor.beta_2, triangle_compressor.alpha_2, triangle_compressor.vm, compressor)

    return outlet_rotor


def get_stator_outlet(stage_inlet, triangle_compressor, compressor, atm) :

    vu_1 = stage_inlet.vu
    wu_2 = triangle_compressor.vm * np.tan(triangle_compressor.beta_1)
    vu_2 = wu_2 + triangle_compressor.u

    v1   =  + np.sqrt(vu_1**2 + triangle_compressor.vm**2)
    v2   =  + np.sqrt(vu_2**2 + triangle_compressor.vm**2)
    alpha_1 = triangle_compressor.alpha_2
    alpha_2 = triangle_compressor.alpha_1
    if v2/v1 < haller_criterion :
        raise ValueError(f'The Haller criterion is not respected in the stator {v2/v1} < 0.7')
    
    DF = 1 - v2/v1 + (vu_1 - vu_2)/(2 * compressor.solidity * v1)
    if DF > 0.6 :   
        raise ValueError(f'The diffusoin factor is not respected in the stator {DF} > 0.6')
    moment_thickness = di.get_interpolation_diffusion(DF)

    loss_coeff       = 2 * moment_thickness * (compressor.solidity /(np.cos(alpha_2))) * (np.cos(alpha_1) /np.cos(alpha_2))**2
    delta_p_tot      = 1/2 * loss_coeff * stage_inlet.rho * v1**2


    T_tot_2  = stage_inlet.T_tot
    p_tot_2  = stage_inlet.p_tot - delta_p_tot

    T_stat_2 = T_tot_2 - (triangle_compressor.vm**2 + vu_2**2)/(2*atm.Cp)
    p_stat_2 = p_tot_2 * (T_stat_2/T_tot_2)**(atm.gamma/(atm.gamma-1))
    rho_2 = p_stat_2 / (atm.R*T_stat_2)

    a = np.sqrt(atm.gamma*atm.R*T_stat_2)
    w = np.sqrt((wu_2**2 + triangle_compressor.vm**2))
    mac = v2/a
    if mac > 0.8:
        raise ValueError(f'Mach number is too high in the stator {mac} > 0.8')
    outlet_stator = cl.Stage(stage_inlet.stage_number + 1, T_stat_2, p_stat_2, T_tot_2, p_tot_2, rho_2, wu_2, vu_2,triangle_compressor.beta_2, triangle_compressor.alpha_2, triangle_compressor.vm, compressor)
    return outlet_stator

def get_OGV_outlet(stage_inlet, triangle_compressor, compressor, atm) :

    wu_2 = - triangle_compressor.u
    vu_1 = stage_inlet.vu
    
    vu_2    =  0
    alpha_2 = 0
    alpha_1 = triangle_compressor.alpha_2
    v1      =  + np.sqrt(vu_1**2 + triangle_compressor.vm**2)
    v2      =  triangle_compressor.vm

    if v2/v1 <= haller_criterion :
        raise ValueError(f'The Haller criterion is not respected in the OGV {v2/v1} < 0.7')

    DF = 1 - v2/v1 + (vu_1 - vu_2)/(2 * compressor.solidity * v1)

    if DF > 0.6 :   
        raise ValueError(f'The diffusoin factor is not respected for the OGV {DF} > 0.6')
    
    moment_thickness = di.get_interpolation_diffusion(DF)
    loss_coeff       = 2 * moment_thickness * (compressor.solidity /(np.cos(alpha_2)))  * (np.cos(alpha_1) /np.cos(alpha_2))**2
    delta_p_tot      =  1/2 * loss_coeff * stage_inlet.rho * v1**2

    T_tot_2  = stage_inlet.T_tot
    p_tot_2  = stage_inlet.p_tot - delta_p_tot

    T_stat_2 = T_tot_2 - (triangle_compressor.vm**2 + vu_2**2)/(2*atm.Cp)
    p_stat_2 = p_tot_2 * (T_stat_2/T_tot_2)**(atm.gamma/(atm.gamma-1))

    rho_2 = p_stat_2 / (atm.R*T_stat_2)
    a = np.sqrt(atm.gamma*atm.R*T_stat_2)
    mac = v2/a
    if mac > 0.8:
        raise ValueError(f'Mach number is too high in the OGV {mac} > 0.8')

    outlet_OGV = cl.Stage(stage_inlet.stage_number + 1, T_stat_2, p_stat_2, T_tot_2, p_tot_2, rho_2, wu_2, vu_2, triangle_compressor.beta_2, triangle_compressor.alpha_2, triangle_compressor.vm, compressor)
    return outlet_OGV