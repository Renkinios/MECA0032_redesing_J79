def get_stator_off(compressor, off_design_stage, flow_area, blade_geometry):
    Ta = 288.15
    Pa = 101325
    gamma = 1.4
    R = 287.1
    Cp = R*gamma/(gamma - 1)
    rho = Pa / (R*Ta)
    
    beta_exit_ref = np.deg2rad(-48.96)
    
    solidity_sensitivity = get_sensitivity_solidity_1_5(off_design_stage.alpha)

    angle_of_attack = off_design_stage.alpha - blade_geometry.stagger_angle_stator
    deviation_angle = angle_of_attack - blade_geometry.design_aoa_stator
    
    alpha_exit = 0.17
    
    # Tolérance
    beta_tolerance = 0.01 * beta_exit_ref  # 1% de tolérance
    
    # Pas d'itération
    step_size = 0.001  # Incrément pour ajuster alpha_exit
    
    max_iterations = 1000
    iteration_count = 0
    
    # Boucle d'itération
    while iteration_count < max_iterations:
        # Calcul des angles et vérifications de décrochage
        if np.abs(off_design_stage.alpha) >= blade_geometry.stall_stator_pos or \
           np.abs(off_design_stage.alpha) <= blade_geometry.stall_stator_neg:
            raise ValueError(f'Alpha stator is in the stall region: {np.rad2deg(alpha_exit)} degrees')
    
        deviation_magnitude = np.abs(angle_of_attack - blade_geometry.design_aoa_stator)
    
        diffusion_factor = (np.cos(alpha_exit) / np.cos(off_design_stage.alpha)) * (
            1.12 + 0.0117 * (deviation_magnitude)**1.43 +
            0.61 * (np.cos(off_design_stage.alpha)**2 / compressor.solidity) *
            (np.tan(off_design_stage.alpha) - np.tan(alpha_exit))
        )
    
        momentum_thickness = get_momentum_thickness_off(diffusion_factor)
        loss_coefficient = (2 * momentum_thickness * (compressor.solidity / np.cos(alpha_exit)) *
                            (np.cos(off_design_stage.alpha) / np.cos(alpha_exit))**2)
    
        velocity_1 = np.sqrt(off_design_stage.vu**2 + off_design_stage.vm**2)
        pressure_drop = 0.5 * loss_coefficient * off_design_stage.rho * velocity_1**2
    
        T_total_exit = off_design_stage.T_tot
        p_total_exit = off_design_stage.p_tot - pressure_drop
    
        mach_number = get_mac(flow_area, T_total_exit, p_total_exit, alpha_exit, compressor.m_dot)
        T_static_exit = T_total_exit / get_f_mac(mach_number)
        p_static_exit = (p_total_exit * (T_static_exit / T_total_exit)**(gamma / (gamma - 1)))
    
        density_exit = p_static_exit / (R * T_static_exit)
        speed_of_sound = np.sqrt(gamma * R * T_static_exit)
        velocity_exit = mach_number * speed_of_sound
    
        tangential_velocity_exit = velocity_exit * np.sin(alpha_exit)
        axial_velocity_exit = velocity_exit * np.cos(alpha_exit)
        relative_tangential_velocity_exit = tangential_velocity_exit - off_design_stage.u
        beta_exit = np.arctan(relative_tangential_velocity_exit / axial_velocity_exit)
    
        # Vérification de la tolérance
        if np.abs(beta_exit - beta_exit_ref) <= beta_tolerance:
            break
    
        # Ajustement de alpha_exit
        if beta_exit < beta_exit_ref:
            alpha_exit += step_size
        else:
            alpha_exit -= step_size
            
        iteration_count+=1   
        
    T_total = T_static_exit + (axial_velocity_exit**2 + tangential_velocity_exit**2) / (2 * Cp)
    p_total = p_static_exit * (T_total / T_static_exit)**(gamma / (gamma - 1)) 

    stator_stage = cl.Stage(off_design_stage.stage_number + 1, T_static_exit, p_static_exit, 
                            T_total, p_total, density_exit, axial_velocity_exit, compressor, 
                            relative_tangential_velocity_exit, tangential_velocity_exit, beta_exit, alpha_exit, mach_number)

    return stator_stage
