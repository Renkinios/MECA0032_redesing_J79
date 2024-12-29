import get_class as cl
import VizData as viz
import design.pitchine_design as pd
import design.difusion_interpolate as di
import numpy as np
import off_design.off_design as od
import matplotlib.pyplot as plt


# Good values for the engine engine have  original, R_hub ~ 0.1 m et Omega = 7385 RPM
m_dot = 77
R_tip = 0.4
reaction_degree = 0.6
solidity = 1.5
poly_eff = 0.85

Work_Coeff =  0.5
RPM = 7103
Nbr_stage = 14
R_hub = 0.18888888888888888


atm = cl.atmosphere()
compressor = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM, atm, Nbr_stage, reaction_degree, Work_Coeff)
print("############## Rotor adimensionalisation ##############")

layout_stage, triangle_compressor = pd.get_pitching_design(compressor, atm) 
# j = 0
# for i in layout_stage:
#     print(f"A_{j} = {i.area}")
#     j +=1

di.viz_interpolation()

print(compressor)
OPR =  layout_stage[-1].p_tot/layout_stage[0].p_tot 
print("############## Rotor Triangle ##############")
print(triangle_compressor)
print("############## OPR ##############")
print(OPR)

viz.viz_pressure_stage(layout_stage)
viz.viz_temperature_stage(layout_stage)
viz.viz_density_stage(layout_stage)
viz.viz_radius(layout_stage, compressor.R_mean)


# ---------------- Chossen blade -------------
stragger_rotor  = np.deg2rad(34)
stall_pos_rotor = np.deg2rad(58)
stall_neg_rotor = np.deg2rad(41)

stragger_stator  = np.deg2rad(24)
stall_pos_stator = np.deg2rad(51)
stall_neg_stator = np.deg2rad(36)

stragger_OGV    = np.deg2rad(15)
stall_pos_OGV   = np.deg2rad(48)
stall_neg_OGV   = np.deg2rad(30)


blade_cascade = cl.blade_cascade(triangle_compressor, stragger_rotor, stragger_stator, stragger_OGV, stall_neg_rotor, stall_pos_rotor, 
                                 stall_neg_stator, stall_pos_stator, stall_neg_OGV, stall_pos_OGV)


# ---------------- OFF DESIGN ----------------

print("############## Off Design ##############")

m_dot_mat  = np.linspace(70, 87, 100)

valid_m_dot = []
compressors_mat = []
layout_stage_off_mat = []

for m_d in m_dot_mat:   
    try:
        compressor = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_d, RPM, atm, Nbr_stage, reaction_degree, Work_Coeff)
        layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)
        valid_m_dot.append(m_d)
        compressors_mat.append(compressor)
        layout_stage_off_mat.append(layout_stage_off)
    except Exception as e:
        continue

compressor       = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM, atm, Nbr_stage, reaction_degree, Work_Coeff)
layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)


path = "../figures/off_design/design_n"
viz.viz_vm_stall_operaing_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_mac_stall_operaing_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_p_tot_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_T_tot_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_p_stat_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_T_stat_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_alpha_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_beta_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], layout_stage_off)
viz.viz_alpha_stator_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade, layout_stage_off)
viz.viz_beta_rotor_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade, layout_stage_off)


m_dot_mat  = np.linspace(80, 87, 100)

valid_m_dot = []
compressors_mat = []
layout_stage_off_mat = []

for m_d in m_dot_mat:   
    try:
        compressor = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_d, RPM*1.05, atm, Nbr_stage, reaction_degree, Work_Coeff)
        layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)
        valid_m_dot.append(m_d)
        compressors_mat.append(compressor)
        layout_stage_off_mat.append(layout_stage_off)
    except Exception as e:
        continue

path = "../figures/off_design/design_n_1_05"
viz.viz_vm_stall_operaing_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_mac_stall_operaing_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_p_tot_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_T_tot_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_p_stat_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_T_stat_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_alpha_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_beta_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_alpha_stator_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade)
viz.viz_beta_rotor_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade)



m_dot_mat  = np.linspace(60, 65, 100)

valid_m_dot = []
compressors_mat = []
layout_stage_off_mat = []

for m_d in m_dot_mat:   
    try:
        compressor = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_d, RPM*0.9, atm, Nbr_stage, reaction_degree, Work_Coeff)
        layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)
        valid_m_dot.append(m_d)
        compressors_mat.append(compressor)
        layout_stage_off_mat.append(layout_stage_off)
    except Exception as e:
        continue


path = "../figures/off_design/design_n_0_9"
viz.viz_vm_stall_operaing_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_mac_stall_operaing_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_p_tot_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_T_tot_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_p_stat_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_T_stat_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_alpha_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_beta_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1])
viz.viz_alpha_stator_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade)
viz.viz_beta_rotor_stall_operating_point(path, layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade)

# compressor_10off    = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM *1.05, atm, Nbr_stage, reaction_degree, Work_Coeff)
# layout_stage_off_10 = od.get_pitching_OFF_design(compressor_10off, atm, layout_stage, blade_cascade)

# compreessor_10down   = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM*0.9, atm, Nbr_stage, reaction_degree, Work_Coeff)
# layout_stage_down_10 = od.get_pitching_OFF_design(compreessor_10down, atm, layout_stage, blade_cascade)

# viz.viz_off_design_severl_vm(layout_stage_off, layout_stage_off_10, layout_stage_down_10)
# viz.viz_off_design_vm(layout_stage_off)
# viz.viz_off_design_p_tot(layout_stage_down_10)
# viz.viz_off_design_T_tot(layout_stage_off)



# viz.viz_data_comp_off_and_on_p_tot(layout_stage_off, layout_stage)
# viz.viz_data_comp_off_and_on_T_tot(layout_stage_off, layout_stage)
# viz.viz_off_design_several_p_tot(layout_stage, layout_stage_off_10, layout_stage_down_10)
dic_operationel_different_n = {}
m_dot_array = np.linspace(60, 90, 500)
n_array = np.linspace(0.9, 1.05, 7)
for j in n_array: 
    dic_cst_n = {'m_dot': [], 'layout_stage_off': []}
    for m_d in (m_dot_array) :
        try:
            compressor       = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_d, RPM*j, atm, Nbr_stage, reaction_degree, Work_Coeff)
            layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)
            dic_cst_n['layout_stage_off'].append(layout_stage_off)
            dic_cst_n['m_dot'].append(m_d)
            viz.viz_off_design_several_p_tot(layout_stage, layout_stage_off, m_d, j)
        except Exception as e:
            pass
    dic_operationel_different_n[j] = dic_cst_n

def eff_OPR_compressor(dic_operationel_different_n) :
    eff_dico_n = {}
    OPR_dico_n = {}
    for key in dic_operationel_different_n.keys():
        eff_dico = []
        OPR_dico = []
        stage_cst_n  = dic_operationel_different_n[key]['layout_stage_off']
        for comp in dic_operationel_different_n[key]['layout_stage_off']:
            OPR = (comp[-1].p_tot/comp[0].p_tot)
            OPR_dico.append(OPR)
            T_out_iso = comp[0].T_tot * OPR**((atm.gamma - 1)/atm.gamma)
            eff_comp_iso = (T_out_iso - comp[0].T_tot)/(comp[-1].T_tot - comp[0].T_tot) *100
            eff_dico.append(eff_comp_iso)
        eff_dico_n[key] = eff_dico
        OPR_dico_n[key] = OPR_dico

    return eff_dico_n, OPR_dico_n

eff_dico_n, OPR_dico_n = eff_OPR_compressor(dic_operationel_different_n)

viz.viz_off_design_OPR(dic_operationel_different_n, OPR_dico_n)
viz.viz_off_design_eff(dic_operationel_different_n, eff_dico_n)