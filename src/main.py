import get_class as cl
import VizData as viz
import design.pitchine_design as pd
import design.difusion_interpolate as di
import numpy as np
import off_design.off_design as od
import matplotlib.pyplot as plt
from VSV.VSV_off_design import VSV_get_pitching_OFF_design, eff_OPR_VSV_compressor, max_min_m_dot_VSV_compressor


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
design_compressor = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM, atm, Nbr_stage, reaction_degree, Work_Coeff)
print("############## Rotor adimensionalisation ##############")

layout_stage, triangle_compressor = pd.get_pitching_design(design_compressor, atm) 
# j = 0
# for i in layout_stage:
#     print(f"A_{j} = {i.area}")
#     j +=1

di.viz_interpolation()

print(design_compressor)
OPR =  layout_stage[-1].p_tot/layout_stage[0].p_tot 
print("############## Rotor Triangle ##############")
print(triangle_compressor)
print("############## OPR ##############")
print(OPR)

viz.viz_pressure_stage(layout_stage)
viz.viz_temperature_stage(layout_stage)
viz.viz_density_stage(layout_stage)
viz.viz_radius(layout_stage, design_compressor.R_mean)


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


layout_stage_off_design = od.get_pitching_OFF_design(design_compressor, atm, layout_stage, blade_cascade)

m_dot_mat  = np.linspace(70, 87, 100)
layout_stage_off_mat = od.max_min_m_dot_off_design_compressor(m_dot_mat, 1, design_compressor, atm, layout_stage, blade_cascade)
path = "../figures/off_design/design_n"
viz.all_plot_out_design_stage(path,layout_stage_off_mat[0], layout_stage_off_mat[-1],layout_stage_off_design)


m_dot_mat  = np.linspace(80, 87, 100)



layout_stage_off_mat = od.max_min_m_dot_off_design_compressor(m_dot_mat, 1.05, design_compressor, atm, layout_stage, blade_cascade)
path = "../figures/off_design/design_n_1_05"
viz.all_plot_out_design_stage(path,layout_stage_off_mat[0], layout_stage_off_mat[-1])




m_dot_mat  = np.linspace(60, 65, 100)

layout_stage_off_mat = od.max_min_m_dot_off_design_compressor(m_dot_mat, 0.9, design_compressor, atm, layout_stage, blade_cascade)

path = "../figures/off_design/design_n_0_9"
viz.all_plot_out_design_stage(path,layout_stage_off_mat[0], layout_stage_off_mat[-1])



m_dot_array = np.linspace(60, 90, 500)
n_array = np.linspace(0.9, 1.05, 7)
eff_dico_n, OPR_dico_n = od.eff_OPR_off_design_compressor(m_dot_array, n_array, design_compressor, atm, layout_stage, blade_cascade)




# ---------------- VSV ----------------
print("############## VSV ##############")
m_dot_mat  = np.linspace(70, 87, 100)

layout_stage_off_mat = max_min_m_dot_VSV_compressor(m_dot_mat, 1, design_compressor, atm, layout_stage, blade_cascade)



path = "../figures/VSV/design_n"
viz.all_plot_VSV_stage(path,layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade,layout_stage_off_design)



m_dot_mat  = np.linspace(55, 80, 100)

layout_stage_off_mat = max_min_m_dot_VSV_compressor(m_dot_mat, 0.95, design_compressor, atm, layout_stage, blade_cascade)

path = "../figures/VSV/design_n_0_9"
viz.all_plot_VSV_stage(path,layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade,layout_stage_off_design)



m_dot_mat  = np.linspace(70, 98, 100)


layout_stage_off_mat = max_min_m_dot_VSV_compressor(m_dot_mat, 1.05, design_compressor, atm, layout_stage, blade_cascade)
path = "../figures/VSV/design_n_1_05"
viz.all_plot_VSV_stage(path,layout_stage_off_mat[0], layout_stage_off_mat[-1], blade_cascade, layout_stage_off_design)




dic_operationel_different_n = {}
m_dot_array = np.linspace(60, 90, 500)
n_array = np.linspace(0.9, 1.05, 7)

eff_dico_n, OPR_dico_n = eff_OPR_VSV_compressor(m_dot_array, n_array, design_compressor, atm, layout_stage, blade_cascade)
viz.viz_off_design_OPR(dic_operationel_different_n, OPR_dico_n, f"../figures/VSV/OPR.pdf")
viz.viz_off_design_eff(dic_operationel_different_n, eff_dico_n, f"../figures/VSV/eff.pdf")