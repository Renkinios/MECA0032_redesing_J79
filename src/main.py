import get_class as cl
import VizData as viz
import design.pitchine_design as pd
import design.difusion_interpolate as di
import numpy as np
import off_design.off_design as od


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

print("############## Sensitivity solidity ##############")

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

# m_dot  = np.linspace(60, 90, 100)

# valid_m_dot = []
# compressors_mat = []
# layout_stage_off_mat = []

# for m_d in m_dot:   
#     try:
#         print(f"m_dot = {m_d}")
#         compressor = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_d, RPM, atm, Nbr_stage, reaction_degree, Work_Coeff)
#         layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)
#         valid_m_dot.append(m_d)
#         compressors_mat.append(compressor)
#         layout_stage_off_mat.append(layout_stage_off)
#     except Exception as e:
#         print(f"Skipping m_dot = {m_d} due to error: {e}")
#         continue

compressor       = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM, atm, Nbr_stage, reaction_degree, Work_Coeff)
layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)

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