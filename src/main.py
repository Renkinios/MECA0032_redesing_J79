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
# OPR =  layout_stage[-1].p_tot/layout_stage[0].p_tot 
print("############## Rotor Triangle ##############")
print(triangle_compressor)
# print("############## OPR ##############")
# print(OPR)

viz.viz_pressure_stage(layout_stage)
viz.viz_temperature_stage(layout_stage)
viz.viz_density_stage(layout_stage)
viz.viz_radius(layout_stage, compressor.R_mean)

print("############## Sensitivity solidity ##############")

# ---------------- Chossen blade -------------
stragger_rotor  = 32 * np.pi/180
stragger_stator = 22 * np.pi/180
stragger_OGV    = 15 * np.pi/180

blade_cascade = cl.blade_cascade(triangle_compressor, stragger_rotor, stragger_stator, stragger_OGV)


# ---------------- OFF DESIGN ----------------

print("############## Off Design ##############")
layout_stage_off = od.get_pitching_OFF_design(compressor, atm, layout_stage, blade_cascade)


compressor_10off    = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM *1.05, atm, Nbr_stage, reaction_degree, Work_Coeff)
layout_stage_off_10 = od.get_pitching_OFF_design(compressor_10off, atm, layout_stage, blade_cascade)

compreessor_10down  = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM*0.9, atm, Nbr_stage, reaction_degree, Work_Coeff)
layout_stage_down_10 = od.get_pitching_OFF_design(compreessor_10down, atm, layout_stage, blade_cascade)

viz.viz_off_design_severl_vm(layout_stage_off, layout_stage_off_10, layout_stage_down_10)
# viz.viz_off_design_vm(layout_stage_off)
# viz.viz_off_design_p_tot(layout_stage_off)
# viz.viz_off_design_T_tot(layout_stage_off)
# viz.viz_data_comp_off_and_on_p_tot(layout_stage_off, layout_stage)
# viz.viz_data_comp_off_and_on_T_tot(layout_stage_off, layout_stage)
viz.viz_off_design_several_p_tot(layout_stage, layout_stage_off_10, layout_stage_down_10)





