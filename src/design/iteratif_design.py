import get_class as cl
import VizData as viz
import pitchine_design as pd
import design.difusion_interpolate as di
import numpy as np

# Contraintes et paramètres initiaux
R_hub = 0.1
R_tip = 0.4
m_dot = 77
reaction_degree = 0.5
solidity = 1.2
poly_eff = 0.85
atm = cl.atmosphere()


# Limites des paramètres
work_coeff_range = (0.5, 1)
nbr_stage_range = (10,17)
rpm_range = (5000, 7385)
R_hub = (0.1,0.2)

# Variables à modifier
work_coeff_values = np.linspace(work_coeff_range[0], work_coeff_range[1], 10)
nbr_stage_values = range(nbr_stage_range[0], nbr_stage_range[1] + 1)
rpm_values = np.linspace(rpm_range[0], rpm_range[1], 10)
R_hub_values = np.linspace(R_hub[0], R_hub[1], 10)

# Fonction pour tester une configuration
def test_compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM, atm, nbr_stage, reaction_degree, work_coeff):
    try:
        # Création du compresseur
        compressor = cl.Compressor(solidity, poly_eff, R_hub, R_tip, m_dot, RPM, atm, nbr_stage, reaction_degree, work_coeff)
        # Calcul des étages et des triangles de vitesse
        layout_stage, triangle_compressor = pd.get_pitching_design(compressor, atm)

        # Calcul de l'OPR
        OPR = layout_stage[-1].p_tot / layout_stage[0].p_tot
        # Vérification de l'OPR
        if (13 <= OPR <= 14) and compressor.flow_coeff < 0.8:
            return compressor, layout_stage, triangle_compressor, OPR
        
        else:
            return None
    except Exception as e:
        print(f"Error encountered: {e}")
        return None

# Liste pour stocker toutes les configurations valides
valid_configurations = []

# Boucle pour trouver toutes les configurations valides
for nbr_stage in nbr_stage_values:
    for work_coeff in work_coeff_values:
        for RPM in rpm_values:
            for solidity in [1,1.5] :
                for R_hub in R_hub_values:

                    # print(f"Testing: Work_Coeff={work_coeff}, Nbr_stage={nbr_stage}, RPM={RPM}")
                    result = test_compressor(
                        solidity, poly_eff, R_hub, R_tip, m_dot, RPM, atm, nbr_stage, reaction_degree, work_coeff
                    )

                    if result is not None:
                        compressor, layout_stage, triangle_compressor, OPR = result
                        print(f"Valid configuration found! OPR={OPR:.2f}")
                        valid_configurations.append({
                            "Work_Coeff": work_coeff,
                            "Flow_coeff": compressor.flow_coeff,
                            "Nbr_stage": nbr_stage,
                            "RPM": RPM,
                            "OPR": OPR,
                            "Solidity": solidity,
                            "R_hub": R_hub,
                        })

# Affichage des résultats
if valid_configurations:
    print("All valid configurations:")
    for config in valid_configurations:
        print(config)
else:
    print("No valid configurations found within the given constraints.")
