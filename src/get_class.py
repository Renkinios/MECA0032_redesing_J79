import numpy as np
import matplotlib.pyplot as plt

class atmosphere :
    def __init__(self):
        self.Ta = 288.15  # K
        self.Pa = 101325     #Pa
        self.gamma = 1.4
        self.R = 287.1
        self.Cp = self.R * self.gamma / (self.gamma - 1)
        self.rho = 1.225          # sea level


class station_compressor:
    def __init__(self, stage_number, T_stat, p_stat, T_tot, p_tot, rho, vm, compressor, w_u=None, v_u=None, beta=None, alpha=None, mac=None,stagger=None):
 
        self.stage_number = stage_number
        self.T_stat = T_stat
        self.p_stat = p_stat
        self.rho = rho
        self.wu = w_u
        self.vu = v_u
        self.vm = vm
        self.u = compressor.RPM * 2 * np.pi / 60 * compressor.R_mean
        self.T_tot = T_tot
        self.p_tot = p_tot
        self.area = compressor.m_dot/(self.rho * vm)
        self.rtip = compressor.R_mean + self.area/(4 * np.pi * compressor.R_mean)
        self.rhub = compressor.R_mean - self.area/(4 * np.pi * compressor.R_mean)
        self.beta = beta
        self.alpha = alpha
        self.mac = mac
        self.stagger = stagger

    def __str__(self):
        """
        Provide a detailed string representation of the stage.
        """
        return (
            f"############ Stage {self.stage_number} Summary ############\n"
            f"  Static Temperature (T_stat): {self.T_stat:.2f} K\n"
            f"  Static Pressure (p_stat): {self.p_stat:.2f} Pa\n"
            f"  Total Temperature (T_tot): {self.T_tot:.2f} K\n"
            f"  Total Pressure (p_tot): {self.p_tot:.2f} Pa\n"
            f"  Axial Velocity (v_u): {self.vu:.2f} m/s\n" if self.vu is not None else "  Axial Velocity (v_u): None\n"
            f"  Tangential Velocity (w_u): {self.wu:.2f} m/s\n" if self.wu is not None else "  Tangential Velocity (w_u): None\n"
            f"  Meriodonial velocity (vm): {self.vm:.2f} m/s\n" 
            f"  Density (rho): {self.rho:.2f} kg/m^3\n"
            f"  Tip Radius (r_tip): {self.rtip:.4f} m\n"
            f"  Hub Radius (r_hub): {self.rhub:.4f} m\n"
            f"  Area: {self.area:.4f} m^2\n" 
            f"  Beta: {np.degrees(self.beta):.2f}°\n" if self.beta is not None else "  Beta: None\n"
            f"  Alpha: {np.degrees(self.alpha):.2f}°" if self.alpha is not None else "  Alpha: None"
        )



class Compressor:
    """
    Class to define the compressor geometry and operating conditions.
    """
    def __init__(self, solidity, eff_poly, R_hub, R_tip, mass_flux, RPM, atm, nb_stage, degree_reaction, work_coeff):
        """
        Initialize the compressor with its physical parameters.

        Args:
            solidity (float): Blade solidity ratio.
            eff_poly (float): Polytropic efficiency.
            R_hub (float): Radius of the hub (m).
            R_tip (float): Radius of the tip (m).
            mass_flux (float): Mass flow rate (kg/s).
            RPM (float): Rotational speed (revolutions per minute).
            atm: Atmospheric conditions (must have a `rho` attribute for density).
            nb_stage (int): Number of compressor stages.
            degree_reaction (float): Degree of reaction (dimensionless, typically between 0 and 1).
            work_coeff (float): Work coefficient (dimensionless).
        """
        # Validate input parameters
        if R_hub <= 0 or R_tip <= 0:
            raise ValueError("R_hub and R_tip must be positive.")
        if R_hub >= R_tip:
            raise ValueError("R_hub must be smaller than R_tip.")
        if not (0 <= degree_reaction <= 1):
            raise ValueError("Degree of reaction must be between 0 and 1.")

        # Initialize compressor attributes
        self.solidity = solidity
        self.eff_poly = eff_poly
        self.R_mean = (R_hub + R_tip) / 2
        self.m_dot = mass_flux
        self.RPM = RPM
        self.number_stage = nb_stage
        self.R = degree_reaction
        self.work_coeff = work_coeff
        self.R_tip = R_tip
        self.R_hub = R_hub

        # Calculate derived parameters
        self.area_entrance = np.pi * (R_tip**2 - R_hub**2)  # Annular area (m^2)


                # Calculate velocities
        v_m = self.m_dot / (self.area_entrance * atm.rho)  # Axial velocity (m/s)
        u   = (self.RPM * 2 * np.pi) / 60 * self.R_mean  # Tip speed (m/s)
        self.flow_coeff = v_m / u  # Flow coefficient

        if self.flow_coeff == 0:
            raise ValueError("Flow coefficient cannot be zero.")

    def __str__(self):
        """
        Provide a detailed string representation of the compressor.
        """
        return (
            f"########### Compressor Summary ############\n"
            f"  Mean Radius (R_mean): {self.R_mean:.3f} m\n"
            f"  Work Coefficient (psi): {self.work_coeff:.3f}\n"
            f"  Degree of Reaction (R): {self.R:.3f}\n"
            f"  Flow Coefficient (phi): {self.flow_coeff:.3f}\n"
            f"  Mass Flux (m_dot): {self.m_dot:.2f} kg/s\n"
            f"  Rotational Speed (RPM): {self.RPM}\n"
            f"  Number of Stages: {self.number_stage}\n"
            f"  Solidity: {self.solidity}\n"
            f"  Polytropic Efficiency: {self.eff_poly:.2f}"
        )


class CompressorTriangle:
    """
    Class to define the velocity triangle of a compressor rotor.
    """
    def __init__(self, compressor, atm):
        """
        Initialize the velocity triangle based on compressor parameters.

        Args:
            compressor: An instance of the `Compressor` class containing
                        required geometric and operational parameters.
        """
        self.vm = compressor.m_dot / (compressor.area_entrance * atm.rho)
        self.u = (compressor.RPM * 2 * np.pi) / 60 * compressor.R_mean

        # Velocity Triangles
        v_u_1 = self.vm * (((1 - compressor.R) - (compressor.work_coeff / 2)) / compressor.flow_coeff)
        w_u_1 = self.u - v_u_1
        v_u_2 = (1 - compressor.R + compressor.work_coeff / 2) * self.u
        w_u_2 = self.u - v_u_2

        self.beta_1 = -np.arctan((self.u - v_u_1) / self.vm)
        self.beta_2 = -np.arctan((self.u - v_u_2) / self.vm)
        self.alpha_1 = np.arctan(v_u_1 / self.vm)
        self.alpha_2 = np.arctan(v_u_2 / self.vm)

        # Set beta_outlet_OGV as None by default
        self.beta_outlet_OGV = None
        self.alpha_outlet_OGV = 0

    def __str__(self, degree=True):
        if degree:
            alpha_1 = np.degrees(self.alpha_1)
            alpha_2 = np.degrees(self.alpha_2)
            beta_1 = np.degrees(self.beta_1)
            beta_2 = np.degrees(self.beta_2)
            unit = "°"
            beta_outlet_OGV = (
                f"{np.degrees(self.beta_outlet_OGV):.2f}{unit}"
                if self.beta_outlet_OGV is not None
                else "No value selected"
            )
        else:
            alpha_1 = self.alpha_1
            alpha_2 = self.alpha_2
            beta_1 = self.beta_1
            beta_2 = self.beta_2
            unit = "rad"
            beta_outlet_OGV = (
                f"{self.beta_outlet_OGV:.2f}{unit}"
                if self.beta_outlet_OGV is not None
                else "No value selected"
            )

        return (
            "########### Compressor Triangle Summary ############\n"
            f"alpha_1: {alpha_1:.2f}{unit}, alpha_2: {alpha_2:.2f}{unit}\n"
            f"beta_1: {beta_1:.2f}{unit}, beta_2: {beta_2:.2f}{unit}\n"
            f"Outlet OGV angle: beta: {beta_outlet_OGV}, alpha: {self.alpha_outlet_OGV:.2f}{unit}\n"
            f"vm: {self.vm:.2f} m/s, u: {self.u:.2f} m/s"
        )


    def viz_triangle_velocity(v_m, w_u, v_u, inlet =True):
        nbr  = 1 if inlet else 2
        sign = 1 if inlet else -1
        

        plt.figure(figsize=(8, 6))

        # Ligne verte (projection verticale de V)
        if inlet:
            plt.quiver(v_m, -w_u, 0,np.abs(v_u + w_u), angles='xy', scale_units='xy', scale=1, color='green')
            plt.text(v_m -  5, -w_u + (v_u + w_u)/2 , fr'$u_{{{nbr}}}$', color='green', fontsize=25, ha='center')
        else:
            plt.quiver(v_m, v_u, 0, np.abs(v_u + w_u), angles='xy', scale_units='xy', scale=1, color='green')
            plt.text(v_m -  5, -w_u + (v_u + w_u)/2 , fr'$u_{{{nbr}}}$', color='green', fontsize=25, ha='center')
        # Ligne bleue (V - vitesse résultante)
        plt.quiver(0, 0, v_m, v_u, angles='xy', scale_units='xy', scale=1, color='navy')
        plt.text(v_m / 2 , v_u / 2 +12 *sign, fr'$v_{{{nbr}}}$', color='navy', fontsize=25, ha='center')

        # Ligne rouge (W - vitesse de glissement)
        plt.quiver(0, 0, v_m, -w_u, angles='xy', scale_units='xy', scale=1, color='red')
        plt.text(v_m / 2, -w_u / 2 - 20 *sign, fr'$w_{{{nbr}}}$', color='red', fontsize=25, ha='center')

        # Configuration des axes
        plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
        plt.axvline(0, color='black', linewidth=0.5, linestyle='--')
        plt.quiver(0, 0, v_m, 0, angles='xy', scale_units='xy', scale=1, color='black')
        plt.text(v_m / 2 , +6, fr'$v_{{m,{nbr}}}$', color='black', fontsize=25, ha='center')

        # Ajustement des limites du graphique
        max_x = max(abs(v_m), abs(v_m)) + 1
        max_y = max(abs(v_u), abs(-w_u)) + 1
        plt.xlim(-1, max_x)
        plt.ylim(-max_y, max_y + 5)

        # Légendes et titres
        plt.xlabel(r"Axe horizontal [m/s]")
        plt.ylabel(r"Axe vertical [m/s]")
        file_type = "inlet" if  inlet else "outlet"
        # Sauvegarde du fichier avec l'indice 'nbr'
        plt.savefig(f"../figures/velocity_triangle_{file_type}.pdf", dpi=300, bbox_inches='tight')
        plt.close()


class blade_cascade :
    def __init__(self, triangle, stragger_stator, stragger_rotor, stragger_OGV, stall_negative_rotor, stall_positive_rotor,
                  stall_negative_stator, stall_positive_stator, stall_negative_OGV, stall_positive_OGV):
        
        self.stragger_stator   = stragger_stator
        self.stragger_rotor    = stragger_rotor
        self.stragger_OGV      = stragger_OGV
        self.aoa_design_stator = triangle.alpha_2 - stragger_stator
        self.aoa_design_rotor  = -triangle.beta_1 - stragger_rotor
        self.aoa_design_OGV    = triangle.alpha_2 - stragger_OGV

        self.deflection_design_stator = np.abs(triangle.alpha_1 - triangle.alpha_2)  # inverse due to the fact that the stator is at the second position
        self.deflection_design_rotor  = np.abs(triangle.beta_2 - triangle.beta_1)
        self.deflection_design_OGV    = np.abs(triangle.alpha_2)
        self.beta1_design      = triangle.beta_1
        self.alpha1_design     = triangle.alpha_1
        self.beta2_design      = triangle.beta_2
        self.alpha2_design     = triangle.alpha_2

        self.stall_rotor_neg  = stall_negative_rotor
        self.stall_rotor_pos  = stall_positive_rotor
        self.stall_stator_neg = stall_negative_stator
        self.stall_stator_pos = stall_positive_stator
        self.stall_OGV_neg    = stall_negative_OGV
        self.stall_OGV_pos    = stall_positive_OGV

        self.u = triangle.u
        self.vm = triangle.vm
