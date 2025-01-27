import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines

# Définition globale des paramètres de police et de taille pour tous les graphiques
plt.rc('font', family='serif')  # Police avec empattements, comme Times
plt.rc('text', usetex=True)  # Utiliser LaTeX pour le texte dans les figures
plt.rcParams.update({
    'font.size': 18,         # Taille générale de la police
    'legend.fontsize': 35,   # Taille de la police des légendes
    'axes.labelsize': 50,    # Taille de la police des étiquettes des axes
    'axes.titlesize': 50,    # Taille de la police pour les titres des axes
    'xtick.labelsize': 20,   # Taille de la police pour les graduations de l'axe des X
    'ytick.labelsize': 20,   # Taille de la police pour les graduations de l'axe des Y
    'text.usetex': True,     # Confirmer l'utilisation de LaTeX pour tout le texte
    'figure.titlesize': 50   # Taille de la police pour les titres des figures
})
color_list = [
    "#007070", "#f07f3c", "#5b57a2", "#7db928", "#e62d31",
    "#005ca9", "#00843b", "#f8aa00", "#5b257d", "#8c8b82"
]
ligne_or_text = "#DCE6EA"

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
    plt.savefig(f"../figures/design/velocity_triangle_{file_type}.pdf", dpi=300, bbox_inches='tight')
    plt.close()

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rcParams.update({
    'font.size': 18,         # Taille générale de la police
    'legend.fontsize': 35,   # Taille de la police des légendes
    'axes.labelsize': 50,    # Taille de la police des étiquettes des axes
    'axes.titlesize': 50,    # Taille de la police pour les titres des axes
    'xtick.labelsize': 20,   # Taille de la police pour les graduations de l'axe des X
    'ytick.labelsize': 20,   # Taille de la police pour les graduations de l'axe des Y
    'text.usetex': True,     # Confirmer l'utilisation de LaTeX pour tout le texte
    'figure.titlesize': 50   # Taille de la police pour les titres des figures
})
def viz_pressure_stage(layout_stage) :
    plt.plot([stage.p_stat*1e-6 for stage in layout_stage], label=r"$p_{stat}$", color=color_list[0], linestyle='-.')
    plt.plot([stage.p_tot*1e-6 for stage in layout_stage], label=r"$p_{tot}$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Presure [MPa]")
    plt.legend()
    plt.xlim(0, len(layout_stage) - 1)
    plt.savefig(f"../figures/design/pressure_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_temperature_stage(layout_stage) :
        # Assurez-vous que les paramètres sont bien appliqués
    plt.plot([stage.T_stat for stage in layout_stage], label=r"$T_{stat}$", color=color_list[0], linestyle='-.')
    plt.plot([stage.T_tot for stage in layout_stage], label=r"$T_{tot}$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Temperature [K]")
    plt.xlim(0, len(layout_stage) - 1)
    plt.legend()
    plt.savefig(f"../figures/design/temperature_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_density_stage(layout_stage) :
    plt.plot([stage.rho for stage in layout_stage], color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Density [kg/m$^3$]")
    plt.xlim(0, len(layout_stage) - 1)
    plt.savefig(f"../figures/design/density_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_radius(layout_stage, R_mean) :
    plt.plot([stage.rtip for stage in layout_stage], label=r"$r_{tip}$", color=color_list[0])
    plt.plot([stage.rhub for stage in layout_stage], label=r"$r_{hub}$", color=color_list[1])
    plt.hlines(R_mean, 0, len(layout_stage), colors='black', linestyle='-.')
    plt.text(5 ,R_mean +0.005, r"Mean radius", color='black', fontsize=22, ha='center')
    plt.xlabel(r"Number of stage [-]")
    plt.xlim(0, len(layout_stage) - 1)
    plt.ylabel(r"Radius [m]")
    plt.legend()
    plt.savefig(f"../figures/design/radius_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()


def viz_off_design_vm(layout_stage_off) :
    plt.plot([stage.vm for stage in layout_stage_off], label=r"$v_m$", color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Mean velocity [m/s]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    plt.savefig(f"../figures/off_design/vm_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_off_design_T_tot(layout_stage_off) :
    plt.plot([stage.T_tot for stage in layout_stage_off], label=r"$T_{tot}$", color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total temperature [K]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/T_tot_stage.pdf", dpi=300, bbox_inches='tight')

    plt.close()

def viz_off_design_p_tot(layout_stage_off) :
    plt.plot([stage.p_tot for stage in layout_stage_off], label=r"$p_{tot}$", color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total pressure [Pa]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()


def viz_data_comp_off_and_on_p_tot(layout_stage_off, layout_stage) :
    plt.plot([stage.p_tot for stage in layout_stage_off], label=r"$p_{tot} off$", color=color_list[0])
    plt.plot([stage.p_tot for stage in layout_stage], label=r"$p_{tot} on$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total pressure [Pa]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_data_comp_off_and_on_T_tot(layout_stage_off, layout_stage) :
    plt.plot([stage.T_tot for stage in layout_stage_off], label=r"$T_{tot} off$", color=color_list[0])
    plt.plot([stage.T_tot for stage in layout_stage], label=r"$T_{tot} on$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total temperature [K]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    plt.savefig(f"../figures/off_design/T_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_off_design_several_p_tot(layout_stage, layout_stage_off, layout_stage_down) :
    plt.plot([stage.p_tot for stage in layout_stage], label=r"$n_{design}$", color=color_list[0])
    plt.plot([stage.p_tot for stage in layout_stage_off], label=r"$ 1.05 \cdot n_{design}$", color=color_list[1])
    plt.plot([stage.p_tot for stage in layout_stage_down], label=r"$0.9 \cdot n_{design}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total pressure [Pa]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')

def viz_off_design_severl_vm(layout_stage, layout_stage_off, layout_stage_down) :
    plt.plot([stage.vm for stage in layout_stage], label=r"$n_{design}$", color=color_list[0])
    plt.plot([stage.vm for stage in layout_stage_off], label=r"$ 1.05 \cdot n_{design}$", color=color_list[1])
    plt.plot([stage.vm for stage in layout_stage_down], label=r"$0.9 \cdot n_{design}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"$V_m [m/s]$")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_off_design_OPR(dic_operationel_different_n, OPR_dico, path):
    fig, ax = plt.subplots(figsize=(10, 6))
    lines = []  # Stocker les lignes pour l'étiquetage

    for key in dic_operationel_different_n.keys():
        # Assigner un label pour labellines et utiliser '_nolegend_' pour exclure de la légende
        label = fr'${key:.2f} n^*$'
        if len(dic_operationel_different_n[key]['m_dot']) == 0:
            continue
        line, = ax.plot(dic_operationel_different_n[key]['m_dot'], 
                        OPR_dico[key], 
                        color='black', 
                        label='_nolegend_')  # Exclure de la légende
        line.set_label(label)  # Définir le label pour labellines
        lines.append(line)

        if key == 1:
            m_dot_ope_n = dic_operationel_different_n[key]['m_dot']
            closest_m_dot = min(m_dot_ope_n, key=lambda x: abs(x - 77))
            ax.scatter(closest_m_dot, 
                       OPR_dico[key][m_dot_ope_n.index(closest_m_dot)], 
                       facecolors='none',  # Pas de remplissage
                       edgecolors='#800020',  # Couleur Bordeaux pour le contour
                       marker='o', 
                       zorder=5)

    # Placer les labels à la fin de chaque courbe
    labelLines(lines, zorder=2.5, align=False, xvals=[line.get_xdata()[-1] for line in lines])

    ax.set_xlabel(r"Mass flow rate $\dot{m}$ [kg/s]")
    ax.set_ylabel(r"Pressure ratio $\Pi$ [-]")

    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def viz_off_design_eff(dic_operationel_different_n, eff_dico_n, path):
    fig, ax = plt.subplots(figsize=(10, 6))
    lines = []  # Stocker les lignes pour l'étiquetage

    for key in dic_operationel_different_n.keys():
        # Assigner un label pour labellines et utiliser '_nolegend_' pour exclure de la légende
        label = fr'${key:.2f} n^*$'
        if len(dic_operationel_different_n[key]['m_dot']) == 0:
            continue
        line, = ax.plot(dic_operationel_different_n[key]['m_dot'], 
                        eff_dico_n[key], 
                        color='black', 
                        label='_nolegend_')  # Exclure de la légende
        line.set_label(label)  # Définir le label pour labellines
        lines.append(line)

        if key == 1:
            m_dot_ope_n = dic_operationel_different_n[key]['m_dot']
            closest_m_dot = min(m_dot_ope_n, key=lambda x: abs(x - 77))
            ax.scatter(closest_m_dot, 
                       eff_dico_n[key][m_dot_ope_n.index(closest_m_dot)], 
                       facecolors='none',  # Pas de remplissage
                       edgecolors='#800020',  # Couleur Bordeaux pour le contour
                       marker='o', 
                       zorder=5)

    # Placer les labels à la fin de chaque courbe
    labelLines(lines, zorder=2.5, align=False, xvals=[line.get_xdata()[-1] for line in lines])

    ax.set_xlabel(r"Mass flow rate $\dot{m}$ [kg/s]")
    ax.set_ylabel(r"Efficiency $\eta$ [\%]")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def viz_vm_stall_operaing_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([stage.vm for stage in layout_stag_design], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([stage.vm for stage in layout_stage_off_neg_stall], label=r"$\dot{m}_{min}$" , color=color_list[1])
    plt.plot([stage.vm for stage in layout_stage_off_pos_stall], label=r"$\dot{m}_{max}$" , color=color_list[2])

    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Mean velocity [m/s]")
    plt.xlim(0, len(layout_stage_off_neg_stall) - 1)
    plt.legend()
    plt.savefig(f"{path}/vm_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_mac_stall_operaing_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([stage.mac for stage in layout_stag_design[1:]], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([stage.mac for stage in layout_stage_off_neg_stall[1:]], label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot([stage.mac for stage in layout_stage_off_pos_stall[1:]], label=r"$\dot{m}_{max}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Mach number [-]")
    plt.xlim(1, len(layout_stage_off_neg_stall) - 2)
    plt.legend()
    plt.savefig(f"{path}/mac_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_p_tot_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([stage.p_tot for stage in layout_stag_design], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([stage.p_tot for stage in layout_stage_off_neg_stall], label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot([stage.p_tot for stage in layout_stage_off_pos_stall], label=r"$\dot{m}_{max}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total pressure [Pa]")
    plt.xlim(0, len(layout_stage_off_neg_stall) - 1)
    plt.legend()
    plt.savefig(f"{path}/p_tot_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_T_tot_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([stage.T_tot for stage in layout_stag_design], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([stage.T_tot for stage in layout_stage_off_neg_stall], label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot([stage.T_tot for stage in layout_stage_off_pos_stall], label=r"$\dot{m}_{max}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total temperature [K]")
    plt.xlim(0, len(layout_stage_off_neg_stall) - 1)
    plt.legend()
    plt.savefig(f"{path}/T_tot_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_p_stat_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([stage.p_stat for stage in layout_stag_design], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([stage.p_stat for stage in layout_stage_off_neg_stall], label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot([stage.p_stat for stage in layout_stage_off_pos_stall], label=r"$\dot{m}_{max}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Static pressure [Pa]")
    plt.xlim(0, len(layout_stage_off_neg_stall) - 1)
    plt.legend()
    plt.savefig(f"{path}/p_stat_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_T_stat_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([stage.T_stat for stage in layout_stag_design], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([stage.T_stat for stage in layout_stage_off_neg_stall], label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot([stage.T_stat for stage in layout_stage_off_pos_stall], label=r"$\dot{m}_{max}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Static temperature [K]")
    plt.xlim(0, len(layout_stage_off_neg_stall) - 1)
    plt.legend()
    plt.savefig(f"{path}/T_stat_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_alpha_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([np.degrees(stage.alpha) for stage in layout_stag_design[1:]], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([np.degrees(stage.alpha) for stage in layout_stage_off_neg_stall[1:]], label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot([np.degrees(stage.alpha) for stage in layout_stage_off_pos_stall[1:]], label=r"$\dot{m}_{max}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r" $\alpha [^\circ]$")
    plt.xlim(0, len(layout_stage_off_neg_stall) - 1)
    plt.legend()
    plt.savefig(f"{path}/alpha_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_alpha_stator_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, blade_cascade, layout_stag_design = []):
    # Calculer les angles alpha
    if len(layout_stag_design) > 0:
        alpha_design = [np.degrees(stage.alpha) for stage in layout_stag_design[2::2]]
    alpha_neg_stall = [np.degrees(stage.alpha) for stage in layout_stage_off_neg_stall[2::2]]
    alpha_pos_stall = [np.degrees(stage.alpha) for stage in layout_stage_off_pos_stall[2::2]]

    x_labels = np.arange(1, len(alpha_neg_stall) + 1)

        # Ajouter les lignes de position de décrochage
    plt.hlines(np.degrees(blade_cascade.stall_stator_neg), x_labels[0], x_labels[-1], colors='black', linestyle='-.',color = color_list[3])
    plt.hlines(np.degrees(blade_cascade.stall_stator_pos), x_labels[0], x_labels[-1], colors='black', linestyle='-.',color = color_list[3], label = "Stall position")
    # Tracer les courbes
    if len(layout_stag_design) > 0:
        plt.plot(x_labels, alpha_design, label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot(x_labels, alpha_neg_stall, label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot(x_labels, alpha_pos_stall, label=r"$\dot{m}_{max}$", color=color_list[2])


    # Configurer les axes
    plt.xlabel(r"Number of stator [-]")
    plt.ylabel(r"$\alpha [^\circ]$")
    plt.xticks(x_labels)  # S'assurer que les ticks sont alignés avec 1, 2, 3, ...
    plt.xlim(x_labels[0], x_labels[-1])
    plt.legend()

    # Sauvegarder et fermer
    plt.savefig(f"{path}/alpha_stall_stator.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_beta_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design = []):
    if len(layout_stag_design) > 0:
        plt.plot([-np.degrees(stage.beta) for stage in layout_stag_design[1:]], label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot([-np.degrees(stage.beta) for stage in layout_stage_off_neg_stall[1:]], label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot([-np.degrees(stage.beta) for stage in layout_stage_off_pos_stall[1:]], label=r"$\dot{m}_{max}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"$- \beta [^\circ]$")
    plt.xlim(0, len(layout_stage_off_neg_stall) - 1)
    plt.legend()
    plt.savefig(f"{path}/beta_stall_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_beta_rotor_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, blade_cascade, layout_stag_design = []):

    # Calculer les angles beta
    if len(layout_stag_design) > 0:
        beta_design = [-np.degrees(stage.beta) for stage in layout_stag_design[1:len(layout_stag_design)-1:2]]
    beta_neg_stall = [-np.degrees(stage.beta) for stage in layout_stage_off_neg_stall[1:len(layout_stage_off_neg_stall)-1:2]]
    beta_pos_stall = [-np.degrees(stage.beta) for stage in layout_stage_off_pos_stall[1:len(layout_stage_off_pos_stall)-1:2]]

    x_labels = np.arange(1, len(beta_neg_stall) + 1)
    # Ajouter les lignes de position de décrochage
    plt.hlines(np.degrees(blade_cascade.stall_rotor_neg), x_labels[0], x_labels[-1], colors='black', linestyle='-.',color = color_list[3])
    plt.hlines(np.degrees(blade_cascade.stall_rotor_pos), x_labels[0], x_labels[-1], colors='black', linestyle='-.',color = color_list[3], label = "Stall position")
    # Tracer les courbes
    if len(layout_stag_design) > 0:
        plt.plot(x_labels, beta_design, label=r"$\dot{m}_{design}$", color=color_list[0])
    plt.plot(x_labels, beta_neg_stall, label=r"$\dot{m}_{min}$", color=color_list[1])
    plt.plot(x_labels, beta_pos_stall, label=r"$\dot{m}_{max}$", color=color_list[2])


    # Configurer les axes
    plt.xlabel(r"Number of rotor [-]")
    plt.ylabel(r"-$\beta [^\circ]$")
    plt.xticks(x_labels)  # S'assurer que les ticks sont alignés avec 1, 2, 3, ...
    plt.xlim(x_labels[0], x_labels[-1])
    plt.legend()

    # Sauvegarder et fermer
    plt.savefig(f"{path}/beta_stall_rotor.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_VSV_stragger_rotor(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,  blade_cascade) :

    
    stragger_neg_stall = [np.degrees(stage.stagger) for stage in layout_stage_off_neg_stall[3:-1:2]]
    stragger_pos_stall = [np.degrees(stage.stagger) for stage in layout_stage_off_pos_stall[3:-1:2]]
    x_labels = np.arange(1, len(stragger_neg_stall) + 1)
    plt.hlines(np.degrees(blade_cascade.stragger_stator), x_labels[0], x_labels[-1], colors=color_list[2], linestyle='-.', label=r"design stragger ")
    plt.plot(x_labels, stragger_neg_stall, label=r"$\dot{m}_{min}$", color=color_list[0])
    plt.plot(x_labels, stragger_pos_stall, label=r"$\dot{m}_{max}$", color=color_list[1])
    plt.xlim(x_labels[0], x_labels[-1])
    plt.xlabel(r"Number of rotor [-]")
    plt.ylabel(r"Stagger angle [°]")
    plt.legend()
    plt.savefig(f"{path}/VSV_stagger_rotor.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def all_plot_VSV_stage(path,layout_stage_off_neg_stall, layout_stage_off_pos_stall,blade_cascade,layout_stag_design = []):
    viz_vm_stall_operaing_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_mac_stall_operaing_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_p_tot_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_T_tot_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_p_stat_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_T_stat_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_alpha_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_beta_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design)
    viz_alpha_stator_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, blade_cascade,layout_stag_design)
    viz_beta_rotor_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, blade_cascade,layout_stag_design)


def all_plot_out_design_stage(path,layout_stage_off_neg_stall, layout_stage_off_pos_stall,layout_stag_design = []):
    viz_vm_stall_operaing_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
    viz_mac_stall_operaing_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
    viz_p_tot_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
    viz_T_tot_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
    viz_p_stat_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
    viz_T_stat_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
    viz_alpha_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
    viz_beta_stall_operating_point(path, layout_stage_off_neg_stall, layout_stage_off_pos_stall, layout_stag_design)
