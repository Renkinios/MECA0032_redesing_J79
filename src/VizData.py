import numpy as np
import matplotlib.pyplot as plt


# Définition globale des paramètres de police et de taille pour tous les graphiques
plt.rc('font', family='serif')  # Police avec empattements, comme Times
plt.rc('text', usetex=True)  # Utiliser LaTeX pour le texte dans les figures
plt.rcParams.update({
    'font.size': 14,       # Taille de police générale
    'legend.fontsize': 14, # Taille de police pour les légendes
    'axes.labelsize': 17,  # Taille de police pour les étiquettes des axes
})

color_list = [
    "#1b4f72",  # Bleu foncé (élégant et classique)
    "#c0392b",  # Rouge foncé (puissant et contrasté)
    "#196f3d",  # Vert foncé (nature et stabilité)
    "#7d3c98",  # Violet foncé (subtil et élégant)
    "#b9770e",  # Orange foncé (chaleureux et dynamique)
    "#5d6d7e",  # Gris-bleu foncé (neutre et moderne)
    "#34495e",  # Bleu nuit (professionnel et sobre)
    "#7e5109",  # Marron foncé (neutre et sérieux)
    "#117864",  # Vert émeraude foncé (frais et équilibré)
    "#6c3483",  # Violet profond (raffiné et audacieux)
    "#884ea0",  # Violet moyen (doux et épuré)
    "#2874a6",  # Bleu profond (épuré et professionnel)
    "#d35400",  # Orange brûlé (dynamique et chaud)
    "#cb4335",  # Rouge brique (contrasté et intense)
    "#1d8348",  # Vert sombre (nature et stabilité)
]

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


def viz_pressure_stage(layout_stage) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.p_stat*1e-6 for stage in layout_stage], label=r"$p_{stat}$", color=color_list[0], linestyle='-.')
    plt.plot([stage.p_tot*1e-6 for stage in layout_stage], label=r"$p_{tot}$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Presure [MPa]")
    plt.legend()
    plt.xlim(0, len(layout_stage) - 1)
    plt.savefig(f"../figures/design/pressure_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_temperature_stage(layout_stage) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.T_stat for stage in layout_stage], label=r"$T_{stat}$", color=color_list[0], linestyle='-.')
    plt.plot([stage.T_tot for stage in layout_stage], label=r"$T_{tot}$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Temperature [K]")
    plt.xlim(0, len(layout_stage) - 1)
    plt.legend()
    plt.savefig(f"../figures/design/temperature_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_density_stage(layout_stage) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.rho for stage in layout_stage], color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Density [kg/m$^3$]")
    plt.xlim(0, len(layout_stage) - 1)
    plt.savefig(f"../figures/design/density_stage.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def viz_radius(layout_stage, R_mean) :
    plt.figure(figsize=(10, 6))
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
    plt.figure(figsize=(10, 6))
    plt.plot([stage.vm for stage in layout_stage_off], label=r"$v_m$", color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Mean velocity [m/s]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/vm_stage.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def viz_off_design_T_tot(layout_stage_off) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.T_tot for stage in layout_stage_off], label=r"$T_{tot}$", color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total temperature [K]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/T_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def viz_off_design_p_tot(layout_stage_off) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.p_tot for stage in layout_stage_off], label=r"$p_{tot}$", color=color_list[0])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total pressure [Pa]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def viz_data_comp_off_and_on_p_tot(layout_stage_off, layout_stage) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.p_tot for stage in layout_stage_off], label=r"$p_{tot} off$", color=color_list[0])
    plt.plot([stage.p_tot for stage in layout_stage], label=r"$p_{tot} on$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total pressure [Pa]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def viz_data_comp_off_and_on_T_tot(layout_stage_off, layout_stage) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.T_tot for stage in layout_stage_off], label=r"$T_{tot} off$", color=color_list[0])
    plt.plot([stage.T_tot for stage in layout_stage], label=r"$T_{tot} on$", color=color_list[1])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total temperature [K]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/T_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def viz_off_design_several_p_tot(layout_stage, layout_stage_off, layout_stage_down) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.p_tot for stage in layout_stage], label=r"$n_{design}$", color=color_list[0])
    plt.plot([stage.p_tot for stage in layout_stage_off], label=r"$ 1.05 \cdot n_{design}$", color=color_list[1])
    plt.plot([stage.p_tot for stage in layout_stage_down], label=r"$0.9 \cdot n_{design}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"Total pressure [Pa]")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.show()

def viz_off_design_severl_vm(layout_stage, layout_stage_off, layout_stage_down) :
    plt.figure(figsize=(10, 6))
    plt.plot([stage.vm for stage in layout_stage], label=r"$n_{design}$", color=color_list[0])
    plt.plot([stage.vm for stage in layout_stage_off], label=r"$ 1.05 \cdot n_{design}$", color=color_list[1])
    plt.plot([stage.vm for stage in layout_stage_down], label=r"$0.9 \cdot n_{design}$", color=color_list[2])
    plt.xlabel(r"Number of stage [-]")
    plt.ylabel(r"$V_m [m/s]$")
    plt.xlim(0, len(layout_stage_off) - 1)
    plt.legend()
    # plt.savefig(f"../figures/off_design/p_tot_stage.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()