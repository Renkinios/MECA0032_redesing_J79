# This code appears to be the property of Koen Hillen, Professor of Turbomachinery at the University of Liège.
# Extracted from the lecture notes.



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

# ------------------------------------------------------------------------------
# specifications
# ------------------------------------------------------------------------------
## --- operating parameters

rHub = 0.18888888888888888
rTip = 0.4
omega = 7103 * 2 * np.pi / 60
reaction = 0.6
psiMid = 0.5
phiMid = 0.735

alpha_1 = np.arctan(((1 - reaction) - (psiMid / 2)) / phiMid)
alpha_2 = np.arctan(((1 - reaction) + (psiMid / 2)) / phiMid)
beta_1  = alpha_2 
beta_2  = alpha_1



# --- derived parameters
rMid = (rHub+rTip)/2

uHub = omega*rHub
uTip = omega*rTip
uMid = omega*rMid

psiTip = psiMid*(rMid*rMid)/(rTip*rTip)
phiTip = phiMid*(rMid/rTip)

psiHub = psiTip*(rTip*rTip)/(rHub*rHub)
phiHub = phiTip* rTip /rHub

vmGlb = phiTip * (uTip)
#-------------------------------------------------------------------------------
# whirl velocity is defined as vu = a + b r + c/r
#-------------------------------------------------------------------------------
def computeWhirlVelocity(a,b,c,r):
    return a + b*r + c/r


#-------------------------------------------------------------------------------
# Integrate the ISRE equation using Gauss-Legendre
# d/dr(vm^2) = -1/r^2 d/dr ((rvu)^2)
#-------------------------------------------------------------------------------
# Sa represente le delta vm en faite, donc ici 
def computeMeridionalVelocitySquared(a,b,c,r,vmHub):
    quadPnt = [-7.745966692414834e-01,0.000000000000000e+00,7.745966692414834e-01]
    quadWgt = [5.555555555555552e-01, 8.888888888888888e-01,5.555555555555552e-01]
    vmSq = r.copy()
    vmSq[0] = vmHub*vmHub
    for i in range(len(r)-1):
        vmSq[i+1] = vmSq[i]
        dr = r[i+1] - r[i]
        for j in range(len(quadPnt)):
            rQP = r[i] + dr*(1+quadPnt[j])/2
            dVmSq = - (4.*b*b*rQP + 6.*a*b + (2.*a*a+4*b*c)/rQP + 2.*a*c/rQP/rQP)
            vmSq[i+1] += dVmSq * quadWgt[j]*dr/2
    
    return vmSq

#-------------------------------------------------------------------------------
# compute the average meridional velocity using trapezoidal rule
#-------------------------------------------------------------------------------
def computeAverageMeridionalVelocity(vmSq,vmHub,r):
    massFlow = 0
    vmHub2 = vmHub*vmHub
    area = 0
    for i in range(len(r)-1):
        dr = r[i+1] - r[i]
        ra = (r[i]+r[i+1])/2
        vm1 = np.sqrt(max(0,vmSq[i] + vmHub2)) # permet que pas negatif mais devrai pas arriver
        vm2 = np.sqrt(max(0,vmSq[i+1] + vmHub2))
        area += ra*dr
        massFlow += (vm1+vm2)/2 *ra *dr
    rh = r[0]
    rt = r[len(r)-1]
    return massFlow/area

#-------------------------------------------------------------------------------
# compute the distributions for a given distribution a + br + c/r
# respecting
# - hub/shroud radii
# - global work and flow coefficient
# - degree of reaction specified at r=rc
# use a regula falsi method to determine hub meridional velocity
#-------------------------------------------------------------------------------

def computeDistributions(a,b,c,r,u,vmTgt,vm,vu,wu):
    vmSq = computeMeridionalVelocitySquared(a,b,c,r,0)
    vmHMin = 0
    vmMin = computeAverageMeridionalVelocity(vmSq,vmHMin,r)
    vmHMax = 2*vmTgt
    vmMax = computeAverageMeridionalVelocity(vmSq,vmHMax,r)
    vmHAve = (vmHMin+vmHMax)/2
    while ((vmHMax-vmHMin) > 0.001*vmTgt):
        vmHAve = (vmHMax + vmHMin)/2
        vmAve = computeAverageMeridionalVelocity(vmSq,vmHAve,r)
        if (vmAve >= vmTgt):
            vmMax = vmAve
            vmHMax = vmHAve
        if (vmAve <= vmTgt):
            vmMin = vmAve
            vmHMin = vmHAve
        vu[:] = computeWhirlVelocity(a,b,c,r)

        vm[:] = np.sqrt(vmSq + vmHAve*vmHAve)
        wu[:] = vu - u

# ------------------------------------------------------------------------------
# plot all relevant distributions
# - flow angles
# - relative flow velocity upstream of rotor and stator
# - axial velocity upstream of the rotor
# ------------------------------------------------------------------------------

def plot_angle(r, u, vm1, vu1, wu1, vm2, vu2, wu2, name):
    # --- Plot angles and reaction in the same figure ---
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("$r$ [m]")
    ax1.set_ylabel("$\\alpha, \\beta~[^\\circ$]")
    ax1.set_xlim(rHub, rTip)
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("$R$ [-]", color=color_list[2])  # Synchronize ylabel color with R
    ax2.tick_params(axis='y', labelcolor=color_list[2])  # Sync tick labels with color_list[2]
    ax2.spines['right'].set_color(color_list[2])  # Set the color of the right spine
    ax2.set_xlim(rHub, rTip)
    ax2.axhline(y=0.5, linestyle='--', color='k')
    ax2.axvline(x=rMid, linestyle='--', color='k')

    # Plot alpha and beta angles
    lns1 = ax1.plot(r, np.arctan(vu1 / vm1) * 180 / np.pi, linestyle='-', color=color_list[0], label="$~\\alpha_1$")
    lns2 = ax1.plot(r, np.arctan(vu2 / vm2) * 180 / np.pi, '-.', color=color_list[0], label="$~\\alpha_2$")
    lns3 = ax1.plot(r, -np.arctan(wu1 / vm1) * 180 / np.pi, '-', color=color_list[1], label="$-\\beta_1$")
    lns4 = ax1.plot(r, -np.arctan(wu2 / vm2) * 180 / np.pi, '-.', color=color_list[1], label="$-\\beta_2$")

    # Plot reaction
    reaction = 1 - (wu1 + wu2) / (2 * u)
    lns5 = ax2.plot(r, reaction, color=color_list[2], label="$~R$")

    # Combine legends
    lns = lns1 + lns2 + lns3 + lns4 + lns5
    labs = [l.get_label() for l in lns]

    # Add legend above the figure
    fig.legend(lns, labs, loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=3, frameon=False)  # Legend closer to the graph

    # Adjust layout to prevent overlap
    fig.tight_layout(rect=[0, 0, 1, 0.85])  # Reduce reserved space above

    # Save figure
    plt.savefig("../figures/ISR/isre_{}_angles.pdf".format(name), format="pdf", transparent=True, bbox_inches='tight')
    plt.close()


def plot_velocities(r, u, vm1, vu1, wu1, vm2, vu2, wu2, name):
    # --- Plot velocities ---
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("$r$ [m]")
    ax1.set_ylabel("$v, w$ [m/s]")
    ax1.set_xlim(rHub, rTip)

    ax2 = ax1.twinx()
    ax2.set_ylabel("$v_m$ [m/s]", color=color_list[2])  # Synchronize ylabel color with vm1/vm2
    ax2.tick_params(axis='y', labelcolor=color_list[2])  # Sync tick labels with color_list[2]
    ax2.spines['right'].set_color(color_list[2])  # Set the color of the right spine
    ax2.set_xlim(rHub, rTip)

    # Plot lines
    lns1 = ax1.plot(r, np.sqrt(vu1**2 + vm1**2), color=color_list[0], label="$v_1$")
    lns2 = ax1.plot(r, np.sqrt(vu2**2 + vm2**2), '-.', color=color_list[0], label="$v_2$")
    lns3 = ax1.plot(r, np.sqrt(wu1**2 + vm1**2), color=color_list[1], label="$w_1$")
    lns4 = ax1.plot(r, np.sqrt(wu2**2 + vm2**2), '-.', color=color_list[1], label="$w_2$")
    lns5 = ax2.plot(r, vm1, color=color_list[2], label="$v_{m1}$")
    lns6 = ax2.plot(r, vm2, '-.', color=color_list[2], label="$v_{m2}$")

    # Vertical line
    ax2.axvline(x=rMid, linestyle='--', color='k')

    # Combine all lines for the legend
    lns = lns1 + lns2 + lns3 + lns4 + lns5 + lns6

    # Add legend explicitly at a controlled location
    fig.legend(
        handles=lns,
        labels=[l.get_label() for l in lns],
        loc='upper center',  # Position: 'upper center'
        bbox_to_anchor=(0.5, 0.95),  # Lower the legend slightly
        ncol=3,  # Number of columns in the legend
        frameon=False  # No frame around the legend
    )

    # Adjust layout and display
    fig.tight_layout(rect=[0, 0, 1, 0.85])  # Adjust layout to avoid overlap with the legend
    plt.savefig("../figures/ISR/isre_{}_speed.pdf".format(name), format="pdf", transparent=True, bbox_inches='tight')




    
def plot_coefficients(r, u, vm1, vu1, wu1, vm2, vu2, wu2, name):
    # --- Plot stage coefficients ---
    coefficients = plt.figure()
    fig, ax1 = plt.subplots()
    
    # Axes configurations
    ax1.set_xlabel("$r$ [m]")
    ax1.set_ylabel("$\\phi$ [-]")
    ax1.set_xlim(rHub, rTip)
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("$\\psi$ [-]", color=color_list[2])  # Sync ylabel color with psi
    ax2.tick_params(axis='y', labelcolor=color_list[2])  # Sync tick labels with psi
    ax2.spines['right'].set_color(color_list[2])  # Set right spine color to match psi
    ax2.set_xlim(rHub, rTip)
    ax2.axvline(x=rMid, linestyle='--', color='k')
    
    # Plot curves
    lns1 = ax1.plot(r, vm1 / u, color=color_list[0], label="$\\phi_1$")
    lns2 = ax1.plot(r, vm2 / u, '-.', color=color_list[1], label="$\\phi_2$")
    lns4 = ax2.plot(r, (vu2 - vu1) / u, color=color_list[2], label="$\\psi$")
    
    # Combine legends
    lns = lns1 + lns2 + lns4
    labs = [l.get_label() for l in lns]
    
    # Add legend
    ax1.legend(lns, labs, loc='best', frameon=False)
    
    # Save the figure
    plt.savefig("../figures/ISR/isre_{}_coefficients.pdf".format(name), format="pdf", transparent=True, bbox_inches='tight')
    plt.close()


# --------------------------------------------------------------------------
# Compute radii
# --------------------------------------------------------------------------

nbPoints = 50
drT = (rTip - rHub) / nbPoints
r = np.arange(rHub, rTip + drT, drT)
u = omega * r
vm1 = r.copy()
vu1 = r.copy()
wu1 = r.copy()
vm2 = r.copy()
vu2 = r.copy()
wu2 = r.copy()



# --------------------------------------------------------------------------
# Free vortex, reaction at mid radius
# --------------------------------------------------------------------------
a1, b1, a2, b2 = 0, 0, 0, 0
c1 = (1 - (2 * reaction + psiMid) / 2) * uMid**2 / omega
c2 = (1 - (2 * reaction - psiMid) / 2) * uMid**2 / omega
print("Free vortex, reaction at midspan: vu1 = {}/r, vu2(r) = {}/r".format(c1, c2))

computeDistributions(a1, b1, c1, r, u, vmGlb, vm1, vu1, wu1)
computeDistributions(a2, b2, c2, r, u, vmGlb, vm2, vu2, wu2)
plot_angle(r, u, vm1, vu1, wu1, vm2, vu2, wu2, "freeVortex_mid")
plot_coefficients(r, u, vm1, vu1, wu1, vm2, vu2, wu2, "freeVortex_mid")
plot_velocities(r, u, vm1, vu1, wu1, vm2, vu2, wu2, "freeVortex_mid")
