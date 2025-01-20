import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Définition globale des paramètres de police et de taille pour tous les graphiques
plt.rc('font', family='serif')  # Police avec empattements, comme Times
plt.rc('text', usetex=True)  # Utiliser LaTeX pour le texte dans les figures
plt.rcParams.update({
    'font.size': 14,       # Taille de police générale
    'legend.fontsize': 14, # Taille de police pour les légendes
    'axes.labelsize': 16,  # Taille de police pour les étiquettes des axes
})

color_list = [
    "#007070", "#f07f3c", "#5b57a2", "#7db928", "#e62d31",
    "#005ca9", "#00843b", "#f8aa00", "#5b257d", "#8c8b82"
]
ligne_or_text = "#DCE6EA"

def sensitivity_15(angle, plot = False):
    Deviation_sensitivity = [0.9285106468770559, 0.9421276651281425, 0.9565957537159332, 
                             0.9659574495024572, 0.9719149029007754, 0.9812766035571258]
    Angle_inlet = [70, 60, 50, 40, 30, 0]
    
    coefficients = np.polyfit(Angle_inlet, Deviation_sensitivity, 2)
    
    polynomial = np.poly1d(coefficients)
    
    new_angles = np.linspace(0, 70, 100)
    
    interpolated_values = polynomial(new_angles)

    if plot == True:
        plt.plot(Angle_inlet, Deviation_sensitivity, 'o', label="Données originales")
        plt.plot(new_angles, interpolated_values, '-', label="Interpolation polynomiale de degré 2")
        plt.xlabel("Angle d'entrée")
        plt.ylabel("Sensibilité à la déviation")
        plt.legend()
        plt.show()
        plt.close()
    
    return polynomial(angle)


def get_sensitivity_solidity_1_5(inlet_angle):
    """
    Interpolates the deviation sensitivity for a given inlet angle using quadratic interpolation.
    
    Parameters:
        inlet_angle (float or array-like): The inlet angle(s) in rad for which deviation sensitivity is to be calculated.
    
    Returns:
        float or ndarray: The interpolated deviation sensitivity value(s).
    """
    # Data
    angle_inlet_data = np.array([70,60,50,40,30,0]) * np.pi / 180
    deviation_sensitivity = np.array([0.9453781512605042, 0.9565826330532212, 
                                       0.9677871148459384, 0.9733893557422969, 
                                       0.9859943977591036, 0.9971988795518206])
    coefficients = np.polyfit(angle_inlet_data, deviation_sensitivity, 2)
    
    polynomial = np.poly1d(coefficients)
    
    
    return polynomial(inlet_angle)


def plot_sensitivity_solidity_1_5():
    """
    Plots the quadratic interpolation of deviation sensitivity for a range of inlet angles.
    """
    # Data
    angle_inlet_data =  np.array([70,60,50,40,30,0]) * np.pi / 180
    deviation_sensitivity = np.array([0.9453781512605042, 0.9565826330532212, 
                                       0.9677871148459384, 0.9733893557422969, 
                                       0.9859943977591036, 0.9971988795518206])
    
    coefficients = np.polyfit(angle_inlet_data, deviation_sensitivity, 2)
    
    polynomial = np.poly1d(coefficients)

    # Generate fine angle grid for smooth curve
    fine_angles = np.linspace(0, 70, 500)
    interpolated_values = polynomial(fine_angles)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(fine_angles, interpolated_values, label='Quadratic Interpolation', color='blue')
    plt.scatter(angle_inlet_data, deviation_sensitivity, color='red', label='Original Data', zorder=5)
    plt.title('Quadratic Interpolation of Deviation Sensitivity')
    plt.xlabel('Inlet Angle (degrees)')
    plt.ylabel('Deviation Sensitivity')
    plt.legend()
    plt.grid()
    plt.show()
    plt.close()


def calculate_polynomial_fit(x, y, degree=10):
    """
    Calculate the polynomial fit for the given data points.

    Parameters:
    - x: list of float, x-coordinates of the data points.
    - y: list of float, y-coordinates of the data points.
    - degree: int, the degree of the polynomial fit.

    Returns:
    - polynomial: np.poly1d, the fitted polynomial function.
    """
    coefficients = np.polyfit(x, y, degree)
    polynomial = np.poly1d(coefficients)
    return polynomial

def get_momentum_thickness_off(Deq_value, plot=False):
    """
    Calculate the momentum thickness for a given Deq value and optionally plot the curve.

    Parameters:
    - Deq_value: float, the Deq value to calculate the momentum thickness for.
    - plot: bool, if True, plot the curve and data points.

    Returns:
    - float, the momentum thickness for the given Deq value.
    """
    # Data points for fitting
    x = [
        1.0037106539853968, 1.0526902063044665, 1.1016698095805253, 1.1539889301036566,
        1.2029684824227262, 1.2541744475158296, 1.3031539998348993, 1.3521336540679474,
        1.4022264127740338, 1.4512059650931035, 1.502411930186207, 1.553617793365332,
        1.60259744759838, 1.654916517164522, 1.7027829650105533, 1.7517625173296227,
        1.8029684824227261, 1.8519480347417958, 1.9020407934478825, 1.951020345766952,
        2.0022263108600553, 2.053432275953159, 2.1012986218852117, 2.153617793365332,
        2.2014841392973845, 2.254916415250544, 2.3038960694835917, 2.3595547582107845,
        2.4074211041428373, 2.4519479328278173
    ]
    y = [
        0.0038764089731951126, 0.003932590453314538, 0.004101129750261678, 0.004269669047208803,
        0.004494384680864197, 0.004775286938050179, 0.0051685418686527075, 0.005617983422785832,
        0.0062359591304549065, 0.006797758501415729, 0.007528092025912503, 0.008314612173939882,
        0.009044945698436665, 0.010000005143411168, 0.010898877964855074, 0.011910113746537849,
        0.012921349528220613, 0.014157306086969947, 0.01539326264571927, 0.01662921406105744,
        0.01797752843663446, 0.019831460703052868, 0.02162921663276302, 0.02415730608696994,
        0.026797753358004568, 0.030505617890841386, 0.035, 0.041516853672524154,
        0.04876404344125067, 0.05960673799720713
    ]

    # Fit a polynomial to the data
    polynomial = calculate_polynomial_fit(x, y)

    # Generate values for plotting
    Deq = np.linspace(min(x), max(x), 100)
    Momentum_thickness = polynomial(Deq)

    # Plot the data and fit if requested
    if plot:
        plt.figure(figsize=(8, 5))
        # plt.plot(x, y, 'o', label="Data points")
        plt.plot(Deq, Momentum_thickness, '-', color = color_list[0])
        plt.xlabel(r"Equivalent diffusion ratio $D_{eq}$")
        plt.ylabel(r"Momentum thickness $\delta_2/c$")
        # plt.legend()
        plt.xlim(np.min(Deq), np.max(Deq))
        plt.savefig("../figures/momentum_thickness_app.pdf", dpi=300, bbox_inches='tight')
        plt.close()

    # Return the value of the polynomial for the given Deq_value
    return polynomial(Deq_value)

get_momentum_thickness_off(1.5, plot=True)