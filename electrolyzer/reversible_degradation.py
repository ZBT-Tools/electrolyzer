import math
import numpy
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import R, physical_constants
from scipy.integrate import solve_ivp

delta_theta_OER_0 = 1.183  # V, Nernst Potential, calculated
delta_theta_0_Ir = 1.379  # V, Nernst Oxide formation bulk equilibrium potential, fitting
v_ox = 373  # 1/s, Forward oxide formation rate constant, fitting
v_red = 0.0066  # 377 # 1/s, Backward oxide formation rate constant, fitting
H_ox = 39720  # J/mol, Partial molar oxide formation activation enthalpy, fitting
lambd = 22863  # J/mol, Oxide dependent kinetic barrier constant, fitting
omega = 59828  # J/mol, Oxide-oxide interaction energy, fitting

i_0_Ir_III = 1.32e-7  # A/cm2 Apparent exchange current density for Ir(III), fitting
i_0_Ir_V = 6.97e-9  # A/cm2 Apparent exchange current density for Ir(V), fitting

a_ox = .5  # Charge transfer coefficient for oxidation assumed
a_red = .5  # Assumed
n_Ir = 2  # Electrons transferred during transition between Ir(III) and Ir(V)

T_anode = 273.15 + 80
# anode charge transfer coefficient
alpha_a = 2
F, _, _ = physical_constants["Faraday constant"]  # Faraday's constant [C/mol]


def reversible_performance(i):
    """
    Krenz, Tobias; Rex, Alexander; Helmers, Lennard; Trinke, Patrick; Bensmann, Boris;
    Hanke-Rauschenbach, Richard (2024):
    Reversible Degradation Phenomenon in PEMWE Cells: An Experimental and Modeling Study.
    In: J. Electrochem. Soc. 171 (12), S. 124501. DOI: 10.1149/1945-7111/ad96e4.

    Model describes the oxidation and reduction of iridium,
    leading to reversible performance losses

    Assumptions: Constant ECSA


    A_ECSA: Active Cell Area
    A_Ir_V: Area of higher oxidation state (A_ECSA = A_Ir_V + A_Ir_III)
    theta = A_Ir_V / A_ECSA

    i = i_Ir_V + i_Ir_III

    eta_Ir: the overpotential of the oxidation of iridium, positive values are the driving force for oxidation, negative values drive for reduction.

    d/dt theta(t, eta_Ir) = r_ox(theta(t,eta_Ir)) - r_red(theta(t,eta_Ir))

    vox is the rate constant of the oxidation,
    vred is the rate constant of the reduction
    Hox is the partial molar oxide formation activation enthalpy
    λ is the oxide dependent kinetic barrier constant,
    α_ox is the charge transfer coefficient of the iridium oxidation and
    α_red is the charge transfer coefficient of the iridium reduction
    Moreover we assume: αox + αred = 1.

    reaction rate of the oxidation
    r_ox = (1-theta) * v_ox e**(-(H_ox + lambda*theta)/(RT)) e**(alpha_ox F / R / T * eta_Ir)

    reaction rate of the reduction
    r_red = theta * v_red e**(-(H_ox + lambda*theta)/(RT)) e**(alpha_red F / R / T * eta_Ir)


    overpotential of the oxidation of iridium eta_ir

    The anode potential (Δtheta_a)
    equilibrium potential for the oxidation/reduction of iridium  (Δtheda_Ir)

    eta_Ir(theta(t), i) = delta_theta_a - delta_theda_Ir(theta(t))

    delta_theta_a = delta_theta_OER_0 + eta_act,OER(i,theta(t))
    delta_theta_Ir = delta_theta_0_Ir + omega * theta / (n_Ir * F)


    """

    return None  # _act_a


def calc_V_act_a(theta, i):
    """
    
    """
    V_act_a = ((R * T_anode) / (alpha_a * F)) * np.arcsinh(
        i / (2 * ((1 - theta) * i_0_Ir_III + theta * i_0_Ir_V)))
    return V_act_a

# def calc_V_act_a_alt(theta, i):
#     b = 46.5/1000
#     V_act_a = b * np.log(i / ((1 - theta) * i_0_Ir_III + theta * i_0_Ir_V))
#
#     return V_act_a

def calc_eta_Ir(V_act_a, theta):
    # The anode potential(Δtheta_a) [13]
    delta_theta_a = delta_theta_OER_0 + V_act_a
    # equilibrium potential for the oxidation/reduction of iridium
    delta_theta_Ir = delta_theta_0_Ir + (omega * theta) / (n_Ir * F)

    # overpotential of the oxidation of iridium eta_ir
    eta_Ir = delta_theta_a - delta_theta_Ir

    return eta_Ir

def calc_r_ox(theta, eta_Ir):
    # reaction rate of the oxidation
    r_ox = ((1 - theta) * v_ox * math.exp(-(H_ox + lambd * theta) / (R * T_anode)) *
            math.exp(a_ox * F / R / T_anode * eta_Ir))

    return r_ox

def calc_r_red(theta, eta_Ir):
    # reaction rate of the reduction
    r_red = (theta * v_red * math.exp(-(H_ox + lambd * theta) / (R * T_anode)) *
             math.exp(- a_red * F / R / T_anode * eta_Ir))
    return r_red


def calc_rate(t, theta, i):
    """
    Calculation of Ir Oxidation state change rate
    """
    if isinstance(theta,numpy.ndarray):
        if len(theta)>1:
            raise Exception
        else:
            theta = theta[0]
    V_act_a = calc_V_act_a(theta, i)
    # V_act_a =calc_V_act_a_alt(theta, i)
    eta_Ir = calc_eta_Ir(V_act_a, theta)
    r_ox = calc_r_ox(theta,eta_Ir)
    r_red = calc_r_red(theta,eta_Ir)

    return r_ox - r_red


# def ode_test(i):

if __name__ == "__main__":

    # # Different starting thetas for constant i
    # # ----------------------------------------------------------------------------------------------
    current_density = 2 # A/cm2
    fig, axs = plt.subplots(2)
    fig.suptitle(
        f'Different Ir-Oxidation states for constant current density {current_density} A/cm2')
    for th in np.arange(0, 1, 0.2):
        result = solve_ivp(calc_rate, (0, 60 * 60 * 24), [th], args=(current_density,))
        V_act_a = calc_V_act_a(theta=result.y[0, :], i=current_density)
        axs[0].plot(result.t, result.y[0, :], label=f"th={th}")
        axs[1].plot(result.t, V_act_a, label=f"th={th}")
    axs[0].set_title('Theta')
    axs[1].set_title('V_act_a [V]')
    #plt.legend()
    # plt.show()

    #
    # Manual integration
    # ----------------------------------------------------------------------------------------------
    t=[0]
    theta=[0]
    for ts in range(60 * 60 * 24):
        t.append(ts+1)
        rate = calc_rate(t=1, theta=theta[-1],i=current_density)
        theta.append(theta[-1]+rate)
    axs[0].plot(t, theta, label=f"manual, th=0")
    axs[0].legend()
    plt.show()


    # # Manual integration with time dependent i
    # # ----------------------------------------------------------------------------------------------
    # fig, axs = plt.subplots(4)
    # # current_density = [2] * 60 * 60 * 24 * 30
    # #                    + \
    # #                    [2] * 60 * 60 * 1 +  \
    # #                     [2] *60 * 60 * 24)
    # #
    # current_density = [2] * 60 * 60 * 24   + \
    #                    [0] * 60 * 5  +  \
    #                     [2] *60 * 60* 24
    #
    #
    # t=[0]
    # theta=[0.8]
    # rate = [0]
    # r_ox=[0]
    # r_red=[0]
    # V_act_a = [0]
    # for cd in current_density:
    #     t.append(t[-1]+1)
    #     r_ox.append(calc_r_ox(theta=theta[-1],
    #                           eta_Ir=calc_eta_Ir(calc_V_act_a(theta=theta[-1], i=cd), theta[-1])))
    #     r_red.append(calc_r_red(theta=theta[-1],
    #                             eta_Ir=calc_eta_Ir(calc_V_act_a(theta=theta[-1], i=cd), theta[-1])))
    #     local_rate = calc_rate(t=1, theta=theta[-1],i=cd)
    #     rate.append(local_rate)
    #     theta.append(theta[-1]+local_rate)
    #     V_act_a.append(calc_V_act_a(theta=theta[-1], i=cd))
    #
    # current_density.insert(0,2)
    # axs[0].plot(t, current_density, label=f"manual")
    # axs[1].plot(t, theta, label=f"theta")
    # axs[2].plot(t, r_ox, label=f"r_ox")
    # axs[2].plot(t, r_red, label=f"r_red")
    # axs[2].plot(t, rate, label=f"Rate")
    # axs[3].plot(t, V_act_a, label=f"V_act_a")
    # axs[0].set_title('Current density')
    # axs[1].set_title('Theta')
    # axs[2].set_title('Red and ox rates')
    # axs[3].set_title('V_act_a')
    # axs[2].legend()
    # plt.show()
