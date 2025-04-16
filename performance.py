import numpy as np
import matplotlib.pyplot as plt

def vel_induce_ACT(disc_load, V_flight = 0, density = 1.225): # 
    '''
    Actuatro disc theory induced velocity, only suitable for alpha=0
    disc_load: [N/m2]
    V_flight: Array of flight velocity
    density: air density
    alpha: AoA of disc [Rad]
    '''
    V_ind_hov = np.sqrt(disc_load/(2*density)) # Hover induced velocity
    # Nondimensionalizing input velocity
    V_flight_nd = V_flight / V_ind_hov
    V_ind_fw_nd = []
    for V_nd in V_flight_nd:
        V_ind_fw_nd.append(np.sqrts(-0.5*V_nd**2+np.sqrt(0.25*V_nd**4+1)))
    return V_ind_fw_nd, V_ind_hov

def power_hov():
    ''' 
    Ideal Hover power, ACT and BEM
    '''



test = True
if test:
    V_test = np.linspace(0,50,100)
    V_ind, V_hov = vel_induce_ACT(disc_load=100, V_flight=V_test)
    plt.plot(V_test/V_hov, V_ind)
    plt.show()