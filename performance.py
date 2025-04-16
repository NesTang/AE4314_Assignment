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
        V_ind_fw_nd.append(np.sqrt(-0.5*V_nd**2+np.sqrt(0.25*V_nd**4+1)))
    return V_ind_fw_nd, V_ind_hov

def power_ACT(thrust, V_ind, FM = 0.7, ):
    ''' 
    Power Calculation With ACT
    thrust: hover thrust
    V_ind: ACT induced hover
    FM: Figure of Merit
    '''
    # Hover
    P_id = V_ind*thrust
    P_ACT_hov = P_id*FM
    # Forward - No vertical speed
    # Not asked
    return P_id, P_ACT_hov

def linebreak(text):
    print(f'-----------------{text}-----------------')
    
    
    
class BEM():
    '''
    Blade Element Method
    '''
    def __init__(self, omega, R, blade_num = 4):

        self.omega = omega # RPM of blade
        self.R = R # Radius of the rotor
        self.chord = 0 # chord of a single blade
        self.pitch = 0 # pitch of blade
        self.V_ind = 0 # induced velocity
        self.V = 0 # forward velocity
        self.blade_num = blade_num # blade count
        self.dlda = 1 # Slope of CL-alpha
        
    def twist(station):    
        if station < 0.1252:
            return 0
        elif station < 0.2143: 
            return 13.75/(0.2143-0.1252) * (station - 0.1252)
        elif station <= 1:
            return (13.75-7.46)/(1-0.2143) * (station - 0.2143) + 13.75
        
        
    def calc_blade_element(self, station, V_flap, az, dr, rho = 1.225):
        '''
        Calculation for one blade element
        station: goes from 0 to 1, indicate location of the station, lower max station value to account for tip loss
        V_flap: flapping velocity of the element
        az: Azimuth of element
        dr: width of the element
        rho: density
        '''
        Up = self.V_ind + V_flap  
        Ut = self.V*np.sin(az) + self.omega*self.R*station
        U = np.sqrt(Up**2+Ut**2)
        AR = self.V/(self.omega*self.R*station) # Advance ratio
        
        
        psi = np.arctan2(Up, Ut) # Angle of inflow
        alpha = self.pitch - psi # angle of attack
        
        cl = self.dlda*alpha # element lift coeff
        cd = cl**2/(np.pi) # element drag coeff # TODO: finish this

        dL = 0.5 * rho * U**2 * cl * self.chord * dr # element lift
        dD = 0.5 * rho * U**2 * cd * self.chord * dr # element drag
        
        dT = dL*np.cos(psi) - dD*np.sin(psi) # element thrust
        dFx = dD*np.cos(psi) - dL * np.sin(psi) # element drag
        dQ = dFx * self.R*station # element torque
        
        return [dL, dD], [dT, dFx], [dQ]

    def calc_tailrotor(self):
        l_tr = 10 #TODO Find actaul tail rotor distance
        T = self.P_hov / self.omega / l_tr
    
    def calc_rotor(self, dT, dQ):
        '''
        Calculations that integrate the elements
        dT: Element thrust, all blades
        dQ: Element torque, all blades
        '''
        solidity = self.blade_num * self.chord / (np.pi*self.R) # Rotor solidity

        self.T = np.sum(dT) # Small angle approx is not applied here. should also equal to weight
        self.Q = np.sum(dQ) # Total torque on the rotor

    def calc_hover(self):
        '''
        calculate total power required by summing up induced power and the profile power
        '''

        P_ind = power_ACT(self.T, self.V_ind, FM=1)[0] # ideal induced power
        P_pro = self.Q
        P_total = P_ind + P_pro
        self.P_hov = P_total
        return	P_total
      

test = True
if test:
    V_test = np.linspace(0,50,100)
    V_ind, V_hov = vel_induce_ACT(disc_load=100, V_flight=V_test)
    blade = BEM(omega=400, R=12.8/2)
    stations = np.linspace(0,1,50)
    for i in range(len(stations)):
        print(stations[i])
        plt.plot(stations, blade.twist(stations[i]))
    # plt.plot(V_test/V_hov, V_ind)5
    plt.show()