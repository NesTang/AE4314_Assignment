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
    
def range_solver(x, y, epsilon=0.01, init_guess = 1000):
    '''
    Take in x and y, uses bisection method to find the slope
    '''
    def calc_error():
        error = np.min(y - (slope_hi+slope_lo)/2 * x)
        return error
    def checker():
        if abs(error) > epsilon:
            return True
        else:
            return False
        
    slope_hi = init_guess
    slope_lo = 0
    error = calc_error()
    
    while checker():
        if error > 0:
            slope_lo = (slope_hi+slope_lo)/2
        else: 
            slope_hi = (slope_hi+slope_lo)/2
        error = calc_error()
    return (slope_hi+slope_lo)/2
    
def min_locator(x):
    min = np.min(x)
    indx = list(x).index(min)
    return min, indx
    
class BEM():
    '''
    Blade Element Method
    '''
    def __init__(self, V_ind, airspeed=0):

        self.omega = 326/60*2*np.pi # RPM of blade [rad/s]
        self.thrust = 46730 # total thrust in N
        self.R = 12.8/2 # Radius of the rotor
        self.solidity = 0.077 
        self.blade_num = 4 # blade count
        self.chord = 0.387 # chord of a single blade
        self.flat = 1.028 # equivalent flat plat area of Lynx
        
        self.V_ind = V_ind # induced velocity
        self.V = airspeed # forward velocity

        self.dlda = 0.768/5.5 # Slope of CL-alpha [-/deg]

        self.CT = self.calc_CT()
        #TODO
        # calculate each blade's thrust and power in hover
        # becasue its in hover, all blade is equal, multipling result by 4
        sample = 1000
        stations = np.linspace(0,1,sample)
        Thrust = []
        Torque = []
        Power = []
        for i in stations:
            dummy, Thrust_temp, Torque_temp = self.calc_blade_element(station=i, dr=self.R/sample, V_ind=self.V_ind[0])
            Thrust.append(Thrust_temp[0])
            Torque.append(Torque_temp[0])
            Power.append(Torque_temp[1])
        # sum up everything
        self.calc_rotor(dT=Thrust, dQ=Torque, dP=Power)

        # applying velocities
        self.P_ind = self.V_ind * self.T # induced power
        self.P_pro = self.P * (1+(self.V/(self.omega*self.R))**2) # profile power
        self.P_par = self.flat * 0.5 * 1.225 * self.V**3 # parasitic power
        self.P_tro = self.calc_tailrotor() # tail rotor power
        self.P_tot = self.P_ind + self.P_pro + self.P_par + self.P_tro

    def calc_CT(self):
        CT = self.thrust / (1.225 * (self.omega*self.R)**2 * (np.pi * self.R**2))
        CL_med = 6.6*CT/self.solidity
        return 
        
    def twist(self, station): # WORKING
        # if station < 0.1252:
        #     return 0
        # elif station < 0.2143: 
        #     return np.deg2rad(13.75/(0.2143-0.1252) * (station - 0.1252))
        # elif station <= 1:
        #     return np.deg2rad((7.46-13.75)/(1-0.2143) * (station - 0.2143) + 13.75)
        if station < 0.1252:
            return 0
        elif station < 0.2143: 
            return np.deg2rad(13.75)
        elif station <= 1:
            return np.deg2rad(7.46)
        
    def calc_blade_element(self, station, dr, V_ind, rho = 1.225):
        '''
        Calculation for one blade element
        station: goes from 0 to 1, indicate location of the station, lower max station value to account for tip loss
        V_flap: flapping velocity of the element
        az: Azimuth of element
        dr: width of the element
        rho: density
        '''
        Up = V_ind
        Ut = self.omega*self.R*station
        U = np.sqrt(Up**2+Ut**2)
        
        
        psi = np.arctan2(Up, Ut) # Angle of inflow
        alpha = self.twist(station) - psi # angle of attack
        
        cl = np.clip(self.dlda*np.rad2deg(alpha),a_min=0,a_max=10) # element lift coeff
        cd = 0.011 # element drag coeff # TODO: finish this

        dL = 0.5 * rho * U**2 * cl * self.chord * dr # element lift
        dD = 0.5 * rho * U**2 * cd * self.chord * dr # element drag
        
        dT = dL*np.cos(psi) - dD*np.sin(psi) # element thrust
        dFx = dD*np.cos(psi) - dL * np.sin(psi) # element drag
        dQ = -1*dFx * self.R*station # element torque
        dP = dQ * self.omega # element power

        return [dL, dD], [dT, dFx], [dQ, dP]

    def calc_tailrotor(self):
        l_tr = 7.76
        T = self.P_pro / self.omega / l_tr
        V_tr = np.sqrt((T/(np.pi*1.1049**2))/(2*1.225))
        k_tr = 1.3
        P_tr = 1.1 * k_tr * T * V_tr
        return P_tr
        
        
    
    def calc_rotor(self, dT, dP, dQ):
        '''
        Calculations that integrate the elements
        dT: Element thrust, all blades
        dQ: Element torque, all blades
        '''
        solidity = self.blade_num * self.chord / (np.pi*self.R) # Rotor solidity

        self.T = self.blade_num*np.sum(dT) # Small angle approx is not applied here. should also equal to weight
        self.Q = self.blade_num*np.sum(dQ) #
        self.P = self.blade_num*np.sum(dP) # Total torque on the rotor
        
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
    V_test = np.linspace(0,100,200)
    V_ind, V_hov = vel_induce_ACT(disc_load=363.2, V_flight=V_test)
    print('V_hov: ', V_hov)

    Lynx_BEM = BEM(V_ind = np.array(V_ind)*V_hov, airspeed = V_test)
    
    
    plt.plot(V_test, Lynx_BEM.P_tot/1000, label='Total', linewidth='2')
    plt.plot(V_test, Lynx_BEM.P_ind/1000, label='Induced', ls='--')
    plt.plot(V_test, Lynx_BEM.P_pro/1000, label='Profile', ls='--')
    plt.plot(V_test, Lynx_BEM.P_par/1000, label='Parasitic', ls='--')
    plt.plot(V_test, Lynx_BEM.P_tro/1000, label='Tail', ls='--')
    plt.grid(True)
    plt.title('Power vs Airspeed')
    plt.xlabel('Airspeed [m/s]')
    plt.ylabel('Power [kW]')
    plt.legend()
    plt.show()
    plt.clf()
    
    plt.plot(V_test, Lynx_BEM.P_tot/1000, label='Total',linewidth='2')
    best_range_slope = range_solver(V_test, Lynx_BEM.P_tot/1000)
    plt.plot(V_test, best_range_slope*V_test, ls ='--', color='gray')
    best_range = min_locator((Lynx_BEM.P_tot/1000)-(best_range_slope*V_test))
    best_endur = min_locator(Lynx_BEM.P_tot/1000)
    plt.scatter(V_test[best_range[1]], Lynx_BEM.P_tot[best_range[1]]/1000,
                label = f'Best Range: V={V_test[best_range[1]]: .2f}m/s, P={Lynx_BEM.P_tot[best_range[1]]/1000: .2f}kW')
    plt.scatter(V_test[best_endur[1]], Lynx_BEM.P_tot[best_endur[1]]/1000,
                label = f'Best Endurance: V={V_test[best_endur[1]]: .2f}m/s, P={Lynx_BEM.P_tot[best_endur[1]]/1000: .2f}kW')
    plt.grid()
    plt.title('Power vs Airspeed')
    plt.xlabel('Airspeed [m/s]')
    plt.ylabel('Power [kW]')
    plt.legend()
    plt.show()