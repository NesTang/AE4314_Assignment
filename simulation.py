import numpy as np
import matplotlib.pyplot as plt


'''3DoF Sim'''

class Sim3DoF():
    def __init__(self, run_sim=False, trim_solve=True, vel_range = [0]):
        '''Parameters'''
        self.My = 0 # TODO
        self.Iy = 2000 # TODO
        self.CdS = 1.028 # equivalent flat plat area of Lynx
        self.Ct = 0.0062100841852850384 # Thrust Coefficent, from performance.py
        self.omega = 326/60*2*np.pi # rotor RPM
        self.R = 12.8/2 # rotor radius
        self.weight = 46730
        self.h = 2.1854 # CG to control plane
        
        '''Position'''
        self.X = 0
        self.Z = 0
        
        '''Controls'''
        self.collective = 0 # blade pitch angle (collective)
        self.control = 0 # control plane angle
        
        '''States'''
        self.u = 0.01 # Horizontal vel
        self.w = 0 # Vertical vel
        self.q = 0 # Pitch rate
        self.theta = 0 # Pitch
        
        '''Extra'''
        self.inflow_i = None
        self.solidity = 0.077 
        self.dlda = 0.768/np.deg2rad(5.5)
        
        '''Logs'''
        self.log_u = [[self.u, 0]]
        self.log_w = [[self.w, 0]]
        self.log_theta = [[self.theta, self.q]]
        self.log_time = [0]
        
        '''Solving for Trim'''
        if trim_solve:
            self.trims = []
            for vel in vel_range:
                self.u = vel
                self.trims.append(self.TrimSolver())
            self.trims = np.array(self.trims).T
                
        '''Start Sim'''
        if run_sim:
            self.MainLoop()
            

        
    def MainLoop(self, dt = 0.01, runtime=60):
        current_sim_time = 0
        while current_sim_time < runtime:
            '''Main Simulation Loop'''
            self.Physics()
            self.InflowSolver()
            du, dw, dq, dtheta = self.EoM()
            '''Update physics'''
            self.u += du*dt
            self.w += dw*dt
            self.q += dq*dt
            self.theta += dtheta*dt
            current_sim_time += dt
            '''Append Logs'''
            self.log_u.append([self.u, du])
            self.log_w.append([self.w, dw])
            self.log_theta.append([self.theta, self.q])
            self.log_time.append(current_sim_time)
        print('Sim Done')
        '''Clean up'''
        self.log_u = np.array(self.log_u).T
        self.log_w = np.array(self.log_w).T
        self.log_theta = np.array(self.log_theta).T 
        self.log_time = np.array(self.log_time)       
        
    def Physics(self):
        self.V = np.sqrt(self.u**2+self.w**2)    
        self.contorl_alpha = self.control - np.arctan2(self.w, self.u)
        self.V_tip = self.omega * self.R
        self.mu = self.V/self.V_tip * np.cos(self.contorl_alpha)
        self.inflow_con = self.V*np.sin(self.contorl_alpha) / self.V_tip
        
    def EoM(self):
        '''Equations of Motion'''
        T = self.Ct * 1.225 * self.V_tip**2 * np.pi * self.R**2
        D = self.CdS * 1/2 * 1.225 * self.V**2
        Fx = -self.weight*self.w/self.V - D*self.u/self.V + T*np.sin(self.control-self.a1)
        Fz =  self.weight*self.u/self.V - D*self.w/self.V - T*np.cos(self.control-self.a1)
        My = -T*self.h*np.sin(self.control-self.a1)
        
        '''Final equations'''
        du = -self.q*self.w + Fx/(self.weight/9.81) - 9.81*np.sin(self.theta)
        dw = self.q*self.u + Fz/(self.weight/9.81)
        dq = My/self.Iy
        dtheta = self.q
        return du, dw, dq, dtheta
        
    def InflowSolver(self):
        def calc_a1():
            # lock number (gamma) is assumed to be 10
            top = 8/3 * self.mu * self.collective - 2 * self.mu*(self.inflow_con+self.inflow_i)-16/10 * self.q/self.omega
            bottom = 1 - 0.5*self.mu**2
            return top/bottom
        def Ct_BEM():
            solidity = 0.077 
            dlda = 0.768/np.deg2rad(5.5)
            
            part1 = 0.25 * dlda *solidity
            part2 = 2/3 * self.collective * (1+3/2*self.mu**2) - (self.inflow_con+self.inflow_i)
            return part1*part2
        def Ct_glau():
            a1 = calc_a1()
            part1 = self.V/self.V_tip*np.cos(self.contorl_alpha-a1)
            part2 = self.V/self.V_tip*np.sin(self.contorl_alpha-a1)+self.inflow_i
            return 2*self.inflow_i*np.sqrt(part1**2+part2**2)
        
        def diff_calc(epsilon=0.001):
            diff = Ct_BEM() - Ct_glau()
            if abs(diff) < epsilon:
                return False, diff
            else: 
                return True, diff
        

        # setting up bisection
        guess_hi = 30
        guess_lo = -25
        self.inflow_i = (guess_hi + guess_lo)/2
        running, diff = diff_calc()
        while running:
            if diff < 0:
                guess_hi = self.inflow_i * 1
            else: 
                guess_lo = self.inflow_i * 1
            self.inflow_i = (guess_hi + guess_lo)/2
            running, diff = diff_calc()
        
        self.Ct = Ct_BEM()
        self.a1 = calc_a1()
    
    def TrimSolver(self):
        '''Trim solver for a single velocity'''
        self.Physics()
        self.Ct = self.weight / (1.225 * (self.omega*self.R)**2 * (np.pi * self.R**2))
        q = 1/2 * 1.225 * self.V
        drag = self.CdS * q
        # bisection again
        def Ct_glau():
            part1 = self.V/self.V_tip*np.cos(drag/self.weight)
            part2 = self.V/self.V_tip*np.sin(drag/self.weight)+self.inflow_i
            return 2*self.inflow_i*np.sqrt(part1**2+part2**2)
        
        def diff_calc(epsilon=0.0001):
            diff = self.Ct - Ct_glau()
            if abs(diff) < epsilon:
                return False, diff
            else: 
                return True, diff
        guess_hi = 5
        guess_lo = -1
        self.inflow_i = (guess_hi + guess_lo)/2
        running, diff = diff_calc()
        while running:
            if diff < 0:
                guess_hi = self.inflow_i * 1
            else: 
                guess_lo = self.inflow_i * 1
            self.inflow_i = (guess_hi + guess_lo)/2
            running, diff = diff_calc()
        
        '''Setting up matrix to solve'''
        matA = np.array([[1+(3/2)*self.mu**2, -8/3*self.mu],
                         [-self.mu, 2/3+self.mu**2]])
        matB = np.array([[-2*self.mu**2 * drag/self.weight - 2*self.mu*self.inflow_i],
                         [4/self.solidity*self.Ct/self.dlda+self.mu*drag/self.weight+self.inflow_i]])
        trim = np.linalg.solve(matA,matB) # [a1, collective]
        self.a1 = trim[0][0]
        
        return trim[0][0], trim[1][0]
        
        
test = True
if test:
    vel = np.linspace(0.1,100,100)
    sim = Sim3DoF(vel_range=vel)
    plt.plot(vel, np.rad2deg(sim.trims[0]), label='cyclic')
    plt.plot(vel, np.rad2deg(sim.trims[1]), label='collective')
    plt.grid(True)
    plt.title('Westland Lynx Trim')
    plt.xlabel('Airspeed [m/s]')
    plt.ylabel('Trim Pitch [deg]')
    plt.legend()
    plt.show()
    plt.clf()
    
else:
    print('simulation.py Testing mode is OFF')