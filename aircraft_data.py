import numpy as np

'''Westland Lynx Data'''
# From Jane's
# NASA Lynx XY170
# https://rotorcraft.arc.nasa.gov/Publications/files/Lau1993TM104000.pdf
class Lynx():
    def __init__(self):
        self.mtow = 5330 # kg
        self.EmptyWeight = 3277 #kg
        self.V_max = 90 # m/s
        self.range = 528e3 # meter
        self.endur = 19200 # second
        self.fuel = 1344.96 # kg
        self.rotor_d = 12.8 # meter
        self.power = 746e3*2 # watt
        self.blade_num = 4

        self.payload = 1361 # kg

        self.disk_area = (self.rotor_d/2)**2 * np.pi
        self.disk_load = self.mtow/self.disk_area # kg/m2
        