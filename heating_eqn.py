##-------------------------------------------------------------------------
##created by : Tristan_Parece
##created date: 04/01/2022
##version = '1.0'
##license = MIT
##-------------------------------------------------------------------------
## Heating Equation Calculator 
"""Assumptions 
-cross sectional area in constant
"""
from math import ceil
import os
import imageio
from gettext import npgettext
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pprint import pprint
import sqlite3
import unittest



class Heat_Eqn:
    def __init__(self, Material,Number_Of_Steps, Shape, Shape_Lengths, Length, Time, Left_Boundary, Right_Boundary, Gif = False, Graph = False, Hat_Function = False, Temperature_0 = None, Location_of_Temp = None ,Length_Of_Hat_Function = None):
        """Creates a Rod of material and initial heat conditions to then apply to the system to calculate relevant 
        qualities of thermal energy evolution.

        Args:
            Material (Text): Material composition of the rod 
            Number_Of_Steps (Integer): The number of iterations for the function to loop over
            Shape (Text): Cross sectional shape of the rod (Rectangle (to use square side lengths the same), Triangle (Equilateral),
                                                            Hexagon, Ellipse (for circle set internal legths the same))
            Shape_Lengths (Tuple of floats): Include the important shape lengths in meters
                                            if Rectangle (a,b) where a and b are the side lengths
                                            if Triangle (a,) where a is the side length
                                            if Hexagon (a,) where a is the side length
                                            if Ellipse (a,b) where a, b are the radiuses of the minor and majour axises
            Length (Float): The Length of the Rod
            Time (Float): Time duration of the calculation in seconds
            Left_Boundary (Float): Temperature in Celsius on the left side boundary 
            Right_Boundary (Float): Temperature in Celsius on the right side boundary

            Hat_Function (Bool): True if hat function applied to initial conditions 
                                False if no hat function exists in initial condition
            Temperature_0 (Float): Temperature in degrees Celsius at start
            Location_of_Temp (Float): Distance from the left end of the rod where the initial temperature conditions occur
            Length_Of_Hat_Function (Float): Distance the plateau of the hat function occurs, if an odd number is used it will have 1 added to it 
        """
        if Hat_Function == True and (Temperature_0 == None or Location_of_Temp == None or Length_Of_Hat_Function == None) :#Catches if the heat source is missing its location
            print('Missing a parameter in relation to the hat function')
        if len(list(Shape_Lengths)) == 1:
            a = Shape_Lengths
            Shape_Lengths = (a,0)
        self.Material = Material
        self.Number_Of_Steps = Number_Of_Steps
        self.Shape = Shape
        self.Shape_Lengths = Shape_Lengths
        self.Length = Length
        self.Time = Time
        self.Left_Boundary = Left_Boundary
        self.Right_Boundary = Right_Boundary
        self.Gif = Gif
        self.Graph = Graph
        self.Hat_Function = Hat_Function
        self.Temperature_0 = Temperature_0
        self.Location_of_Temp = Location_of_Temp
        self.Length_Of_Hat_Function = Length_Of_Hat_Function
    def area(self):
        """ Takes the Shape Profile and the Shape Lengths to determine the Cross sectional area
        Note: This function is not necessary to the heating equation, but could be a useful additional function in the future

        Returns:
            Cross Sectional area (Float): the cross sectional area of the given shape profile
        """
        shape = self.Shape
        a,b = self.Shape_Lengths
        if shape.lower() == 'ellipse':
            return np.pi*a*b
        elif shape.lower() == 'rectangle':
            return a*b
        elif shape.lower() == 'triangle':
            return ((3**(1/2))/4)*(a**2)
        elif shape.lower() == 'hexagon':
            return ((3*(3**(1/2)))/2)*(a**2)
        else: # Catch for when the shape is not supported
            print("Looks like you choose a shape that is not currently supported, review docstring of Heat_Eqn.")
            return None
    def mat_constants(self):
        """
        Takes the inputed material and searched the created Thermal database 
        Note: This Thermal database is incomplete and should be updated in the future
        Returns:
            Thermal Conductivity (Float): Units (W/m K)
            Thermal Diifusivity (Float): Units (mm^2/s)
            Specific Heat Capacity (Float): Units (J/kg K)
            Thermal Effusivity (Float): Units (W s^0.5 / m^2 K)
            Material Density (Float): Units (kg/m^3)
        """
        con = sqlite3.connect("Thermal.db")
        cur =con.cursor()
        mat = self.Material.lower()
        cur.execute(f"SELECT * FROM thermal_prop WHERE lower(Material) = '{mat}'")
        lis = cur.fetchall()
        if len(lis) == 0: #Catch for when material is incorrectly inputed or when material does not exist in database
            print("Looks like the material you selected either is not in the database or it is misspelled.")
            print("please use one of the following materials")
            cur.execute("SELECT Material from thermal_prop")
            pprint(cur.fetchall())
            con.close()
            return None
        else:
            M,K,alpha,c,e,rho = lis[0]
            con.close()
            return K,alpha,c,e,rho
    def heateqn(self):
        """Takes the inputs of the class, and calculates the heat transfer equation over the elloted time, 
        it can then either send the information to create a graph, gif or the 2 dimensional array

        Returns:
            Either nothing or the 2 dimensional temperature array U[:,:]
        """
        k,alpha,c,e,rho = self.mat_constants()
        if k == None:#confirms the Material consants have succsessfully been collected from the database
            print("An Error Occured please, double check your inputs and try again.")
            return None
        m = self.Number_Of_Steps
        alphaprime = k / (c*rho)
        xf = self.Length
        dx = xf/m
        tf = self.Time
        dt = (dx**2)/(2*alphaprime) #Von Neumann Stability Analysis
        t = np.arange(0,tf,dt)
        x = np.empty(m)
        U = np.ones((len(t),len(x)))*0
        if self.Hat_Function == True:
            deltafunction = int(ceil(self.Length_Of_Hat_Function/2))
            T_0 = self.Location_of_Temp
            startingcondition = int(((T_0/xf)*m)-1)
            U[0,startingcondition-deltafunction:startingcondition+deltafunction] = self.Temperature_0
        bl = self.Left_Boundary
        br = self.Right_Boundary
        for j in range(0,len(t)-1):
            for i in range(1,m-1):
                x[i] = (alphaprime)*((U[j,i+1]-(2*U[j,i])+U[j,i-1])/(dx**2)) #Central 2nd Order Difference for all points besides the first and last indicy
            x[0] = (alphaprime)*((U[j,1]-(2*U[j,0])+bl)/(dx**2)) #Central 2nd Order Difference for first indicy
            x[m-1] = (alphaprime)*((br-(2*U[j,m-1])+U[j,m-2])/(dx**2)) #Central 2nd Order Difference for last indicy
            U[j+1,:]= (U[j,:]+(x[:]*dt)) #forward 1st order difference to slove left hand side
        xn = np.linspace(dx/2, xf - dx/2, m)
        if self.Temperature_0 == None:
            if br >= bl:
                Th = br
            else:
                Th = bl
        else:
            if self.Temperature_0 > bl and self.Temperature_0 > br: #to find Highest Temperature
                Th = self.Temperature_0
            elif br >= bl:
                Th = br
            else:
                Th = bl
        if self.Gif == True:
            self.gif_maker(xn,U,xf,Th,t)
        elif self.Graph == True:
            self.graph_sketcher(xn,U,xf,Th,t)
        else:
            return U
    def graph_sketcher(self,xn,U,xf,Th,t):
        """
        Creates a plot of the temperature function in real time that is updating itself
        Note: this is slow but I find I can run many more points
        """
        if self.Hat_Function == True:
            plt.title(f'{self.Material} with boundary conditions {self.Left_Boundary} C and {self.Right_Boundary} C and {self.Temperature_0} C hat function')
        else:
            plt.title(f'{self.Material} with boundary conditions {self.Left_Boundary} C and {self.Right_Boundary} C')
        plt.axis([0,xf,0,Th])
        plt.xlabel('Position (m)')
        plt.ylabel('Temperature(C)')
        for p in range(len(t)):
            plt.clf()
            plt.figure(1)
            if self.Hat_Function == True:
                plt.title(f'{self.Material} with boundary conditions {self.Left_Boundary} C and {self.Right_Boundary} C and {self.Temperature_0} C hat function')
            else:
                plt.title(f'{self.Material} with boundary conditions {self.Left_Boundary} C and {self.Right_Boundary} C')
            plt.axis([0,xf,0,Th])
            plt.xlabel('Position (m)')
            plt.ylabel('Temperature(C)')
            plt.plot(xn,U[p,:])
            plt.pause(0.0005)
        plt.show()
    def gif_maker(self,xn,U,xf,Th,t):
        """
        Gif writing function passes on the Max values of the temperature, length, the x axis line space
        and the Termperature array with the first dimension being the time and the second being the temperature.
        Note: This seems function seems to crash on me aytime there is more than 121 Frames, unsure if its code related or 
        computer related.
        """
        giflists=[]
        for p in range(len(t)):
            plt.figure(1)
            plt.plot(xn,U[p,:])
            if self.Hat_Function == True:
                plt.title(f'{self.Material} with boundary conditions {self.Left_Boundary} C and {self.Right_Boundary} C and {self.Temperature_0} C hat function')
            else:
                plt.title(f'{self.Material} with boundary conditions {self.Left_Boundary} C and {self.Right_Boundary} C')
            plt.axis([0,xf,0,Th])
            plt.xlabel('Position (m)')
            plt.ylabel('Temperature(C)')
            filename = f'{p}.png'
            giflists.append(filename)
            plt.savefig(filename)
            plt.close()
        if self.Hat_Function == True:
            with imageio.get_writer(f'{self.Material}_{self.Left_Boundary}_{self.Right_Boundary}_HAT.gif', mode ='I') as writer:
                for giflist in giflists:
                    image = imageio.imread(giflist)
                    writer.append_data(image)
                for giflist in giflists:
                        os.remove(giflist)
        else:
            with imageio.get_writer(f'{self.Material}_{self.Left_Boundary}_{self.Right_Boundary}.gif', mode ='I') as writer:
                for giflist in giflists:
                    image = imageio.imread(giflist)
                    writer.append_data(image)
                for giflist in giflists:
                    os.remove(giflist)
class TestHeatEquation(unittest.TestCase):
    def testHeatEqnMatCon(self):
        """Tests to make sure mat_constants pulls correct information"""
        expected = 225.94
        heat = Heat_Eqn('Aluminum',1000, 'ellipse', (1.0,1.0),1,1000,100,100)
        K,alpha,c,e,rho = heat.mat_constants()
        actual = K
        self.assertEqual(expected,actual, "Material information extracted properly")
    def testHeatEqnMatConFail(self):
        """Tests to make sure the mat_constats handle failure properly"""
        expected = None
        heat = Heat_Eqn('Aluminium',1000, 'ellipse', (1.0,1.0),1,1000,100,100)#like the british
        K= heat.mat_constants()
        actual = K
        self.assertEqual(expected,actual, "Material information extraction Fails as expected")
    def testHeatEqnArea (self):
        """Tests the area function works properly"""
        expected = np.pi *2*3
        heat = Heat_Eqn('Aluminum',1000, 'ellipse', (2.0,3.0),1,1000,100,100)
        actual= heat.area()
        self.assertEqual(expected,actual,"Area properly delt with")
    def testHeatEqnAreazeroSL(self):
        """tests to see if input a zero/missing for side length is delt with as expected"""
        expected = 0
        heat = Heat_Eqn('Aluminum', 1000,'ellipse', (2.0,0.0),1,1000,100,100)
        actual= heat.area()
        self.assertEqual(expected,actual,"Area fails as expected")
    def testHeatEqnAreaPentagon(self):
        """tests to make sure area deals with unsupported type as expected"""
        expected = None
        heat = Heat_Eqn('Aluminum',1000, 'Pentagon', (2.0,1.0),1,1000,100,100)
        actual= heat.area()
        self.assertEqual(expected,actual,"Area fails as expected")
    def testHeatEqnAreaCapital(self):
        """tests to make sure area can handle capital leters"""
        expected = np.pi *2*3
        heat = Heat_Eqn('Aluminum', 1000, 'Ellipse', (2.0,3.0),1,1000,100,100)
        actual= heat.area()
        self.assertEqual(expected,actual,"Area calculates as expected")
if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    con = sqlite3.connect("Thermal.db")
    cur =con.cursor()
#    cur.execute("CREATE TABLE thermal_prop (Material TEXT, Conductivity REAL, Diffusivity REAL, Specific_Heat REAL, Effusivity REAL, Density REAL)")
    data = [("Air", 0.025, 19.4, 1004,6,1.29),
            ('Aluminum', 225.94, 91,921,23688,2698),
            ('Copper', 397.48,116,385,36983,8940),
            ('Gold',317.98,129,128,28027,19300),
            ('Iron',71.965,20.4,448,15924,7870),
            ('Magnesium',150.62,86.2,1004,16221,1740),
            ('Nitrogen',0.026,19.6,1042,6,1.251),
            ('Platnium',69.036,24.1,134,14065,21400),
            ('Plutonium',8.201,3.19,134,4592,19200),
            ('High Density Ployethylene',0.502,0.23,2301,1048,950),
            ('Medium Density Ployethylene',0.414,0.19,2301,944,935),
            ('Low Density Ployethylene',0.331,0.17,2092,798,920),
            ('Polycarbonate',0.192,0.09,1674,660,1350),
            ('Rock',1.757,0.81,837,1955,2600),
            ('Silver',426.77,172,236,32520,10500),
            ('Soil',0.837,0.62,1046,1067,1300),
            ('Steam',0.023,19.9,1966,5,0.598),
            ('Stainless Steel',22.928,6.56,460,8955,7600),
            ('Mild Steel', 41.84,10.8,502,12760,7750),
            ('Carbon Steel', 71.128,19.7,460,16040,7860),
            ('Titanium', 20.92,8.89,523,7017,4500),
            ('Tungsten', 196.65,76.1,134,22543,19300),
            ('Uranium',26.443,11.8,117,7694,19100),
            ('Water',0.603,0.14,4184,1588,1000),
            ('Ice',2.092,0.54,4217,2844,917),
            ('Wood Spruce With Grain',0.23,0.45,1255,344,410),
            ('Wood Spruce Across Grain',0.126,0.24,1255,254,410),
            ('Wood Teak With Grain',0.172,0.12,2301,503,640),
            ('Wood Teak Across Grain',0.142,0.12,2301,420,540),
            ('Wood Fir Across Grain',0.109,0.11,2301,336,450),
            ('Wood Pine Acrosse Grain',0.13,0.11,2301,398,530),
            ('Wood Maple Across Grain',0.176,0.11,2301,536,710)
            ]
#    for row in data:
#        cur.execute("INSERT INTO thermal_prop VALUES (?, ?, ?, ?, ?, ?)",row)
#    con.commit()
#    cur.execute("SELECT * from Thermal_Prop") # Pulls all the data from the database
#    header = list(map(lambda x: x[0], cur.description)) #pulls the Headers from the table Thermal_Prop
#    data = cur.fetchall()
#    con.close()
#    print(header)
#    pprint(data)
#    heat = Heat_Eqn('Aluminum',100, 'ellipse', (1.0,1.0),1,6000,100,100,False,True)#set up for graph
    heat = Heat_Eqn('Aluminum',10, 'ellipse', (1.0,1.0),1,6000,100,100,True)#set up for a quick gif
#    print(heat.mat_constants())
    heat.heateqn()
