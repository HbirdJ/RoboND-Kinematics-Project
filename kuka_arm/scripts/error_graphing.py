from mpmath import *
from sympy import *
from sympy import pi
import matplotlib as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pickle   
from time import strftime

q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

class kuka_path(object):

    def __init__(self):
        self.time = strftime("%Y-%m-%d %H:%M:%S")
        
        self.Px = []
        self.Py = []
        self.Pz = []
        
        self.poses = 0
        
        self.error = []
        
        self.theta1 = []
        self.theta2 = []
        self.theta3 = []
        self.theta3DH = []
        self.theta4 = [] 
        self.theta5 = [] 
        self.theta6 = []
        
        
    def test_data(self):
        self.Px = [2.153, 1.63336344034018, 0.678332525582959]
        self.Py = [0, 0.969143676640062, 1.3388888170805]
        self.Pz = [1.946, 2.3751247189682, 2.74260627437288]
        
        self.poses = 3
        
        self.error = [.25, .5, 1]
        
        self.theta1 = [0, pi/6, pi/3] 
        self.theta2 = [0,-pi/20,-pi/10] 
        self.theta3DH = [0,-pi/20,-pi/10] 
        self.theta4 = [0,pi/16,pi/8] 
        self.theta5 = [0,pi/8,pi/4] 
        self.theta6 = [0,pi/20,pi/10]
    
    def FK(self):
        if self.poses == 0:
            print("No poses in object")
            return
        
        s = { alpha0:     0, a0:     0, d1:  .75, q1: q1,
              alpha1: -pi/2, a1:   .35, d2:    0, q2: q2-pi/2,
              alpha2:     0, a2:  1.25, d3:    0, q3: q3,
              alpha3: -pi/2, a3: -.054, d4:  1.5, q4: q4,
              alpha4:  pi/2, a4:     0, d1:    0, q5: q5,
              alpha5: -pi/2, a5:     0, d1:    0, q6: q6,
              alpha6:     0, a6:  1.25, d7: .303, q7: 0,}

        t0_1 = Matrix([[0], [0], [d1 - 0.42]])
        t0_2 = Matrix([[a1*cos(q1)], [a1*sin(q1)], [d1]])
        t0_3 = Matrix([[(a1 + a2*sin(q2))*cos(q1)], [(a1 + a2*sin(q2))*sin(q1)], [a2*cos(q2) + d1]])
        t0_4 = Matrix([[(a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*cos(q1) - 0.54*cos(q1)*cos(q2 + q3)], [(a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*sin(q1) - 0.54*sin(q1)*cos(q2 + q3)], [a2*cos(q2) + a3*cos(q2 + q3) + d1 - d4*sin(q2 + q3) + 0.54*sin(q2 + q3)]])
        t0_5 = Matrix([[(a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*cos(q1)], [(a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*sin(q1)], [a2*cos(q2) + a3*cos(q2 + q3) + d1 - d4*sin(q2 + q3)]])
        t0_6 = Matrix([[-0.193*(sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) + (a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*cos(q1) + 0.193*cos(q1)*cos(q5)*cos(q2 + q3)], [-0.193*(sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) + (a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*sin(q1) + 0.193*sin(q1)*cos(q5)*cos(q2 + q3)], [a2*cos(q2) + a3*cos(q2 + q3) + d1 - d4*sin(q2 + q3) - 0.193*sin(q5)*cos(q4)*cos(q2 + q3) - 0.193*sin(q2 + q3)*cos(q5)]])
        t0_G = Matrix([[-d7*((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) - cos(q1)*cos(q5)*cos(q2 + q3)) + (a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*cos(q1)], [-d7*((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) - sin(q1)*cos(q5)*cos(q2 + q3)) + (a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*sin(q1)], [a2*cos(q2) + a3*cos(q2 + q3) + d1 - d4*sin(q2 + q3) - d7*(sin(q5)*cos(q4)*cos(q2 + q3) + sin(q2 + q3)*cos(q5))]])
        T0_G = Matrix([
            [-(sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) + cos(q1)*cos(q5)*cos(q2 + q3), ((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*cos(q5) + sin(q5)*cos(q1)*cos(q2 + q3))*sin(q6) - (sin(q1)*cos(q4) - sin(q4)*sin(q2 + q3)*cos(q1))*cos(q6), ((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*cos(q5) + sin(q5)*cos(q1)*cos(q2 + q3))*cos(q6) + (sin(q1)*cos(q4) - sin(q4)*sin(q2 + q3)*cos(q1))*sin(q6), -d7*((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) - cos(q1)*cos(q5)*cos(q2 + q3)) + (a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*cos(q1)],
            [-(sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) + sin(q1)*cos(q5)*cos(q2 + q3), ((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*cos(q5) + sin(q1)*sin(q5)*cos(q2 + q3))*sin(q6) + (sin(q1)*sin(q4)*sin(q2 + q3) + cos(q1)*cos(q4))*cos(q6), ((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*cos(q5) + sin(q1)*sin(q5)*cos(q2 + q3))*cos(q6) - (sin(q1)*sin(q4)*sin(q2 + q3) + cos(q1)*cos(q4))*sin(q6), -d7*((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) - sin(q1)*cos(q5)*cos(q2 + q3)) + (a1 + a2*sin(q2) + a3*sin(q2 + q3) + d4*cos(q2 + q3))*sin(q1)],
            [                                    -sin(q5)*cos(q4)*cos(q2 + q3) - sin(q2 + q3)*cos(q5),                                                                -(sin(q5)*sin(q2 + q3) - cos(q4)*cos(q5)*cos(q2 + q3))*sin(q6) + sin(q4)*cos(q6)*cos(q2 + q3),                                                                -(sin(q5)*sin(q2 + q3) - cos(q4)*cos(q5)*cos(q2 + q3))*cos(q6) - sin(q4)*sin(q6)*cos(q2 + q3),                                                a2*cos(q2) + a3*cos(q2 + q3) + d1 - d4*sin(q2 + q3) - d7*(sin(q5)*cos(q4)*cos(q2 + q3) + sin(q2 + q3)*cos(q5))],
            [                                                                                       0,                                                                                                                                                            0,                                                                                                                                                            0,                                                                                                                                                             1]])

        transforms = [t0_1, t0_2, t0_3, t0_4, t0_5, t0_6, t0_G]

        sol_start = {q1: path.theta1[0], q2: path.theta2[0], q3: path.theta3DH[0], q4: path.theta4[0], q5: path.theta5[0], q6: path.theta6[0], 
                                     d1:  .75, a1:   .35, a2:  1.25, a3: -.054, d4:  1.5, d7: .303}

        length = self.poses-1

        sol_end = {q1: path.theta1[length], q2: path.theta2[length], q3: path.theta3DH[length], q4: path.theta4[length], q5: path.theta5[length], q6: path.theta6[length], 
                                     d1:  .75, a1:   .35, a2:  1.25, a3: -.054, d4:  1.5, d7: .303}

        self.kuka_start_x = [0]
        self.kuka_start_y = [0]
        self.kuka_start_z = [0]
        self.kuka_end_x = [0]
        self.kuka_end_y = [0]
        self.kuka_end_z = [0]

        for trans in transforms:
            self.kuka_start = trans.evalf(subs=sol_start)
            self.kuka_start_x.append(self.kuka_start[0])
            self.kuka_start_y.append(self.kuka_start[1])
            self.kuka_start_z.append(self.kuka_start[2])
            self.kuka_end = trans.evalf(subs=sol_end)
            self.kuka_end_x.append(self.kuka_end[0])
            self.kuka_end_y.append(self.kuka_end[1])
            self.kuka_end_z.append(self.kuka_end[2])
    
    def datadump(self):
        pickle.dump(self, open("datadump.dat", "ab"))
    
    def plot(self):

        fig = plt.figure(figsize=(20, 14))
        mpl.rcParams['font.size'] = 16
        ax = fig.add_subplot(111, projection='3d', xlim=[-1,3.5], ylim=[-1,3.5], zlim=[0,3.5]) 
        ax.view_init(elev=25, azim=210)
        
        marks = [1,2,3,4,5,6]
        
        ax.plot(self.kuka_start_x, self.kuka_start_y, self.kuka_start_z, c='grey', linestyle='--',
                marker='*',  markevery=marks, label='start position')
        ax.plot(self.kuka_end_x, self.kuka_end_y, self.kuka_end_z, c='black', linestyle='--', 
                marker='*', markevery=marks, label='end position')
        p = ax.scatter(self.Px, self.Py, self.Pz, c=self.error, cmap='Reds', marker='^',
                   depthshade=False, label='calculated gripper positions')
        ax.legend()
    
        fig.colorbar(p, label='FK vs IK position error')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
    

        plt.savefig('test.png', bbox_inches='tight')
        plt.close()
        
        return