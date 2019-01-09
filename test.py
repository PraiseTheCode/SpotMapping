import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Sector:
    def __init__(self, theta_upper, theta_lower, phi_left, phi_right):
        self.theta_upper = theta_upper
        self.theta_lower = theta_lower
        self.phi_left = phi_left
        self.phi_right = phi_right
        
        self.phi_cent = phi_left + (phi_right-phi_left)/2
        self.theta_cent = theta_upper + (theta_lower - theta_upper)/2
        self.S = (phi_right-phi_left)*(np.cos(self.theta_upper)-np.cos(self.theta_lower))
    def cos_gamma(self, phase, inc):
        return np.cos(inc)*np.cos(self.theta_cent) + np.sin(inc)*np.sin(self.theta_cent)*np.cos(phase + self.phi_cent)
    def A(self, cosg):
        if (cosg > 0):
            return self.S*cosg
        else :
            return 0
    def LimbDarkening(self, cosg, epsilon):
        return 1-epsilon*(1-cosg)
        
class StarSurface:
    def __init__(self, epsilon, inc):
        self.epsilon = epsilon
        self.inc = inc
    
    def Split(self, N, M):
        self.subpartition = np.empty(shape=(N,M), dtype=object)
        step_th = np.pi/N
        step_phi = 2*np.pi/M
        self.phases = np.arange(0,1,0.01)
        
        self.normal_intens_average = np.ones(shape=(N*M,1))
        self.normal_intens = np.ones(shape=(N*M,1))
        self.G = np.zeros(shape=(len(self.phases),N*M))
        self.N = N
        self.M = M
        
        count = 0
        for i in range(N):
            for j in range(M):
                self.subpartition[i,j] = Sector(step_th*i,step_th*(i+1),step_phi*j,step_phi*(j+1))
                self.subpartition[i,j].index = count
                count += 1
                
        for p in range(len(self.phases)):
            for i in range(N):
                for j in range(M):
                    index = self.subpartition[i,j].index
                    cosg = self.subpartition[i,j].cos_gamma(self.phases[p]*2*np.pi,self.inc)
                    A = self.subpartition[i,j].A(cosg)
                    L = self.subpartition[i,j].LimbDarkening(cosg,self.epsilon)
                    self.G[p,index] = A*L
                        
    def SetSpot(self, center_theta, center_phi, rad):
        step_theta = np.pi/self.N
        step_phi = 2*np.pi/self.M
        iin = 0
        iim = 0
        for i in range(self.N):
            for j in range(self.M):
                if ((np.abs(self.subpartition[i,j].theta_cent - center_theta)<=(step_theta/2))):
                    if ((np.abs(self.subpartition[i,j].phi_cent - center_phi)<=((step_phi)/2+0.001))):
                        #print('!!!')
                        #print(self.subpartition[i,j].theta_cent*180/np.pi)
                        #print(self.subpartition[i,j].phi_cent*180/np.pi)
                        iin = i
                        iim = j
        for i in range(self.N):
            for j in range(self.M):
                if ((np.sqrt((iin-i)**2 + (iim-j)**2)<rad) or
                (np.sqrt((iin-i)**2 + (iim-j+2*np.pi/step_phi)**2)<rad)) :
                    index = self.subpartition[i,j].index
                    self.normal_intens[index] = 0
                        
                    
class ColouredSphere():
    def __init__(self, star):
        self.colors = np.zeros(shape=(N,M,4))
        count = 0
        maxint = np.max(star.normal_intens)
        for i in range(star.N):
            for j in range(star.M):
                aaa = star.normal_intens[count,0]/maxint
                if(aaa < 0): aaa = 0
                self.colors[i][j][0] = aaa
                self.colors[i][j][1] = aaa
                self.colors[i][j][2] = aaa
                self.colors[i][j][3] = 1
                count += 1
        self.colors = self.colors.transpose(1,0,2)
        
        stheta, sphi = np.linspace(0, np.pi, star.N), np.linspace(0, 2*np.pi, star.M)
        STHETA, SPHI = np.meshgrid(stheta, sphi)
        self.X = np.sin(STHETA) * np.cos(SPHI)
        self.Y = np.sin(STHETA) * np.sin(SPHI)
        self.Z = np.cos(STHETA)
        
    def DrawGrid(self, inc, nphases, rows, columns):
            ps = 1.0/nphases
            fig = plt.figure()
            for ii in range(nphases):
                phase = ps*ii
                ax = fig.add_subplot(rows,columns,ii+1, projection='3d', xmargin=0, ymargin=0, xlim=(-0.65,0.65), ylim=(-0.65,0.65), zlim=(-0.65,0.65))
                ax.set_aspect('equal')
                ax.set_axis_off()
                plot = ax.plot_surface(
                        self.X, self.Y, self.Z, rstride=1, cstride=1, 
                        facecolors = self.colors,
                        linewidth=0, antialiased=False, alpha=1)
                ax.view_init(elev=(90-inc), azim=-360*phase)
                ax.set_title(round(phase,3), fontsize=10)
            plt.show()
    
#######

N = 36
M = 72
inc = 90
LD = 0.0

test = StarSurface(LD, inc)
test.Split(N,M)

spot_theta = 80*np.pi/180
spot_phi = 0*np.pi/180
rad_sect = 5
test.SetSpot(spot_theta, spot_phi, rad_sect)

'''spot_theta = 80*np.pi/180
spot_phi = 45*np.pi/180
rad_sect = 5
test.SetSpot(spot_theta, spot_phi, rad_sect)

spot_theta = 95*np.pi/180
spot_phi = 25*np.pi/180
rad_sect = 5
test.SetSpot(spot_theta, spot_phi, rad_sect)

fluxes = test.G.dot(test.normal_intens)

plt.plot(test.phases, fluxes,'ok')
plt.show()
with open('D:/heart.dat',"w") as out:
    out.write('1'+'\n')
    out.write('2'+'\n')
    out.write('3'+'\n')
    out.write('4'+'\n')
    for i in range(len(fluxes)):
        row = str(test.phases[i]) + '\t' + str(-2.5*np.log10(fluxes[i,0]/3.14)) + '\t' + str(fluxes[i,0]) + '\n'
        out.write(row)'''


testPlot = ColouredSphere(test)
testPlot.DrawGrid(inc, 10, 2, 5)