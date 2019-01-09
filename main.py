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
    def __init__(self, epsilon, inc, filename, nf):
        self.epsilon = epsilon
        self.inc = inc
        self.observed = np.loadtxt(filename, skiprows=4)
        self.nf = nf
        self.F0 = np.mean(self.observed[:,nf])
    
    def Split(self, N, M):
        self.subpartition = np.empty(shape=(N,M), dtype=object)
        step_th = np.pi/N
        step_phi = 2*np.pi/M
        
        self.normal_intens_average = np.zeros(shape=(N*M,1))
        self.normal_intens_average += self.F0/(np.pi*(1-self.epsilon/3))
        self.normal_intens = np.zeros(shape=(N*M,1))
        self.X = np.zeros(shape=(N*M,1))
        self.fluxes = self.observed[:,self.nf]
        self.fluxes.shape = (len(self.fluxes),1)
        self.G = np.zeros(shape=(len(self.fluxes),N*M))
        self.E = np.eye(N*M)
        self.N = N
        self.M = M
        
        count = 0
        for i in range(N):
            for j in range(M):
                self.subpartition[i,j] = Sector(step_th*i,step_th*(i+1),step_phi*j,step_phi*(j+1))
                self.subpartition[i,j].index = count
                count += 1
                
        for p in range(len(self.observed[:,0])):
            for i in range(N):
                for j in range(M):
                    index = self.subpartition[i,j].index
                    cosg = self.subpartition[i,j].cos_gamma(self.observed[p,0]*2*np.pi,self.inc)
                    A = self.subpartition[i,j].A(cosg)
                    L = self.subpartition[i,j].LimbDarkening(cosg,self.epsilon)
                    self.G[p,index] = A*L
        
    def SolveTikhonov(self, alpha):
        GT = self.G.transpose()
        mat1 = np.linalg.inv(GT.dot(self.G) + alpha*self.E)
        self.X = (mat1.dot(GT)).dot(self.fluxes-self.G.dot(self.normal_intens_average))
        self.normal_intens = self.X + self.normal_intens_average
        
    def GetTheorCurve(self, anyphases):
        tf = np.zeros(shape=(len(anyphases)))
        count = 0
        for ap in anyphases:
            for i in range(self.N):
                for j in range(self.M):
                    index = self.subpartition[i,j].index
                    cosg = self.subpartition[i,j].cos_gamma(ap*2*np.pi,self.inc)
                    A = self.subpartition[i,j].A(cosg)
                    L = self.subpartition[i,j].LimbDarkening(cosg,self.epsilon)
                    tf[count] += A*L*self.normal_intens[index,0] 
            count += 1
        return tf

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
                ax.set_title(np.round(phase,3), fontsize=10)
            plt.show()
    
####### собсно XX Tri

N = 36
M = 72
inc = 50
LD = 0.6
regpar = 0.0003
filename = 'LightCurve.dat'

XXTri = StarSurface(LD, inc*np.pi/180, filename, 2)
XXTri.Split(N,M)
XXTri.SolveTikhonov(regpar)

plt.plot(XXTri.observed[:,0],XXTri.observed[:,XXTri.nf],'ok')
tap = np.arange(0,1,0.01)
plt.plot(tap,XXTri.GetTheorCurve(tap),'-r')
plt.plot(XXTri.observed[:,0],XXTri.G.dot(XXTri.normal_intens), 'og')
plt.show()

####### отрисовка поверхности на разные фазы 

XXTriPlot = ColouredSphere(XXTri)
XXTriPlot.DrawGrid(inc, 8, 2, 4)

'''####### сердешко

N = 36
M = 72
inc = 90
LD = 0.0
regpar = 0.000002
filename = 'heart.dat'

heart = StarSurface(LD, inc*np.pi/180, filename, 2)
heart.Split(N,M)
heart.SolveTikhonov(regpar)

plt.plot(heart.observed[:,0],heart.observed[:,heart.nf],'ok')
tap = np.arange(0,1,0.01)
plt.plot(tap,heart.GetTheorCurve(tap),'-r')
plt.plot(heart.observed[:,0],heart.G.dot(heart.normal_intens), 'og')
plt.show()

####### отрисовка поверхности на разные фазы

heartPlot = ColouredSphere(heart)
heartPlot.DrawGrid(inc, 10, 2, 5)'''