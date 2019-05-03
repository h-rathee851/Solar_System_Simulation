#I think I need 2 classes, one that creates the planet objects and calculates an array for position
#velocity, acceleration, Energy etc (methods), then a simulation class to take those values and 
#visually simulate them #need to build simulation on computer at lab
#methods will be 
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
class MakePlanet():
    G = 6.67*(10**-11)
    dt = 10000000
    msun = 1.989*(10**30)
    def __init__(self, xpos, ypos, mass, vx, vy):
        self.xpos = xpos
        self.ypos = ypos
        self.radius = np.array([self.xpos, self.ypos])
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.v = np.array([self.vx, self.vy])
        self.a = np.array([0.0,0.0])
        self.a_prev = self.a
        self.a_new = self.a
        #self.f = self.init_force(MakePlanet.G,MakePlanet.msun)
    def init_accel(self,G):
        #self.f = G*msun*self.mass/(self.radius)
       # f = np.array([0.0,0.0])
        with open('Compsimfile.txt','r') as f:
            data = np.loadtxt('Compsimfile.txt', delimiter = " ", skiprows=2, usecols=[1,2,3])
            for line in data:
                for i in range(len(line)):
                    if self.mass != data[i][2]: 
                        self.a += G*(data[i][2])*(np.array([self.xpos, self.ypos])-(np.array([data[i][0], 0])))/(np.linalg.norm(np.array([self.xpos, self.ypos])-(np.array([data[i][0], 0])))**3)
                    else:
                        self.a += 0
           # print(self.a)

    def calc_v(self):
        self.v 
        #print(self.v)
    def k_energy(self):
        self.ke = .5*self.mass*(norm(self.v))**2 
    def potential_energy(self):
        self.pe = 0
        #want to use update acceleration method in here
        with open('Compsimfile.txt','r') as f:
            data = np.loadtxt('Compsimfile.txt', delimiter = " ", skiprows=2, usecols=[1,2,3])
            for line in data:
                for i in range(len(line)):
                    if self.mass != data[i][2]:
                        self.pe -= data[i][2]*(np.linalg.norm(self.radius)**2)*self.a/(np.array([self.xpos, self.ypos]-np.array([data[i][0], 0])))
                    else:
                        self.pe += 0
    def tot_energy(self):
        self.energy_array = ([])
        self.te = self.pe + self.ke
        self.energy_array.append(self.te)
        with open("Energy.txt", 'a') as f:
            f.write(self.te)
        
    def update_r(self,dt):
        self.radius = (self.radius + self.v*dt + ((1/6)*((4*self.a-self.a_prev)*(dt**2))))
    
    def update_v(self,dt):
        #MakeSimulation.calc_acc
        self.v = (self.v + ((1/6)*(2*self.a_new + 5*self.a - self.a_prev)*dt))
        self.a_prev = self.a
        self.a = self.a_new
        #print("BodyAcc =",self.a,self.a_prev,self.a_new)
        
        
class MakeSimulation():
    #calc force method is in simulation
    #do i need to call body class here?
    G = 6.67E-11
    def __init__(self):
        #manually create objects
        self.timesteps = 1000
        self.dt = 100000
        filename = 'Compsimfile.txt'
        self.read_data(filename)
        self.lis = []
        self.xmax = 2*np.pi
        #self.lis.append(MakePlanet(self.data[0,0], self.data[0,1], self.data[0,2]))
        #print(self.data[0,1])
        for i in range(0, 5):
            self.lis.append(MakePlanet(self.data[i,0], self.data[i,1], self.data[i,2], self.data[i,3], self.data[i,4]))
    
    def read_data(self,filename):
            with open(filename,'r') as f:
                self.data = np.loadtxt(filename, delimiter = " ", skiprows=2, usecols=[1,2,3,4,5])


    def run(self,G):
        #print(self.lis[0].mass)
        for p in range(0, len(self.lis)):
            self.lis[p].calc_v()
            #self.lis[p].init_accel(G)
            #print(self.lis[p].calc_v()
            self.calc_acc(G, p)
            self.lis[p].tot_energy
            self.lis[p].update_r(self.dt)
            self.lis[p].update_v(self.dt)
            #print(self.lis[2].radius)
                
#use norm and change grav equation write equivalency for r vector 
    def calc_acc(self,G, i):
        #j.accel = 0
            a = np.array([0.0,0.0])
        #print(G*self.data[1][2]*np.array([self.lis[1].xpos, self.lis[1].ypos])/(np.linalg.norm(np.array([self.lis[1].xpos,self.lis[1].ypos])-np.array([self.lis[2].xpos, self.lis[2].ypos]))))
        #print(np.linalg.norm(np.array([self.lis[1].xpos,self.lis[1].ypos])-np.array([self.lis[2].xpos, self.lis[2].ypos])))
        # i in range(0, len(self.lis)):
            for j in range(0, len(self.lis)):  
                if self.lis[i].mass != self.lis[j].mass:
                    #print("Math:",(np.linalg.norm(np.array([self.lis[i].radius])-np.array([self.lis[j].radius]))**3))
                    #print("Array",((self.lis[i].radius)-(self.lis[j].radius)))
                    #self.lis[i].a -= G*self.data[j][2]/((np.linalg.norm(self.lis[i].radius-self.lis[j].radius))**3)
                    a -= (G*self.data[j][2]*((self.lis[i].radius)-(self.lis[j].radius))/(np.linalg.norm(np.array([self.lis[i].radius])-np.array([self.lis[j].radius]))**3))
            self.lis[i].a_new = a
            #print("SimAcc =",a,self.lis[i].a_prev,self.lis[i].a_new)
    def calc_orb(self):
        self.timesteps = 0
        for i in range(0, len(self.lis)):
            if self.lis[i].ypos != 0:
                self.timesteps += 1
            else:
                self.period = (2*self.timesteps)/(self.lis[3].timesteps)
            print(self.lis[i], "period is", self.period)
    #method to graph total energy?
    def init(self):
        return self.patches    
    
    '''def runSim(self):
        theta = np.linspace(0, 2*np.pi, 500)
        self.xpos = np.cos(theta)
        self.ypos = np.sin(theta)
'''
    def animate(self,i):
        G = 6.67E-11
        self.run(G)
        for i in range(0, len(self.patches)):
            #self.patches[i].center = (self.lis[i].xpos,self.lis[i].ypos)
            self.patches[i].center = (self.lis[i].radius)
            #print((self.lis[i].xpos,self.lis[i].ypos))
        return self.patches
            
    def run_visual(self):
        fig = plt.figure()
        ax = plt.axes()
        self.patches = []
        color = ['r','b','g']
        for i in range(0, len(self.lis)):
            self.patches.append(plt.Circle(((self.lis[i].xpos,self.lis[i].ypos)), 10**10, color = 'b', animated = True))
        for i in range(0, len(self.patches)):
            ax.add_patch(self.patches[i])
        ax.axis('scaled')
        ax.set_xlim(-3E11, 3E11)
        ax.set_ylim(-3E11, 3E11)
        numFrames = 1000
        anim = FuncAnimation(fig, self.animate, numFrames, repeat=False, interval = 20, blit=True)
        plt.show()


        

    
            
   # def animate()
   # def display()





def main():
    #read in data to create planets
    '''
    names = ['Sun','Mercury','Venus','Earth','Mars']
    for i in range(0,len(names)):
        test_planet = MakePlanet(10, 0, 5E20, 10, 10)
    '''
    sim = MakeSimulation()
    #sim.calc_acc(6.67E-11)
    sim.run_visual()
    #sim.init()
    #sim.animate()
    #sim.run_visual()
    
main()
    
