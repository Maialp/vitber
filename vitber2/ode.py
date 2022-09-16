import numpy as np
from math import factorial
import matplotlib.pyplot as plt

#Constant parametes 
sigma_0 = 1000              #water density [kg/m^2] 
sigma = 500                 #Ship density [kg/m^2] 
R = 10                      #Ship radius [m]
A_s = 1/2 * np.pi * R**2    #Ship cross section 
m = A_s * sigma      #Ship's mass per length unit (in z direcion)
A_0 = sigma*np.pi*R**2 / (2*sigma_0)
I_c = 1/2 * m * R**2 * (1-(32)/(9*np.pi**2)) #Ships moment of intertia      #/(9*np.pi**2) Legg til for riktig I_c
g = -981/100          #Gravitational acceleration9*np.pi**2))   
h = 0.42* R                 #Distanve from center of mass to physical senter
w_0 = np.sqrt(m*g*h/I_c)

#Oppgave 1b)
def newton(f, df, x0, m = 0, tol=1.e-8, max_iter=30):
    ''' Solve f(x)=0 by Newtons method
        The output of each iteration is printed
        Input:
        f, df:   The function f and its derivate f'.
        x0:  Initial values
        tol: The tolerance
      Output:
        The root and the number of iterations
    '''
    x = x0
    #print(f"k ={0:3d}, x = {x:18.15f}, f(x) = {f(x):10.3e}")
    for k in range(max_iter):
        fx = f(x,m)
        if abs(fx) < tol:           # Aksepterer løsningne 
            break 
        x = x - fx/df(x)            # Newton-iterasjon
        #print(f"k ={k+1:3d}, x = {x:18.15f}, f(x) = {f(x):10.3e}")
    return x, k+1

def f3(x, m=0):
    sigma_m = m/(1/2 *np.pi*R**2)
    return x - np.sin(x) -np.pi*((sigma+sigma_m)/sigma_0)

def df3(x):
    return 1-np.cos(x)

def get_beta(m=0): #Hvor m=load mass. Brukes til oppgave 2
    beta, its = newton(f3, df3, 1 , m, tol=1.e-8, max_iter=30)
    return beta

beta = get_beta()

def Y_MB(gamma=beta):
    '''
    Input: Gamma (sector angle)
    Output: 
    '''
    return R * (4*(np.sin(gamma/2)**3))/(3*(gamma-np.sin(gamma)))

def RK4_step(f, t_i, w_i, h):
    ''' 
    One step using Runge_kutta
    Intput: 
        w_i: [angle theta, frequency ohmega], 
        t_i: time
        h: step length
        f: function to be solved
    Return: next function value
    '''
    k1 = f(t_i,w_i)
    k2 = f(t_i + h/2, w_i + h*k1/2)
    k3 = f(t_i + h/2, w_i + h*k2/2)
    k4 = f(t_i+ h, w_i + h*k3)

    return t_i + h,  w_i + h/6 *(k1+ 2*k2 + 2*k3 + k4)

y_M0 = R * np.cos(beta / 2)
y_C0 = y_M0 - h
y_B0 = y_M0 - Y_MB()
y_D0 = y_M0 - R

y_M0 = R * np.cos(beta / 2)
y_C0 = y_M0 - h
y_B0 = y_M0 - Y_MB()
y_D0 = y_M0 - R

class Enviroments:
    def __init__(self, k_f=0.8, mu_last=0.3, F_0 = 0, w_w =0):
        self.k_f = k_f #water friction coef.
        self.mu_last = mu_last #load friction coef.
        self.F_0 = F_0 #wind force amplitude
        self.w_w = w_w #wind force frequency


class BoatDynamics:
    def __init__(self, fence, loadmass, theta_area, env):
        #enviroments. Such as wind, friction coefficients etc.¨
        self.env = env

        #simulation parameters
        self.fence = fence #bolean. Fences on boat sides.
        self.loadmass = loadmass #mass of load beeing caried
        self.theta_area = theta_area #sier om vi skal regne på areal som funksjon av theta

        #some events:
        self.load_at_rest = False #tells the load to rest.
        self.load_of_board = False  #load has fallen ofboard.
        self.boat_tipped = False #boat has tipped

    def f(self, t, state):

        #derivert av state:
        state_dot = np.zeros(len(state))
        
        #pakker ut state vektoren
        x_dot, y_dot = state[0], state[1]
        x, y = state[2], state[3]
        theta_dot = state[4]
        theta = state[5]
        s_Ldot = state[6]
        s_L = state[7]

        beta = get_beta(self.loadmass)
        A0 = 0.5 * R**2 * (beta - np.sin(beta))

        delta_y = y - (R * np.cos(beta / 2) - h) # y - y_C0, tyngde punktets relative bevegelse.

        #sektor vinkel:
        gamma = 2 * np.arccos(np.cos(beta/2) - 4 / (3 * np.pi) * (1 - np.cos(theta)) + delta_y / R)
 
        #krefter på lasten + litt logikk for bevegelse med gjerde
        if(self.loadmass != 0):
            friction_x =  self.env.mu_last * self.loadmass * abs(g) * np.cos(theta) #friksjon på lasten
            gravity_x =  self.loadmass * abs(g) * np.sin(theta) #tyngdekraft på lasten
            
            F_sx = - np.sign(theta) * abs(friction_x) + np.sign(theta) * abs(gravity_x)
            F_sy = - self.loadmass * abs(g) * np.cos(theta) #y-komponent
            F_Ly = - self.loadmass * abs(g) * (np.cos(theta))**2 #kontakt kraft fra last i x retning
            F_Lx = self.loadmass * abs(g) * np.cos(theta) * np.sin(theta) # -=- i y retning.


            if(abs(s_L) >= R and self.fence):
                if(abs(theta) < 0.1):
                    self.load_at_rest = False
                
                else:
                    self.load_at_rest = True



            elif(abs(s_L) >= R and not self.fence):
                self.load_of_board = True

        else:
            s_Ldot, F_sx, F_sy, F_Ly, F_Lx = 0, 0, 0, 0, 0 #Ingen last til å begynne med
        
        if(self.load_of_board):
            self.loadmass = 0
            s_Ldotdot, F_sx, F_sy, F_Ly, F_Lx = 0, 0, 0, 0, 0
        
        #areal
        if(self.theta_area): A = 0.5 * R**2 * (gamma - np.sin(gamma)) #fortrengt vann
        else: A = 0.5 * R**2 * (beta - np.sin(beta)) #ikke fortrengt vann

        #krefter på båten:
        F_B = A * sigma_0 * abs(g) #boyuansy force. Fancy french
        F_G = - (m + self.loadmass) * abs(g) #gravity
        f_small = - self.env.k_f * R * gamma * theta_dot #friksjons kraft fra vannet på båten
        F_w = self.env.F_0 * np.cos(self.env.w_w * t) #kraft fra vind

        #dreie momenter:
        tau_B = - F_B * h * np.sin(theta) #dreiemoment fra bouyancy kraft
        tau_f = f_small * (y - (R * (np.cos(gamma/2)) - 1))
        tau_w = F_w * y
        tau_L = F_sy * s_L

        #regner ut summen av kreftene på båten: 
        F_sum_x = f_small + F_w + F_Lx #newtons
        F_sum_y = F_G + F_B + F_Ly #newtons

        #regner ut summen av dreiemomenter på båten:
        tau_sum = tau_B + tau_f + tau_w + tau_L

        #lager ny state:
        state_dot[0] = F_sum_x / m  #akselerasjon i x retning
        state_dot[1] = F_sum_y / m  #akselerasjon i y retning
        state_dot[2] = x_dot #x -> xdot
        state_dot[3] = y_dot #y -> ydot
        state_dot[4] = tau_sum / I_c  #theta_dot -> theta_dotdot
        state_dot[5] = theta_dot #theta -> theta_dot

        if(self.loadmass != 0): state_dot[6] = F_sx / self.loadmass
        else: state_dot[6] = 0

        state_dot[7] = s_Ldot

        #sjekker for kantring:
        if(abs(theta) > (np.pi - gamma) / 2):
            self.boat_tipped = True

        return state_dot

def solve_boat_dynamics(x0, xend, y0, mass_load = 0.01*m, theta_area=True, fence=True, h=0.01, method=RK4_step, env=Enviroments()):
    if(len(y0) != 8): #om init state er mindre enn 8 legger vi til 0 på resten.
        for _ in range(8 - len(y0)):
            y0.append(0)
        
    tipped = False #har disse her slik at vi kun printer at båten har tippet en gang.
    load_ofboard = False

    dynamics = BoatDynamics(fence, mass_load, theta_area, env) #lager dynamikk objektet gitt om vi skal ha gjerde, last eller ta hensyn til at 
    #arealet endrer seg som funksjon av theta. Slik oppgaven spurte om i 2b.

    #initialiserer
    x_num = np.array([x0])
    y_num = np.array([y0])
    xn = x0
    yn = y0
    
    while xn <= xend:
        if(not dynamics.boat_tipped): #om båten ikke har flipper gjør vi en ny numerisk beregning:
            xn, yn = method(dynamics.f, xn, yn, h) #Do one step

        #sjekker flagg:
        if(dynamics.load_at_rest):
            yn[6] = 0 #setter lastens hastigheten til 0 om ode sier at lasten hviler.
            
        if(dynamics.load_of_board):
            yn[6], yn[7] = 0, 0 #setter posisjon og hastighet til 0. Massen blir satt til 0 i ode.
            if(not load_ofboard): #printer til konsoll
                load_ofboard = True
                print("Load fell of board at", round(xn, 3), "s")

        if(dynamics.boat_tipped):
            yn =  np.array([0, 0, yn[2], yn[3], 0, np.sign(yn[5])*np.pi / 2, 0, 0]) #lager ny state som er "flipped"
            xn += h
            if(not tipped): #printer til konsoll
                print("Boat tipped at", round(xn, 3), "s")
                tipped = True

        y_num = np.concatenate((y_num, np.array([yn]))) #legger til ny state i container.
        x_num = np.append(x_num,xn)
    
    if(not tipped): #printer viktige events:
        print("Boat did not tip during the simulation.")
    if(not load_ofboard and dynamics.loadmass != 0):
        print("Load mass did not fall over board.")

    return  x_num, y_num