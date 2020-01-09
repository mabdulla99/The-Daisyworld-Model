import numpy as np 
import matplotlib.pyplot as plt

plt.rc('xtick',labelsize = 15) 
plt.rc('ytick',labelsize = 15) # for fontsize of plot axes ticks

############################## defining initial conditions

init_g = 1 
init_w = 0 
init_b = 0
S = 917 # W/m^2 
sigma = 5.67e-8 # W/K^4 
albedo_b = 0.25 
albedo_w = 0.75 
albedo_g = 0.5 
q = 20 
gamma = 0.3

tmax = 200 
dt = 0.01 
N = int(tmax/dt) 
time = np.linspace(0, tmax, N+1)

############################## defining functions

def L(t): 
  return 0.6 + 1.2*(t/tmax)

def albedo_p(t, alpha_w, alpha_b, alpha_g): 
  return alpha_w*albedo_w + alpha_b*albedo_b + alpha_g*albedo_g
  
def T_p(t, alpha_w, alpha_b, alpha_g): 
  return ((S*L(t)*(1-albedo_p(t, alpha_w, alpha_b, alpha_g)))/(sigma))**(0.25) - 273

def T_w(t, alpha_w, alpha_b, alpha_g): 
  return q*(albedo_p(t, alpha_w, alpha_b, alpha_g) - albedo_w) + T_p(t, alpha_w, alpha_b, alpha_g)
  
def T_b(t, alpha_w, alpha_b, alpha_g): 
  return q*(albedo_p(t, alpha_w, alpha_b, alpha_g) - albedo_b) + T_p(t, alpha_w, alpha_b, alpha_g)
  
def G_w(t, alpha_w, alpha_b, alpha_g): 
  return 1 - ((1/17.5)**2)*(22.5 - T_w(t, alpha_w, alpha_b, alpha_g))**2
  
def G_b(t, alpha_w, alpha_b, alpha_g): 
  return 1 - ((1/17.5)**2)*(22.5 - T_b(t, alpha_w, alpha_b, alpha_g))**2
  
def alphaw_dot(t, alpha_w, alpha_b, alpha_g): 
  return alpha_w*(alpha_g*G_w(t, alpha_w, alpha_b, alpha_g) - gamma) + 0.001
  
def alphab_dot(t, alpha_w, alpha_b, alpha_g): 
  return alpha_b*(alpha_g*G_b(t, alpha_w, alpha_b, alpha_g) - gamma) + 0.001
  
def alphag_dot(t, alpha_w, alpha_b, alpha_g): 
  return alphaw_dot(t, alpha_w, alpha_b, alpha_g) - alphab_dot(t, alpha_w, alpha_b, alpha_g)
