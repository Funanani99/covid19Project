# -*- coding: utf-8 -*-
"""
Created on Thu Feb  29 01:01:19 2024

@author: Funanani.Raphulu
"""

from __future__ import division
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


#Epidemiological Model parameters
sigma= 0.60; # number of people recovering
beta_A = 0.00031; # contact rate
beta_S = 0.0003; # contact rate
beta_E = 0.020; # number of people in contact
pi_H = 0.9; #number of people recovering
Lambda_H = 0.0145; #birth rate
alpha_S= 0.1;    #Natural recovery rateand death rate
alpha_A = 0.4;    #Natural recovery rate and death rate
alpha_L=0.054; #imfection rate
alpha_E = 0.04;    #Natural recovery rate
V_0 = 100;
p = 0.15;
mu_H= 0.0161;

#intervations
alpha_V = 0.3; #Health Education (Quarantined)
#tr = 0; # Health treatments
phi_V = 0.5; # Vaccination
# Initial value conditions
N =12000;
L_H0 = 14;
I_S0 = 1;
I_A0 = 4;
S_H0 = N;
D_H0 = 0;
V_E0 = 0;
V_H0 = 0;
def diff_eqns(x,t):
	dxdt = np.zeros((7),'d')
	S_H = x[0]; L_H = x[1]; I_S = x[2]; I_A = x[3]; D_H = x[4]; V_E = x[5]; V_H = x[6];
    
	dxdt[0] = Lambda_H -(mu_H + phi_V + beta_S*I_S + beta_A*I_A +(beta_E*V_E)/(V_0 +  V_E))*S_H + alpha_V*V_H +pi_H*alpha_A*I_A + sigma*alpha_S*I_S;
	dxdt[1] = (beta_S*I_S + beta_A*I_A +(beta_E*V_E)/(V_0 +  V_E))*S_H - (mu_H + alpha_L)*L_H;
	dxdt[2] =p*alpha_L*L_H - (mu_H + alpha_S)*I_S;	      	
	dxdt[3] = (1-p)*alpha_L*L_H - (mu_H + alpha_A)*I_A;
	dxdt[4] = (1- pi_H)*alpha_A*I_A + (1- sigma)*alpha_S*I_S;   
	dxdt[5] = alpha_S*I_S +alpha_A*I_A - alpha_E*V_E;
	dxdt[6] = phi_V*S_H - (alpha_V + mu_H )*V_H;
	
	return dxdt

Initials = (S_H0, L_H0, I_S0, I_A0, D_H0, V_E0, V_H0)
t = np.arange(0,999,0.1)
M = odeint(diff_eqns,Initials,t)

# Plotting the population dynamics of humans



tre =[0.3, 0.5,0.90]
str_tau = "" #For legend
for phi_V in tre:
    M = odeint(diff_eqns,Initials,t)
    S_H = M[:,0]; L_H = M[:, 1]; I_S = M[:, 2]; I_A = M[:,3]; D = M[:, 4]; V_E = M[:, 5];V_H = M[:, 6];
    inf1 = S_H # Infected human population
    inf2 = L_H
    inf3 = I_S # Infected human population
    inf4 = I_A
    inf5 = D
    inf6 = V_E
    inf7 = V_H
    
    plt.subplot(221)
    plt.plot(t, inf3, lw=3)
    plt.grid()
    plt.legend(("$\\phi_V ="+str(tre[0])+"$","$\\phi_V="+str(tre[1])+"$","$\\phi_V ="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel("Infected symptomatic individuals")
    plt.yticks()
    
    plt.subplot(223)
    plt.plot(t, inf4, lw=3)
    plt.grid()
    plt.legend(("$\\phi_V ="+str(tre[0])+"$","$\\phi_V="+str(tre[1])+"$","$\\phi_V ="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel("Infected asymptomatic individuals")
    plt.yticks()
    
    plt.subplot(222)
    plt.plot(t, inf6, lw=3)
    plt.grid()
    plt.legend(("$\\phi_V ="+str(tre[0])+"$","$\\phi_V="+str(tre[1])+"$","$\\phi_V="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel('Viral load')
    plt.yticks()
    
    plt.subplot(224)
    plt.plot(t, inf5, lw=3)
    plt.grid()
    plt.legend(("$\\phi_V ="+str(tre[0])+"$","$\\phi_V="+str(tre[1])+"$","$\\phi_V ="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel('number of Death')
    plt.yticks()
plt.show()

tre =[0.3, 0.5,0.90]
str_tau = "" #For legend
for alpha_V in tre:
    M = odeint(diff_eqns,Initials,t)
    S_H = M[:,0]; L_H = M[:, 1]; I_S = M[:, 2]; I_A = M[:,3]; D = M[:, 4]; V_E = M[:, 5];V_H = M[:, 6];
    inf1 = S_H # Infected human population
    inf2 = L_H
    inf3 = I_S # Infected human population
    inf4 = I_A
    inf5 = D
    inf6 = V_E
    inf7 = V_H
    
    plt.subplot(221)
    plt.plot(t, inf3, lw=3)
    plt.grid()
    plt.legend(("$\\alpha_V ="+str(tre[0])+"$","$\\alpha_V="+str(tre[1])+"$","$\\alpha_V ="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel("Infected symptomatic individuals")
    plt.yticks()
    
    plt.subplot(223)
    plt.plot(t, inf4, lw=3)
    plt.grid()
    plt.legend(("$\\alpha_V ="+str(tre[0])+"$","$\\alpha_V="+str(tre[1])+"$","$\\alpha_V ="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel("Infected asymptomatic individuals")
    plt.yticks()
    
    plt.subplot(222)
    plt.plot(t, inf6, lw=3)
    plt.grid()
    plt.legend(("$\\alpha_V ="+str(tre[0])+"$","$\\alpha_V="+str(tre[1])+"$","$\\alpha_V="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel('Viral load')
    plt.yticks()
    
    plt.subplot(224)
    plt.plot(t, inf5, lw=3)
    plt.grid()
    plt.legend(("$\\alpha_V ="+str(tre[0])+"$","$\\alpha_V="+str(tre[1])+"$","$\\alpha_V ="+str(tre[2])+"$"), loc="best")
    plt.xlabel('Time(Days)')
    plt.ylabel('number of Death')
    plt.yticks()
plt.show()
