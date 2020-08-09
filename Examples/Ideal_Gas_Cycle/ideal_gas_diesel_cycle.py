# -*- coding: utf-8 -*-
from thermodynamics import ideal_gas_cycle

k,Ra,cp = 1.4,0.27714,1
T_atm = 340
P_atm = 0.95
V1 = 33e-4

P1 = P_atm
T1 = T_atm   
h1 = h_atm = cp*T_atm
mass = (P1*V1*100)/(Ra*T1)
s1 = 0
v1=V1/mass

epsilon = 17
phi = 2.362
r=1.5

processes = [{'1--2':'dS'},
             {'2--3':'dP'},
             {'3--4':'dS'},
             {'4--1':'dV'}]

final_vals = [{'Pf' :P_atm*(epsilon**k) }, #2
             {'vf' : v1*(phi/epsilon) }, #3
             {'vf' : v1}, #4
             {'Pf' : P_atm}, #1
             {'Pf' :P_atm*(epsilon**k) }] #2
         
             
diesel = ideal_gas_cycle(k,Ra,cp,s1,P1,T1,v1,h1,mass,processes,final_vals)

states = diesel.states_values()
diesel.states_values_to_csv('diesel')

processes = diesel.processes_values()
diesel.processes_values_to_csv('diesel')

h_list = diesel.h_values()
T_list = diesel.T_values()
q_list = diesel.q_values()
w_list = diesel.w_values()
wt_list = diesel.wt_values()

diesel.plot_T_s_diagram('diesel')
diesel.plot_P_v_diagram('diesel')
    



