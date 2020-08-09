# -*- coding: utf-8 -*-
from thermodynamics import ideal_gas_cycle

k,Ra,cp = 1.4,8.314/30,1
T_atm = 308
P_atm = 1.0
h_atm = cp*T_atm
V1 = 6e-4
T4 = 800

epsilon = 9.5

s1 = 0
P1 = P_atm
T1 = h1 = T_atm 
mass = (P1*V1*100)/(Ra*T1)
v1=V1/mass

processes = [{'1--2':'dS'},
             {'2--3':'dV'},
             {'3--4':'dS'},
             {'4--1':'dV'}]

final_vals =[{'vf' :v1*(1/epsilon) }, #2
             {'Tf' : T4*(epsilon**(k-1)) }, #3
             {'Tf' : T4}, #4
             {'Pf' : P1}, #1
             {'vf' :v1*(1/epsilon) }] #2
         

otto = ideal_gas_cycle(k,Ra,cp,s1,P1,T1,v1,h1,mass,processes,final_vals)

states = otto.states_values()
otto.states_values_to_csv('otto')

processes = otto.processes_values()
otto.processes_values_to_csv('otto')

h_list = otto.h_values()
T_list = otto.T_values()
q_list = otto.q_values()
w_list = otto.w_values()
wt_list = otto.wt_values()


otto.plot_T_s_diagram('otto')
otto.plot_P_v_diagram('otto')


