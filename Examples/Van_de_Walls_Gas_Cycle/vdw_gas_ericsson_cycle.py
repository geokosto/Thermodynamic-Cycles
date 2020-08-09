# -*- coding: utf-8 -*-
from thermodynamics import VDW_fluid_cycle

R = 8.314
k = 1.3 
cv_0 = R/(k-1)
Tc = 419.7 # K
Pc = 40*1e2 #kPa
Zc = 3/8

P1 = P4 = 0.5*Pc # kPa
P2 = P3 = 0.8*Pc # kPa
T1 = T2 = 0.96*Tc
T3 = T4 = 3*Tc
s1 = h1 = u1 = 0

processes = [{'1--2':'dT'},
             {'2--3':'dP'},
             {'3--4':'dT'},
             {'4--1':'dP'}]

final_vals =[{'Pf' : P2}, #2
             {'Tf' : T3}, #3
             {'Pf' : P4}, #4
             {'Tf' :T1}, #1
             {'Pf' :P2}] #2


ericsson = VDW_fluid_cycle(R,cv_0,Zc,Tc,Pc,P1,T1,s1,h1,u1,processes,final_vals)

a = ericsson.a()
b = ericsson.b()
vc = ericsson.vc()
h_list = ericsson.h_values()
T_list = ericsson.T_values()
q_list = ericsson.q_values()
w_list = ericsson.w_values()
wt_list = ericsson.wt_values()

states = ericsson.states_values()
ericsson.states_values_to_csv('ericsson')
processes = ericsson.processes_values()
ericsson.processes_values_to_csv('ericsson')

ericsson.plot_T_s_diagram('ericsson')
ericsson.plot_P_v_diagram('ericsson')



