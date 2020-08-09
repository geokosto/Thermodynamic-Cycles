# -*- coding: utf-8 -*-
from thermodynamics import two_phase_cycle
from CoolProp.CoolProp import PropsSI


fluid = 'IF97::Water'
PH = 30 *1e5
PM = 10 *1e5
PL = 1.5 * 1e5
Τ_max = PropsSI('T','P',PH,'S',PropsSI('S','P',PM,'Q',1,fluid),fluid) 
s1 = PropsSI('S','P',PL,'Q',0,fluid)
P1 = PL
T1 = PropsSI('T','P',PL,'Q',0,fluid) 
x1 = 0
h1 = PropsSI('H','P',PL,'Q',0,fluid) 

####    SIMPLE RANKINE CYCLE 
eta_isentropic = [1.0,0.9]
processes = [{'1--2':'dq'},
             {'2--3':'dP'},
             {'3--4':'dP'},
             {'4--5':'dP'},
             {'5--6':'dq'},
             {'6--1':'dP'}]
final_vals =[{'Pf' :PH }, #2
             {'xf' : 0 }, #3
             {'xf' : 1}, #4
             {'Tf' : Τ_max}, #5     
             {'Pf' : PL }, #6
             {'xf' :0 }, #1
             {'Pf' :PH }] #2

simple = two_phase_cycle(fluid,eta_isentropic,s1,P1,T1,x1,h1,processes,final_vals)

simple.states_values_to_csv('simple_rankine')
simple.processes_values_to_csv('simple_rankine')
simple.plot_T_s_diagram('simple_rankine')
simple.plot_h_s_diagram('simple_rankine') 


####    REHEATED RANKINE CYCLE 
eta_isentropic = [1.0,0.9,0.87]
processes = [{'1--2':'dq'},
             {'2--3':'dP'},
             {'3--4':'dP'},
             {'4--5':'dP'},
             {'5--6':'dq'},
             {'6--7':'dP'},
             {'7--8':'dq'},
             {'8--1':'dP'}]
final_vals =[{'Pf' :PH }, #2
             {'xf' : 0 }, #3
             {'xf' : 1}, #4
             {'Tf' : Τ_max}, #5 
             {'Pf' : PM }, #6
             {'Tf' : Τ_max }, #7
             {'Pf' : PL }, #8
             {'xf' :0 }, #1
             {'Pf' :PH }] #2

reheated = two_phase_cycle(fluid,eta_isentropic,s1,P1,T1,x1,h1,processes,final_vals)

reheated.states_values_to_csv('reheated_rankine')
reheated.processes_values_to_csv('reheated_rankine')
reheated.plot_T_s_diagram('reheated_rankine')
reheated.plot_h_s_diagram('reheated_rankine') 


####    REGENERATIVE RANKINE CYCLE 
eta_isentropic = [1,1,0.9,0.87]
processes = [{'1--a':'dq'},
             {'a--b':'dP'},
             {'b--c':'dq'},
             {'c--3':'dP'},
             {'3--4':'dP'},
             {'4--5':'dP'},
             {'5--6':'dq'},
             {'6--7':'dP'},
             {'7--8':'dq'},
             {'8--1':'dP'}]
final_vals =[{'Pf' : PM}, #a
             {'xf' : 0}, #b
             {'Pf' : PH}, #c
             {'xf' : 0}, #3 
             {'xf' : 1}, #4
             {'Tf' : Τ_max}, #5
             {'Pf' : PM}, #6
             {'Tf' : Τ_max}, #7
             {'Pf' : PL}, #8
             {'xf' : 0}, #1
             {'Pf' : PM}] #a
regenerative = two_phase_cycle(fluid,eta_isentropic,s1,P1,T1,x1,h1,processes,final_vals,
                               regenerative=True,processes_with_m_extruded_index=[0,7,8,9],mixed_state=1 )

m_extruded = regenerative.m_extruded()

regenerative.states_values_to_csv('regenerative_rankine')
regenerative.processes_values_to_csv('regenerative_rankine')
regenerative.plot_T_s_diagram('regenerative_rankine')
regenerative.plot_h_s_diagram('regenerative_rankine') 

