# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy import optimize

class two_phase_cycle:
    
    def __init__(self,fluid,eta_isentropic,s1,P1,T1,x1,h1,processes,next_state_values,mass=1,regenerative=False,processes_with_m_extruded_index=[],mixed_state=0):
        self.fluid = fluid
        self.eta_isentropic = eta_isentropic
        self.s1 = s1
        self.P1 = P1
        self.T1 = T1
        self.x1 = x1
        self.h1 = h1
        self.mass = mass
        self.processes = processes
        self.next_state_values = next_state_values
        self.regenerative = regenerative
        self.processes_with_m_extruded_index = processes_with_m_extruded_index
        self.mixed_state = mixed_state
        
    def quantities(self):
        def process(self,eta_isentropic,process_type,**kwargs):
          if 'Ti' in kwargs:
            Ti = kwargs['Ti']
          else:
            Ti = '-'
          if 'Tf' in kwargs:
            Tf = kwargs['Tf']
          else:
            Tf = '-'
          if 'Pi' in kwargs:
            Pi = kwargs['Pi']
          else:
            Pi = '-'
          if 'Pf' in kwargs:
            Pf = kwargs['Pf']
          else:
            Pf = '-' 
          if 'si' in kwargs:
            si = kwargs['si']
          else:
            si = '-' 
          if 'sf' in kwargs:
            sf = kwargs['sf']
          else:
            sf = '-' 
          if 'vi' in kwargs:
            vi = kwargs['vi']
          else:
            vi = '-' 
          if 'vf' in kwargs:
            vf = kwargs['vf']
          else:
            vf = '-' 
          if 'hi' in kwargs:
            hi = kwargs['hi']
          else:
            hi = '-' 
          if 'hf' in kwargs:
            hf = kwargs['hf']
          else:
            hf = '-'
          if 'xi' in kwargs:
            xi = kwargs['xi']
          else:
            xi = '-' 
          if 'xf' in kwargs:
            xf = kwargs['xf']
          else:
            xf = '-' 
          if process_type == 'dP':
              Pf=Pi
              if Tf =='-':
                  if xf!='-'and xf>=0 and xf<=1:
                      Tf = PropsSI('T','P',Pf,'Q',xf,self.fluid)
                  elif sf!='-' :
                      Tf = PropsSI('T','P',Pf,'S',sf,self.fluid)
              if hf=='-':
                  if xf!='-':
                      hf = PropsSI('H','P',Pf,'Q',xf,self.fluid)
                  else:
                      hf = PropsSI('H','P',Pf,'T',Tf,self.fluid)
              if sf=='-' :
                  if xf!='-':
                      sf = PropsSI('S','P',Pf,'Q',xf,self.fluid)
                  else:
                      sf = PropsSI('S','P',Pf,'T',Tf,self.fluid)
              
              if xf=='-':
                  xf = (hf - PropsSI('H','P',Pf,'Q',0,self.fluid))/(PropsSI('H','P',Pf,'Q',1,self.fluid) - PropsSI('H','P',Pf,'Q',0,self.fluid))
              wt = 0
              q=hf - hi
              w=-1
          elif process_type == 'dV':
              pass
          elif process_type == 'dT':
              pass
          elif process_type == 'dq':        
              sf_is = si
              #if hf=='-':
                #hf_is = PropsSI('H','P',Pf,'S',sf_is,self.fluid)
              if eta_isentropic==1 or eta_isentropic==1.0:
                  #hf = hf_is
                  hf = PropsSI('H','P',Pf,'S',sf_is,self.fluid)
                  sf = sf_is
              else:
                  x_is = (sf_is - PropsSI('S','P',Pf,'Q',0,self.fluid))/(PropsSI('S','P',Pf,'Q',1,self.fluid) - PropsSI('S','P',Pf,'Q',0,self.fluid))
                  hf_is = x_is*PropsSI('H','P',Pf,'Q',1,self.fluid)+(1-x_is)*PropsSI('H','P',Pf,'Q',0,self.fluid)
                  hf = hi+eta_isentropic*(hf_is - hi)
                  sf = PropsSI('S','P',Pf,'H',hf,self.fluid)
              if Tf=='-':
                  Tf = PropsSI('T','P',Pf,'S',sf,self.fluid)
              if xf=='-':
                  xf = (hf - PropsSI('H','P',Pf,'Q',0,self.fluid))/(PropsSI('H','P',Pf,'Q',1,self.fluid) - PropsSI('H','P',Pf,'Q',0,self.fluid))
              wt =hi - hf
              q=0
              w=-1
          else:
            return print('process_type type is not valid!')
          result = {'Pf':round(Pf,2),'Tf':round(Tf,2),'xf':xf,'sf':sf,'hf':hf,'q':q,'wt':wt,'w':w}
          return result
        j = 0
        state1 = {'si': self.s1,'Pi' : self.P1,'Ti' : self.T1  ,'xi':self.x1,'hi':self.h1}
        processes_values=[]
        state1.update(self.next_state_values[0])
        states_values = [[state1.get('si'),state1.get('Pi'),state1.get('Ti'),state1.get('xi'),state1.get('hi')]]
        state_i = state1
        states = range(1,len(self.processes)+1)
        for state in states :
            process_type = list(self.processes[state-1].values())[0]
            if process_type=='dq':
                state_i = process(self,self.eta_isentropic[j],process_type,**state_i)
                j+=1
            else:
                state_i = process(self,1,process_type,**state_i)
            processes_values.append([state_i.get('q'),state_i.get('wt'),state_i.get('w')])
            state_i = dict(zip(['Pi','Ti','xi','si','hi'],state_i.values()))
            state_i.update(self.next_state_values[state])
            states_values.append([state_i.get('si'),state_i.get('Pi'),state_i.get('Ti'),state_i.get('xi'),state_i.get('hi')])
        return states_values,processes_values  
    
    def states_values(self):
        return self.quantities()[0]
    
    def h_values(self):
        return  np.array(self.states_values()).T[4][0:-1]
    
    def T_values(self):
        return  np.array(self.states_values()).T[2][0:-1]
    
    def m_extruded(self):
        return (self.h_values()[2]-self.h_values()[1])/(self.h_values()[7]-self.h_values()[1])

    def processes_values(self):
        result = self.quantities()[1]
        if self.regenerative == False:
            return result
        else:
            for col in range(len(result[0])):
                result[self.mixed_state][col] = 0
            for row in self.processes_with_m_extruded_index:
                for col in range(len(result[0])):
                    result[row][col] = result[row][col]*(1.0-self.m_extruded())
            return result
                
    def q_values(self):
        return  np.array(self.processes_values()).T[0]
    
    def wt_values(self):
        return  np.array(self.processes_values()).T[1]
    
    def w_values(self):
        return  np.array(self.processes_values()).T[2]

    def states_values_to_csv(self,name):
        pr_keys=[list(i.keys())[0] for i in self.processes]
        states_values = self.states_values()[:-1]
        states_values = np.array(states_values).T  
        df= pd.DataFrame()
        df['States'] = [pr_keys[i][0] for i in range(len(pr_keys))]
        df['s [kJ/kgK]'] = [round(i*1e-3,3) for i in states_values[0]]   
        df['P [bar]'] = [round(i*1e-5,3) for i in states_values[1]]
        df['T [K]'] = [round(i,3) for i in states_values[2]]
        df['$\\theta\ [ ^\circ C ]$'] = [round(i-273.15,3) for i in states_values[2]]
        x=[]
        for i in states_values[3]:
            if i<0:
                x.append('{} : Subcooled Liquid'.format(round(i,3)))
            elif i>1:
                x.append('{} : Superheated Steam'.format(round(i,3)))
            elif i == 0 or i ==0.0:
                x.append('{} : Saturated Liquid'.format(round(i,3)))
            elif i == 1 or i ==1.0:
                x.append('{} : Saturated Steam'.format(round(i,3)))
            else:
                x.append('{} : 2 Phase'.format(round(i,3)))
        df['x [--]'] = x
        df['h [kJ]'] = [round(i*1e-3,3) for i in states_values[4]]       
        df.to_csv('states_results_'+name+'.csv', index=False)
        
        
    def processes_values_to_csv(self,name):
        pr_keys=[list(i.keys())[0] for i in self.processes]
        pr_type=[list(i.values())[0] for i in self.processes]
        for i,element in enumerate(pr_type):
            if element=='dq':
                pr_type[i]='Adiabatic'
            elif element=='dP':
                pr_type[i]='Isobaric'
            elif element=='dV':
                pr_type[i]='Isochoric'
            elif element=='dT':
                pr_type[i]='Isothermal'
        processes_values = np.array(self.processes_values()).T  
        df= pd.DataFrame()
        df['Process'] = pr_keys
        df['Type'] = pr_type
        df['q [kJ/kg]'] = ['{:.5e}'.format(i*1e-3) for i in processes_values[0]]    
        #df['Q [kJ]'] = ['{:.2e}'.format(i*mass) for i in processes_values[0]]    
        df['wt [kJ/kg]'] = ['{:.5e}'.format(i*1e-3) for i in processes_values[1]]
        #df['Wt [kJ]'] = ['{:.2e}'.format(i*mass) for i in processes_values[1]]
        df['w [kJ/kg]'] = ['{:.5e}'.format(i*1e-3) for i in processes_values[2]]
        #df['W [kJ]'] = ['{:.2e}'.format(i*mass) for i in processes_values[2]]
        df.to_csv('processes_results_'+name+'.csv', index=False)
        
    
    def plot_T_s_diagram(self,name):       
        pr_keys=[list(i.keys())[0] for i in self.processes]
        plt.figure()
        x,y = [],[]
        states_values = self.states_values()
        Pmin = np.min(np.array(states_values).T[1])
        Pmin = self.P1
        Pmax = PropsSI('Pcrit',self.fluid)
        for P in np.linspace(0.5*Pmin,Pmax,1000):
            x.append(PropsSI('S','P',P,'Q',0,self.fluid))
            y.append(PropsSI('T','P',P,'Q',0,self.fluid))  
        for P in np.linspace(Pmax,0.5*Pmin,1000):
            x.append(PropsSI('S','P',P,'Q',1,self.fluid))
            y.append(PropsSI('T','P',P,'Q',1,self.fluid))
        plt.plot(x,y,'k--',label='Saturation Curve')
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][0],states_values[i][2],color='black')
            plt.text(states_values[i][0]*1.0005,states_values[i][2]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0]=='dq':
                plt.plot( (states_values[i][0],states_values[i+1][0]),
                             (states_values[i][2],states_values[i+1][2]),'b-')
            elif list(self.processes[i].values())[0]=='dP' :
                x,y=[],[]
                h_range = np.linspace(states_values[i][4],states_values[i+1][4],1000,endpoint=True)
                for j in range(len(h_range)):
                    x.append(PropsSI('S','P',states_values[i][1],'H',h_range[j],self.fluid))
                    y.append(PropsSI('T','P',states_values[i][1],'H',h_range[j],self.fluid))
                plt.plot(x,y,'b-')
            elif list(self.processes[i].values())[0]=='dV':
                pass
        if self.regenerative==True:
            x,y=[],[]
            h_range = np.linspace(states_values[2][4],states_values[7][4],1000,endpoint=True)
            for j in range(len(h_range)):
                x.append(PropsSI('S','P',states_values[7][1],'H',h_range[j],self.fluid))
                y.append(PropsSI('T','P',states_values[7][1],'H',h_range[j],self.fluid))
            plt.plot(x,y,'b-')
        plt.xlabel('Specific Entropy'+' $s\ [\dfrac{J}{kgK}]$')
        plt.ylabel('Temperature'+' $T\ [K]$')
        plt.grid()
        plt.legend()
        plt.show() 
        plt.savefig('T-s_diagram_'+name+'.png', dpi=300,bbox_inches='tight')
           
    
    def plot_h_s_diagram(self,name):
        pr_keys=[list(i.keys())[0] for i in self.processes]
        plt.figure()
        x,y = [],[]
        states_values = self.states_values()
        Pmin = np.min(np.array(states_values).T[1])
        Pmax = PropsSI('Pcrit',self.fluid)
        for P in np.linspace(0.5*Pmin,Pmax,1000):
            x.append(PropsSI('S','P',P,'Q',0,self.fluid))
            y.append(PropsSI('H','P',P,'Q',0,self.fluid))  
        for P in np.linspace(Pmax,0.5*Pmin,1000):
            x.append(PropsSI('S','P',P,'Q',1,self.fluid))
            y.append(PropsSI('H','P',P,'Q',1,self.fluid))
        plt.plot(x,y,'k--',label='Saturation Curve')
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][0],states_values[i][4],color='black')
            plt.text(states_values[i][0]*1.0005,states_values[i][4]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0]=='dq':
                plt.plot( (states_values[i][0],states_values[i+1][0]),
                             (states_values[i][4],states_values[i+1][4]),'b-')
            elif list(self.processes[i].values())[0]=='dP' :
                x,y=[],[]
                h_range = np.linspace(states_values[i][4],states_values[i+1][4],1000,endpoint=True)
                for j in range(len(h_range)):
                    x.append(PropsSI('S','P',states_values[i][1],'H',h_range[j],self.fluid))
                plt.plot(x,h_range,'b-')
            elif list(self.processes[i].values())[0]=='dV':
                pass
        if self.regenerative==True:
            x,y=[],[]
            h_range = np.linspace(states_values[2][4],states_values[7][4],1000,endpoint=True)
            for j in range(len(h_range)):
                x.append(PropsSI('S','P',states_values[7][1],'H',h_range[j],self.fluid))
            plt.plot(x,h_range,'b-')
            
        plt.xlabel('Specific Entropy'+' $s\ [\dfrac{J}{kgK}]$')
        plt.ylabel('Specific Enthalpy'+' $h\ [\dfrac{J}{kg}]$')
        plt.legend()
        plt.grid()
        plt.show()
        plt.savefig('h-s_diagram_'+name+'.png', dpi=300,bbox_inches='tight')
        
    
class ideal_gas_cycle:
    def __init__(self,k,Ra,cp,s1,P1,T1,v1,h1,mass,processes,next_state_values):
        self.k = k
        self.cp = cp
        self.Ra = Ra
        self.s1 = s1
        self.P1 = P1
        self.T1 = T1
        self.v1 = v1
        self.h1 = h1
        self.mass = mass
        self.processes = processes
        self.next_state_values = next_state_values
    def quantities(self):
        def process(self,process_type,**kwargs):
          if 'Ti' in kwargs:
            Ti = kwargs['Ti']
          else:
            Ti = '-'
          if 'Tf' in kwargs:
            Tf = kwargs['Tf']
          else:
            Tf = '-'
          if 'Pi' in kwargs:
            Pi = kwargs['Pi']
          else:
            Pi = '-'
          if 'Pf' in kwargs:
            Pf = kwargs['Pf']
          else:
            Pf = '-' 
          if 'si' in kwargs:
            si = kwargs['si']
          else:
            si = '-' 
          if 'sf' in kwargs:
            sf = kwargs['sf']
          else:
            sf = '-' 
          if 'vi' in kwargs:
            vi = kwargs['vi']
          else:
            vi = '-' 
          if 'vf' in kwargs:
            vf = kwargs['vf']
          else:
            vf = '-' 
          if 'hi' in kwargs:
            hi = kwargs['hi']
          else:
            hi = '-' 
          if 'hf' in kwargs:
            hf = kwargs['hf']
          else:
            hf = '-' 
          if process_type == 'dP':
              ## Must know vf or Tf to find Tf or vf
              # pressure
              Pf=Pi
              # temperature
              if Tf =='-':
                if Ti!='-' and vi!='-' and vf!='-': 
                  Tf = Ti*(vf/vi)
              # volume
              if vf == '-':
                if Ti!='-' and Tf!='-' and vi!='-': 
                  vf = vi*(Tf/Ti)
              # enthalpy
              if Tf!='-':
                hf = self.cp*Tf
              # entropy
              if si!='-' and Ti!='-' and vi!='-' :
                if vf!='-':
                  sf = si + self.cp*np.log(vf/vi)
                elif Tf!='-':
                  sf = si + self.cp*np.log(Tf/Ti)
              wt = 0
              q=self.cp*(Tf-Ti)
              ###w=Pf*1e5*(vf-vi)
              w=self.Ra*(Tf-Ti)
          elif process_type == 'dV':
              ## Must know Pf or Tf to find Tf or Pf
              # volume
              vf=vi
              # temperature
              if Tf =='-':
                if Ti!='-' and Pi!='-' and Pf!='-': 
                  Tf = Ti*(Pf/Pi)
              # Pressure
              if Pf == '-':
                if Ti!='-' and Tf!='-' and Pi!='-': 
                  Pf = Pi*(Tf/Ti)
              # enthalpy
              if Tf!='-':
                hf = self.cp*Tf
              # entropy
              if si!='-' and Ti!='-' and Pi!='-' :
                if Pf!='-':
                  sf = si + (self.cp/self.k)*np.log(Pf/Pi)
                elif Tf!='-':
                  sf = si + (self.cp/self.k)*np.log(Tf/Ti)
              wt =-vf*(Pf-Pi)
              q=(self.cp/self.k)*(Tf-Ti)
              w=0
          elif process_type == 'dT':
              ##p1V1 = p2V2
              pass
        
          elif process_type == 'dS':
            ### Must have Tf or Pf to find Pf or Tf
              sf = si
              if Pf=='-':
                if Pi!='-'and Ti!='-'and Tf!='-':
                  Pf = Pi*(Tf/Ti)**(self.k/(self.k-1))
                elif Pi!='-' and vi!='-' and vf !='-':
                  Pf = Pi*(vi/vf)**(self.k)  
              if Tf=='-':
                if Pf !='-':
                  Tf = Ti*(Pf/Pi)**((self.k-1)/self.k)
              if vf=='-':
                if vi != '-' and Pi != '-' and Pf != '-':
                  vf = vi*(Pi/Pf)**(1/self.k)
              if Tf!='-':
                hf = self.cp*Tf
              wt =self.cp*(Ti-Tf)
              q=0
              w=(self.cp/self.k)*(Ti-Tf)
          else:
            return print('Process type is not valid!')
          result = {'Pf':Pf,'Tf':Tf,'vf':vf,'sf':sf,'hf':hf,'q':q,'wt':wt,'w':w}
          return result
        state1 = {'si': self.s1,'Pi' : self.P1,'Ti' : self.T1  ,'vi':self.v1,'hi':self.h1}
        processes_values=[]
        state1.update(self.next_state_values[0])
        states_values = [[state1.get('si'),state1.get('Pi'),state1.get('Ti'),state1.get('vi'),state1.get('hi')]]
        state_i = state1
        states = range(1,len(self.processes)+1)
        for state in states :
            process_type = list(self.processes[state-1].values())[0]
            state_i = process(self,process_type,**state_i)
            processes_values.append([state_i.get('q'),state_i.get('wt'),state_i.get('w')])
            state_i = dict(zip(['Pi','Ti','vi','si','hi'],state_i.values()))
            state_i.update(self.next_state_values[state])
            states_values.append([state_i.get('si'),state_i.get('Pi'),state_i.get('Ti'),state_i.get('vi'),state_i.get('hi')])
        return states_values,processes_values  
    
    def states_values(self):
        return self.quantities()[0]
    def processes_values(self):
        return self.quantities()[1]
    
    def h_values(self):
        return  np.array(self.states_values()).T[4][0:-1]
    
    def T_values(self):
        return  np.array(self.states_values()).T[2][0:-1]
    
    def q_values(self):
        return  np.array(self.processes_values()).T[0]
    
    def wt_values(self):
        return  np.array(self.processes_values()).T[1]
    
    def w_values(self):
        return  np.array(self.processes_values()).T[2]
    
    def states_values_to_csv(self,name):
        states_values = self.states_values()[:-1]
        states_values = np.array(states_values).T  
        df= pd.DataFrame()
        df['States'] = range(1,len(states_values.T)+1)
        df['s [kJ/kgK]'] = [round(i,2) for i in states_values[0]]   
        df['S [kJ/K]'] = ['{:.2e}'.format(i*self.mass) for i in states_values[0]]
        df['P [bar]'] = [round(i,2) for i in states_values[1]]
        df['T [K]'] = [round(i,2) for i in states_values[2]]
        df['v [m3/kg]'] = ['{:.2e}'.format(i) for i in states_values[3]]
        df['V [m3]'] = ['{:.2e}'.format(i*self.mass) for i in states_values[3]]
        df['h [kJ]'] = [round(i,2) for i in states_values[4]]
        df.to_csv('states_results_'+name+'.csv', index=False)     
        
    def processes_values_to_csv(self,name):
        pr_keys=[list(i.keys())[0] for i in self.processes]
        pr_type=[list(i.values())[0] for i in self.processes]
        for i,element in enumerate(pr_type):
            if element=='dS':
                pr_type[i]='Isentropic'
            elif element=='dP':
                pr_type[i]='Isobaric'
            elif element=='dV':
                pr_type[i]='Isochoric'
            elif element=='dT':
                pr_type[i]='Isothermal'
        processes_values = np.array(self.processes_values()).T  
        df= pd.DataFrame()
        df['Process'] = pr_keys
        df['Type'] = pr_type
        df['q [kJ/kg]'] = ['{:.2e}'.format(i) for i in processes_values[0]]    
        df['Q [kJ]'] = ['{:.2e}'.format(i*self.mass) for i in processes_values[0]]    
        df['wt [kJ/kg]'] = ['{:.2e}'.format(i) for i in processes_values[1]]
        df['Wt [kJ]'] = ['{:.2e}'.format(i*self.mass) for i in processes_values[1]]
        df['w [kJ/kg]'] = ['{:.2e}'.format(i) for i in processes_values[2]]
        df['W [kJ]'] = ['{:.2e}'.format(i*self.mass) for i in processes_values[2]]
        df.to_csv('processes_results_'+name+'.csv', index=False)
        
    
    def plot_T_s_diagram(self,name):
        isobaric_x = lambda Sa,Sb : np.linspace(Sa,Sb)
        isobaric_y = lambda Sa,Sb,Ta: Ta*np.exp((np.linspace(Sa,Sb)-Sa)/self.cp)
        isochoric_x = lambda Sa,Sb : np.linspace(Sa,Sb)
        isochoric_y = lambda Sa,Sb,Ta: Ta*np.exp((np.linspace(Sa,Sb)-Sa)/(self.cp/self.k))
        pr_keys=[list(i.keys())[0] for i in self.processes]
        states_values = self.states_values()
        plt.figure()
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][0],states_values[i][2],color='black')
            plt.text(states_values[i][0]*1.0005,states_values[i][2]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0]=='dS':
                plt.plot( (states_values[i][0],states_values[i+1][0]),
                             (states_values[i][2],states_values[i+1][2]),'b-')
            elif list(self.processes[i].values())[0]=='dP' :
                plt.plot( isobaric_x(states_values[i][0],states_values[i+1][0]),
                         isobaric_y(states_values[i][0],states_values[i+1][0],states_values[i][2]),'b-')
            elif list(self.processes[i].values())[0]=='dV':
                plt.plot( isochoric_x(states_values[i][0],states_values[i+1][0]),
                         isochoric_y(states_values[i][0],states_values[i+1][0],states_values[i][2]),'b-')
        plt.xlabel('Specific Entropy'+' $s\ [\dfrac{kJ}{kgK}]$')
        plt.ylabel('Temperature'+' $T\ [K]$')
        plt.grid()
        plt.show() 
        plt.savefig('T-s_diagram_'+name+'.png', dpi=300,bbox_inches='tight')
           
    def plot_P_v_diagram(self,name):
        isentropic_x = lambda va,vb : np.linspace(va,vb)
        isentropic_y = lambda va,vb,Pa: Pa*(np.linspace(va,vb)/va)**(-self.k)
        pr_keys=[list(i.keys())[0] for i in self.processes]
        states_values = self.states_values()
        plt.figure()
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][3],states_values[i][1],color='black')
            plt.text(states_values[i][3]*1.005,states_values[i][1]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0]=='dP' or list(self.processes[i].values())[0]=='dV':
                plt.plot( (states_values[i][3],states_values[i+1][3]),
                             (states_values[i][1],states_values[i+1][1]),'b-')
            elif list(self.processes[i].values())[0]=='dS'  :
                plt.plot( isentropic_x(states_values[i][3],states_values[i+1][3]),
                         isentropic_y(states_values[i][3],states_values[i+1][3],states_values[i][1]),'b-')
        plt.xlabel('Specific Volume'+' $v\ [\dfrac{m^{3}}{kg}]$')
        plt.ylabel('Pressure'+' $P\ [bar]$')
        plt.grid()
        plt.show()
        plt.savefig('P-v_diagram_'+name+'.png', dpi=300,bbox_inches='tight')

class VDW_fluid_cycle:

    def __init__(self,R,cv_0,Zc,Tc,Pc,P1,T1,s1,h1,u1,processes,next_state_values):
        self.R = R
        self.cv_0 = cv_0
        self.Zc = Zc
        self.Tc = Tc
        self.Pc = Pc
        self.s1 = s1
        self.h1 = h1
        self.u1 = u1
        self.P1 = P1
        self.T1 = T1
        self.processes = processes
        self.next_state_values = next_state_values
    
    def vc(self):
        return  self.Zc*(self.R*self.Tc/self.Pc)
    
    def a(self):
        return (9/8.)* self.R*self.Tc*self.vc()
       
    def b(self):
        return self.vc()/3.
     
    def quantities(self):
        VDW = lambda v,R,a,b,P,T: ((R*T)/(v-b)) - (a/(v**2)) - P
        VDWprime = lambda v,R,a,b,P,T :((2*a)/(v**3)) - ((R*T)/((v-b)**2))
        def process(self,process_type,**kwargs):
          if 'Ti' in kwargs:
            Ti = kwargs['Ti']
          else:
            Ti = '-'
          if 'Tf' in kwargs:
            Tf = kwargs['Tf']
          else:
            Tf = '-'
          if 'Pi' in kwargs:
            Pi = kwargs['Pi']
          else:
            Pi = '-'
          if 'Pf' in kwargs:
            Pf = kwargs['Pf']
          else:
            Pf = '-' 
          if 'si' in kwargs:
            si = kwargs['si']
          else:
            si = '-' 
          if 'sf' in kwargs:
            sf = kwargs['sf']
          else:
            sf = '-' 
          if 'vi' in kwargs:
            vi = kwargs['vi']
          else:
            vi = '-' 
          if 'vf' in kwargs:
            vf = kwargs['vf']
          else:
            vf = '-' 
          if 'hi' in kwargs:
            hi = kwargs['hi']
          else:
            hi = '-' 
          if 'hf' in kwargs:
            hf = kwargs['hf']
          else:
            hf = '-'
          if 'ui' in kwargs:
            ui = kwargs['ui']
          else:
            ui = '-' 
          if 'hf' in kwargs:
            uf = kwargs['uf']
          else:
            uf = '-'

          if process_type == 'dP':
              Pf=Pi
              sol = optimize.newton(VDW, x0=(self.R*Tf/(Pf)),fprime=VDWprime, args=(self.R,self.a(),self.b(),Pf,Tf),full_output=True)
              vf = sol[0]
              v_error = np.abs(VDW(vf,self.R,self.a(),self.b(),Pf,Tf))
              if v_error>1e-10:
                  print('Big Error in calculation of v (error>1e-10)')
                            
              sf = si + self.cv_0*np.log(Tf/Ti) + self.R*np.log((vf - self.b())/(vi - self.b()))
              uf = ui + self.cv_0*(Tf-Ti) - self.a()*((1/vf) - (1/vi))
              Delta_u = uf - ui
              Delta_s =sf - si
              hf = hi + Delta_u + (Pf*vf - Pi*vi )
              Delta_h = hf - hi
                            
              wt = 0        
              w=Pf*(vf-vi)
              q = Delta_u +w
              
          elif process_type == 'dV':
              pass
          elif process_type == 'dT':
              Tf=Ti
              sol = optimize.newton(VDW, x0=(self.R*Tf/(Pf)),fprime=VDWprime, args=(self.R,self.a(),self.b(),Pf,Tf),full_output=True)
              vf = sol[0]
              v_error = np.abs(VDW(vf,self.R,self.a(),self.b(),Pf,Tf))
              if v_error>1e-10:
                  print('Big Error in calculation of v (error>1e-10)')
              
              sf = si + self.R*np.log((vf - self.b())/(vi - self.b()))
              uf = ui - self.a()*((1/vf) - (1/vi))
              Delta_s = sf - si
              Delta_u = uf - ui
              hf = hi + Delta_u + (Pf*vf - Pi*vi )
              Delta_h = hf - hi
              
              w=self.R*Ti*np.log((vf - self.b())/(vi - self.b())) + self.a()*((1/vf) - (1/vi))
              q = Delta_u +w
              wt = q - Delta_h
          elif process_type == 'dq': 
              pass
          else:
            return print('Process type is not valid!')
          return {'Pf':Pf,'Tf':Tf,'vf':vf,'sf':sf,'hf':hf,'uf':uf,'q':q,'wt':wt,'w':w,'Delta_s':Delta_s,'Delta_h':Delta_h,'Delta_u':Delta_u}
        sol = optimize.newton(VDW, x0=(self.R*self.T1/(self.P1)),fprime=VDWprime, args=(self.R,self.a(),self.b(),self.P1,self.T1),full_output=True)
        v1 = sol[0]
        state1 = {'si': self.s1,'Pi' : self.P1,'Ti' : self.T1  ,'vi':v1,'ui':self.u1,'hi':self.h1}
        processes_values=[]
        state1.update(self.next_state_values[0])
        states_values = [[state1.get('si'),state1.get('Pi'),state1.get('Ti'),state1.get('vi'),state1.get('ui'),state1.get('hi')]]
        state_i = state1
        states = range(1,len(self.processes)+1)
        for state in states :
            process_type = list(self.processes[state-1].values())[0]
            state_i = process(self,process_type,**state_i)
            processes_values.append([state_i.get('q'),state_i.get('wt'),state_i.get('w'),state_i.get('Delta_s'),state_i.get('Delta_h'),state_i.get('Delta_u')])
            state_i = dict(zip(['Pi','Ti','vi','si','hi','ui'],state_i.values()))
            state_i.update(self.next_state_values[state])
            states_values.append([state_i.get('si'),state_i.get('Pi'),state_i.get('Ti'),state_i.get('vi'),state_i.get('ui'),state_i.get('hi')])
        return states_values,processes_values  
    
    def states_values(self):
        return self.quantities()[0]
    
    def h_values(self):
        return  np.array(self.states_values()).T[4][0:-1]
    
    def T_values(self):
        return  np.array(self.states_values()).T[2][0:-1]
    
    def processes_values(self):
        return self.quantities()[1]

    def q_values(self):
        return  np.array(self.processes_values()).T[0]
    
    def wt_values(self):
        return  np.array(self.processes_values()).T[1]
    
    def w_values(self):
        return  np.array(self.processes_values()).T[2]
    
    def states_values_to_csv(self,name):
        pr_keys=[list(i.keys())[0] for i in self.processes]
        states_values = self.states_values()[:-1]
        states_values = np.array(states_values).T  
        df= pd.DataFrame()
        df['State'] = [pr_keys[i][0] for i in range(len(pr_keys))]
        df['$s\ [kJ/kmolK]$'] = [round(i,3) for i in states_values[0]]   
        df['$P\ [bar]$'] = [round(i*1e-2,3) for i in states_values[1]]
        df['$T\ [K]$'] = [round(i,3) for i in states_values[2]]
        df['$v\ [m^ 3/kmol]$'] = [round(i,3) for i in states_values[3]]  
        df['$h\ [kJ/kmol]$'] = [round(i,3) for i in states_values[4]]
        df['$u\ [kJ/kmol]$'] = [round(i,3) for i in states_values[5]] 
        df.to_csv('states_results_'+name+'.csv', index=False)
                
    def processes_values_to_csv(self,name):
        pr_keys=[list(i.keys())[0] for i in self.processes]
        pr_type=[list(i.values())[0] for i in self.processes]
        for i,element in enumerate(pr_type):
            if element=='dq':
                pr_type[i]='Adiabatic'
            elif element=='dP':
                pr_type[i]='Isobaric'
            elif element=='dV':
                pr_type[i]='Isochoric'
            elif element=='dT':
                pr_type[i]='Isothermal'
        processes_values = np.array(self.processes_values()).T  
        df= pd.DataFrame()
        df['Process'] = pr_keys
        df['Type'] = pr_type
        df['$q\ [kJ/kmol]$'] = ['{:.3f}'.format(i) for i in processes_values[0]]    
        df['$w_t\ [kJ/kmol]$'] = ['{:.3f}'.format(i) for i in processes_values[1]]
        df['$w\ [kJ/kmol]$'] = ['{:.3f}'.format(i) for i in processes_values[2]]
        df['$\\Delta s\ [kJ/kmolK]$'] = ['{:.3f}'.format(i) for i in processes_values[3]]
        df['$\\Delta h\ [kJ/kmol]$'] = ['{:.3f}'.format(i) for i in processes_values[4]]
        df['$\\Delta u\ [kJ/kmol]$'] = ['{:.3f}'.format(i) for i in processes_values[5]]
        df.to_csv('processes_results_'+name+'.csv', index=False)
        
    def plot_T_s_diagram(self,name): 
        VDW = lambda T,R,a,b,P,v: ((R*T)/(v-b)) - (a/(v**2)) - P
        pr_keys=[list(i.keys())[0] for i in self.processes]
        plt.figure()
        states_values = self.states_values()
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][0],states_values[i][2],color='black')
            plt.text(states_values[i][0]*1.0005,states_values[i][2]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0]=='dT':
                plt.plot( (states_values[i][0],states_values[i+1][0]),
                             (states_values[i][2],states_values[i+1][2]),'b-')
            elif list(self.processes[i].values())[0]=='dP' :
                s = np.linspace(states_values[i][0],states_values[i+1][0],1000,endpoint=True)
                P = np.ones(len(s))*states_values[i][1]
                v = np.linspace(states_values[i][3],states_values[i+1][3],1000,endpoint=True)
                T=[]
                for i in range(len(s)):
                    sol = optimize.newton(VDW, x0=(P[i]*v[i]/(self.R)), args=(self.R,self.a(),self.b(),P[i],v[i]),full_output=True)
                    Ti = sol[0]
                    T_error = np.abs(VDW(Ti,self.R,self.a(),self.b(),P[i],v[i]))
                    if T_error>1e-10:
                        print('Big Error in calculation of v (error>1e-10)')
                    T.append(Ti)
                plt.plot(s,T,'b-')
            elif list(self.processes[i].values())[0]=='dV':
                pass
        plt.xlabel('Molecular Entropy'+' $s\ [\dfrac{kJ}{kmolK}]$')
        plt.ylabel('Temperature'+' $T\ [K]$')
        plt.grid()
        plt.show() 
        plt.savefig('T-s_diagram_'+name+'.png', dpi=300,bbox_inches='tight')
           
    
    def plot_P_v_diagram(self,name):
        VDW = lambda P,R,a,b,v,T: ((R*T)/(v-b)) - (a/(v**2)) - P
        pr_keys=[list(i.keys())[0] for i in self.processes]
        plt.figure()
        states_values = self.states_values()
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][3],states_values[i][1],color='black')
            plt.text(states_values[i][3]*1.0005,states_values[i][1]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0]=='dP':
                plt.plot( (states_values[i][3],states_values[i+1][3]),
                             (states_values[i][1],states_values[i+1][1]),'b-')
            elif list(self.processes[i].values())[0]=='dT' :
                s = np.linspace(states_values[i][0],states_values[i+1][0],1000,endpoint=True)
                T = np.ones(len(s))*states_values[i][2]
                v = np.linspace(states_values[i][3],states_values[i+1][3],1000,endpoint=True)
                P=[]
                for i in range(len(s)):
                    sol = optimize.newton(VDW, x0=(self.R*T[i]/(v[i])), args=(self.R,self.a(),self.b(),v[i],T[i]),full_output=True)
                    Pi = sol[0]
                    P_error = np.abs(VDW(Pi,self.R,self.a(),self.b(),v[i],T[i]))
                    if P_error>1e-10:
                        print('Big Error in calculation of P (error>1e-10)')
                    P.append(Pi)
                plt.plot(v,P,'b-')
            elif list(self.processes[i].values())[0]=='dV':
                pass          
        plt.xlabel('Molecular Volume'+' $v\ [\dfrac{m^3}{kmol}]$')
        plt.ylabel('Pressure'+' $P\ [kPa]$')
        plt.grid()
        plt.show()
        plt.savefig('h-s_diagram_'+name+'.png', dpi=300,bbox_inches='tight')
        