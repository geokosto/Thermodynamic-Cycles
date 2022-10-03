import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize


class VdwFluidCycle:

    def __init__(self, R, cv_0, Zc, Tc, Pc, P1, T1, s1, h1, u1, processes, next_state_values):
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
        return self.Zc*(self.R*self.Tc/self.Pc)

    def a(self):
        return (9/8.) * self.R*self.Tc*self.vc()

    def b(self):
        return self.vc()/3.

    def quantities(self):
        def VDW(v, R, a, b, P, T): return ((R*T)/(v-b)) - (a/(v**2)) - P

        def VDWprime(v, R, a, b, P, T): return (
            (2*a)/(v**3)) - ((R*T)/((v-b)**2))

        def process(self, process_type, **kwargs):
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
                Pf = Pi
                sol = optimize.newton(VDW, x0=(self.R*Tf/(Pf)), fprime=VDWprime,
                                      args=(self.R, self.a(), self.b(), Pf, Tf), full_output=True)
                vf = sol[0]
                v_error = np.abs(VDW(vf, self.R, self.a(), self.b(), Pf, Tf))
                if v_error > 1e-10:
                    print('Big Error in calculation of v (error>1e-10)')

                sf = si + self.cv_0 * \
                    np.log(Tf/Ti) + self.R * \
                    np.log((vf - self.b())/(vi - self.b()))
                uf = ui + self.cv_0*(Tf-Ti) - self.a()*((1/vf) - (1/vi))
                Delta_u = uf - ui
                Delta_s = sf - si
                hf = hi + Delta_u + (Pf*vf - Pi*vi)
                Delta_h = hf - hi

                wt = 0
                w = Pf*(vf-vi)
                q = Delta_u + w

            elif process_type == 'dV':
                pass
            elif process_type == 'dT':
                Tf = Ti
                sol = optimize.newton(VDW, x0=(self.R*Tf/(Pf)), fprime=VDWprime,
                                      args=(self.R, self.a(), self.b(), Pf, Tf), full_output=True)
                vf = sol[0]
                v_error = np.abs(VDW(vf, self.R, self.a(), self.b(), Pf, Tf))
                if v_error > 1e-10:
                    print('Big Error in calculation of v (error>1e-10)')

                sf = si + self.R*np.log((vf - self.b())/(vi - self.b()))
                uf = ui - self.a()*((1/vf) - (1/vi))
                Delta_s = sf - si
                Delta_u = uf - ui
                hf = hi + Delta_u + (Pf*vf - Pi*vi)
                Delta_h = hf - hi

                w = self.R*Ti * \
                    np.log((vf - self.b())/(vi - self.b())) + \
                    self.a()*((1/vf) - (1/vi))
                q = Delta_u + w
                wt = q - Delta_h
            elif process_type == 'dq':
                pass
            else:
                return print('Process type is not valid!')
            return {'Pf': Pf, 'Tf': Tf, 'vf': vf, 'sf': sf, 'hf': hf, 'uf': uf, 'q': q, 'wt': wt, 'w': w, 'Delta_s': Delta_s, 'Delta_h': Delta_h, 'Delta_u': Delta_u}
        sol = optimize.newton(VDW, x0=(self.R*self.T1/(self.P1)), fprime=VDWprime,
                              args=(self.R, self.a(), self.b(), self.P1, self.T1), full_output=True)
        v1 = sol[0]
        state1 = {'si': self.s1, 'Pi': self.P1, 'Ti': self.T1,
                  'vi': v1, 'ui': self.u1, 'hi': self.h1}
        processes_values = []
        state1.update(self.next_state_values[0])
        states_values = [[state1.get('si'), state1.get('Pi'), state1.get(
            'Ti'), state1.get('vi'), state1.get('ui'), state1.get('hi')]]
        state_i = state1
        states = range(1, len(self.processes)+1)
        for state in states:
            process_type = list(self.processes[state-1].values())[0]
            state_i = process(self, process_type, **state_i)
            processes_values.append([state_i.get('q'), state_i.get('wt'), state_i.get(
                'w'), state_i.get('Delta_s'), state_i.get('Delta_h'), state_i.get('Delta_u')])
            state_i = dict(
                zip(['Pi', 'Ti', 'vi', 'si', 'hi', 'ui'], state_i.values()))
            state_i.update(self.next_state_values[state])
            states_values.append([state_i.get('si'), state_i.get('Pi'), state_i.get(
                'Ti'), state_i.get('vi'), state_i.get('ui'), state_i.get('hi')])
        return states_values, processes_values

    def states_values(self):
        return self.quantities()[0]

    def h_values(self):
        return np.array(self.states_values()).T[4][0:-1]

    def T_values(self):
        return np.array(self.states_values()).T[2][0:-1]

    def processes_values(self):
        return self.quantities()[1]

    def q_values(self):
        return np.array(self.processes_values()).T[0]

    def wt_values(self):
        return np.array(self.processes_values()).T[1]

    def w_values(self):
        return np.array(self.processes_values()).T[2]

    def states_values_to_csv(self, name):
        pr_keys = [list(i.keys())[0] for i in self.processes]
        states_values = self.states_values()[:-1]
        states_values = np.array(states_values).T
        df = pd.DataFrame()
        df['State'] = [pr_keys[i][0] for i in range(len(pr_keys))]
        df['$s\ [kJ/kmolK]$'] = [round(i, 3) for i in states_values[0]]
        df['$P\ [bar]$'] = [round(i*1e-2, 3) for i in states_values[1]]
        df['$T\ [K]$'] = [round(i, 3) for i in states_values[2]]
        df['$v\ [m^ 3/kmol]$'] = [round(i, 3) for i in states_values[3]]
        df['$h\ [kJ/kmol]$'] = [round(i, 3) for i in states_values[4]]
        df['$u\ [kJ/kmol]$'] = [round(i, 3) for i in states_values[5]]
        df.to_csv('states_results_'+name+'.csv', index=False)

    def processes_values_to_csv(self, name):
        pr_keys = [list(i.keys())[0] for i in self.processes]
        pr_type = [list(i.values())[0] for i in self.processes]
        for i, element in enumerate(pr_type):
            if element == 'dq':
                pr_type[i] = 'Adiabatic'
            elif element == 'dP':
                pr_type[i] = 'Isobaric'
            elif element == 'dV':
                pr_type[i] = 'Isochoric'
            elif element == 'dT':
                pr_type[i] = 'Isothermal'
        processes_values = np.array(self.processes_values()).T
        df = pd.DataFrame()
        df['Process'] = pr_keys
        df['Type'] = pr_type
        df['$q\ [kJ/kmol]$'] = ['{:.3f}'.format(i)
                                for i in processes_values[0]]
        df['$w_t\ [kJ/kmol]$'] = ['{:.3f}'.format(i)
                                  for i in processes_values[1]]
        df['$w\ [kJ/kmol]$'] = ['{:.3f}'.format(i)
                                for i in processes_values[2]]
        df['$\\Delta s\ [kJ/kmolK]$'] = ['{:.3f}'.format(i)
                                         for i in processes_values[3]]
        df['$\\Delta h\ [kJ/kmol]$'] = ['{:.3f}'.format(i)
                                        for i in processes_values[4]]
        df['$\\Delta u\ [kJ/kmol]$'] = ['{:.3f}'.format(i)
                                        for i in processes_values[5]]
        df.to_csv('processes_results_'+name+'.csv', index=False)

    def plot_T_s_diagram(self, name):
        def VDW(T, R, a, b, P, v): return ((R*T)/(v-b)) - (a/(v**2)) - P
        pr_keys = [list(i.keys())[0] for i in self.processes]
        plt.figure()
        states_values = self.states_values()
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][0],
                        states_values[i][2], color='black')
            plt.text(states_values[i][0]*1.0005, states_values[i][2]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0] == 'dT':
                plt.plot((states_values[i][0], states_values[i+1][0]),
                         (states_values[i][2], states_values[i+1][2]), 'b-')
            elif list(self.processes[i].values())[0] == 'dP':
                s = np.linspace(
                    states_values[i][0], states_values[i+1][0], 1000, endpoint=True)
                P = np.ones(len(s))*states_values[i][1]
                v = np.linspace(
                    states_values[i][3], states_values[i+1][3], 1000, endpoint=True)
                T = []
                for i in range(len(s)):
                    sol = optimize.newton(VDW, x0=(
                        P[i]*v[i]/(self.R)), args=(self.R, self.a(), self.b(), P[i], v[i]), full_output=True)
                    Ti = sol[0]
                    T_error = np.abs(
                        VDW(Ti, self.R, self.a(), self.b(), P[i], v[i]))
                    if T_error > 1e-10:
                        print('Big Error in calculation of v (error>1e-10)')
                    T.append(Ti)
                plt.plot(s, T, 'b-')
            elif list(self.processes[i].values())[0] == 'dV':
                pass
        plt.xlabel('Molecular Entropy'+' $s\ [\dfrac{kJ}{kmolK}]$')
        plt.ylabel('Temperature'+' $T\ [K]$')
        plt.grid()
        plt.show()
        plt.savefig('T-s_diagram_'+name+'.png', dpi=300, bbox_inches='tight')

    def plot_P_v_diagram(self, name):
        def VDW(P, R, a, b, v, T): return ((R*T)/(v-b)) - (a/(v**2)) - P
        pr_keys = [list(i.keys())[0] for i in self.processes]
        plt.figure()
        states_values = self.states_values()
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][3],
                        states_values[i][1], color='black')
            plt.text(states_values[i][3]*1.0005, states_values[i][1]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0] == 'dP':
                plt.plot((states_values[i][3], states_values[i+1][3]),
                         (states_values[i][1], states_values[i+1][1]), 'b-')
            elif list(self.processes[i].values())[0] == 'dT':
                s = np.linspace(
                    states_values[i][0], states_values[i+1][0], 1000, endpoint=True)
                T = np.ones(len(s))*states_values[i][2]
                v = np.linspace(
                    states_values[i][3], states_values[i+1][3], 1000, endpoint=True)
                P = []
                for i in range(len(s)):
                    sol = optimize.newton(VDW, x0=(
                        self.R*T[i]/(v[i])), args=(self.R, self.a(), self.b(), v[i], T[i]), full_output=True)
                    Pi = sol[0]
                    P_error = np.abs(
                        VDW(Pi, self.R, self.a(), self.b(), v[i], T[i]))
                    if P_error > 1e-10:
                        print('Big Error in calculation of P (error>1e-10)')
                    P.append(Pi)
                plt.plot(v, P, 'b-')
            elif list(self.processes[i].values())[0] == 'dV':
                pass
        plt.xlabel('Molecular Volume'+' $v\ [\dfrac{m^3}{kmol}]$')
        plt.ylabel('Pressure'+' $P\ [kPa]$')
        plt.grid()
        plt.show()
        plt.savefig('h-s_diagram_'+name+'.png', dpi=300, bbox_inches='tight')
