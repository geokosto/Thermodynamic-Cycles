import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from .fluid import Fluid
from .thermocycle import ThermoCycle
from .solve_cycle import solve
from .solve_process_factory import TwoPhaseFluidSolveProcess


class TwoPhaseCycle(ThermoCycle):

    def __init__(
        self,
        fluid: Fluid,
        initial_state: dict,
        processes: list[dict],
        # eta_isentropic: list,
        # mass=1,
        # regenerative=False,
        # processes_with_m_extruded_index=[],
        # mixed_state=0
    ):
        self.fluid = fluid
        self.initial_state = initial_state
        self.processes = processes
        # self.eta_isentropic = eta_isentropic

        # self.mass = mass
        # self.regenerative = regenerative
        # self.processes_with_m_extruded_index = processes_with_m_extruded_index
        # self.mixed_state = mixed_state

        self.solved_processes: list[dict] = solve(
            self.fluid,
            self.initial_state,
            self.processes,
            TwoPhaseFluidSolveProcess,
            # self.eta_isentropic
        )

    # def m_extruded(self):
    #     return (self.h_values()[2]-self.h_values()[1])/(self.h_values()[7]-self.h_values()[1])

    def plot_t_s_diagram(self, name):
        pr_keys = [list(i.keys())[0] for i in self.processes]
        plt.figure()
        x, y = [], []
        states_values = self.states_values()
        Pmin = np.min(np.array(states_values).T[1])
        Pmin = self.P1
        Pmax = PropsSI('Pcrit', self.fluid)
        for P in np.linspace(0.5*Pmin, Pmax, 1000):
            x.append(PropsSI('S', 'P', P, 'Q', 0, self.fluid))
            y.append(PropsSI('T', 'P', P, 'Q', 0, self.fluid))
        for P in np.linspace(Pmax, 0.5*Pmin, 1000):
            x.append(PropsSI('S', 'P', P, 'Q', 1, self.fluid))
            y.append(PropsSI('T', 'P', P, 'Q', 1, self.fluid))
        plt.plot(x, y, 'k--', label='Saturation Curve')
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][0],
                        states_values[i][2], color='black')
            plt.text(states_values[i][0]*1.0005, states_values[i][2]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0] == 'dq':
                plt.plot((states_values[i][0], states_values[i+1][0]),
                         (states_values[i][2], states_values[i+1][2]), 'b-')
            elif list(self.processes[i].values())[0] == 'dP':
                x, y = [], []
                h_range = np.linspace(
                    states_values[i][4], states_values[i+1][4], 1000, endpoint=True)
                for j in range(len(h_range)):
                    x.append(
                        PropsSI('S', 'P', states_values[i][1], 'H', h_range[j], self.fluid))
                    y.append(
                        PropsSI('T', 'P', states_values[i][1], 'H', h_range[j], self.fluid))
                plt.plot(x, y, 'b-')
            elif list(self.processes[i].values())[0] == 'dV':
                pass
        if self.regenerative == True:
            x, y = [], []
            h_range = np.linspace(
                states_values[2][4], states_values[7][4], 1000, endpoint=True)
            for j in range(len(h_range)):
                x.append(
                    PropsSI('S', 'P', states_values[7][1], 'H', h_range[j], self.fluid))
                y.append(
                    PropsSI('T', 'P', states_values[7][1], 'H', h_range[j], self.fluid))
            plt.plot(x, y, 'b-')
        plt.xlabel('Specific Entropy'+' $s\ [\dfrac{J}{kgK}]$')
        plt.ylabel('Temperature'+' $T\ [K]$')
        plt.grid()
        plt.legend()
        plt.show()
        plt.savefig('T-s_diagram_'+name+'.png', dpi=300, bbox_inches='tight')

    def plot_h_s_diagram(self, name):
        pr_keys = [list(i.keys())[0] for i in self.processes]
        plt.figure()
        x, y = [], []
        states_values = self.states_values()
        Pmin = np.min(np.array(states_values).T[1])
        Pmax = PropsSI('Pcrit', self.fluid)
        for P in np.linspace(0.5*Pmin, Pmax, 1000):
            x.append(PropsSI('S', 'P', P, 'Q', 0, self.fluid))
            y.append(PropsSI('H', 'P', P, 'Q', 0, self.fluid))
        for P in np.linspace(Pmax, 0.5*Pmin, 1000):
            x.append(PropsSI('S', 'P', P, 'Q', 1, self.fluid))
            y.append(PropsSI('H', 'P', P, 'Q', 1, self.fluid))
        plt.plot(x, y, 'k--', label='Saturation Curve')
        for i in range(len(states_values)-1):
            plt.scatter(states_values[i][0],
                        states_values[i][4], color='black')
            plt.text(states_values[i][0]*1.0005, states_values[i][4]*1.01,
                     pr_keys[i][0])
            if list(self.processes[i].values())[0] == 'dq':
                plt.plot((states_values[i][0], states_values[i+1][0]),
                         (states_values[i][4], states_values[i+1][4]), 'b-')
            elif list(self.processes[i].values())[0] == 'dP':
                x, y = [], []
                h_range = np.linspace(
                    states_values[i][4], states_values[i+1][4], 1000, endpoint=True)
                for j in range(len(h_range)):
                    x.append(
                        PropsSI('S', 'P', states_values[i][1], 'H', h_range[j], self.fluid))
                plt.plot(x, h_range, 'b-')
            elif list(self.processes[i].values())[0] == 'dV':
                pass
        if self.regenerative == True:
            x, y = [], []
            h_range = np.linspace(
                states_values[2][4], states_values[7][4], 1000, endpoint=True)
            for j in range(len(h_range)):
                x.append(
                    PropsSI('S', 'P', states_values[7][1], 'H', h_range[j], self.fluid))
            plt.plot(x, h_range, 'b-')

        plt.xlabel('Specific Entropy'+' $s\ [\dfrac{J}{kgK}]$')
        plt.ylabel('Specific Enthalpy'+' $h\ [\dfrac{J}{kg}]$')
        plt.legend()
        plt.grid()
        plt.show()
        plt.savefig('h-s_diagram_'+name+'.png', dpi=300, bbox_inches='tight')
