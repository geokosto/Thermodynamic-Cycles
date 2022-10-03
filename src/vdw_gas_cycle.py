# import numpy as np
from .fluid import Fluid
# import matplotlib.pyplot as plt
from .thermocycle import ThermoCycle
from .solve_cycle import solve
from .solve_process_factory import VdwGasSolveProcess


class VdwGasCycle(ThermoCycle):
    def __init__(
        self, fluid: Fluid,
        initial_state: dict,
        processes: list[dict],
    ):
        self.fluid = fluid
        self.initial_state = initial_state
        self.processes = processes
        self.solved_processes: list[dict] = solve(
            self.fluid,
            self.initial_state,
            self.processes,
            VdwGasSolveProcess
        )

    # def plot_t_s_diagram(self,name):
    #     isobaric_x = lambda sa,sb : np.linspace(sa,sb)
    #     isobaric_y = lambda sa,sb,ta: ta*np.exp((np.linspace(sa,sb)-sa)/self.cp)
    #     isochoric_x = lambda sa,sb : np.linspace(sa,sb)
    #     isochoric_y = lambda sa,sb,ta: ta*np.exp((np.linspace(sa,sb)-sa)/(self.cp/self.k))
    #     pr_keys=[list(i.keys())[0] for i in self.processes]
    #     states_values = self.states_values()
    #     plt.figure()
    #     for i in range(len(states_values)-1):
    #         plt.scatter(states_values[i][0],states_values[i][2],color='black')
    #         plt.text(states_values[i][0]*1.0005,states_values[i][2]*1.01,
    #                  pr_keys[i][0])
    #         if list(self.processes[i].values())[0]=='ds':
    #             plt.plot( (states_values[i][0],states_values[i+1][0]),
    #                          (states_values[i][2],states_values[i+1][2]),'b-')
    #         elif list(self.processes[i].values())[0]=='dp' :
    #             plt.plot( isobaric_x(states_values[i][0],states_values[i+1][0]),
    #                      isobaric_y(states_values[i][0],states_values[i+1][0],states_values[i][2]),'b-')
    #         elif list(self.processes[i].values())[0]=='dv':
    #             plt.plot( isochoric_x(states_values[i][0],states_values[i+1][0]),
    #                      isochoric_y(states_values[i][0],states_values[i+1][0],states_values[i][2]),'b-')
    #     plt.xlabel('specific Entropy'+' $s\ [\dfrac{kJ}{kgK}]$')
    #     plt.ylabel('temperature'+' $t\ [K]$')
    #     plt.grid()
    #     plt.show()
    #     plt.savefig('t-s_diagram_'+name+'.png', dpi=300,bbox_inches='tight')

    # def plot_p_v_diagram(self,name):
    #     isentropic_x = lambda va,vb : np.linspace(va,vb)
    #     isentropic_y = lambda va,vb,pa: pa*(np.linspace(va,vb)/va)**(-self.k)
    #     pr_keys=[list(i.keys())[0] for i in self.processes]
    #     states_values = self.states_values()
    #     plt.figure()
    #     for i in range(len(states_values)-1):
    #         plt.scatter(states_values[i][3],states_values[i][1],color='black')
    #         plt.text(states_values[i][3]*1.005,states_values[i][1]*1.01,
    #                  pr_keys[i][0])
    #         if list(self.processes[i].values())[0]=='dp' or list(self.processes[i].values())[0]=='dv':
    #             plt.plot( (states_values[i][3],states_values[i+1][3]),
    #                          (states_values[i][1],states_values[i+1][1]),'b-')
    #         elif list(self.processes[i].values())[0]=='ds'  :
    #             plt.plot( isentropic_x(states_values[i][3],states_values[i+1][3]),
    #                      isentropic_y(states_values[i][3],states_values[i+1][3],states_values[i][1]),'b-')
    #     plt.xlabel('specific volume'+' $v\ [\dfrac{m^{3}}{kg}]$')
    #     plt.ylabel('pressure'+' $p\ [bar]$')
    #     plt.grid()
    #     plt.show()
    #     plt.savefig('p-v_diagram_'+name+'.png', dpi=300,bbox_inches='tight')
