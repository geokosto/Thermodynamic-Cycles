# -*- coding: utf-8 -*-
from src.two_phase_fluid_cycle import TwoPhaseCycle
from CoolProp.CoolProp import PropsSI
from src.fluid import Fluid


def main():

    # general paramemeters
    fluid = Fluid("IF97::Water")
    ph = 30 * 1e5
    pm = 10 * 1e5
    pl = 1.5 * 1e5
    t_max = PropsSI(
        'T', 'P', ph, 'S', PropsSI(
            'S', 'P', pm, 'Q', 1, fluid.name
        ), fluid.name
    )
    s1 = PropsSI('S', 'P', pl, 'Q', 0, fluid.name)
    p1 = pl
    t1 = PropsSI('T', 'P', pl, 'Q', 0, fluid.name)
    x1 = 0
    h1 = PropsSI('H', 'P', pl, 'Q', 0, fluid.name)

    initial_state = {
        "pi": p1,
        "si": s1,
        "ti": t1,
        "xi": x1,
        "hi": h1,
    }

    # SIMPLE RANKINE CYCLE
    processes_rankine_simple = [
        {"prcs_name": '1--2', "prcs_type": "dq",
         "pf": ph, "eta_isentropic": 1.0},
        {"prcs_name": '2--3', "prcs_type": "dp",
         "xf": 0, "eta_isentropic": 0},
        {"prcs_name": '3--4', "prcs_type": "dp",
         "xf": 1, "eta_isentropic": 0},
        {"prcs_name": '4--5', "prcs_type": "dp",
         "tf": t_max, "eta_isentropic": 0},
        {"prcs_name": '5--6', "prcs_type": "dq",
         "pf": pl, "eta_isentropic": 0.9},
        {"prcs_name": '6--1', "prcs_type": "dp",
         "xf": 0, "eta_isentropic": 0}
    ]

    rankine_simple = TwoPhaseCycle(
        fluid,
        initial_state,
        processes_rankine_simple,
    )

    [print(process, "\n") for process in rankine_simple.solved_processes]
    rankine_simple.states_values_to_csv('rankine_simple')
    rankine_simple.processes_values_to_csv('rankine_simple')
    # simple.plot_T_s_diagram('simple_rankine')
    # simple.plot_h_s_diagram('simple_rankine')

    # REHEATED RANKINE CYCLE
    processes_rankine_reheated = [
        {"prcs_name": '1--2', "prcs_type": "dq",
         "pf": ph, "eta_isentropic": 1.0},
        {"prcs_name": '2--3', "prcs_type": "dp",
         "xf": 0, "eta_isentropic": 0},
        {"prcs_name": '3--4', "prcs_type": "dp",
         "xf": 1, "eta_isentropic": 0},
        {"prcs_name": '4--5', "prcs_type": "dp",
         "tf": t_max, "eta_isentropic": 0},
        {"prcs_name": '5--6', "prcs_type": "dq",
         "pf": pm, "eta_isentropic": 0.9},
        {"prcs_name": '6--7', "prcs_type": "dp",
         "tf": t_max, "eta_isentropic": 0},
        {"prcs_name": '7--8', "prcs_type": "dq",
         "pf": pl, "eta_isentropic": 0.87},
        {"prcs_name": '8--1', "prcs_type": "dp",
         "xf": 0, "eta_isentropic": 0}
    ]

    rankine_reheated = TwoPhaseCycle(
        fluid,
        initial_state,
        processes_rankine_reheated,
    )
    [print(process, "\n") for process in rankine_reheated.solved_processes]
    rankine_reheated.states_values_to_csv('rankine_reheated')
    rankine_reheated.processes_values_to_csv('rankine_reheated')
    # rankine_reheated.plot_T_s_diagram('reheated_rankine')
    # rankine_reheated.plot_h_s_diagram('reheated_rankine')

    # REGENERATIVE AND REHEATED RANKINE CYCLE
    processes_rankine_regenerative = [
        {"prcs_name": '1--a', "prcs_type": "dq",
         "pf": pm, "eta_isentropic": 1.0},
        {"prcs_name": 'a--b', "prcs_type": "dp",
         "xf": 0, "eta_isentropic": None},
        {"prcs_name": 'b--c', "prcs_type": "dq",
         "pf": ph, "eta_isentropic": 1.0},
        {"prcs_name": 'c--3', "prcs_type": "dp",
         "xf": 0, "eta_isentropic": None},
        {"prcs_name": '3--4', "prcs_type": "dp",
         "xf": 1, "eta_isentropic": None},
        {"prcs_name": '4--5', "prcs_type": "dp",
         "tf": t_max, "eta_isentropic": None},
        {"prcs_name": '5--6', "prcs_type": "dq",
         "pf": pm, "eta_isentropic": 0.9},
        {"prcs_name": '6--7', "prcs_type": "dp",
         "tf": t_max, "eta_isentropic": None},
        {"prcs_name": '7--8', "prcs_type": "dq",
         "pf": pl, "eta_isentropic": 0.87},
        {"prcs_name": '8--1', "prcs_type": "dp",
         "xf": 0, "eta_isentropic": None}
    ]
    rankine_regenerative = TwoPhaseCycle(
        fluid,
        initial_state,
        processes_rankine_regenerative,
    )
    [print(process, "\n") for process in rankine_regenerative.solved_processes]
    # m_extruded = rankine_regenerative.m_extruded()
    rankine_regenerative.states_values_to_csv('rankine_regenerative')
    rankine_regenerative.processes_values_to_csv('rankine_regenerative')
    # rankine_regenerative.plot_T_s_diagram('regenerative_rankine')
    # rankine_regenerative.plot_h_s_diagram('regenerative_rankine')


if __name__ == "__main__":
    main()
