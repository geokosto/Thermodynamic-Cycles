from src.vdw_gas_cycle import VdwGasCycle
from src.fluid import Fluid


def main():
    r = 8.314
    k = 1.3
    cv_0 = r/(k-1)
    tc = 419.7  # K
    pc = 40*1e2  # kPa
    zc = 3/8

    fluid = Fluid("air", r=r, cv_0=cv_0, tc=tc, pc=pc, zc=zc)

    p1 = p4 = 0.5*pc  # kPa
    p2 = 0.8*pc  # kPa
    t1 = 0.96*pc
    t3 = 3*tc
    s1 = h1 = u1 = 0
    v1 = 1.345

    initial_state = {
        "pi": p1,
        "si": s1,
        "ti": t1,
        "hi": h1,
        "ui": u1,
        "vi": v1
    }

    # Ericsson cycle
    processes = [
        {"prcs_name": '1--2', "prcs_type": "dt", "pf": p2},
        {"prcs_name": '2--3', "prcs_type": "dp", "tf": t3},
        {"prcs_name": '3--4', "prcs_type": "dt", "pf": p4},
        {"prcs_name": '4--1', "prcs_type": "dp", "tf": t1}
    ]

    ericsson = VdwGasCycle(fluid, initial_state, processes)

    [print(_, "\n") for _ in ericsson.solved_processes]

    ericsson.states_values_to_csv('ericsson')
    ericsson.processes_values_to_csv('ericsson')


if __name__ == "__main__":
    main()
