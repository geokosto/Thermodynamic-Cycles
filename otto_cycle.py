from src.ideal_gas_cycle import IdealGasCycle
from src.fluid import Fluid


def main():
    # general parameters
    k, ra, cp = 1.4, 8.314/30, 1
    t_atm = 308
    p_atm = 1.0
    h_atm = cp*t_atm
    v1 = 6e-4
    t4 = 800
    epsilon = 9.5
    s1 = 0
    p1 = p_atm
    mass = (p_atm*v1*100)/(ra*t_atm)
    v1 = v1/mass

    ideal_gas_fluid = Fluid("air", ra, cp, k)

    initial_state = {
        "pi": p_atm,
        "si": s1,
        "ti": t_atm,
        "vi": v1,
        "hi": h_atm,
    }

    # otto cycle
    processes = [
        {"prcs_name": '1--2', "prcs_type": "ds", "vf": v1*(1/epsilon)},
        {"prcs_name": '2--3', "prcs_type": "dv", "tf": t4*(epsilon**(k-1))},
        {"prcs_name": '3--4', "prcs_type": "ds", "tf": t4},
        {"prcs_name": '4--1', "prcs_type": "dv", "pf": p1}
    ]

    otto = IdealGasCycle(ideal_gas_fluid, initial_state, processes)

    [print(_, "\n") for _ in otto.solved_processes]

    otto.states_values_to_csv('otto')
    otto.processes_values_to_csv('otto')
    otto.create_diagram("otto", ("v", "p"))


if __name__ == "__main__":
    main()
