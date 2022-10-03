# from typing import Callable
from .fluid import Fluid
from .solve_process_factory import SolveProcessFactory
# from .utils import process_type_factory_dict


def solve(
    fluid: Fluid,
    initial_state: dict,
    processes: list[dict],
    factory: SolveProcessFactory,
    # eta_isentropic: list
) -> list[dict]:

    process_itter = initial_state.copy()
    solved_processes = []
    for idx, process in enumerate(processes):
        if idx > 0:
            process_itter = {
                f"{key[0]}i": process_itter[f"{key[0]}f"] for key in initial_state
            }
        process_itter.update(process)
        process_itter = getattr(
            factory(
                fluid, process_itter
            ),
            process.get("prcs_type")
        )()
        solved_processes.append(process_itter.copy())
    return solved_processes
