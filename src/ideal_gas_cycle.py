from .fluid import Fluid
from .thermocycle import ThermoCycle
from .solve_cycle import solve
from .solve_process_factory import IdealGasSolveProcess
from .diagram_factory import DiagramFactory, IdealGasDiagram


class IdealGasCycle(ThermoCycle):
    def __init__(
        self, fluid: Fluid,
        initial_state: dict,
        processes: list[dict],
    ):
        self.fluid = fluid
        self.initial_state = initial_state
        self.processes = processes
        self.diagram_factory: DiagramFactory = IdealGasDiagram
        self.solved_processes: list[dict] = solve(
            self.fluid,
            self.initial_state,
            self.processes,
            IdealGasSolveProcess
        )
