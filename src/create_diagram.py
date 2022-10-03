from .diagram_factory import DiagramFactory
from .fluid import Fluid
import matplotlib.pyplot as plt


def create_diagram(
    name: str,
    diagram_type: tuple,
    fluid: Fluid,
    processes: dict,
    factory: DiagramFactory

):
    plt.figure()
    for idx, process in enumerate(processes):
        process_values = getattr(
            factory(fluid, process),
            process.get("prcs_type")
        )
        x_vals = process_values.get(diagram_type[0])
        y_vals = process_values.get(diagram_type[1])

        # plots state points
        plt.scatter(x_vals[0], y_vals[0], color='black')
        # plot state texts
        plt.text(
            x_vals[0]*1.0005, y_vals[0]*1.01,
            process.get("prcs_name")
        )
        # plot state texts
        plt.plot(x_vals, y_vals, 'b-')

    # plt.xlabel('specific Entropy'+' $s\ [\dfrac{kJ}{kgK}]$')
    # plt.ylabel('temperature'+' $t\ [K]$')
    plt.xlabel(diagram_type[0])
    plt.ylabel(diagram_type[1])
    plt.grid()
    plt.show()
    plt.savefig(
        f"{diagram_type[0]}-{diagram_type[1]}_diagram_{name}.png",
        dpi=300,
        bbox_inches='tight'
    )
