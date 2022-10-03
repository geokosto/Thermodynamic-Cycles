from .utils import prcs_type_map_dict, map_moisture
import pandas as pd
import pathlib
import matplotlib.pyplot as plt


class ThermoCycle:

    def processes_values(self, variable: str) -> list:
        return [process.get(f"{variable}") for process in self.solved_processes]

    def states_values_to_csv(self, name: str) -> pd.DataFrame:
        df = pd.DataFrame()
        df['states'] = [x[0] for x in self.processes_values("prcs_name")]
        df['s [kJ/kgK]'] = self.processes_values("si")
        df['p [bar]'] = self.processes_values("pi")
        df['t [K]'] = self.processes_values("ti")
        df['theta [ oC ]'] = [x-273.15 for x in df['t [K]']]
        df['v [m3/kg]'] = self.processes_values("vi")
        df['h [kJ/kg]'] = self.processes_values("hi")
        if self.processes_values("xi")[0] is not None:
            df['x [-]'] = self.processes_values("xi")
            df['fluid type'] = [map_moisture(x)
                                for x in self.processes_values("xi")]
        path = pathlib.Path(__file__).parent.resolve()
        df.to_csv(
            path.joinpath(
                f'..//data//states_results_{name}.csv'), index=False
        )
        return df

    def processes_values_to_csv(self, name: str) -> pd.DataFrame:
        df = pd.DataFrame()
        df['process'] = self.processes_values("prcs_name")
        df['type'] = [prcs_type_map_dict.get(
            x) for x in self.processes_values("prcs_type")]
        df['q [kJ/kg]'] = self.processes_values("q")
        df['wt [kJ/kg]'] = self.processes_values("wt")
        df['w [kJ/kg]'] = self.processes_values("w")
        path = pathlib.Path(__file__).parent.resolve()
        df.to_csv(
            path.joinpath(
                f'..//data//processes_results_{name}.csv'), index=False
        )
        return df

    def create_diagram(self, name: str, diagram_type: tuple):
        plt.figure()
        for process in self.solved_processes:
            process_values = getattr(
                self.diagram_factory(self.fluid, process),
                process.get("prcs_type")
            )()
            # print(process_values)
            x_vals = process_values.get(diagram_type[0])
            y_vals = process_values.get(diagram_type[1])

            # plots state points
            plt.scatter(x_vals[0], y_vals[0], color='black')
            # plot state texts
            plt.text(
                x_vals[0]*1.0005, y_vals[0]*1.01,
                process.get("prcs_name")[0]
            )
            # plot state texts
            plt.plot(x_vals, y_vals, 'b-')

        # plt.xlabel('specific Entropy'+' $s\ [\dfrac{kJ}{kgK}]$')
        # plt.ylabel('temperature'+' $t\ [K]$')
        plt.xlabel(diagram_type[0])
        plt.ylabel(diagram_type[1])
        plt.grid()
        # plt.show()
        path = pathlib.Path(__file__).parent.resolve()
        plt.savefig(
            path.joinpath(
                f"..//data//{diagram_type[0]}-{diagram_type[1]}_diagram_{name}.png"
            ),
            dpi=300,
            bbox_inches='tight'
        )
