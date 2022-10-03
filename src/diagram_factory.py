from numpy import linspace, exp
# import matplotlib.pyplot as plt
# from CoolProp.CoolProp import PropsSI
from abc import ABC, abstractmethod
from .fluid import Fluid
# from scipy import optimize
# from .utils import vdw, vdwprime


class DiagramFactory(ABC):

    @abstractmethod
    def dp(self):
        pass

    @abstractmethod
    def dv(self):
        pass

    @abstractmethod
    def ds(self):
        pass

    @abstractmethod
    def dt(self):
        pass

    @abstractmethod
    def dq(self):
        pass


class IdealGasDiagram(DiagramFactory):
    def __init__(
        self,
        fluid: Fluid,
        process: dict,
    ) -> None:
        self.fluid = fluid
        self.process = process

    def dp(self):
        results = {
            "p": linspace(self.process["pi"], self.process["pf"]),
            "s": linspace(self.process["si"], self.process["sf"]),
            "t": self.process["ti"]*exp(
                (
                    linspace(
                        self.process["si"], self.process["sf"]
                    )-self.process["si"]
                )/self.fluid.cp
            )
        }
        return results

    def dv(self):
        results = {
            "p": linspace(self.process["pi"], self.process["pf"]),
            "v": linspace(self.process["vi"], self.process["vf"]),
            "s": linspace(self.process["si"], self.process["sf"]),
            "t": self.process["ti"] * exp(
                (
                    linspace(
                        self.process["si"], self.process["sf"]
                    )-self.process["si"]
                ) / (
                    self.fluid.cp/self.fluid.k
                )
            )
        }
        return results

    def ds(self):
        results = {
            "s": linspace(self.process["si"], self.process["sf"]),
            "t": linspace(self.process["ti"], self.process["tf"]),
            "v": linspace(self.process["vi"], self.process["vf"]),
            "p": self.process["pi"]*(
                linspace(
                    self.process["vi"], self.process["vf"]
                )/self.process["vi"]
            )**(-self.fluid.k)
        }
        return results

    def dt(self):
        raise ValueError(
            'dq process type for ideal gas diagram is not yet implemented!'
        )

    def dq(self):
        raise ValueError(
            'dq process type for ideal gas diagram is not yet implemented!'
        )


# plt.xlabel('specific volume'+' $v\ [\dfrac{m^{3}}{kg}]$')
# plt.ylabel('pressure'+' $p\ [bar]$')
