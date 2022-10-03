from numpy import log, abs
from CoolProp.CoolProp import PropsSI
from abc import ABC, abstractmethod
from .fluid import Fluid
from scipy import optimize
from .utils import vdw, vdwprime


class SolveProcessFactory(ABC):

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


class IdealGasSolveProcess(SolveProcessFactory):
    def __init__(
        self,
        fluid: Fluid,
        process: dict,
    ) -> None:
        self.fluid = fluid
        self.process = process

    def dp(self) -> dict:
        if "pi" in self.process:
            self.process["pf"] = self.process["pi"]
        if "tf" not in self.process:
            if "ti" in self.process and "vi" in self.process and "vf" in self.process:
                self.process["tf"] = self.process["ti"] * (
                    self.process["vf"]/self.process["vi"]
                )
        if "vf" not in self.process:
            if "ti" in self.process and "tf" in self.process and "vi" in self.process:
                self.process["vf"] = self.process["vi"] * (
                    self.process["tf"]/self.process["ti"]
                )
        if "tf" in self.process:
            self.process["hf"] = self.fluid.cp*self.process["tf"]
        if "si" in self.process and "ti" in self.process and "vi" in self.process:
            if "vf" in self.process:
                self.process["sf"] = self.process["si"] + self.fluid.cp * log(
                    self.process["vf"]/self.process["vi"]
                )
            elif "tf" in self.process:
                self.process["sf"] = self.process["si"] + self.fluid.cp * log(
                    self.process["tf"]/self.process["ti"]
                )
        # self.process["xf"] = self.process["xi"]

        self.process["wt"] = 0
        self.process["q"] = self.fluid.cp * (
            self.process["tf"]-self.process["ti"])
        self.process["w"] = self.fluid.ra * (
            self.process["tf"]-self.process["ti"]
        )
        return self.process

    def dv(self) -> dict:
        # if "vi" in self.process:
        self.process["vf"] = self.process["vi"]
        if "tf" not in self.process:
            if "ti" in self.process and "pi" in self.process and "pf" in self.process:
                self.process["tf"] = self.process["ti"] * (
                    self.process["pf"]/self.process["pi"]
                )
        if "pf" not in self.process:
            if "ti" in self.process and "tf" in self.process and "pi" in self.process:
                self.process["pf"] = self.process["pi"] * (
                    self.process["tf"]/self.process["ti"]
                )
        if "tf" in self.process:
            self.process["hf"] = self.fluid.cp*self.process["tf"]
        if "si" in self.process and "ti" in self.process and "pi" in self.process:
            if "pf" in self.process:
                self.process["sf"] = self.process["si"] + (
                    self.fluid.cp/self.fluid.k
                )*log(
                    self.process["pf"]/self.process["pi"]
                )
            elif "tf" in self.process:
                self.process["sf"] = self.process["si"] + (
                    self.fluid.cp/self.fluid.k
                )*log(
                    self.process["tf"]/self.process["ti"]
                )
        # self.process["xf"] = self.process["xi"]

        self.process["wt"] = -self.process["vf"] * (
            self.process["pf"]-self.process["pi"]
        )
        self.process["q"] = (self.fluid.cp/self.fluid.k) * (
            self.process["tf"]-self.process["ti"]
        )
        self.process["w"] = 0
        return self.process

    def ds(self) -> dict:
        self.process["sf"] = self.process["si"]
        if "pf" not in self.process:
            if "pi" in self.process and "ti" in self.process and "tf" in self.process:
                self.process["pf"] = self.process["pi"] * (
                    self.process["tf"]/self.process["ti"]
                )**(
                    self.fluid.k/(self.fluid.k-1)
                )
            elif "pi" in self.process and "vi" in self.process and "vf" in self.process:
                self.process["pf"] = self.process["pi"]*(
                    self.process["vi"]/self.process["vf"]
                )**(
                    self.fluid.k
                )
        if "tf" not in self.process:
            if "pf" in self.process:
                self.process["tf"] = self.process["ti"] * (
                    self.process["pf"]/self.process["pi"]
                )**(
                    (self.fluid.k-1)/self.fluid.k
                )
        if "vf" not in self.process:
            if "vi" in self.process and "pi" in self.process and "pf" in self.process:
                self.process["vf"] = self.process["vi"] * (
                    self.process["pi"]/self.process["pf"]
                )**(
                    1/self.fluid.k
                )
        if "tf" in self.process:
            self.process["hf"] = self.fluid.cp*self.process["tf"]
        # self.process["xf"] = self.process["xi"]

        self.process["wt"] = self.fluid.cp * (
            self.process["ti"]-self.process["tf"]
        )
        self.process["q"] = 0
        self.process["w"] = (self.fluid.cp/self.fluid.k) * (
            self.process["ti"]-self.process["tf"]
        )
        return self.process

    def dt(self):
        raise ValueError(
            'dt process type for ideal gas is not yet implemented!'
        )

    def dq(self):
        raise ValueError(
            'dq process type for ideal gas is not yet implemented!'
        )


class VdwGasSolveProcess(SolveProcessFactory):
    def __init__(
        self,
        fluid: Fluid,
        process: dict,
    ) -> None:
        self.fluid = fluid
        self.process = process
        self.vc = fluid.zc*(fluid.r*fluid.tc/fluid.pc)
        self.a = (9/8.) * fluid.r*fluid.tc*self.vc
        self.b = self.vc/3.

    def _vdw(self, v, p, t):
        # return ((self.fluid.r*t)/(v-self.b)) - (self.a/(v**2)) - p
        return vdw(v, self.fluid.r, self.a, self.b, p, t)

    def _vdwprime(self, v, p, t):
        return vdwprime(v, self.fluid.r, self.a, self.b, p, t)

    def dp(self) -> dict:
        self.process["pf"] = self.process["pi"]
        sol = optimize.newton(
            self._vdw, x0=(
                self.fluid.r*self.process["tf"]/(self.process["pf"])
            ),
            fprime=self._vdwprime,
            args=(
                self.process["pf"], self.process["tf"]
            ), full_output=True
        )
        self.process["vf"] = sol[0]
        v_error = abs(
            self._vdw(
                sol[0], self.process["pf"], self.process["tf"]
            )
        )
        if v_error > 1e-10:
            raise ValueError('Big Error in calculation of v (error>1e-10)')
        self.process["sf"] = self.process["si"] + self.fluid.cv_0 * log(
            self.process["tf"]/self.process["ti"]
        ) + self.fluid.r * log(
            (self.process["vf"] - self.b)/(self.process["vi"] - self.b)
        )
        self.process["uf"] = self.process["ui"] + self.fluid.cv_0*(
            self.process["tf"]-self.process["ti"]
        ) - self.a*(
            (1/self.process["vf"]) - (1/self.process["vi"])
        )
        self.process["hf"] = self.process["hi"] + (
            self.process["uf"] - self.process["ui"]
        ) + (
            self.process["pf"]*self.process["vf"] -
            self.process["pi"]*self.process["vi"]
        )

        self.process["wt"] = 0
        self.process["w"] = self.process["pf"]*(
            self.process["vf"]-self.process["vi"]
        )
        self.process["q"] = self.process["uf"] - \
            self.process["ui"] + self.process["w"]
        return self.process

    def dv(self) -> dict:
        raise ValueError(
            'dv process type for van der waals gas is not yet implemented!'
        )

    def ds(self) -> dict:
        raise ValueError(
            'ds process type for van der waals gas is not yet implemented!'
        )

    def dt(self):
        self.process["tf"] = self.process["ti"]
        sol = optimize.newton(
            self._vdw, x0=(
                self.fluid.r*self.process["tf"]/(self.process["pf"])
            ), fprime=self._vdwprime,
            args=(
                self.process["pf"], self.process["tf"]
            ), full_output=True)
        self.process["vf"] = sol[0]
        v_error = abs(
            self._vdw(sol[0], self.process["pf"], self.process["tf"]))
        if v_error > 1e-10:
            raise ValueError('Big Error in calculation of v (error>1e-10)')
        self.process["sf"] = self.process["si"] + self.fluid.r*log(
            (self.process["vf"] - self.b)/(self.process["vi"] - self.b)
        )
        self.process["uf"] = self.process["ui"] - self.a*(
            (1/self.process["vf"]) - (1/self.process["vi"])
        )
        self.process["hf"] = self.process["hi"] + (
            self.process["uf"] - self.process["ui"]
        ) + (
            self.process["pf"]*self.process["vf"] -
            self.process["pi"]*self.process["vi"]
        )

        self.process["w"] = self.fluid.r*self.process["ti"] * log(
            (self.process["vf"] - self.b)/(self.process["vi"] - self.b)
        ) + self.a*(
            (1/self.process["vf"]) - (1/self.process["vi"])
        )
        self.process["q"] = (self.process["uf"] -
                             self.process["ui"]) + self.process["w"]
        self.process["wt"] = self.process["q"] - (
            self.process["hf"] - self.process["hi"]
        )
        return self.process

    def dq(self):
        raise ValueError(
            'dq process type for van der waals gas is not yet implemented!'
        )


class TwoPhaseFluidSolveProcess(SolveProcessFactory):

    def __init__(
        self,
        fluid: Fluid,
        process: dict,
    ) -> None:
        self.fluid = fluid
        self.process = process
        self.eta_is = process.get("eta_isentropic")

    def dp(self) -> dict:
        self.process["pf"] = self.process["pi"]
        if "tf" not in self.process:
            if "xf" in self.process and self.process["xf"] >= 0 and self.process["xf"] <= 1:
                self.process["tf"] = PropsSI(
                    'T', 'P', self.process["pf"],
                    'Q', self.process["xf"], self.fluid.name
                )
            elif "sf" in self.process:
                self.process["tf"] = PropsSI(
                    'T', 'P', self.process["pf"],
                    'S', self.process["sf"], self.fluid.name
                )
        if "hf" not in self.process:
            if "xf" in self.process:
                self.process["hf"] = PropsSI(
                    'H', 'P', self.process["pf"],
                    'Q', self.process["xf"], self.fluid.name
                )
            else:
                self.process["hf"] = PropsSI(
                    'H', 'P', self.process["pf"],
                    'T', self.process["tf"], self.fluid.name
                )
        if "sf" not in self.process:
            if "xf" in self.process:
                self.process["sf"] = PropsSI(
                    'S', 'P', self.process["pf"],
                    'Q', self.process["xf"], self.fluid.name
                )
            else:
                self.process["sf"] = PropsSI(
                    'S', 'P', self.process["pf"],
                    'T', self.process["tf"], self.fluid.name
                )
        if "xf" not in self.process:
            self.process["xf"] = (
                self.process["hf"] - PropsSI(
                    'H', 'P', self.process["pf"], 'Q', 0, self.fluid.name
                )
            )/(
                PropsSI(
                    'H', 'P', self.process["pf"], 'Q', 1, self.fluid.name
                ) - PropsSI(
                    'H', 'P', self.process["pf"], 'Q', 0, self.fluid.name
                )
            )
        # self.process["vf"] = self.process["vi"]

        self.process["wt"] = 0
        self.process["q"] = self.process["hf"] - self.process["hi"]
        self.process["w"] = -1
        return self.process

    def dv(self) -> dict:
        raise ValueError(
            'dv process type for two phase fluid is not yet implemented!'
        )

    def ds(self) -> dict:
        raise ValueError(
            'ds process type for two phase fluid is not yet implemented!'
        )

    def dt(self) -> dict:
        raise ValueError(
            'dt process type for two phase fluid is not yet implemented!'
        )

    def dq(self) -> dict:
        self.process["sf_is"] = self.process["si"]
        if self.eta_is > 0.999:
            self.process["hf"] = PropsSI(
                'H', 'P', self.process["pf"],
                'S', self.process["sf_is"], self.fluid.name
            )
            self.process["sf"] = self.process["sf_is"]
        else:
            self.process["x_is"] = (
                self.process["sf_is"] - PropsSI(
                    'S', 'P', self.process["pf"], 'Q', 0, self.fluid.name
                )
            )/(
                PropsSI(
                    'S', 'P', self.process["pf"], 'Q', 1, self.fluid.name
                ) - PropsSI(
                    'S', 'P', self.process["pf"], 'Q', 0, self.fluid.name
                )
            )
            self.process["hf_is"] = self.process["x_is"]*PropsSI(
                'H', 'P', self.process["pf"], 'Q', 1, self.fluid.name
            )+(
                1-self.process["x_is"]
            )*PropsSI(
                'H', 'P', self.process["pf"], 'Q', 0, self.fluid.name
            )
            self.process["hf"] = self.process["hi"]+self.eta_is * (
                self.process["hf_is"] - self.process["hi"]
            )
            self.process["sf"] = PropsSI(
                'S', 'P', self.process["pf"],
                'H', self.process["hf"], self.fluid.name
            )
        if "tf" not in self.process:
            self.process["tf"] = PropsSI(
                'T', 'P', self.process["pf"],
                'S', self.process["sf"], self.fluid.name
            )
        if "xf" not in self.process:
            self.process["xf"] = (
                self.process["hf"] - PropsSI(
                    'H', 'P', self.process["pf"], 'Q', 0, self.fluid.name
                )
            )/(
                PropsSI(
                    'H', 'P', self.process["pf"], 'Q', 1, self.fluid.name
                ) - PropsSI(
                    'H', 'P', self.process["pf"], 'Q', 0, self.fluid.name
                )
            )
        # self.process["vf"] = self.process["vi"]

        self.process["wt"] = self.process["hi"] - self.process["hf"]
        self.process["q"] = 0
        self.process["w"] = -1
        return self.process
