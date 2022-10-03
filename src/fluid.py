from dataclasses import dataclass


@ dataclass
class Fluid:
    name: str = None
    ra: float = None
    cp: float = None
    k: float = None
    r: float = None
    cv_0: float = None
    tc: float = None
    pc: float = None
    zc: float = None
