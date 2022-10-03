prcs_type_map_dict = {
    "ds": "isentropic",
    "dp": "isobaric",
    "dv": "isochoric",
    "dt": "isothermal",
    "dq": "adiabatic"
}


def map_moisture(x: float) -> str:
    if x < 0:
        result = "subcooled liquid"
    elif x > 1:
        result = "superheated steam"
    elif x == 0 or x == 0.0:
        result = "saturated liquid"
    elif x == 1 or x == 1.0:
        result = "saturated stream"
    else:
        result = "two phase"
    return result


def vdw(v, r, a, b, p, t):
    return ((r*t)/(v-b)) - (a/(v**2)) - p


def vdwprime(v, r, a, b, p, t):
    return (
        (2*a)/(v**3)) - ((r*t)/((v-b)**2)
                         )
