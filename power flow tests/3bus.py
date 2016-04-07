import pp
pp.define_constants()
from numpy import array

ppc = {"version": '2'}

##-----  Power Flow Data  -----##
## system MVA base
ppc["baseMVA"] = 100.0

## bus data
# bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
ppc["bus"] = array([
    [0, 3, 0, 0, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.9],
    [1, 2, 0, 0, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.9],
    [2, 1, 100, 0, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.9],
])

## generator data
# bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
# Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
ppc["gen"] = array([
    [0,   0, 0, 100, -100, 1.02, 100, 1, 200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 100, 0, 100, -100, 1.02, 100, 1, 200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
])

## branch data
#fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
ppc["branch"] = array([
    [0, 2, 0.00744, 0.0372, 0.0775, 250, 250, 250, 0, 0, 1, -360, 360],
    [1, 2, 0.00644, 0.0372, 0.0775, 250, 250, 250, 0, 0, 1, -360, 360],
])

## generator cost data
# 1 startup shutdown n x1 y1 ... xn yn
# 2 startup shutdown n c(n-1) ... c0
ppc["gencost"] = array([
    [2, 0, 0, 2, 100, 0],
    [2, 0, 0, 2, 100, 0],
])

results=pp.runpf(ppc)[0]
results["gen"][:,(PG, QG)]
results["gen"][:,VG]
results["bus"][:,(VM, VA)]
