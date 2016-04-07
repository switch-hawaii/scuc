from numpy import array
import pp
pp.define_constants()

ppc = {"version": '2'}

##-----  Power Flow Data  -----##
## system MVA base
ppc["baseMVA"] = 100.0

t=17.625

## bus data
# bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
ppc["bus"] = array([
    [0, 3, 0, 0, 0, 0, 1, 0.9703, 0, 345, 1, 1.05, 0.95],
    [1, 1, 20*t, -20.3121*t, 0, 0, 1, 1.03, -64.3755, 345, 1, 1.03, 0.95],
])

## generator data
# bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
# Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
ppc["gen"] = array([
    [0, 447.67, 153.55, 900, -900, 0.9703, 100, 1, 900, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
])

## branch data
#fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
ppc["branch"] = array([
    [0, 1, 0.04, 0.20, 0.0, 900, 900, 900, 0, 0, 1, -360, 360],
])

## generator cost data
# 1 startup shutdown n x1 y1 ... xn yn
# 2 startup shutdown n c(n-1) ... c0
ppc["gencost"] = array([
    [2, 1500, 0, 3, 0, 1, 0],
])

results=pp.runpf(ppc)[0]
results["bus"][:,VM]
