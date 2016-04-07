from collections import defaultdict
import numpy as np
import pp
pp.define_constants()
c = pp.case24_ieee_rts()

# drop the reactor on bus 6 and most of the capacitance on 6-10
bus_6 = (c["bus"][:,BUS_I]==6).nonzero()[0]
branch_6_10 = (c["branch"][:,(F_BUS, T_BUS)]==[6, 10]).all(1).nonzero()[0]
c["bus"][bus_6, BS]=0
print c["bus"][bus_6, [PD, QD, BS]]
c["bus"][bus_6, [PD, QD, BS]] = [90, 20, 0]
print c["bus"][bus_6, [PD, QD, BS]]

c["branch"][branch_6_10, BR_B] -= 2
print c["branch"][branch_6_10, BR_B]

# take branch 15-24 out of service to put more pressure on 14-16
branch_15_24 = (c["branch"][:,(F_BUS, T_BUS)]==[15, 24]).all(1).nonzero()[0]
c["branch"][branch_15_24, BR_STATUS] = 0

# reduce all loads to 91%
c["bus"][:,[PD, QD]] *= 0.91

# return only digits from a string
digits = lambda p: filter(lambda x: x.isdigit(), p)

# convenience function for calculating total generation on each bus
def gen_by_bus(gen_levels):
    total = defaultdict(float)
    for g in gen_levels:
        bus = int(digits(g[0].split("_")[1]))
        total[bus] += g[1]
    return total

# set generation levels
gen_levels = [
    ('U76_2c', 76.0), ('U76_2d', 76.0), ('U350_23c', 350.0), ('U20_2b', 0.0), 
    ('U155_23b', 155.0), ('U20_2a', 0.0), ('U20_1b', 0.0), ('U50_22b', 50.0), 
    ('U155_15f', 155.0), ('U155_23a', 155.0), ('U400_18', 400.0), ('U100_7c', 100.0), 
    ('U100_7b', 0.0), ('U100_7a', 62.048748556526377), ('U50_22e', 50.0), ('U50_22f', 50.0), 
    ('U155_16', 155.0), ('U197_13c', 0.0), ('U197_13b', 0.0), ('U197_13a', 69.435341456343849), 
    ('U76_1d', 76.0), ('SyncCond_14', 0.0), ('U76_1c', 76.0), ('U12_15d', 0.0), ('U12_15e', 0.0), 
    ('U20_1a', 0.0), ('U50_22c', 50.0), ('U50_22d', 50.0), ('U12_15a', 0.0), ('U12_15b', 0.0),
    ('U12_15c', 0.0), ('U400_21', 400.0), ('U50_22a', 50.0)
]
# calculate average production for each unit size on each bus
total_gen = defaultdict(float)
count = defaultdict(int)
for g in gen_levels:
    if g[0].startswith('U'):    # filter out the synchronous condenser
        size, bus = tuple(int(digits(p)) for p in g[0].split('_'))
        output = g[1]
        total_gen[bus, size] += output
        count[bus, size] += 1

# assign average production to each unit size on each bus
print c["gen"][:,[BUS_I, PMAX, PG]]
for g in c["gen"]:
    if g[PMAX] > 0:     # filter out synchronous condenser
        g[PG] = total_gen[g[GEN_BUS], g[PMAX]] / count[g[GEN_BUS], g[PMAX]]

# solve the power flow
results = pp.runpf(c)[0]
