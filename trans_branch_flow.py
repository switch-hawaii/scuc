"""This is a nodal alternative to trans_dispatch. 

Note: this would work better with a single all-zone energy balance, like the CAISO, and then
it would calculate flows and avoid congestion only on congested lines. However, load_zones
are widely used in SWITCH, so we can't eliminate them entirely.

For now, we deactivate the zonal Energy_Balance constraint from load_zones.py and replace it
with an all-zone energy balance here (which sums the components of
LZ_Energy_Components_Produce and LZ_Energy_Components_Consume across all zones, and then 
adds in transmission losses). 

An alternative option would be to retain the zonal Energy_Balance constraint and calculate 
inter-zonal transfer factors based on PTDFs of all inter-zonal branches. Then we would also
need to calculate factors showing expected transmission losses for each zone based on
injections and withdrawals at all generation projects and load nodes. But this would be getting
close to a complete DC power flow, which wouldn't offer any advantages over the iterative
AC power flow (with minimal DC components) that we are currently using.

Key things done by this module:

- replace zonal load balancing constraint with system-level load balancing constraint
- calculate induced branch flows and constrain within branch limits (list of congested 
  branches and their shift factors is built up iteratively)
- calculate energy losses and add them to the system load balancing constraint 
  (collection of loss models is built up incrementally)

"""  

import os, numpy as np, math
from collections import defaultdict
from pyomo.environ import *
try:
    import pypower
except ImportError:
    print "======================================================================"
    print "The {} module requires the pypower package.".format(__name__)
    print "Please install pypower (pip install pypower) before using this module."
    print "======================================================================"
    raise
from pypower.makePTDF import makePTDF
from pypower.runpf import runpf
from pypower.runopf import runopf
from pypower.ppoption import ppoption
from pypower.idx_bus import *
from pypower.idx_brch import *
from pypower.idx_gen import *
from pypower.idx_cost import *

def define_components(m):
    # create bus-level and system-level energy balance lists
    m.Bus_Energy_Components_Produce = []
    m.Bus_Energy_Components_Consume = []
    m.System_Energy_Components_Produce = []
    m.System_Energy_Components_Consume = []

    # base MVA to use for power flows (not very important)
    m.base_mva = Param(default = 100.0)
    
    # extra properties of generator technologies, needed for power flow
    m.g_q_max_mvar = Param(m.GENERATION_TECHNOLOGIES, default=0.0)
    m.g_q_min_mvar = Param(m.GENERATION_TECHNOLOGIES, default=0.0)
    
    # don't allow any power production from synchronous condenser
    # TODO: separate unit size from p_max, to allow for devices like this 
    # (non-zero unit size but zero max real power)?
    m.SYNCHRONOUS_CONDENSERS = Set(initialize = lambda m:
        [p for p in m.PROJECTS if m.proj_gen_tech[p] == 'SyncCond']
    )
    m.Synchronous_Condenser_No_Real_Power = Constraint(m.SYNCHRONOUS_CONDENSERS, m.TIMEPOINTS, 
        rule=lambda m, p, t: 
            m.DispatchProj[p, t] == 0 if (p, t) in m.PROJ_DISPATCH_POINTS
            else Constraint.Skip
    )
    
    # list of all buses in the system
    m.BUSES = Set(dimen=1, ordered=True)
    
    # bus data (mostly passed through to pypower; see pypower.idx_bus for descriptions)
    m.bus_load_zone = Param(m.BUSES, within=m.LOAD_ZONES) # not used
    m.bus_v_setpoint_pu = Param(m.BUSES)
    m.bus_type = Param(m.BUSES)
    m.bus_g_mw = Param(m.BUSES)
    m.bus_b_mvar = Param(m.BUSES)
    m.bus_area = Param(m.BUSES)
    m.bus_base_kv = Param(m.BUSES)
    m.bus_loss_zone = Param(m.BUSES)
    m.bus_v_max_pu = Param(m.BUSES)
    m.bus_v_min_pu = Param(m.BUSES)

    # bus-level load data
    m.bus_demand_mw = Param(m.BUSES, m.TIMEPOINTS, default=0.0)
    m.bus_demand_mvar = Param(m.BUSES, m.TIMEPOINTS, default=0.0)
    # add loads to the bus-level energy components list and remove from load-zone level
    m.Bus_Energy_Components_Consume.append('bus_demand_mw')
    m.LZ_Energy_Components_Consume.remove('lz_demand_mw')
    
    # TODO: add bus-level EV loads to the consumption side (in a different module)

    # branch data; most of these are just passed through to pypower (see pypower.idx_brch for definitions)

    # list of all branches (transmission lines and transformers) in the system
    m.BRANCHES = Set(dimen=1, ordered=True)

    # should this branch be treated as a contingency?
    m.branch_contingency = Param(m.BRANCHES, within=Binary)

    # list of contingencies (branches that could be taken off-line, or None for base case)
    m.CONTINGENCIES = Set(ordered=True, initialize=lambda m: 
        [None]+[b for b in m.BRANCHES if m.branch_contingency[b]]
    )
    
    # start and end bus for each branch
    m.branch_from = Param(m.BRANCHES, within=m.BUSES)
    m.branch_to = Param(m.BRANCHES, within=m.BUSES)

    # impedance of each branch
    m.branch_r_pu = Param(m.BRANCHES)	
    m.branch_x_pu = Param(m.BRANCHES)
    m.branch_b_pu = Param(m.BRANCHES)

    # branch ratings (note: for now we use the same rating for normal and emergency operation)
    m.branch_rating_mw_long_term = Param(m.BRANCHES)
    m.branch_rating_mw_short_term = Param(m.BRANCHES)
    m.branch_rating_mw_emergency = Param(m.BRANCHES)
    
    # voltage tap ratio (for transformers); will be 0.0 for transmission lines
    m.branch_tap_ratio = Param(m.BRANCHES)
    
    # transformer phase shift angle (degrees)
    m.branch_phase_shift = Param(m.BRANCHES)
        
    # is the branch in service?
    m.branch_status = Param(m.BRANCHES, within=Binary)
    
    # create an ordered list of generators. This is the same as PROJECTS, except that
    # the list can be referenced in a consistent order, so it can be matched to 
    # arrays used in the power flow analysis.
    m.GENERATORS = Set(initialize=lambda m: [p for p in m.PROJECTS], ordered=True)

    # generator-bus linkages
    m.proj_bus = Param(m.GENERATORS)
    def BUS_GENERATORS_init(m, b):
        if not hasattr(m, "BUS_GENERATORS_dict"):
            m.BUS_GENERATORS_dict = defaultdict(list)
            for g in m.GENERATORS:
                m.BUS_GENERATORS_dict[m.proj_bus[g]].append(g)
        return m.BUS_GENERATORS_dict[b]
    m.BUS_GENERATORS = Set(m.BUSES, initialize=BUS_GENERATORS_init)
    
    # temporary rule to force all gens on at bus 15
    m.All_on_15 = Constraint(m.PROJ_DISPATCH_POINTS, rule=lambda m, p, tp: 
        m.DispatchProj[p, tp] == m.DispatchUpperLimit[p, tp] if m.proj_bus[p] == 15 else Constraint.Skip
    )
    
    # add generator dispatch to the bus level energy components and remove from the load-zone level
    m.Bus_Gen_Dispatch = Expression(
        m.BUSES, m.TIMEPOINTS,
        rule=lambda m, b, t: sum(
            m.DispatchProj[g, t] 
                for g in m.BUS_GENERATORS[b]
                    if (g, t) in m.PROJ_DISPATCH_POINTS)
        )
    m.Bus_Energy_Components_Produce.append('Bus_Gen_Dispatch')
    m.LZ_Energy_Components_Produce.remove('LZ_NetDispatch')


    # zero-based index of each bus, branch and generator
    # (needed by pypower instead of the actual bus number)
    m.bus_idx = Param(m.BUSES, initialize=lambda m, b: m.BUSES.ord(b) - 1)
    m.branch_idx = Param(m.BRANCHES, initialize=lambda m, b: m.BRANCHES.ord(b) - 1)
    m.gen_idx = Param(m.GENERATORS, initialize=lambda m, g: m.GENERATORS.ord(g) - 1)

    # gather various elements needed for power flow analyses
    # and attach them to the model for later reference
    def Get_Power_Flow_Data_rule(m):
        m.bus_array = np.array([
            (
                m.bus_idx[b],  # should automatically end up in sequence from 0 to len(m.BUSES)-1
                m.bus_type[b],
                0.0,    # real power demand (to fill in later)
                0.0,    # reactive power demand (to fill in later)
                m.bus_g_mw[b],
                m.bus_b_mvar[b],
                m.bus_area[b],
                1.0,    # voltage magnitude (power flow will calculate)
                0.0,    # voltage angle (power flow will calculate)
                m.bus_base_kv[b],
                m.bus_loss_zone[b],
                m.bus_v_max_pu[b],
                m.bus_v_min_pu[b]
            )
            for b in m.BUSES
        ])

        m.branch_array = np.array([
            (
                m.bus_idx[m.branch_from[b]],
                m.bus_idx[m.branch_to[b]],
                m.branch_r_pu[b],
                m.branch_x_pu[b],
                m.branch_b_pu[b],
                m.branch_rating_mw_long_term[b],
                m.branch_rating_mw_short_term[b],
                m.branch_rating_mw_emergency[b],
                # convert all 0 tap settings to 1 to simplify later usage
                1.0 if m.branch_tap_ratio[b] == 0.0 else m.branch_tap_ratio[b],
                m.branch_phase_shift[b],
                m.branch_status[b]
            )
            for b in m.BRANCHES
        ])

        # make a list of the bus indexes for the "from" and "to" end of each branch
        # (useful for retrieving properties of the relevant buses later)
        # m.branch_bus_idx = m.branch_array[:,[0, 1]].astype(int)
        
        # tabulate data for all generators (only the columns needed to define a power flow)
        m.gen_array = np.zeros((len(m.GENERATORS), 21))
        for i, g in enumerate(m.GENERATORS):
            r = m.gen_array[i]
            gt = m.proj_gen_tech[g]
            r[GEN_BUS] = m.bus_idx[m.proj_bus[g]]
            r[QMAX] = m.g_q_max_mvar[gt]
            r[QMIN] = m.g_q_min_mvar[gt]
            r[VG] = m.bus_v_setpoint_pu[m.proj_bus[g]]
            r[GEN_STATUS] = 1

        # dummy cost array to be used for minimizing losses in opf (while maintaining the scheduled generation level)
        m.gencost_array = np.array([
            [
                2,  # polynomial model
                0.0,    # startup cost
                0.0,    # shutdown cost
                2,      # number of cost coefficients
                100.0,  # cost per MWh
                0.0     # constant cost term
            ]
            for g in m.GENERATORS
        ])
        
        # Create a dictionary to hold arrays of power transfer distribution factors (ptdfs) 
        # for each power flow contingency.
        # These will be calculated as needed and then cached in the dictionary.
        # The ptdf array has a coefficient for each bus / branch pair, showing
        # the change in flow expected on each branch due to change in injections at each bus.
        # (This is a numpy array with one row for each branch, one column for each bus.)
        m.contingency_ptdf_dict = dict()

    m.Get_Power_Flow_Data = BuildAction(rule=Get_Power_Flow_Data_rule)

    # calculate bus-level injections and withdrawals
    m.Bus_Total_Production = Expression(m.BUSES, m.TIMEPOINTS, rule=lambda m, b, t:
        sum(
            getattr(m, component)[b, t]
            for component in m.Bus_Energy_Components_Produce
        )
    )
    m.Bus_Total_Load = Expression(m.BUSES, m.TIMEPOINTS, rule=lambda m, b, t:
        sum(
            getattr(m, component)[b, t]
            for component in m.Bus_Energy_Components_Consume
        )
    )
    m.Bus_Net_Injection = Expression(m.BUSES, m.TIMEPOINTS, rule=lambda m, b, t:
        m.Bus_Total_Production[b, t] - m.Bus_Total_Load[b, t]
    )

    # calculate system-wide tranmsission losses
    # This is setup as a decision variable which must exceed the losses that would
    # be forecast based on lost sensitivities found in all prior runs
    m.TRANS_LOSS_CONSTRAINTS = Set(initialize=[], ordered=True)

    # expected change in losses per MW injected at each bus (from power flow analysis and previous solutions)
    # amount of extra transmission losses expected for each MW injected at each bus,
    # based on conditions during a particular timepoint
    # We retain collections of sensitivity factors from previous runs, to build up
    # a "bowl" showing the estimates of transmission losses for states surrounding
    # the state at the end of each prior run.
    m.bus_trans_loss_sensitivity = Param(m.TRANS_LOSS_CONSTRAINTS, m.TIMEPOINTS, m.BUSES, mutable=True)

    # baseline transmission losses (from power flow analysis of previous solutions)
    m.transmission_losses_baseline = Param(m.TRANS_LOSS_CONSTRAINTS, m.TIMEPOINTS, mutable=True)
    
    m.TransmissionLosses = Var(m.TIMEPOINTS)
    
    m.TransmissionLosses_Calculation = Constraint(m.TRANS_LOSS_CONSTRAINTS, m.TIMEPOINTS, 
        rule=lambda m, i, tp:
            m.TransmissionLosses[tp]
            >=
            m.transmission_losses_baseline[i, tp] 
            + sum(m.bus_trans_loss_sensitivity[i, tp, b] * m.Bus_Net_Injection[b, tp] for b in m.BUSES)
    )
    
    # don't allow negative transmission losses in the initial state (before baseline has been defined)
    m.TransmissionLossesInitial_Calculation = Constraint(m.TIMEPOINTS, rule=lambda m, t:
        Constraint.Skip if len(m.TRANS_LOSS_CONSTRAINTS) > 0 else m.TransmissionLosses[t] >= 0
    )

    # add transmission losses to the system energy balance
    m.System_Energy_Components_Consume.append("TransmissionLosses")

    # enforce constraints on branch flows under base case and all contingencies
    
    # list of congested branches in contingencies, indexed by
    # a (contingency, branch) tuple, where contingency identifies
    # a branch that has been removed (or None for the base case),
    # and branch identifies a branch that is congested (at some point) 
    # in that contingency.
    m.CONGESTED_BRANCHES = Set(dimen=2, initialize=lambda m: [])
    # list of active congested branches (in contingencies) for each timepoint
    # (contingency, branch, tp)
    m.CONGESTED_BRANCHES_BY_TP = Set(dimen=3, initialize=lambda m: [])

    # power transfer distribution factors, sensitivities and baseline flows
    # for each bus for congested branches
    m.congested_ptdf = Param(m.CONGESTED_BRANCHES, m.BUSES, mutable=True)
    m.congested_branch_sensitivity = Param(m.CONGESTED_BRANCHES_BY_TP, mutable=True)
    m.congested_branch_base_flow = Param(m.CONGESTED_BRANCHES_BY_TP, mutable=True)

    # penalty factor to allow violation of branch limits
    # (should get cleared to zero once the model reaches its final state)
    m.branch_limit_violation_penalty = Param()

    # choose/calculate the size of the branch limit violation each timepoint
    m.ExceedBranchLimitMW = Var(m.CONGESTED_BRANCHES_BY_TP, within=NonNegativeReals)

    # calculate total penalty for exceeding branch limits
    m.Total_Branch_Limit_Violation_Penalty = Expression(m.TIMEPOINTS, rule=lambda m, tp: 
        sum(
            m.ExceedBranchLimitMW[o, b, t] * m.branch_limit_violation_penalty
                for o, b, t in m.CONGESTED_BRANCHES_BY_TP if t==tp
        )
    )
    m.cost_components_tp.append('Total_Branch_Limit_Violation_Penalty')

    # constrain flow on each branch to be below 99% of the limit 
    # (this leaves a little slack for convergence, which just checks 
    # that they are below 100% of the limit, but based on a more detailed model)
    m.Max_Branch_Flow = Constraint(m.CONGESTED_BRANCHES_BY_TP, rule=lambda m, out_br, br, tp:
        m.congested_branch_base_flow[out_br, br, tp] 
        + 
        m.congested_branch_sensitivity[out_br, br, tp] 
            * sum(m.congested_ptdf[out_br, br, b] * m.Bus_Net_Injection[b, tp] for b in m.BUSES)
        <=
        0.99 * m.branch_rating_mw_long_term[br] 
        + m.ExceedBranchLimitMW[out_br, br, tp]
    )

    # create system-wide energy balance requirement
    # (all bus-level components, any remaining zone-level components, plus all system-level components)
    m.System_Energy_Balance = Constraint(
        m.TIMEPOINTS,
        rule=lambda m, t: 
        sum(m.Bus_Total_Production[b, t] for b in m.BUSES)
        -
        sum(m.Bus_Total_Load[b, t] for b in m.BUSES)
        +
        sum(
            getattr(m, component)[lz, t]
                for component in m.LZ_Energy_Components_Produce
                    for lz in m.LOAD_ZONES
        )
        -
        sum(
            getattr(m, component)[lz, t]
                for component in m.LZ_Energy_Components_Consume
                    for lz in m.LOAD_ZONES
        )
        +
        sum(
            getattr(m, component)[t]
                for component in m.System_Energy_Components_Produce
        )
        -
        sum(
            getattr(m, component)[t]
                for component in m.System_Energy_Components_Consume
        )
        == 0.0
    )

def reconstruct_energy_balances(m):
    # reconstruct expressions which may have been changed, as well as their dependencies
    # (is there a more systematic way to do this?)
    for c in [
        m.Bus_Total_Production, 
        m.Bus_Total_Load, 
        m.Bus_Net_Injection,
        m.TransmissionLosses_Calculation, 
        m.TransmissionLossesInitial_Calculation,
        m.System_Energy_Balance
    ]:
        if c.is_constructed():
            c.reconstruct()

def define_dynamic_components(m):
    # deactivate load zone energy balance requirement
    # (There should be nothing left in it anyway.)
    # This has to be scheduled to be called later (during construction time), 
    # since Energy_Balance.deactivate cannot be called before its index set is constructed. 
    m.Deactivate_LZ_Energy_Balance = BuildAction(rule=lambda m: m.Energy_Balance.deactivate())
    
    # reconstruct expressions which may have been changed in other modules
    # (e.g., demand response may add a term to m.Bus_Energy_Components_Consume)
    # This is scheduled to occur at this point in the construction sequence
    # (as if we had just defined them here, with deferred construction), rather than immediately.
    m.Reconstruct_Energy_Balances = BuildAction(rule=reconstruct_energy_balances)


def contingency_ptdf(m, out_branch):
    if out_branch not in m.contingency_ptdf_dict:
        m.contingency_ptdf_dict[out_branch] = get_ptdf_for_contingency(m, out_branch)
    return m.contingency_ptdf_dict[out_branch]
    
def get_ptdf_for_contingency(m, out_branch):
    """Return an array of power transfer distribution factors (ptdfs) for each 
    bus / branch pair defined in model m, under the specified contingency
    (branch out_branch taken out of service, or base case if out_branch is None)."""
    branch_adj = np.copy(m.branch_array)
    if out_branch is not None:
        # deactivate the line for this contingency
        branch_adj[m.branch_idx[out_branch], BR_STATUS] = 0
    return makePTDF(value(m.base_mva), m.bus_array, branch_adj)

def mag(*args):
    return np.sqrt(sum(np.multiply(a, a) for a in args))

def p2c(mag, angle):
    """Convert a number reported as magnitude and angle (in degrees) into a complex number."""
    return mag * np.exp(1j * np.deg2rad(angle))
    
def r2c(real, imag):
    """Convert a number reported as real and imaginary parts into a complex number."""
    return real + 1j * imag

def post_iterate(m):
    print "objective value:", value(m.Minimize_System_Cost)
    # print "cost_components_tp:", m.cost_components_tp
    print "\ncosts by period and component:"
    costs = [
        (p, tp_cost, value(sum(getattr(m, tp_cost)[t] * m.tp_weight_in_year[t] for t in m.PERIOD_TPS[p])))
            for p in m.PERIODS for tp_cost in m.cost_components_tp
    ]
    print costs
    converged = run_power_flows(m)
    return converged

def run_power_flows(m):
    """Run power flow analyses for every timepoint of the current solution, for every possible
    contingency, and update the model with new base levels and sensitivity factors for losses 
    and branch flows. Return True if all branches are within their limits.
    
    Branch flow parameters are set for each line that is congested (flow above 95% of limit)
    in any contingency. Parameters from previous iterations are retained even if no longer congested,
    to keep solutions stable (as is done by CAISO; not clear if this is necessary).
    
    Loss sensitivities are set only for the base case. Loss sensitivity factors from previous cases
    are retained in order to build up a collection of constraints on the loss calculation, forcing
    increasing accuracy near the final operating state. (Loss sensitivity factors are not shared
    across timepoints, but they could be.)
    """

    # TODO: modify convergence rules to only report non-convergence if 
    # new data comes from the AC power flow. (If the linear model already
    # matches the losses and base flows in the AC power flow, then further
    # iterations will not help. In this case, if some lines are still overloaded
    # it's because the model is choosing to pay the penalty rather than relieve
    # them, which may mean it is impossible to avoid the overload. That should
    # be reported, but shouldn't stop convergence.)
    converged = True
    
    # create an iteration ID and add it to the list of iterations
    if len(m.TRANS_LOSS_CONSTRAINTS) == 0:
        iter = 1
    else:
        iter = m.TRANS_LOSS_CONSTRAINTS.last() + 1
    m.TRANS_LOSS_CONSTRAINTS.add(iter)
        
    slack_buses = m.bus_array[:,BUS_I].compress(m.bus_array[:,BUS_TYPE]==REF).astype(int)
    if len(slack_buses) != 1:
        raise ValueError("The system should have exactly one slack bus.")
    slack_bus = slack_buses[0]
    slack_gens = np.nonzero(m.gen_array[:,GEN_BUS]==slack_bus)[0]

    #import pdb; pdb.set_trace()

    for tp in m.TIMEPOINTS:
    
        case = dict(
            version=2,
            baseMVA=value(m.base_mva),
            bus=np.copy(m.bus_array),
            branch=np.copy(m.branch_array),
            gen=np.copy(m.gen_array),
            gencost=np.copy(m.gencost_array),
            # augment the case dictionary with other elements that aren't used for the power flow
            # but are useful for calculating base levels for losses and branch flows, to make the DC
            # model match the AC results.
            net_injection=np.array([value(m.Bus_Net_Injection[b, tp]) for b in m.BUSES]),
            slack_bus=slack_bus,
            slack_gens=slack_gens
        )

        # report the total load
        case['bus'][:,[PD,QD]] = [
            [value(m.Bus_Total_Load[b, tp]), value(m.bus_demand_mvar[b, tp])] for b in m.BUSES
        ]

        # use current generator output levels
        case['gen'][:,PG] = [
            value(m.DispatchProj[g, tp]) if (g, tp) in m.PROJ_DISPATCH_POINTS else 0.0
                for g in m.GENERATORS
        ]
        case['gen'][:,[PMIN,PMAX]] = case['gen'][:,[PG]]
        # relax constraints on the "slack" generators when running opf
        # import pdb; pdb.set_trace()
        case['gen'][slack_gens[:,np.newaxis],[PMIN,PMAX]] = [0,10000]
        
        print ""
        print "================================================================="
        print "Iteration {}, timepoint {}".format(iter, tp)
        print "generation: {}".format([(m.GENERATORS[i+1], g) for i, g in enumerate(case['gen'][:,PG])])
        print "generation (from LP): {}".format([(b, value(m.Bus_Total_Production[b, tp])) for b in m.BUSES])
        print "loads: {}".format([(m.BUSES[i+1], d) for i, d in enumerate(case['bus'][:,PD])])
        print "================================================================="
        print ""
        # import pdb; pdb.set_trace()
        # calculate expected flow along all branches, for every contingency

        for out_branch in m.CONTINGENCIES:
            # TODO: skip if the line outage distribution factors for this line don't
            # suggest any risk (e.g., all affected lines will move less than 50% 
            # toward overload if this line fails) and this line is not already in the
            # list of managed constraints (may also be able to focus only on the lines
            # that are expected to be congested if this fails)
            print ""
            print "================================================================="
            print "Iteration {}, timepoint {}, testing with contingency: {}".format(iter, tp, out_branch)
            print "================================================================="
            print ""
            if out_branch is not None:
                # running a contingency case
                out_idx = m.branch_idx[out_branch]
                # deactivate the line
                case["branch"][out_idx, BR_STATUS] = 0
                # run ac power flow
                solved_case = calc_power_flow(case)
                # return the line to its original status
                case["branch"][out_idx, BR_STATUS] = m.branch_array[out_idx, BR_STATUS]
            else:
                # run ac power flow and evaluate base-case properties
                solved_case = calc_power_flow(case)
                # record the loss sensitivity factors for each bus in the current state
                # and adjust the converged status
                if not apply_bus_loss_sensitivities(m, iter, tp, solved_case):
                    converged = False
                    print "Iteration {}, timepoint {} failed to converge in apply_bus_loss_sensitivities()".format(iter, tp)
            
            # record sensitivity factors for any congested branches
            for branch_idx, branch in enumerate(solved_case["branch"]):
                branch_id = m.BRANCHES[branch_idx+1]
                # calculate apparent power (proportional to current) at each end of line
                # note: it may be better to normalize power based on bus voltage, but we don't
                S = (mag(branch[PF], branch[QF]), mag(branch[PT], branch[QT]))
                P = (branch[PF], -branch[PT])   # amount flowing in the "from -> to" direction at either end
                # identify the end that's a bigger problem (could do both but we don't)
                congested_end = 0 if S[0] >= S[1] else 1
                # if any line is overloaded, we need to keep iterating
                if S[congested_end] > branch[RATE_A]:
                    converged = False
                    print "Iteration {} failed to converge at timepoint {} with {} out; branch {} is overloaded, {} > {}".format(iter, tp, out_branch, branch_id, S[congested_end], branch[RATE_A])
                    # if iter == 6 and tp == 20 and out_branch == "A_15_24":
                    #     import pdb; pdb.set_trace()
                    if (out_branch, branch_id, tp) not in m.CONGESTED_BRANCHES_BY_TP:
                        print "This is not currently in the linear program."
                    else:
                        print "Linear program total flow is:", value(
                            m.congested_branch_base_flow[out_branch, branch_id, tp] 
                            + 
                            m.congested_branch_sensitivity[out_branch, branch_id, tp] 
                                * sum(m.congested_ptdf[out_branch, branch_id, b] * m.Bus_Net_Injection[b, tp] for b in m.BUSES)
                        )
                        # print "Linear program sensitivities:", [
                        #     (b, value(m.congested_ptdf[out_branch, branch_id, b])) for b in m.BUSES
                        # ]
                    
                    
                # check for congested branches
                # (or reevaluate ones from prior contingencies)
                if  (
                    S[congested_end] > 0.95 * branch[RATE_A]
                    or (out_branch, branch_id, tp) in m.CONGESTED_BRANCHES_BY_TP
                ):
                    # expand the list of congested branches/contingencies we are considering
                    # and store the PTDFs for the congested line (if not already calculated)
                    if (out_branch, branch_id) not in m.CONGESTED_BRANCHES:
                        m.CONGESTED_BRANCHES.add((out_branch, branch_id))
                        for bus in m.BUSES:
                            m.congested_ptdf[out_branch, branch_id, bus] = ( 
                                contingency_ptdf(m, out_branch)[branch_idx, m.bus_idx[bus]]
                            )

                    # record the updated sensitivity and base flow for this contingency/branch/timepoint
                    if (out_branch, branch_id, tp) not in m.CONGESTED_BRANCHES_BY_TP:
                        m.CONGESTED_BRANCHES_BY_TP.add((out_branch, branch_id, tp))
                    
                    # note: |S| = sqrt(P**2 + Q**2), so d|S|/dP = P/|S|; this includes a sign compensation for the current direction of flow
                    congestion_sensitivity = P[congested_end] / S[congested_end]

                    m.congested_branch_sensitivity[out_branch, branch_id, tp] = congestion_sensitivity

                    m.congested_branch_base_flow[out_branch, branch_id, tp] = (
                        S[congested_end] - congestion_sensitivity * np.dot(
                            solved_case["net_injection"], 
                            contingency_ptdf(m, out_branch)[branch_idx,:]
                        )
                    )

                if S[congested_end] > branch[RATE_A]:
                    print "total flow will now be calculated as:", value(
                        m.congested_branch_base_flow[out_branch, branch_id, tp]
                        + congestion_sensitivity * np.dot(
                            solved_case["net_injection"], 
                            contingency_ptdf(m, out_branch)[branch_idx,:]
                        )
                    )

    # reconstruct the constraint that forces calculation of transmission losses
    m.TransmissionLosses_Calculation.reconstruct()
    m.TransmissionLossesInitial_Calculation.reconstruct()
                                            
    # reconstruct model element(s) that depend on the CONGESTED_BRANCHES lists and parameters
    m.ExceedBranchLimitMW.reconstruct()
    m.Max_Branch_Flow.reconstruct()
    m.Total_Branch_Limit_Violation_Penalty.reconstruct()
    m.SystemCostPerPeriod.reconstruct()
    m.Minimize_System_Cost.reconstruct()

    return converged

def calc_power_flow(case):
    # we run an opf, but real power output is fixed everywhere except a single quasi-slack bus,
    # so it just adjusts the voltage setpoints to minimize losses (and hopefully get more 
    # reasonable solutions than a raw power flow)
    #results = runopf(case)
    #results = runpf(case, ppopt=ppoption(ENFORCE_Q_LIMS=1))[0]
    results = runpf(case, ppopt=ppoption(OUT_ALL=0))[0]
    slack_bus = case["slack_bus"]
    slack_gens = case["slack_gens"]
    # add in the extra slack generation introduced by the model, so the results
    # show the operating state accurately (even if it differs from the proposed state)
    results["net_injection"][slack_bus] += (
        np.sum(results["gen"][slack_gens, PG] - case["gen"][slack_gens, PG])
    )
    return results
    

def apply_bus_loss_sensitivities(m, i, tp, case):
    """Calculate the incremental transmission losses for injections on any bus, based on the
    current operating state. Extend the list of contraints on TransmissionLosses to include this state."""

    bls = bus_loss_sensitivities(m, case)

    for b in m.BUSES:
        m.bus_trans_loss_sensitivity[i, tp, b] = bls[m.bus_idx[b]]

    # calculate the baseline loss level that would result in the current total losses 
    # with the current injections
    # note: PT is negative if power is flowing out of the "to" end, 
    # so the first line measures power _lost_ in each branch
    total_trans_losses = np.sum(case["branch"][:,PF] + case["branch"][:,PT])  

    # m.transmission_losses_baseline[i, tp] = total_trans_losses - np.dot(case["net_injection"], bls)
    m.transmission_losses_baseline[i, tp] = total_trans_losses - np.dot(case["net_injection"], bls)

    # check whether the model is estimating transmission losses accurately:
    print "Linearized transmission losses: {}; actual losses: {}".format(value(m.TransmissionLosses[tp]), total_trans_losses)
    print "Loss estimation error: {}".format(value(m.TransmissionLosses[tp]) / total_trans_losses - 1)
    #import pdb; pdb.set_trace()

    # we declare the model converged when its estimate of losses at the current operating point is at least 99% of the true value
    # (note: this is one-sided because due to imperfect sensitivity values, sometimes the model's estimate
    # of losses at the current operating point can slightly exceed the true value)
    converged = (value(m.TransmissionLosses[tp]) >= 0.99 * total_trans_losses)
    return converged

def bus_loss_sensitivities(m, case):
    """calculate bus loss sensitivities (expected change in losses when 1 more MW is injected 
    at each bus) based on power flow reported in the AC power flow analysis.
    """
    br_arr = case["branch"]
    bus_arr = case["bus"]
    base_mva = case["baseMVA"]

    # for branch j: 
    # Loss_j = |I_j|^2 R_j = ( |S_j|^2 / |V_j|^2 ) R_j = (P_j^2 + Q_j^2) * R_j / |V_j|^2
    # so dL_j / dP_j = 2 P_j R_j / |V_j|^2

    # simple approach to calculating sensitivities
    # note: V =~ 1 (esp. in DC model), so I =~ S for each line
    # If we take V_j as 1 (standard assumption for DC power flow), introduce the PTDFs 
    # (dP_j / dP_i for injections at bus i), and sum over all branches, we get the equation below: 
    bls = np.dot(
        (2 * (br_arr[:,PF] / base_mva) * br_arr[:,BR_R]), 
        contingency_ptdf(m, None)
    )

    # NOTE: below is a way to get loss sensitivities that are a few percent more accurate,
    # but it is a lot harder to explain.
    # This code assumes V_j may not be 1.0. Also, since V_j and P_j can vary between 
    # the two ends of branch j, this calculates
    # separate sensitivity factors for each end of the line, and then averages them.
    # (The simpler code above only uses power at the "from" end, and assumes V=1.0 at both ends.)

    # # get complex voltage at "from" and "to" end of each line
    # v_bus = p2c(bus_arr[:,[VM]], bus_arr[:,[VA]])
    # # voltage at each end of each branch
    # v_br = v_bus[m.branch_bus_idx, 0]
    # # account for tap settings (assumed to apply between "from" bus and the line)
    # # (see, e.g., get_losses() in matpower or the matpower manual)
    # v_br[:,1] /= br_arr[:,TAP]
    # # get sensitivities at both ends of the line (may be slightly different)
    # blsft = np.dot(
    #     (2 * (np.c_[br_arr[:,[PF]], -br_arr[:,[PT]]] / base_mva) * br_arr[:,[BR_R]]
    #         / (abs(v_br)**2)).T,
    #     contingency_ptdf(m, None)
    # )
    # bls = np.mean(blsft.T, axis=1)

    return bls

def load_inputs(m, switch_data, inputs_dir):
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'branches.tab'),
        auto_select=True,
        index=m.BRANCHES,
        param=(
            m.branch_contingency,
            m.branch_from,
            m.branch_to,
            m.branch_r_pu,
            m.branch_x_pu,
            m.branch_b_pu,
            m.branch_rating_mw_long_term,
            m.branch_rating_mw_short_term,
            m.branch_rating_mw_emergency,
            m.branch_tap_ratio,
            m.branch_phase_shift,
            m.branch_status,
        ))
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'buses.tab'),
        auto_select=True,
        index=m.BUSES,
        param=(
            m.bus_load_zone,
            m.bus_v_setpoint_pu,
            m.bus_type,
            m.bus_g_mw,
            m.bus_b_mvar,
            m.bus_area,
            m.bus_base_kv,
            m.bus_loss_zone,
            m.bus_v_max_pu,
            m.bus_v_min_pu,
        ))
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'project_info.tab'),
        auto_select=True,
        param=(
            m.proj_bus,
        ))
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'bus_loads.tab'),
        auto_select=True,
        param=(
            m.bus_demand_mw,
            m.bus_demand_mvar,
        ))
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'generator_info.tab'),
        auto_select=True,
        param=(
            m.g_q_max_mvar,
            m.g_q_min_mvar,
        ))
    switch_data.load_aug(filename=os.path.join(inputs_dir, 'power_flow.dat'))
