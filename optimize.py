import pulp as p
from pulp import LpVariable as var
from itertools import product
import os
import sys

EPS = 0.00001
REGIMES = ["min=1", "min=m^3/(R^3 L)", "min=m^4/(L^4 R^2)"]
BOOLS = [True, False]

def compare_le0_ge1(a,b,direction):
    if direction:
        return a >= b
    else:
        return a <= b

def solve_case(large12, large23, regime12, regime23, 
               balance12, balance23, anyhex12, anyhex23):
    for param in [large12, large23, balance12, balance23, anyhex12, anyhex23]:
        assert param in BOOLS
    assert regime12 in REGIMES and regime23 in REGIMES

    LP = p.LpProblem('cappedwalksLP', p.LpMaximize)
    #### VARIABLES ####
    # sizes of the parts
    W1 = var("W1", 0, 1)
    W2 = var("W2", 0, 1)
    W3 = var("W3", 0, 1)
    # degrees
    dstar = var("dstar", 0, 1)
    d1 = var("d1", 0, 1)
    d2 = var("d2", 0, 1)
    d3 = var("d3", 0, 1)
    # number of c6's
    t = var("t", 0)

    ####### OBJECTIVE: [capped 3-walks] ######
    LP += W1+d1+d2+d3

    ####### CONSTRAINTS ######
    # First some basic constraints: 
    # at most m edges in the graph, 
    # enforce that dstar is the max degree
    LP += dstar <= .4+EPS
    for deg in [d1, d2, d3]:
        LP += deg <= dstar
    LP += dstar <= 1-W1
    LP += d2 <= 1-W2
    LP += d3 <= 1-W3

    # compare to find max(W1,W2) max(W2,W3)
    LP += compare_le0_ge1(W1,W2,large12)
    maxW1W2 = W1 if large12 else W2
    minW1W2 = W2 if large12 else W1
    LP += compare_le0_ge1(W2, W3, large23)
    maxW2W3 = W2 if large23 else W3
    minW2W3 = W3 if large23 else W2

    # now we figure out whether we are supposed to award C6s or not
    LP += compare_le0_ge1(2*minW1W2, maxW1W2, balance12)
    edgereq12 = (2/3)*(W1+W2) if balance12 else maxW1W2
    LP += compare_le0_ge1(2*minW2W3, maxW2W3, balance23)
    edgereq23 = (2/3)*(W3+W2) if balance23 else maxW2W3

    LP += compare_le0_ge1(W1+d1, edgereq12, anyhex12)
    LP += compare_le0_ge1(W2+d2, edgereq23, anyhex23)

    # now we compute the number of C6's
    # 1. enforce that we are in the right regime.
    # 2. award C6's corresponding to this regime

    def add_constraints_count_c6s(e, L, R, regime, LP):
        m3R3L = 3*e - 3*R - L
        m4L4R2 = 4*e - 4*L - 2*R
        base = 6*e - 3*L - 3*R

        REGIMES = ["min=1", "min=m^3/(R^3 L)", "min=m^4/(L^4 R^2)"]

        if regime == "min=1":
            LP += 0 <= m3R3L
            LP += 0 <= m4L4R2
            return base
        elif regime == "min=m^3/(R^3 L)":
            LP += m3R3L <= 0
            LP += m3R3L <= m4L4R2
            return base + m3R3L
        elif regime == "min=m^4/(L^4 R^2)":
            LP += m4L4R2 <= 0
            LP += m4L4R2 <= m3R3L
            return base + m4L4R2
        else: 
            assert False

    t12 = 0; t23 = 0
    if anyhex12:
        t12 = add_constraints_count_c6s(d1+W1, minW1W2, maxW1W2, regime12, LP)
    if anyhex23:
        t23 = add_constraints_count_c6s(d2+W2, minW2W3, maxW2W3, regime23, LP)
    
    # Number of C6's must be less than number of capped walks
    LP += W1+d1+d2+d3 >= t12
    LP += W1+d1+d2+d3 >= t23
   
    #### SOLVE THE LP ####
    status = LP.solve(p.PULP_CBC_CMD(msg=0))

    def clean(x):
        try:
            y = x.value()
        except:
            y = x
        return round(y, 4)


    ### Display solution ####
    achieved_alpha = p.value(LP.objective)
    if p.LpStatus[LP.status] == "Optimal":
        print("Extremal example:")
        print("Wsizes:", [clean(W) for W in [W1,W2,W3]])
        print("degrees:", [clean(deg) for deg in [d1,d2,d3]])
        print("C6's:", [clean(t12), clean(t23)])
        print("alpha:", clean(achieved_alpha))
        if achieved_alpha > 1.6+EPS:
            print("Counter Example!")
            assert False
        return achieved_alpha
    else: 
        print("INFEASIBLE")
        return 1

max_alpha = 1
ct = 0
for thiscase in product(BOOLS, BOOLS, REGIMES, REGIMES,
                        BOOLS, BOOLS, BOOLS, BOOLS):
    large12, large23, regime12, regime23,_,_,_,_ = thiscase
    ct += 1
    W1W2 = "W1 >= W2" if large12 else "W1 <= W2"
    W2W3 = "W2 >= W3" if large23 else "W2 <= W3"
    print(f"ct:{ct} \t{W1W2}, {W2W3}, \t{regime12}, {regime23}\n")
    print(thiscase)
    alpha = solve_case(*thiscase)
    if alpha > max_alpha:
        max_alpha = alpha
    print("____"*50)

print("MAXALPHA", max_alpha)

