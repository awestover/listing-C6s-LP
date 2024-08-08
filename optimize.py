import pulp as p
from pulp import LpVariable as var
from itertools import product
import os
import sys

EPS = 0.00001
REGIMES = ["max=m", "max=RL13", "max=LR12"]

def compare_le0_ge1(a,b,direction):
    if direction == 0:
        return a <= b
    elif direction == 1:
        return a >= b
    assert False

def solve_case(large12, large23, regime12, regime23):
    assert large12 in [0,1] and large23 in [0,1]
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

    # now we compute the number of C6's
    # 1. enforce that we are in the right regime.
    # 2. award C6's corresponding to this regime

    def add_constraints_count_c6s(e, L, R, regime, LP):
        RL13 = R+L*(1/3)
        LR12 = L+R*(1/2)

        if regime == "max=m":
            LP += e >= RL13
            LP += e >= LR12
            return 6*e - 3*(L+R)
        elif regime == "max=RL13":
            LP += e <= RL13
            LP += RL13 <= LR12
        elif regime == "max=LR12":
            LP += e <= LR12
            LP += RL13 >= LR12
        else: 
            assert False
        return 0

    t12 = add_constraints_count_c6s(d1+W1, minW1W2, maxW1W2, regime12, LP)
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
for thiscase in product([0,1], [0,1], REGIMES, REGIMES):
    large12, large23, regime12, regime23 = thiscase
    ct += 1
    W1W2 = "W1 >= W2" if large12 else "W1 <= W2"
    W2W3 = "W2 >= W3" if large23 else "W2 <= W3"
    print(f"ct:{ct} \t{W1W2}, {W2W3}, \t{regime12}, {regime23}\n")
    alpha = solve_case(*thiscase)
    if alpha > max_alpha:
        max_alpha = alpha
    print("____"*50)

print("MAXALPHA", max_alpha)

