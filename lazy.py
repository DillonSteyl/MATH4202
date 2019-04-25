from read import *
from gurobipy import *
from math import *
from itertools import *
import matplotlib.pyplot as plt
from pandas import DataFrame, read_excel

#-------------------------
# Helper Functions
#-------------------------

def cost(X0,Y0,X1,Y1):
    '''
    Returns the cost between points (X0,Y0) and (X1,Y1)
    '''
    return floor(sqrt((X1-X0)**2 + (Y1-Y0)**2) + 0.5)

#-------------------------
# Initialization
#-------------------------

# Data String
# Format:
# 'Customers-Products-Vehicles-Periods-Version'
datastr = "30-3-3-3-5"

[customers, products, vehicles, period, filename] = getParams(datastr)

[VehicleCap, X_Coords,Y_Coords,I0,Production,Hold,Dem,Capacity] = \
    readFile(customers, products, vehicles, period, filename)

V = list(range(customers+1))
P = list(range(products))
K = list(range(vehicles))
T = list(range(period+1))

Cost = {(i,j): cost(X_Coords[i],Y_Coords[i],X_Coords[j],Y_Coords[j])
                 for i in V for j in V}

VC = [VehicleCap for k in K]

m = Model()

#-------------------------
# Variables
#-------------------------

X = {}
for i in V:
    for j in V[i:]:
        for k in K:
            for t in T[1:]:
                
                if i == 0:
                    X[i,j,k,t] = m.addVar(vtype=GRB.INTEGER, ub=2)
                else:
                    X[i,j,k,t] = m.addVar(vtype=GRB.BINARY)
                
                X[j,i,k,t] = X[i,j,k,t]

Y = {(i,k,t): m.addVar(vtype=GRB.BINARY) 
        for i in V for k in K for t in T[1:]}
Q = {(i,p,k,t): m.addVar() for i in V for p in P for k in K for t in T[1:]}
I = {(i,p,t): m.addVar() for i in V for p in P for t in T}

#-------------------------
# Objective / Paramaters
#-------------------------

m.setObjective(quicksum(Hold[i,p]*I[i,p,t] for i in V for p in P for t in T) + 
               quicksum(Cost[i,j]*X[i,j,k,t] for i in V for j in V
                            for k in K for t in T[1:] if i < j), GRB.MINIMIZE)

m.setParam('LazyConstraints',1)
m._X = X
m._Y = Y

#-------------------------
# Constraints
#-------------------------

initialInventory = {(i,p): m.addConstr(I[i,p,0] == I0[i,p])
                    for i in V for p in P}

inventoryTransferSupplier = {(p,t): m.addConstr(
        I[0,p,t] == I[0,p,t-1] + Production[p,t] 
            - quicksum(Q[i,p,k,t] for k in K for i in V[1:]))
    for p in P for t in T[1:]}

inventoryTransfer = {(i,p,t): m.addConstr(
        I[i,p,t] == I[i,p,t-1] + quicksum(Q[i,p,k,t] for k in K) - Dem[i,p,t])
    for i in V[1:] for p in P for t in T[1:]}

inventoryCapacity = {(i,t): m.addConstr(
        quicksum(I[i,p,t] for p in P) <= Capacity[i])
    for i in V[1:] for t in T[1:]}

legalDelivery = {(i,t): m.addConstr(
        quicksum(Q[i,p,k,t] for p in P for k in K) <= Capacity[i] - 
        quicksum(I[i,p,t-1] for p in P))
    for i in V[1:] for t in T[1:]}

deliveryCapacity = {(i,p,k,t): m.addConstr(Q[i,p,k,t] <= Capacity[i]*Y[i,k,t]) 
    for i in V[1:] for p in P for k in K for t in T[1:]}

vehicleCapacity = {(k,t): m.addConstr(
        quicksum( Q[i,p,k,t] for i in V[1:] for p in P) <= VC[k]*Y[0,k,t])
    for k in K for t in T[1:]}

nodeConnections = {(i,k,t): m.addConstr(
        (quicksum(X[i,j,k,t] for j in V if i < j) + 
        quicksum(X[j,i,k,t] for j in V if j < i)) == 2*Y[i,k,t])
    for i in V for k in K for t in T[1:]}

### VALID INEQUALITIES ###

vIneq1 = {(i,k,t): m.addConstr(X[i,0,k,t] <= 2 * Y[i,k,t]) 
    for i in V for k in K for t in T[1:]}

vIneq2 = {(i,j,k,t): m.addConstr(X[i,j,k,t] <= Y[i,k,t])
    for i in V[1:] for j in V[1:] if i < j for k in K for t in T[1:]}

vIneq3 = {(i,k,t): m.addConstr(Y[i,k,t] <= Y[0,k,t]) 
    for i in V[1:] for k in K for t in T[1:]}

vIneq4 = {(i,p,t): m.addConstr(
        quicksum(Y[i,k,t] for k in K for l in range(t)) >= 
        ((quicksum(Dem[i,p,t] - I[i,p,0] for k in K for l in range(t-1)))/Capacity[i] ) )
    for i in V[1:] for p in P for t in T[1:]}

vIneq5 = {(k,t): m.addConstr(Y[0,k,t] <= Y[0,k-1,t])
    for k in K[1:] for t in T[1:]}

vIneq6 = {(i,k,t): m.addConstr(Y[i,k,t] <= quicksum( Y[j,k-1,t] for j in V[1:] if j < i))
    for i in V[1:] for k in K[1:] for t in T[1:]}

#-------------------------
# Subtour Elimination Constraints
#-------------------------

'''
Find the shortest subtour in a given set of edges
that doesn't include the supplier node.
'''
def subtour(model, edges, k, t):
    n = 0
    done = [False]*(customers+1)
    cycles = []
    lengths = []
    connected = [[] for i in V]
    
    for i in V:
        inc = model.cbGetSolution(model._Y[i,k,t])
        if inc < 0.5:
            done[i] = True
        else:
            n += 1
    
    for (i,j) in edges:
        connected[i].append(j)
    
    while True:
        current = done.index(False)            
        thiscycle = [current] # start cycle with node that hasn't been investigated
        while True:
            done[current] = True
            neighbors = [x for x in connected[current] if not done[x]]
            if len(neighbors) == 0:
                # no non-done neighbours left -> cycle is complete!
                break
            current = neighbors[0]
            thiscycle.append(current)
        
        cycles.append(thiscycle)
        lengths.append(len(thiscycle))
        if sum(lengths) == n:
            # All nodes visited.
            break
    
    if len(cycles) == 1:
        return cycles[0]
    else:
        nonZeroCycles = []
        nonZeroLengths = []
        for cycle in cycles:
            if 0 not in cycle:
                nonZeroCycles.append(cycle)
                nonZeroLengths.append(len(cycle))
        return nonZeroCycles[nonZeroLengths.index(min(nonZeroLengths))]
                

def subtourelim(model, where):

    if where == GRB.callback.MIPSOL:
        for k in K:
            for t in T[1:]:
                selected = [] # edges that are in this route
                number_nodes = 0 # number of nodes in this route
                for i in V:
                    # Append edges that are in this route
                    sol = model.cbGetSolution([model._X[i,j,k,t] for j in V])
                    selected += [(i,j) for j in V if sol[j] > 0.5] 
                    inc = model.cbGetSolution(model._Y[i,k,t])
                    if inc > 0.5:
                        # If this node is included, increase number_nodes
                        number_nodes += 1
                if len(selected) > 0:
                    # Shortest subtour in the current route.
                    tour = subtour(model, selected, k, t)
                    if len(tour) < number_nodes:
                        # Add a subtour elimination constraint
                        expr = 0
                        for i in range(len(tour)):
                            for j in range(i+1, len(tour)):
                                expr += model._X[tour[i],tour[j],k,t]
                        model.cbLazy(expr <= len(tour)-1)

### OPTIMIZE! ###
m.optimize(subtourelim)

#-------------------------
# Output
#-------------------------

## Plots ##

colour = ['r','g','b','y','m']
labels = [str(i) for i in V]
for t in T[1:]:
    if sum(X[i,j,k,t].x for i in V for j in V for k in K) < 0.5:
        print("No delivery on day", t)
    else:
        plt.figure()
        plt.plot(X_Coords, Y_Coords, 'k', marker='o', linestyle='None')
        for k in K:
            for label, x, y in zip(labels, X_Coords, Y_Coords):
                plt.annotate(label, xy = (x-20, y))
            for i in V:
                for j in V:
                    if X[i,j,k,t].x > 0.9:
                        plt.plot([X_Coords[i], X_Coords[j]],
                            [Y_Coords[i], Y_Coords[j]],
                            colour[k])
        plt.axis('off')
        plt.title("Day " + str(t))
        plt.show()
    
      
## Vehicle Reports ##

for k in K:
    print("\n\n***************************")
    print("Vehicle Report - Vehicle {}".format(k))
    print("Capacity:", VC[k])
    for t in T[1:]:
        print("------------------\n")
        print("Day", t)
        print("Product Delivered: ",
                  round(sum(Q[i,p,k,t].x for i in V[1:] for p in P)))
        for p in P:
            print("\nProduct", p)
            demstr = "Demand:".ljust(15)
            invstr = "Inventory (S):".ljust(15)
            invstr2 = "Inventory (E):".ljust(15)
            recstr = "Received:".ljust(15)
            inroute = "In Route:".ljust(15)
            for i in V:
                routebool = False
                demstr += str(round(Dem[i,p,t])).ljust(7)
                invstr += str(round(I[i,p,t-1].x)).ljust(7)
                invstr2 += str(round(I[i,p,t].x)).ljust(7)
                recstr += str(round(Q[i,p,k,t].x)).ljust(7)
                for j in V:
                    if X[i,j,k,t].x > 0.9:
                        routebool = True
                    if X[j,i,k,t].x > 0.9:
                        routebool = True
                inroute += str(routebool).ljust(7)
            print(invstr)
            print(demstr)
            print(recstr)
            print(invstr2)            
            print(inroute)
            
## Cost Analysis

print("As Calculated:", sum(Cost[i,j]*X[i,j,k,t].x for i in V for j in V
                            for k in K for t in T[1:]))

print("\nHolding Costs:")
for t in T:
    print("Day {}: {}".format(t, sum(Hold[i,p] * I[i,p,t].x for i in V for p in P)))
    
print("\nTotal Cost:")
print("$"+str(round(m.objVal,2))+" (with Day 0 Holding)")
print("$"+str(round(m.objVal-sum(Hold[i,p]*I[i,p,0].x for i in V for p in P),2))
    +" (without Day 0 Holding)")

print("\nSolution Time:")
print(m.Runtime)

# Check Solutions
solutions_df = DataFrame(read_excel("Solutions.xls"))
solutions_df = solutions_df.set_index('Key')
col_list = list(solutions_df.columns)

# Find Row
row_arr = re.split('-', datastr)
row_arr[2] = 'k'
row_str = '-'.join(row_arr)
s = solutions_df.loc[row_str]

# Find Column
index_list = list(s.index)
start_ind = index_list.index('K = {}'.format(str(round(vehicles))))
indexes = list(range(start_ind, start_ind+4))
sol_info = list(s[indexes])

print("\n***********************\nSolution:")
print("LB:", sol_info[0])
print("UB:", sol_info[1])
print("Time (s):", sol_info[2])
print("Gap (%):", sol_info[3])