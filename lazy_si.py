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
datastr = "10-1-1-3-1"

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

m._X = X
m._Y = Y
m._Q = Q
m._I = I
m._customers = customers
m._vehicles = vehicles
m._products = products
m._period = period
m._I0 = I0
m._Hold = Hold
m._Dem = Dem
m._Capacity = Capacity
m._Production = Production
m._Cost = Cost
m._VC = VC
m._Best = 1e10

m.setParam('LazyConstraints',1)
m.setParam('OutputFlag', 1)

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
        ((quicksum(Dem[i,p,t] - I[i,p,0] for k in K for l in range(t-1)))/Capacity[i] ))
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
                

def subtourelim_si(model, where):

    if where == GRB.callback.MIPSOL:
        currentSol = model.cbGet(GRB.callback.MIPSOL_OBJBST)

        if currentSol < m._Best - 0.1:

            # Solve a subproblem...
            mSI = Model()
            SI_V = list(range(model._customers + 1))
            SI_P = list(range(model._products))
            SI_K = list(range(model._vehicles))
            SI_T = list(range(model._period + 1))
            
            SI_X = {(i,j,k,t): model.cbGetSolution(model._X[i,j,k,t])
                for i in SI_V for j in SI_V for k in SI_K for t in SI_T[1:]}

            SI_Y = {(i,k,t): model.cbGetSolution(model._Y[i,k,t])
                for i in SI_V for k in SI_K for t in SI_T[1:]}
            
            # all edges
            edges = [(i,j,k,t) for i in SI_V for j in SI_V for k in SI_K for t in SI_T[1:] if SI_X[i,j,k,t] > 0.5]

            # define neighbours of a node            
            neighbours = {(i,k,t): [j for j in SI_V if (i,j,k,t) in edges]
                for i in SI_V for k in SI_K for t in SI_T[1:]}

            # (i,k,t): routing cost reduction from removal of i
            SI_A = {}
            for i in SI_V[1:]:
                for k in SI_K:
                    for t in SI_T[1:]:
                        if neighbours[i,k,t] == []:
                            SI_A[i,k,t] = 0
                        else:
                            if len(neighbours[i,k,t]) == 1:
                                SI_A[i,k,t] = 2*model._Cost[0,i]
                            else:
                                SI_A[i,k,t] = (
                                    sum(model._Cost[i,j] for j in neighbours[i,k,t]) -
                                    model._Cost[neighbours[i,k,t][0], neighbours[i,k,t][1]])
            
            # insertionChange(ii,i,j,k,t): routing cost change if ii is inserted between i and j
            insertionChange = {}
            for ii in SI_V[1:]:
                for (i,j,k,t) in edges:
                    insertionChange[ii,i,j,k,t] = model._Cost[ii,i] + model._Cost[ii,j] - model._Cost[i,j]

            # SI_B(ii,k,t): best routing cost reduction from insertion
            # insertNodes(ii,k,t): (i,j) such that insertionChange is minimal
            SI_B = {}
            insertNodes = {}
            for ii in SI_V[1:]:
                for k in SI_K:
                    for t in SI_T[1:]:
                        try:
                            thisChange = min(insertionChange[ii,i,j,k,t] for i in SI_V for j in SI_V if (i,j,k,t) in edges)
                            for key, change in insertionChange.items():
                                if change == thisChange:
                                    insertNodes[ii,k,t] = (key[1], key[2])
                            SI_B[ii,k,t] = sum(SI_X[i,j,k,t] * model._Cost[i,j] for i in SI_V for j in SI_V) + thisChange
                        except ValueError:
                            SI_B[ii,k,t] = 2*model._Cost[0,ii]
                            insertNodes[ii,k,t] = tuple([0])

            # (i,k,t): 1 if node i visited by vehicle k in time t
            SI_r = {(i,k,t): model.cbGetSolution(model._Y[i,k,t])
                for i in SI_V for k in SI_K for t in SI_T[1:]}

            # (i,k,t): 1 if customer is removed
            SI_u = {(i,k,t): mSI.addVar(vtype = GRB.BINARY) 
                for i in SI_V for k in SI_K for t in SI_T}
            # (i,k,t): 1 if customer is inserted
            SI_v = {(i,k,t): mSI.addVar(vtype = GRB.BINARY)
                for i in SI_V for k in SI_K for t in SI_T}
            
            SI_Q = {(i,p,k,t): mSI.addVar()
                for i in SI_V for p in SI_P for k in SI_K for t in SI_T[1:]}

            SI_I = {(i,p,t): mSI.addVar()
                for i in SI_V for p in SI_P for t in SI_T}

            mSI.setObjective(quicksum( model._Hold[i,p]*SI_I[i,p,t]
                        for i in SI_V for p in SI_P for t in SI_T) + 
                    quicksum( SI_B[i,k,t]*SI_v[i,k,t] - SI_A[i,k,t]*SI_u[i,k,t] 
                        for i in SI_V[1:] for k in SI_K for t in SI_T[1:]), GRB.MINIMIZE)

            # Constraints
            SI_initialInventory = {(i,p): mSI.addConstr(SI_I[i,p,0] == model._I0[i,p])
                        for i in SI_V for p in SI_P}

            SI_inventoryTransferSupplier = {(p,t): mSI.addConstr(
                    SI_I[0,p,t] == SI_I[0,p,t-1] + model._Production[p,t] 
                        - quicksum(SI_Q[i,p,k,t] for k in SI_K for i in SI_V[1:]))
                for p in SI_P for t in SI_T[1:]}

            SI_inventoryTransfer = {(i,p,t): mSI.addConstr(
                    SI_I[i,p,t] == SI_I[i,p,t-1] + quicksum(SI_Q[i,p,k,t] for k in SI_K) - model._Dem[i,p,t])
                for i in SI_V[1:] for p in SI_P for t in SI_T[1:]}

            SI_inventoryCapacity = {(i,t): mSI.addConstr(
                    quicksum(SI_I[i,p,t] for p in SI_P) <= model._Capacity[i])
                for i in SI_V[1:] for t in SI_T[1:]}

            # quicksum(SI_Q[i,p,k,t] for p in SI_P for k in SI_K) 
            SI_legalDelivery = {(i,t): mSI.addConstr(
                    quicksum(SI_Q[i,p,k,t] for p in SI_P for k in SI_K) <=
                        (model._Capacity[i] - quicksum(SI_I[i,p,t-1] for p in SI_P)))
                for i in SI_V[1:] for t in SI_T[1:]}

            SI_1 = {(i,p,k,t): mSI.addConstr(
                SI_Q[i,p,k,t] <= (SI_r[i,k,t] - SI_u[i,k,t] + SI_v[i,k,t])*model._Capacity[i])
            for i in SI_V[1:] for p in SI_P for k in SI_K for t in SI_T[1:]}

            SI_2 = {(i,k,t): mSI.addConstr(SI_v[i,k,t] <= 1 - SI_r[i,k,t])
                for i in SI_V[1:] for k in SI_K for t in SI_T[1:]}

            SI_3 = {(i,k,t): mSI.addConstr(SI_u[i,k,t] <= SI_r[i,k,t])
                for i in SI_V[1:] for k in SI_K for t in SI_T[1:]}

            # EPSILON
            epsilon = 1
            SI_4 = {(k,t): mSI.addConstr(
                quicksum(SI_u[i,k,t] for i in SI_V[1:]) + quicksum(SI_v[i,k,t] for i in SI_V[1:])
                    <= epsilon) for k in SI_K for t in SI_T[1:]}
            
            SI_5 = {(k,t): mSI.addConstr(
                quicksum(SI_Q[i,p,k,t] for i in SI_V[1:] for p in SI_P) <= int(model._VC[k]))
                for k in SI_K for t in SI_T[1:]}
            
            # Optimize!
            mSI.setParam('OutputFlag', 0)
            mSI.optimize()

            #PLOT BEFORE
            # colour = ['r','g','b','y','m']
            # labels = [str(i) for i in V]
            # for t in T[1:]:
            #     if sum(SI_X[i,j,k,t] for i in V for j in V for k in K) < 0.5:
            #         pass
            #     else:
            #         plt.figure(figsize=(15,10))
            #         plt.plot(X_Coords, Y_Coords, 'k', marker='o', linestyle='None')
            #         for k in K:
            #             #print("Vehicle", k)
            #             #plt.figure()
            #             for label, x, y in zip(labels, X_Coords, Y_Coords):
            #                 plt.annotate(label, xy = (x-20, y))
            #             for i in V:
            #                 #print(i, Y[i,k,t].x)
            #                 for j in V:
            #                     if SI_X[i,j,k,t] > 0.9:
            #                         plt.plot([X_Coords[i], X_Coords[j]],
            #                             [Y_Coords[i], Y_Coords[j]],
            #                             'k')
            #         plt.axis('off')
            #         plt.title("Day " + str(t) + "BEFORE")
            #         plt.show()

            model._newX = SI_X

            for ii in SI_V[1:]:
                for k in SI_K:
                    for t in SI_T[1:]:
                        if SI_u[ii,k,t].x > 0.5:
                            print("REMOVAL:")
                            print("Vertex {}, Vehicle {}, Day {}".format(ii, k, t))
                            if len(neighbours[ii,k,t]) == 2:
                                for j in SI_V:
                                    model._newX[ii,j,k,t] = 0
                                    model._newX[j,ii,k,t] = 0
                                model._newX[neighbours[ii,k,t][0], neighbours[ii,k,t][1], k, t] = 1
                                model._newX[neighbours[ii,k,t][1], neighbours[ii,k,t][0], k, t] = 1
                            else:
                                model._newX[0,ii,k,t] = 0
                                model._newX[ii,0,k,t] = 0

                        if SI_v[ii,k,t].x > 0.5:
                            print("INSERTION:")
                            print("Vertex {}, Vehicle {}, Day {}".format(ii, k, t))
                            if len(insertNodes[ii,k,t]) == 2:
                                (i, j) = insertNodes[ii,k,t]
                                model._newX[i,j,k,t] = 0
                                model._newX[j,i,k,t] = 0

                                model._newX[i,ii,k,t] = 1
                                model._newX[ii,i,k,t] = 1

                                model._newX[j,ii,k,t] = 1
                                model._newX[ii,j,k,t] = 1
                            else:
                                model._newX[0,ii,k,t] = 2
                                model._newX[ii,0,k,t] = 2

            model._newQ = {}
            model._newI = {}
            model._newY = {}
            for i in SI_V:
                for p in SI_P:
                    for k in SI_K:
                        for t in SI_T[1:]:
                            model._newQ[i,p,k,t] = SI_Q[i,p,k,t].x
            
            for i in SI_V:
                for p in SI_P:
                    for t in SI_T:
                        model._newI[i,p,t] = SI_I[i,p,t].x

            for i in SI_V:
                for k in SI_K:
                    for t in SI_T[1:]:
                        model._newY[i,k,t] = SI_Y[i,k,t]
                        if SI_u[i,k,t].x > 0.5:
                            #removed
                            model._newY[i,k,t] = 0
                        if SI_v[i,k,t].x > 0.5:
                            #inserted
                            model._newY[i,k,t] = 1

            # New variables!!!
            model._newVars = []
            for i in V:
                for j in V[i:]:
                    for k in K:
                        for t in T[1:]:
                            model._newVars.append(model._newX[i,j,k,t])

            for i in V:
                for k in K:
                    for t in T[1:]:
                        model._newVars.append(model._newY[i,k,t])

            for i in V:
                for p in P:
                    for k in K:
                        for t in T[1:]:
                            model._newVars.append(model._newQ[i,p,k,t])

            for i in V:
                for p in P:
                    for t in T:
                        model._newVars.append(model._newI[i,p,t])

            routing = sum( 
                model._Cost[i,j] * model._newX[i,j,k,t]
                for i in SI_V for j in SI_V if i < j for k in SI_K for t in SI_T[1:])

            holding = sum(
                model._Hold[i,p]*SI_I[i,p,t].x 
                for i in SI_V for p in SI_P for t in SI_T)

            print("*************** UPDATE:")
            print("Current Sol:", currentSol)
            print("SI Routing:", routing)
            print("SI Holding:", holding)
            print("SI Routing + Holding:", routing+holding)

            # Plot AFTER
            # colour = ['r','g','b','y','m']
            # labels = [str(i) for i in V]
            # for t in T[1:]:
            #     if sum(model._newX[i,j,k,t] for i in V for j in V for k in K) < 0.5:
            #         pass
            #     else:
            #         plt.figure()
            #         plt.plot(X_Coords, Y_Coords, 'k', marker='o', linestyle='None')
            #         for k in K:
            #             #print("Vehicle", k)
            #             #plt.figure()
            #             for label, x, y in zip(labels, X_Coords, Y_Coords):
            #                 plt.annotate(label, xy = (x-20, y))
            #             for i in V:
            #                 #print(i, Y[i,k,t].x)
            #                 for j in V:
            #                     if model._newX[i,j,k,t] > 0.9:
            #                         plt.plot([X_Coords[i], X_Coords[j]],
            #                             [Y_Coords[i], Y_Coords[j]],
            #                             colour[k] + ':')
            #         plt.axis('off')
            #         plt.title("Day " + str(t))
            #         plt.show()

            if routing+holding <= model._Best:
                # update solution vector?
                # print("\n BETTER \n")
                model._Best = routing+holding

        ## SUBTOUR ELIMINATION
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

    if where == GRB.Callback.MIPNODE and \
      model.cbGet(GRB.callback.MIPNODE_STATUS)==GRB.status.OPTIMAL and \
      model._Best < model.cbGet(GRB.Callback.MIPNODE_OBJBST) - 0.1:
        try:
            model.cbSetSolution(model.getVars(), model._newVars)
            model.cbUseSolution()
        except UnboundLocalError:
            pass
        except AttributeError:
            pass

### OPTIMIZE! ###

m.optimize(subtourelim_si)

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