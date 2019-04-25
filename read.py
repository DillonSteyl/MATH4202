import re

'''
These functions can be imported to implementation files so that data files can be read easily.
It is assumed that the data files are in a folder called "data" in the same location as the python file
calling the functions.
'''

def getParams(datString):
    if 'mmirp' in datString:
        datString = (datString.split('mmirp-'))[1] # remove mmirp from start
    
    if '.dat' in datString:
        datString = (datString.split('.dat'))[0] # remove.dat from end
    
    # Read paramaters directly from data string
    params = re.split('[-, .]', datString) # Split String
    customers = int(params[0])
    products = int(params[1])
    vehicles = int(params[2])
    period = int(params[3])
    filename = 'data/mmirp-' + datString + '.dat'
    
    return [customers, products, vehicles, period, filename]

def readFile(customers, products, vehicles, period, filename):    
    dat = open(filename, 'r')
    VC = re.split(' ', dat.readline())[-1] # Ignore the first line 
    
    # Coordinates of supply node
    coords = re.split('[ ]', dat.readline())
    coords = [float(x) for x in coords]
    
    X = [0 for i in range(customers+1)]
    Y = [0 for i in range(customers+1)]
    X[0] = coords[0]
    Y[0] = coords[1]

    T = range(period+1)
    
    I0 = {(i,p): 0 for i in range(customers+1) for p in range(products)}
    Production = {(p,t): 0 for p in range(products) for t in T[1:]}
    Hold = {(i,p): 0 for i in range(customers+1) for p in range(products)}
    Demand = {(i,p,t): 0 for i in range(customers+1) for p in range(products)
                    for t in T[1:]}
    Capacity = [0 for i in range(customers+1)]
    
    # Supplier Node Information
    initial = re.split(' ', dat.readline())
    
    for p in range(products):
        I0[0, p] = float(initial.pop(0))
        
    for p in range(products):
        for t in T[1:]:
            Production[p,t] = float(initial.pop(0))
            
    for p in range(products):
        Hold[0,p] = float(initial.pop(0))
    
    # Customer Information
    for c in range(customers):
        customer = re.split(' ', dat.readline())
        
        i = int(customer.pop(0)) # index
        X[i] = float(customer.pop(0))
        Y[i] = float(customer.pop(0))
        
        for p in range(products):
            I0[i,p] = float(customer.pop(0))
            
        for p in range(products):
            for t in T[1:]:
                Demand[i,p,t] = float(customer.pop(0))
        
        Capacity[i] = float(customer.pop(0))
        
        for p in range(products):
            Hold[i,p] = float(customer.pop(0))
    
    dat.close()
    
    return [VC,X,Y,I0,Production,Hold,Demand,Capacity]