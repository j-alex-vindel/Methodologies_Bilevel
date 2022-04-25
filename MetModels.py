''' Created on 21/04/2022
This module helps in the construction of linear models for the case of bilevel programming. 


-> build_master
   Builds the master problem to a gurobi model. It adds v,y variables.And, Stoichiometric and Knapsack constraints
 > input = object containing the data for the model
-> build_follower
   Builds the follower problem to a gurobi model. It adds v, variables and Stoichiometric constraints
 > input = object contaning the data for the model
-> optimize_follower
    Updates the LB and UB parameters from the inner model and optimizes it returning the inner values
 > input =  model (object), LB and UB parameters, a vector of y's values, the infeasibility and biomas index
-> checkstatus
   Checks the status when a model is infeasible or unbounded
 > input = the model status after it optimization

Author: @j-alex-vindel
'''

import gurobipy as gp

from gurobipy import GRB


def build_master(obj):

    m = gp.Model()
    mv = m.addVars(obj.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS, name='mv')
    my = m.addVars(obj.M,vtype=GRB.BINARY,name='my')

    m.setObjective(1*mv[obj.chemical],GRB.MAXIMIZE)

    m.addConstrs((gp.quicksum(obj.S[i,j]*mv[j] for j in obj.M) == 0 for i in obj.N),name='S1')
    m.addConstr(sum(1-my[j] for j in obj.KO) == obj.k, name='kpc')

    m.update()
    print(m)
    return m

def build_follower(obj):
    
    f = gp.Model()

    fv = f.addVars(obj.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')

    f.setObjective(2000*fv[obj.biomas] + fv[obj.chemical],GRB.MAXIMIZE)
    f.addConstrs((gp.quicksum(obj.S[i,j]*fv[j] for j in obj.M) == 0 for i in obj.N),name='S2')

    f.update()
    print(f)    
    return f

def optimize_follower(model,LB,UB,yhat,M,infeas,biomas):
    
    print(f'>> Optimizing Inner model')
    model.setAttr('LB',model.getVars(),[LB[j]*yhat[j] for j in M])
    model.setAttr('UB',model.getVars(),[UB[j]*yhat[j] for j in M])
    model.Params.OptimalityTol = infeas
    model.Params.IntFeasTol = infeas
    model.Params.FeasibilityTol = infeas
    model.optimize()
    if model.status == GRB.OPTIMAL:
        vij = [model.getVarByName('fv[%d]'%a).x for a in M]
    elif model.status in (GRB.INFEASIBLE, GRB.UNBOUNDED,GRB.INF_OR_UNBD):
        vij = [2000 if i == biomas else yhat[i] for i in M]
    print(f'>> Optimization completed!')
    return vij, model.Runtime

def checkstatus(modelstatus):
    if modelstatus == 4:
        text = 'Model was proven to be either infeasible or unbounded\n'
    elif modelstatus == 3:
        text = 'Model was proven to be infeasible\n'
    elif modelstatus == 5:
        text = 'Model was proven infeasible\n'
    return print(text)

def fba(obj):
    f = gp.Model()

    fv = f.addVars(obj.M,lb=obj.LB,ub=obj.UB,vtype=GRB.CONTINUOUS,name='fv')

    f.setObjective(2000*fv[obj.biomas] + fv[obj.chemical],GRB.MAXIMIZE)
    f.addConstrs((gp.quicksum(obj.S[i,j]*fv[j] for j in obj.M) == 0 for i in obj.N),name='S2')
    print(f'>> Calculating FBA')
    f.Params.OptimalityTol = obj.infeas
    f.Params.IntFeasTol = obj.infeas
    f.Params.FeasibilityTol = obj.infeas
    f.Params.OutputFlag = 0
    f.optimize()
    if f.status == GRB.OPTIMAL:
        vij = [f.getVarByName('fv[%d]'%a).x for a in obj.M]
    elif f.status in (GRB.INFEASIBLE, GRB.UNBOUNDED,GRB.INF_OR_UNBD):
        vij = [2000 if i == obj.biomas else 1 for i in obj.M]
    print(f'>> FBA Completed!')
    return vij
    
    
