''' Created on 21/04/2022
Creates an object with the data to build linear models for bilevel programming in the metabolic engineering.
The object after its creation has three methods to solve the bilevel problem

 > methods :
    * looping - enumeration approach
    * milp - reformulating the problem as a single level milp
    * callbacks - a faster enumeration technique using callbacks

input:
rxn      - a list of reactions [j]
met      - a list of metabolites [i]
S        - a matrix [i][j] with the stoichimetric values
LB       - Lower bounds for the reactions
UB       - Upper bounds for the reactions
biomas   - index of the biomas (cellular growth)
chemical - index of the chemical of interes to optimize
KO       - a list of potential candidates for knockouts (index)
WT       - a list with the values of the FBA 
k        - number of reactions to knockout

Parameters with default value
BM         = 1000 -> big number for the milp reformulation
target     = .5   -> the precentage of requiered production of biomas (celluar growth)
infeas     = 1e-6 -> the parameter for infeasibilities
time_limit = 1000 -> limits the time the solver spends trying to find a solution

Author: @j-alex-vindel
'''

import gurobipy as gp
import MetModels
from gurobipy import GRB 
from itertools import combinations
import math as ma

class MetData:
    def __init__(self,BM=1000,time_limit=1000,infeas=1e-6,target=.5,**kwargs):
        
        self.infeas = infeas
        self.target = target
        self.rxn = kwargs['rxn']
        self.met = kwargs['met']
        self.S   = kwargs['S']
        self.LB  = kwargs['LB']
        self.UB  = kwargs['UB']
        self.k   = kwargs['k']
        self.KO  = kwargs['KO']
        self.biomas = kwargs['biomas']
        self.chemical = kwargs['chemical']

        self.minprod = kwargs['WT'][self.biomas]*self.target
        self.M = [j for j in range(len(self.rxn))]
        self.N = [i for i in range(len(self.met))]
        self.time_limit = time_limit
        self.BM = BM
        self.WT = MetModels.fba(self)
        #try to implement the master and follower here in the __init__ may be easier to have them built already.
        self.master = MetModels.build_master(self)
        self.follower = MetModels.build_follower(self)

    def __str__(self):
        return f"The model has {len(self.N)} metabolites (i) and {len(self.M)} reactions (j). The FBA_biomas: {round(self.WT[self.biomas],5)}.The FBA_chemical: {round(self.WT[self.chemical],5)}"

    def looping(self):
        LB_looping = self.LB.copy()
        
        LB_looping[self.biomas] = self.WT[self.biomas]*self.target
        master = self.master.copy()
        follower = self.follower.copy()
        # follower = MetModels.build_follower(self)
        # master = MetModels.build_master(self)
        UB_control = max(self.UB) + 1
        LB_control = -1*UB_control
        iterations = 0 
        m_time = 0
        f_time = 0
        LB_tracker = {}
        biomas_tracker = {}
        mv = [master.getVarByName('mv[%d]'%a) for a in self.M]
        my = [master.getVarByName('my[%d]'%a) for a in self.M]
        master.addConstrs((LB_looping[j]*my[j] <= mv[j] for j in self.M), name='LB')
        master.addConstrs((mv[j] <= my[j]*self.UB[j] for j in self.M), name='UB')
        master.update()
        k_limit = len(list(combinations(self.KO,self.k)))

        while abs((LB_control - UB_control)/LB_control) >= 1e-4:
            iterations += 1
            if iterations >= k_limit or LB_control>UB_control:
                print(f"-->> TERMINATED <<--")
                break
            master.Params.OptimalityTol = self.infeas
            master.Params.IntFeasTol = self.infeas 
            master.Params.FeasibilityTol = self.infeas
            master.Params.TimeLimit = self.time_limit
            master.optimize()
            m_time += master.Runtime
            if master.status in (GRB.OPTIMAL,GRB.TIME_LIMIT):
                ys = [master.getVarByName('my[%d]'%a).x if a in self.KO else 1 for a in self.M]
                voj = [master.getVarByName('mv[%d]'%a).x for a in self.M]

                vij,ftime = MetModels.optimize_follower(follower,LB_looping,self.UB,ys,self.M,self.infeas,self.biomas)
                f_time += ftime
                y_hat = [i for i,y_value in enumerate(ys) if y_value < 1e-6]
                y_dot = [b for b in self.KO if (abs(vij[b])< 1e-6) and (b not in y_hat)]
                py = list(combinations(y_dot,self.k))
                strats = list(combinations(y_hat,self.k))

                for strat in strats:
                    if vij[self.biomas] != 2000:
                        LB_tracker[strat] = vij[self.chemical]
                        biomas_tracker[strat] = vij[self.biomas]
                if abs(vij[self.biomas] - voj[self.biomas]) >= 1e-6:
                    master.addConstr(sum(my[index] for index in y_hat)>= self.k-1, name='cut_k[%d]'%iterations)

                    for i,comb in enumerate(py):
                        master.addConstr((vij[self.biomas] <= mv[self.biomas] + 
                        (ma.ceil(vij[self.biomas]*10)/10)* (sum(my[i]for i in comb))),name='vij_k_%s[%s]'%(iterations,comb))
                UB_control = voj[self.chemical]
                LB_control = max(LB_control,vij[self.chemical])
                master.update()
                master.reset(1)
                follower.reset(1)

            elif master.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
                MetModels.checkstatus(master.status)
                break
    
        print(f"Best Strategy: ({[self.rxn[index] for index in max(LB_tracker,key=LB_tracker.get)]} -> Max Value: {LB_tracker[max(LB_tracker,key=LB_tracker.get)]}")
        print(f'Biomas value: {biomas_tracker[max(LB_tracker,key=LB_tracker.get)]}')
        print(f"Total time: {m_time + f_time}")
        Total_time = m_time+f_time
        strat = [self.rxn[index] for index in max(LB_tracker,key=LB_tracker.get)]
        value = LB_tracker[max(LB_tracker,key=LB_tracker.get)]
        return strat,value,Total_time

    def milp(self):
        LB_milp = self.LB.copy()
        
        LB_milp[self.biomas] = self.WT[self.biomas]*self.target 
        master = self.master.copy()

        mv = [master.getVarByName('mv[%s]'%a) for a in self.M]
        my = [master.getVarByName('my[%s]'%a) for a in self.M]

        l = master.addVars(self.N,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='lambda')
        a1 = master.addVars(self.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
        b1 = master.addVars(self.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
        a2 = master.addVars(self.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
        b2 = master.addVars(self.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
        a = master.addVars(self.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha')
        b = master.addVars(self.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta')
        master.addConstrs((my[j] == 1 for j in self.M if j not in self.KO), name='y_essentials')
        # Dual objective

        master.addConstr((mv[self.biomas] == (sum(a1[j]*self.UB[j] - b1[j]*LB_milp[j] for j in self.M)
        + sum(a2[j]*self.UB[j] - b2[j]*LB_milp[j] for j in self.M))),name='dual-objective')

        # Dual constraints

        master.addConstrs((gp.quicksum(self.S.transpose()[i,j]*l[j] for j in self.N)
                - b[i]
                + a[i] - b2[i] + a2[i]
                == 0 for i in self.M if i !=self.biomas)
                ,name='S_dual')

        master.addConstr((gp.quicksum(self.S.transpose()[self.biomas,j]*l[j] for j in self.N)
                - b[self.biomas]
                + a[self.biomas]
                - b2[self.biomas] + a2[self.biomas] == 1), name='Sdual_t')

        # linearization

        master.addConstrs((a1[j] <= self.BM*my[j] for j in self.M),name='l1_a1')

        master.addConstrs((a1[j] >= -self.BM*my[j] for j in self.M),name='l2_a1')

        master.addConstrs((a1[j] <= a[j] + self.BM*(1-my[j]) for j in self.M),name='l3_a1')

        master.addConstrs((a1[j] >= a[j] - self.BM*(1-my[j]) for j in self.M),name='l4_a1')

        master.addConstrs((b1[j] <= self.BM*my[j] for j in self.M),name='l1_b1')

        master.addConstrs((b1[j] >= -self.BM*my[j] for j in self.M),name='l2_b1')

        master.addConstrs((b1[j] <= b[j] + self.BM*(1-my[j]) for j in self.M),name='l3_b1')

        master.addConstrs((b1[j] >= b[j] - self.BM*(1-my[j]) for j in self.M),name='l4_b1')

        # Bounds

        master.addConstrs((LB_milp[j]*my[j] <= mv[j] for j in self.M), name='LB')
        master.addConstrs((mv[j] <= self.UB[j]*my[j] for j in self.M), name='UB')

        master.addConstrs((LB_milp[j] <= mv[j] for j in self.M),name='lb')
        master.addConstrs((mv[j] <= self.UB[j] for j in self.M),name='ub')
        
        master.update()
        master.Params.OptimalityTol = self.infeas
        master.Params.IntFeasTol = self.infeas
        master.Params.FeasibilityTol = self.infeas
        master.Params.NodefileStart = 0.5
        master.optimize()
        s = master.Runtime
        del_strat = []
        if master.status == GRB.OPTIMAL:
            chem = master.getObjective().getValue()
            ys = [master.getVarByName('my[%d]'%j).x for j in self.M]
            vs = [master.getVarByName('mv[%d]'%j).x for j in self.M]
        print('*'*4,'SOLUTION','*'*4)
        print('Time (s):',s,sep=' -> ')
        print('Chemical Overproduction:',chem,sep=' -> ')
        print('Biomass production:',vs[self.biomas],sep=' -> ')
        print('**** Deletion Strategy: ****')
        for i in self.M:
            if ys[i] < .5:
                print('*'*2,i,self.rxn[i],sep=' -> ')
                del_strat.append(self.rxn[i])
        # print('*'*3,'Pool Solutions:','*'*3)
        # for e in range(nsolutions):
        #     m.setParam(GRB.Param.SolutionNumber,e)
        #     print('Succinate Overproduction:','%g'%m.PoolObjVal,sep=' -> ')
        #     print('Biomass Production:',v[biomas].Xn,sep=' -> ')
        #     if e <= nsolutions:
        #         m.setParam(GRB.Param.SolutionNumber,e)
        #         for i in M:
        #             if y[i].XN < .5:
        #                 print('**Knockout Strategy:', rxn[i],sep=' -> ')

        if master.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
            MetModels.checkstatus(master.status)
            ys = ['$' for i in self.M]
            vs = ['~' for i in self.M]
            print(f'Chemical: {vs[self.chemical]}')
            print(f'Biomass: {vs[self.biomas]}')
            del_strat = 'all'

        print('*'*4,' FINISHED!!! ','*'*4)

        return  del_strat, vs, s

    def callbacks(self):
        LB_callback = self.LB.copy()
        maxvalue = []
        minprod = self.WT[self.biomas]*self.target
        def inner(imodel,yoj):
            global vij
            imodel.setAttr('LB',imodel.getVars(),[LB_callback[j]*yoj[j] for j in self.M])
            imodel.setAttr('UB',imodel.getVars(),[self.UB[j]*yoj[j] for j in self.M])
            imodel.Params.OptimalityTol = self.infeas
            imodel.Params.IntFeasTol = self.infeas
            imodel.Params.FeasibilityTol = self.infeas
            imodel.optimize()

            if imodel.status == GRB.OPTIMAL:
                vij = [imodel.getVarByName('fv[%s]'%a).x for a in self.M]
            elif imodel.status in (GRB.INFEASIBLE, GRB.UNBOUNDED,GRB.INF_OR_UNBD):
                vij = [2000 if i == self.biomas else yoj[i] for i in self.M]
            return vij
        def lazycall(model,where):
            if where == GRB.Callback.MIPSOL:
                model._voj = model.cbGetSolution(model._vars)
                model._yoj = model.cbGetSolution(model._varsy)

                model._vij = inner(model._inner,model._yoj)

                knockset = [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]
                knockset_inner = [i for i,y in enumerate(model._vij) if abs(model._vij[i]) < 1e-6 and i in self.KO]
                ki = list(combinations(knockset_inner,self.k))
                # print('****Knockset Len****',len(ki))
                if len(knockset) !=self.k:
                    # print('Error knocking out')
                    return
                    #print('***','Begin Lazy Constraints','***')
                maxvalue.append(round(model._vij[self.biomas],5))
                if abs(model._vij[self.biomas] - model._voj[self.biomas]) >= 1e-6:
                    if model._vij[self.biomas] != 2000:
                        for i,comb in enumerate(ki):
                            # print(f'**** Lazy constrain {i}: {comb}*****')
                            model.cbLazy(round(max(maxvalue),5) <= model._vars[self.biomas] +
                                    (ma.ceil(model._vij[self.biomas]*10)/10) * (sum(model._varsy[f] for f in comb)))

                    else:
                        model.cbLazy(sum(model._varsy[j] for j in knockset) >= self.k-1)

                    # print('*** ENd Lazy Constraints ***')

            elif where == GRB.Callback.MIPNODE:
                #print('*** Begin Lazy CTR Callback (MIPNODE) ***')
                model._ryoj = model.cbGetNodeRel(model._varsy)
                for i,y in enumerate(model._ryoj):
                    if model._ryoj[y] >= 0.8:
                        model._ryoj[y] = 1.0
                    elif model._ryoj[y] <= 0.2:
                        model._ryoj[y] = 0.0
                    else:
                        model._ryoj[y] = 1.0
                #print('Rounded Solution', sum(model._ryoj.values()))
                # print('*** Deletion Strategy ***')
                if sum(model._ryoj.values()) == len(model._ryoj)-self.k:
                    model._vij = inner(model._inner,model._ryoj)

                    model.cbSetSolution(model._vars, model._vij)
                    model.cbSetSolution(model._varsy, model._ryoj)
                # print('*** Set Solution Passed ***', model.cbUseSolution())
        m = gp.Model()
        mv = m.addVars(self.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='mv')
        my = m.addVars(self.M,vtype=GRB.BINARY,name='my')

        m.setObjective(1*mv[self.chemical],GRB.MAXIMIZE)

        m.addConstrs((gp.quicksum(self.S[i,j]*mv[j] for j in self.M) == 0 for i in self.N),name='Stoichiometry')

        m.addConstrs((my[j] == 1 for j in self.M if j not in self.KO))

        m.addConstr(sum(1-my[j] for j in self.KO) == self.k, name='knapsack')

        m.addConstr(mv[self.biomas] >= minprod, name='target')

        m.addConstrs((LB_callback[j]*my[j] <= mv[j] for j in self.M),name='LB')
        m.addConstrs((mv[j] <= self.UB[j]*my[j] for j in self.M),name='UB')

        m._vars = mv
        m._varsy = my
        m.Params.lazyConstraints = 1

        imodel = gp.Model()
        fv = imodel.addVars(self.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
        imodel.params.LogToConsole = 0
        imodel.setObjective(2000*fv[self.biomas] + fv[self.chemical], GRB.MAXIMIZE)

        imodel.addConstrs((gp.quicksum(self.S[i,j]*fv[j] for j in self.M) == 0 for i in self.N),name='S2')
        imodel.addConstr(fv[self.biomas] >= minprod, name='target2')

        imodel.update()
        m._inner = imodel.copy()
        m._innerv = fv

        m.Params.OptimalityTol = self.infeas
        m.Params.IntFeasTol = self.infeas
        m.Params.FeasibilityTol = self.infeas
        # m.Params.NodefileStart = 0
        # m.Params.Threads = 4
        m.optimize(lazycall)
        # m.setParam(GRB.Param.PoolSolutions, 10)
        # m.setParam(GRB.Param.PoolSearchMode, 2)
        # m.setParam(GRB.Param.PoolGap, 0.01)
        s = m.Runtime
        # nsolutions = m.SolCount

        if m.status == GRB.OPTIMAL:
            ys = [m.getVarByName('my[%d]'%j).x for j in self.M]
            vs = [m.getVarByName('mv[%d]'%j).x for j in self.M]
            del_strat = [self.rxn[i] for i in self.M if ys[i] < .5]
        elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            ys = ['all' for i in self.M]
            vs = ['~' for i in self.M]
            del_strat = ['all']



        print('*** Best Solution ***')
        print('Biomass outer v:',vs[self.biomas],sep=' -> ')
        print('Biomass inner v:',vij[self.biomas],sep=' -> ')

        print('Chemical Overproduction:',vs[self.chemical],sep=' -> ')
        print('Deletion Strategy:',[self.rxn[i] for i in self.M if ys[i]<.5],sep=' -> ')

        print('Time in seconds: %d'%s,'Time in minutes: %d'%(s/60),sep=' -> ')

        return del_strat, vs, s

