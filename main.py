import time

import networkx as nx
import operations
# from networkx.algorithms.flow import maximum_flow
# from networkx.algorithms.flow import minimum_cut
# import matplotlib.pyplot as plt
from networkx.algorithms.connectivity import cuts
# import time
# from networkx import k_core
from networkx.algorithms.flow.shortestaugmentingpath import shortest_augmenting_path
import networks as ns
import network_trade
import copy
import collections
import MinCut_Sample


def Expand_Single_Node(k, lmd, G, D, Cprime):
    C = set(copy.deepcopy(Cprime))
    W = []
    for v in set(G.nodes()).difference(set(C)):
        if len(set(nx.neighbors(G, v)).intersection(set(C))) >= k:
            W.append(v)
    # print("num of expandable nodes",len(W))
    while len(W) > 0 and len(C) < lmd:
        u = W.pop()
        C.add(u)
        if len(D[len(C)]) == 0 or D[len(C)][0] < k:
            D[len(C)] = (min(k, len(C) - 1), C)
        # D[len(C)].append((min(k,len(C)-1),C))
        if len(W) == 0:
            for v in set(G.nodes()).difference(set(C)):
                if len(set(nx.neighbors(G, v)).intersection(set(C))) >= k:
                    W.append(v)
    return C


def min_node_cut(G, k, s, t, H=None):
    if H is None:
        H = MinCut_Sample.build_auxiliary_node_connectivity(G)
    mapping = H.graph.get("mapping", None)
    if G.has_edge(s, t) or G.has_edge(t, s):
        return None

    edge_cut = MinCut_Sample.minimum_st_edge_cut(H, f"{mapping[s]}B", f"{mapping[t]}A", auxiliary=H, cutoff=k)
    if edge_cut is None:
        cut = None
    else:
        node_cut = {H.nodes[node]["id"] for edge in edge_cut for node in edge}
        cut = node_cut - {s, t}
    if cut is not None and len(cut) < k:
        return cut

def update_D(D, C, k):
    if len(D[len(C)]) == 0 or D[len(C)][0] < k:
        D[len(C)] = (min(k, len(C) - 1), C)


def Remove_bunch_of_nodes(k, lmd, og, D, W, R):
    oldW = list(W)
    resW = set(oldW[len(og) - lmd:])
    oldW = set(W) - resW
    # maxC = set()
    while len(oldW) > 0:
        # start=time.time()
        g = nx.Graph(nx.subgraph(og, resW))
        VCCs = operations.k_vcc_exact(g, k, lmd)
        print("vccs done", len(VCCs))
        if len(VCCs) > 0:
            maxC = VCCs[0]
            for C in VCCs:
                update_D(D, C, k)
                if len(C) <= lmd:
                    extendedC = Expand_Single_Node(k, lmd, og, D, C)
                else:
                    extendedC = C

                print("%d nodes expand, after expand %d" % (len(extendedC) - len(C), len(extendedC)))
                if len(extendedC) > len(maxC): maxC = extendedC

                if len(extendedC) == lmd:
                    R.append(extendedC)
                    print("find R by expanding with residual w=", len(oldW))
                    return R
                update_D(D, extendedC, k)
            if len(maxC) >= lmd:
                resW = maxC
                break
            resW.update(extendedC)
            oldW = set(oldW) - set(extendedC)
        if len(oldW) == 0: break
        cur = oldW.pop()
        resW.add(cur)
    print("original %d after delete %d" % (len(W), len(resW)))
    print("remove single Node", len(resW) - lmd)
    maxg = nx.Graph(nx.subgraph(og, resW))
    gprime = Remove_Single_Node_exact(k, lmd, maxg, D)

    # print("remove Single Node:", end - start)
    update_D(D, set(gprime.nodes()), k)
    # R=Remove_Group_of_Nodes_exact(k,lmd,gprime,D,set(gprime.nodes()),R)
    if len(gprime) <= lmd:
        R.append(set(gprime.nodes()))
    return R


def discard(L, VCCs):
    R = []
    for x in range(len(L)):
        flag = False
        for y in range(x + 1, len(L)):
            if len(set(VCCs[L[x]]) - set(VCCs[L[y]])) == 0:  # x in y
                flag = True
                break
        if not flag: R.append(x)
    return [L[x] for x in range(len(L)) if x in R]


def solver(G, lmd, D, Delta=None):
    n = len(G)
    # print(lmd)
    # to get the upper bound
    upperBound = 1
    # if len(kSize)>0 and max(kSize.values())>=lmd:
    #     for k in range(len(G),0,-1):
    #         if kSize[k]>=lmd:
    #             upperBound=k
    #             break
    # else:
    for i in range(1, n):
        subg = nx.Graph(nx.k_core(G, i))
        if len(subg) < lmd:
            subg.clear()
            upperBound = i - 1
            break
    Delta = 1

    VCCs = []
    g = copy.deepcopy(G)
    if upperBound == 1: return (1, None, Delta)
    for k in range(upperBound + 1, 1, -1):
        print("k=", k)
        if len(kSize[k]) > 0:
            VCCs.append(kSize[k])
        else:
            VCCs.extend(operations.k_vcc_exact(g, k, lmd))
        print("main vccs done", len(VCCs))
        if len(VCCs) > 0:
            maxC = VCCs[0]
            Cprime = None
            for C in VCCs:
                print("VCC size=", len(C))
                update_D(D, C, k)
                # if len(D[len(C)])==0 or D[len(C)][0]<k:
                #     D[len(C)]=(min(k,len(C)-1),C)
                if len(C) > len(maxC):
                    maxC = C
                if len(C) == lmd:
                    Cprime = C
            kSize[k] = set(maxC)
            if len(maxC) >= lmd and Delta == 1:
                Delta = k
            if Cprime:
                print("1st output")
                return (min(k, len(Cprime) - 1), Cprime, min(Delta, len(Cprime) - 1))

            L = [x for x in range(len(VCCs)) if len(VCCs[x]) <= lmd]
            H = [x for x in range(len(VCCs)) if len(VCCs[x]) > lmd]
            L = discard(L, VCCs)

            toDelidx = []
            toAdd = []
            for Cidx in L:
                Cprime = Expand_Single_Node(k, lmd, G, D, VCCs[Cidx])
                print("%d nodes expand, after expand %d" % (len(Cprime) - len(VCCs[Cidx]), len(Cprime)))

                if len(Cprime) == lmd:
                    print("2nd output")
                    return (min(k, len(Cprime) - 1), Cprime, min(Delta, len(Cprime) - 1))
                toDelidx.append(Cidx)
                toAdd.append(Cprime)
                # print(len(Cprime))
                # print("expand done")
            for Cidx in H:
                subg = nx.Graph(nx.subgraph(G, VCCs[Cidx]))
                R = []
                R = Remove_bunch_of_nodes(k, lmd, subg, D, VCCs[Cidx], R)
                print("remove bunch of node done")
                # newg=Remove_Single_Node_exact(k,lmd,subg,D)
                # if len(newg) == lmd:
                #     print("3-1 output")
                #     return (min(k, len(newg) - 1), set(newg.nodes()), min(Delta, len(newg) - 1))
                # R = []
                # R=Remove_Group_of_Nodes_exact(k,lmd,newg,D,VCCs[Cidx],R)
                for C in R:
                    if len(C) == lmd:
                        print("3-2rd output")
                        return (min(k, len(C) - 1), C, min(Delta, len(C) - 1))
                for C in R:
                    if len(C) < lmd:
                        Cprime = Expand_Single_Node(k, lmd, G, D, C)
                        if len(Cprime) == lmd:
                            print("3-3nd output")
                            return (min(k, len(Cprime) - 1), Cprime, min(Delta, len(Cprime) - 1))
                toDelidx.append(Cidx)
                toAdd.extend(R)
            VCCs = [VCCs[x] for x in range(len(VCCs)) if x not in toDelidx]
            VCCs.extend(toAdd)
    print("4th output")
    return (1, None, Delta)
if __name__ == '__main__':

    # geng=nx.karate_club_graph()
    # geng=ns.dolphins()
    # geng=ns.CollegeMsgNetwork()
    # geng=ns.election_Data()
    # geng=ns.MySmall()
    # geng=ns.Email_eucore()
    # geng=ns.facebook_combined()
    # geng=ns.USpowergrid_n4941()
    # geng=ns.USairport500()
    # geng=ns.USairport_2010()
    # geng=ns.celegans_n306()
    # geng=ns.escorts()
    # geng=ns.primarySchool()
    # geng=ns.highschool()
    # geng=ns.advogado()
    geng=ns.zacharyclub()
    # geng=ns.enron()
    # geng=ns.FacebookWOSN()
    # print(len(geng))
    # geng=nx.read_edgelist("dataEpoch//soc-hamsterster.edges.txt")
    # geng = nx.read_edgelist("dataEpoch//M1Anonymized.txt")
    # geng=nx.read_edgelist("dataEpoch//out.subelj_euroroad_euroroad.txt")
    # geng = nx.read_edgelist("dataEpoch//soc-buzznet.txt")
    geng.remove_edges_from(nx.selfloop_edges(geng))


    # geng=ns.onionNetwork1000()
    # geng=ns.richClub1000()d
    # print(geng.nodes())
    # g = nx.Graph()
    g = geng
    # edges = geng.edges()
    # g.add_edges_from([(str(e[0]), str(e[1])) for e in edges])
    # nodes=list(g.nodes())
    # mxg=nx.relabel_nodes(g,mapping)
    # print(len(geng))
    n = len(g)
    print(len(g))
    print(g.number_of_edges())
    # t=15   #2: 17  20    3:
    kSize = collections.defaultdict(set)
    de =  n // 2 // 50
    # de=1
    # before 259350 is 1
    # print(g.nodes())
    for t in range(1,n//2):
        # if t<n//2: continue
        lmd = n - t
        print("=====================================")
        print("size of cur g", len(g))
        print("size", lmd)
        if len(g)<lmd:
            print("skip cur input")
            continue
        D = collections.defaultdict(tuple)
        start = time.time()
        ansk, ansC, ansdelta = solver(g, lmd, D)
        print(ansC)
        end = time.time()
        print("time spend:", end - start)

        print("output: k=%d, optimal<=%d" % (ansk, ansdelta))
        if ansC: print("size of the output is", len(ansC))
        # Delta=ansdelta
        if ansk == ansdelta:
            print("optimal solution")
            # g = nx.Graph(nx.subgraph(g, ansC))
            # print(D[max(D.keys())])
            if ansdelta > 1:
                g = nx.Graph(nx.subgraph(g, D[max(D.keys())][1]))
            else:
                g = nx.Graph(nx.subgraph(g, D[max(D.keys())][1]))

            D.clear()
            continue
        for i in range(lmd, n):
            if len(D[i]) > 0:
                # for j in range(len(D[i])):
                if D[i][0] > ansk and D[i][0] <= ansdelta:
                    print("not optimal, a close size is %d" % i, "with vertex connectivity %d"%D[i][0])
                    ansk = D[i][0]
                if ansk == ansdelta: break


