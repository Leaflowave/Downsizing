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


def Remove_Single_Node(k, lmd, maxg, D):
    W = set(maxg.nodes())

    n = len(W)
    if n <= lmd:
        return maxg
    sorted(W, key=lambda x: nx.degree(maxg, x))
    for _ in range(n):
        flag = False
        S = set()
        for u in set(W.intersection(set(maxg.nodes()))):
            if Single_Removable_exact(maxg, u, k, W):
                S.update(set(nx.neighbors(maxg, u)))
                maxg.remove_node(u)
                if len(D[len(maxg)]) == 0 or D[len(maxg)][0] < k:
                    D[len(maxg)] = (min(k, len(maxg) - 1), set(maxg.nodes()))
                # D[len(maxg)].append((min(k,len(maxg)-1),set(maxg.nodes())))
                flag = True
                if len(maxg) == lmd:
                    return maxg
        W = S
        print("one pass")
        if not flag: break
    return maxg


def Remove_Single_Node_exact(k, lmd, maxg, D):
    W = set(maxg.nodes())
    n = len(W)
    if n <= lmd:
        return maxg
    sorted(W, key=lambda x: nx.degree(maxg, x))
    for _ in range(n):
        flag = False
        S = set()

        while len(W) > 0:
            u = W.pop()
            if Single_Removable_exact(maxg, u, k, W):
                print("remove 1 node")
                S.update(set(nx.neighbors(maxg, u)) - {u})
                maxg.remove_node(u)

                if len(D[len(maxg)]) == 0 or D[len(maxg)][0] < k:
                    D[len(maxg)] = (min(k, len(maxg) - 1), set(maxg.nodes()))
                # D[len(maxg)].append((min(k,len(maxg)-1),set(maxg.nodes())))
                flag = True
                if len(maxg) == lmd:
                    return maxg
            # print("residual #of nodes to delete",len(W))
        W = S.intersection(set(maxg.nodes()))
        print("one pass")
        if not flag: break
    return maxg


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


def Single_Removable(og, u, k):
    """
    Test whether the residual graph after deleting node u is still k-VC.

    we first sort the neighbors according to their degrees (ascending order).   O(nlogn)
    Next we compute a sparse certificate of the residual graph.  O(m+n)
    Finally we check if residual graph is k-VC.

    The checking method utilize Even1975 method with we modify the step 3 by compute only the neighbors of u,
    since if u0 is k-VC to all neighbors of u then finnally u0=u is k-VC to all nodes.

    the total number of calls of FF alg is at most k*(k-1)/2+delta-k.
    Each call cost O(km) time.
    Thus, O(nlogn+m+n+(delta-k+k^2/2)(mk))=O(mkdelta+mk^3)
    :param og: a graph
    :param u: current node u
    :param k: connectivity requirement
    :return: boolean, True if the residual graph after deleing u is still k-VC
    """
    # Nu=list(nx.neighbors(og,u))
    # degDict = dict(og.degree())
    # Sort nodes by degree.
    # nodes = sorted(degrees, key=degrees.get)
    # degNu={x:degDict[x] for x in Nu}
    # Nu=sorted(degNu,key=degNu.get)
    # Nu=[x[1] for x in degNu]
    # n=len(Nu)
    g = copy.deepcopy(og)
    g.remove_node(u)
    sparseg = operations.sparse_certificate(g, k)
    return MinCut_Sample.FastCheck(sparseg, k) is None


def Single_Removable_exact(og, u, k, W):
    Nu = list(nx.neighbors(og, u))
    if len(Nu) < k:
        print("wrong-----------------------------")
        return False
    degNu = dict(nx.degree(og, Nu))
    sorted(Nu, key=lambda x: degNu[x])
    if degNu[Nu[0]] <= k:
        print("rule 1")
        W.difference_update(set(nx.neighbors(og, Nu[0])))
        return False  # rule 1: if a neighbor's degree is k then u undeletable.
    # if degNu[0][0]==
    # if len(set(nx.neighbors(og, Nu[0])).intersection(Nu)) == len(Nu) - 1:
    #     print("rule 2")
    #     return True  # rule 2: if u's neighborhood forms a clique, then deletable
    n = len(Nu)
    g = copy.deepcopy(og)
    g.remove_node(u)
    sparseg, sidegroup = operations.sparse_certificate(g, k)
    # sparseg=g
    # flowg=operations.directed_flow_graph(sparseg)
    flowg = MinCut_Sample.build_auxiliary_node_connectivity(sparseg)

    kneighbors = collections.defaultdict(set)  # rule 3: if x is k-vc to y's k neighbors, then x is k-vc to y.
    for i in range(n):
        kneighbors[Nu[i]].update(set(nx.neighbors(g, Nu[i])))
    for s in range(k - 1):
        for t in range(s + 1, k):
            if g.has_edge(Nu[s], Nu[t]):
                kneighbors[Nu[t]].add(Nu[s])
                kneighbors[Nu[s]].add(Nu[t])
                continue
            if len(kneighbors[Nu[s]].intersection(kneighbors[Nu[t]])) >= k:
                kneighbors[Nu[t]].add(Nu[s])
                kneighbors[Nu[s]].add(Nu[t])
                continue
            cut = min_node_cut(g, k, Nu[s], Nu[t], flowg)
            # if nx.maximum_flow_value(flowg, Nu[s] + "#", Nu[t], flow_func=shortest_augmenting_path, cutoff=k) < k:
            if cut is not None and len(cut) < k:
                print("rule 3")
                # print(cut)
                W.difference_update(set(cut))
                # print("left #of nodes to try",len(W))
                return False
            kneighbors[Nu[t]].add(Nu[s])
            kneighbors[Nu[s]].add(Nu[t])
    u0 = "u0"
    sparsegj = copy.deepcopy(sparseg)
    for x in range(k):
        sparsegj.add_edge(u0, Nu[x])
        kneighbors[u0].add(Nu[x])
    # sparsegj.add_edges_from([(u0, Nu[x]) for x in range(k)])
    for j in range(k - 1, n - 1):
        # sparsegj.add_edges_from([(u0, Nu[x]) for x in range(k,j+1)])
        # kneighbors[u0].update(set([Nu[x] for x in range(j+1)]))
        if len(kneighbors[u0].intersection(kneighbors[Nu[j + 1]])) >= k:
            sparsegj.add_edge(u0, Nu[j + 1])
            kneighbors[u0].add(Nu[j + 1])
            continue
        # flowg = operations.directed_flow_graph(sparsegj)
        cut = min_node_cut(sparsegj, k, u0, Nu[j + 1])
        # if nx.maximum_flow_value(flowg, u0 + "#", Nu[j+1], flow_func=shortest_augmenting_path, cutoff=k) < k:
        if cut is not None and len(cut) < k:
            print("rule 4")
            W.difference_update(set(cut))
            return False
        sparsegj.add_edge(u0, Nu[j + 1])
        kneighbors[u0].add(Nu[j + 1])
    print("rule 5")
    return True


def Remove_Group_of_Nodes_exact(k, lmd, g, D, W: set, R):
    if len(g) <= lmd: return R
    tmp = set(W).intersection(set(g.nodes()))
    while len(tmp) > 0:
        v = list(tmp)[0]
        W.remove(v)
        g.remove_node(v)
        VCCs = operations.k_vcc_exact(g, k, lmd)
        if len(VCCs) > 0:
            maxC = VCCs[0]
            for C in VCCs:
                update_D(D, C, k)
                if len(C) > len(maxC): maxC = C
            if len(maxC) <= lmd:
                R.extend(VCCs)
                return R
            maxg = nx.Graph(nx.subgraph(g, maxC))
            print("remove single node")
            gprime = Remove_Single_Node_exact(k, lmd, maxg, D)
            if len(D[len(gprime)]) == 0 or D[len(gprime)][0] < k:
                D[len(gprime)] = (min(k, len(gprime) - 1), set(gprime.nodes()))
            # D[len(gprime)].append((min(k,len(gprime)-1),set(gprime.nodes())))
            if len(gprime) <= lmd:
                R.append(set(gprime.nodes()))
                return R
            else:
                print("remove group node")
                return Remove_Group_of_Nodes_exact(k, lmd, gprime, D, W, R)
        tmp = set(W).intersection(set(g.nodes()))
    return R


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


def Remove_bunch_of_nodes_random(k, lmd, og, D, W, R):
    oldW = list(W)
    resW = set(oldW[len(og) - lmd:])
    oldW = set(W) - resW
    # maxC = set()
    while len(oldW) > 0:
        # start=time.time()
        g = nx.Graph(nx.subgraph(og, resW))
        # VCCs, _ = operations.k_vcc(g, k, lmd)
        VCCs= operations.k_vcc_exact(g, k, lmd)
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
    gprime = Remove_Single_Node(k, lmd, maxg, D)

    # print("remove Single Node:", end - start)
    update_D(D, set(gprime.nodes()), k)
    # R=Remove_Group_of_Nodes_exact(k,lmd,gprime,D,set(gprime.nodes()),R)
    if len(gprime) <= lmd:
        R.append(set(gprime.nodes()))
    return R


def Remove_Group_of_Nodes(k, lmd, g, D, W: set, R):
    tmp = set(W).intersection(set(g.nodes()))
    while len(tmp) > 0:
        v = list(tmp)[0]
        W.remove(v)
        g.remove_node(v)
        VCCs = operations.k_vcc(g, k, lmd)
        if len(VCCs) > 0:
            maxC = VCCs[0]
            for C in VCCs:
                if len(D[len(C)]) == 0 or D[len(C)][0] < k:
                    D[len(C)] = (min(k, len(C) - 1), C)

                # D[len(C)].append((min(k,len(C)-1),C))
                if len(C) > len(maxC): maxC = C
            if len(maxC) <= lmd:
                R.extend(VCCs)
                return R
            maxg = nx.Graph(nx.subgraph(g, maxC))
            print("remove single node")
            gprime = Remove_Single_Node(k, lmd, maxg, D)
            if len(D[len(gprime)]) == 0 or D[len(gprime)][0] < k:
                D[len(gprime)] = (min(k, len(gprime) - 1), set(gprime.nodes()))

            # D[len(gprime)].append((min(k,len(gprime)-1),set(gprime.nodes())))
            if len(gprime) <= lmd:
                R.append(set(gprime.nodes()))
                return R
            else:
                print("remove group node")
                return Remove_Group_of_Nodes(k, lmd, gprime, D, W, R)
        tmp = set(W).intersection(set(g.nodes()))
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


def solver_random(G, lmd, D, Delta=None):
    n = len(G)
    upperBound = 1
    for i in range(1, n):
        subg = nx.Graph(nx.k_core(G, i))
        if len(subg) < lmd:
            subg.clear()
            upperBound = i - 1
            break
    Delta = 1
    if upperBound < 1:
        update_D(D,set(list(G.nodes())),upperBound)
        return (1, None, Delta)

    # g = copy.deepcopy(G)
    for k in range(upperBound + 1, 1, -1):
        g = operations.sparse_certificate_random(G, k)
        VCCs = []
        print("k=", k)
        if len(kSize[k]) > 0:
            VCCs.append(kSize[k])
        else:
            vcc= operations.k_vcc_exact(g, k, lmd)
            # vcc, canc = operations.k_vcc(g, k, lmd)
            VCCs.extend(vcc)
            # if len(vcc) == 0 and k == 2:
            #     update_D(D, canc, 1)
            #     break
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
                    return (min(k, len(Cprime) - 1), Cprime, min(k, len(Cprime) - 1))
                toDelidx.append(Cidx)
                toAdd.append(Cprime)
                # print(len(Cprime))
                # print("expand done")
            for Cidx in H:
                subg = nx.Graph(nx.subgraph(G, VCCs[Cidx]))
                R = []

                R = Remove_bunch_of_nodes_random(k, lmd, subg, D, VCCs[Cidx], R)
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


