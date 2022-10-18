import copy

import networkx as nx
import random
import math
import sys
import operations

from networkx.algorithms.connectivity.utils import build_auxiliary_node_connectivity, build_auxiliary_edge_connectivity
from networkx.algorithms.flow.shortestaugmentingpath import shortest_augmenting_path

sys.setrecursionlimit(10000)


def minimum_st_node_cut(G, s, t, flow_func=None, auxiliary=None, residual=None, cutoff=None):
    if auxiliary is None:
        H = build_auxiliary_node_connectivity(G)
    else:
        H = auxiliary

    mapping = H.graph.get("mapping", None)
    if mapping is None:
        raise nx.NetworkXError("Invalid auxiliary digraph.")
    kwargs = dict(flow_func=flow_func, residual=residual, auxiliary=H, cutoff=cutoff)

    # The edge cut in the auxiliary digraph corresponds to the node cut in the
    # original graph.
    edge_cut = minimum_st_edge_cut(H, f"{mapping[s]}B", f"{mapping[t]}A", **kwargs)
    if edge_cut is None: return None
    # Each node in the original graph maps to two nodes of the auxiliary graph
    node_cut = {H.nodes[node]["id"] for edge in edge_cut for node in edge}
    return node_cut - {s, t}


def minimum_st_edge_cut(G, s, t, flow_func=None, auxiliary=None, residual=None, cutoff=None):
    if flow_func is None:
        flow_func = shortest_augmenting_path

    if auxiliary is None:
        H = build_auxiliary_edge_connectivity(G)
    else:
        H = auxiliary

    kwargs = dict(capacity="capacity", flow_func=flow_func, residual=residual, cutoff=cutoff)

    cut_value, partition = minimum_cut(H, s, t, **kwargs)
    # print(cut_value)
    if cut_value is None: return None
    reachable, non_reachable = partition
    # Any edge in the original graph linking the two sets in the
    # partition is part of the edge cutset
    cutset = set()
    for u, nbrs in ((n, G[n]) for n in reachable):
        cutset.update((u, v) for v in nbrs if v in non_reachable)

    return cutset


def minimum_cut(flowG, _s, _t, capacity="capacity", flow_func=None, **kwargs):
    if flow_func is None:
        flow_func = shortest_augmenting_path

    if not callable(flow_func):
        raise nx.NetworkXError("flow_func has to be callable.")

    R = flow_func(flowG, _s, _t, capacity=capacity, value_only=True, **kwargs)
    if R.graph["flow_value"] == kwargs.get("cutoff"):
        return (None, None)
    # Remove saturated edges from the residual network
    cutset = [(u, v, d) for u, v, d in R.edges(data=True) if d["flow"] == d["capacity"]]
    R.remove_edges_from(cutset)

    # Then, reachable and non reachable nodes from source in the
    # residual network form the node partition that defines
    # the minimum cut.
    non_reachable = set(dict(nx.shortest_path_length(R, target=_t)))
    partition = (set(flowG) - non_reachable, non_reachable)
    # Finally add again cutset edges to the residual network to make
    # sure that it is reusable.
    if cutset is not None:
        R.add_edges_from(cutset)
    return (R.graph["flow_value"], partition)


def FastCheck(g1, k):
    # G = operations.sparse_certificate_random(g1, k)
    G=g1
    # G=copy.deepcopy(g1)
    # edges=list(G.edges())

    # vertices=list(G.nodes())
    # degrees = dict(nx.degree(G))
    # for v, deg in degrees.items():
    #    if deg < k:
    #       # print("neighbors")
    #       return set(nx.neighbors(G, v))
    if k==2:
        cur=set(g1.nodes())
        for i in cur:
            g = copy.deepcopy(G)
            g.remove_node(i)
            if nx.is_connected(g):
                continue
            else:
                cut={i}
                del g
                return cut
        del g
        return None



    H = build_auxiliary_node_connectivity(G)
    mapping = H.graph.get("mapping", None)
    g = directed_flow_graph(G)

    # print("larger end")
    m = G.number_of_edges()

    for _ in range(28 * k):  # 56k
        e, f = random.sample(list(G.edges()), 2)
        if random.random() < 0.5:
            s = e[1]
        else:
            s = e[0]
        if random.random() < 0.5:
            t = f[1]
        else:
            t = f[0]
        if G.has_edge(s, t) or G.has_edge(t, s):
            continue

        edge_cut = minimum_st_edge_cut(H, f"{mapping[s]}B", f"{mapping[t]}A", auxiliary=H, cutoff=k)
        if edge_cut is None:
            cut = None
        else:
            node_cut = {H.nodes[node]["id"] for edge in edge_cut for node in edge}
            cut = node_cut - {s, t}
        if cut is not None and len(cut) < k:
            # print("28k:",cut)
            del G
            # del degrees
            del H
            del g
            # del adjg
            del mapping
            return cut
        # print(_)
        # if cut:print(len(cut))
    print("large end")
    adjg = nx.to_dict_of_lists(g)
    for i in range(int(math.log2(m / (7 * k)) + 1),0,-1):
        print(i)
        # if i <= 2:
        #     for _ in range(1):
        #         for s in list(G.nodes()):
        #             res1 = LocalEC(adjg, str(s) + '#', math.pow(2, i), k)
        #             if res1:
        #                 del G
        #                 # del degrees
        #                 del H
        #                 del g
        #                 del adjg
        #                 del mapping
        #                 return dir2Vertex(res1)
        # else:
        for _ in range(int(m / math.pow(2, i))):
            # print("inner", _)
            e = random.choice(list(G.edges()))
            s = e[1]
            res1 = LocalEC(adjg, str(s) + '#', math.pow(2, i), k)
            if res1:
                return dir2Vertex(res1)
            s = e[0]
            res1 = LocalEC(adjg, str(s) + '#', math.pow(2, i), k)
            if res1:
                del G
                # del degrees
                del H
                del g
                del adjg
                del mapping
                return dir2Vertex(res1)
    print("small end")
    del G
    # del degrees
    del H
    del g
    del adjg
    del mapping
    return None


def dir2Vertex(visited):
    a = set()
    for i in visited:
        if type(i) is str:
            ele = int(i.strip('#'))
        else:
            ele = i
        a.add(ele)
        if ele in a:
            a.remove(ele)
        else:
            a.add(ele)
    return a


def LocalEC_old(adjg, x, mu, k):
    """
    :param x: x\in V
    :param mu: target volume mu>k
    :param k: target cut-size k>=1
    # :param gamma: slack gamma<=k
    :return:
    """
    marked = set()
    count = 0
    for _ in range(k):
        visited = set()
        path = []
        T = dfs(adjg, x, marked, count, visited, path, mu, k)
        if T is False: return False
        if T is None:
            return dir2Vertex(visited)
        else:
            for e in T[:-1]:
                marked.add((e[1], e[0]))
                adjg[e[1]].append(e[0])
                adjg[e[0]].remove(e[1])
    return False


def LocalEC(adjg, x, mu, k):
    """
    :param x: x\in V
    :param mu: target volume mu>k
    :param k: target cut-size k>=1
    # :param gamma: slack gamma<=k
    :return:
    """
    # marked=set()
    # count=0
    for _ in range(k):
        # visited=set()
        Tree = nx.DiGraph()
        # print("begin dfs")
        dfs_star(adjg, x, Tree, mu, k)
        # if T is False: return False
        if len(Tree) < 8 * mu * k: return dir2Vertex(list(Tree.nodes()))
        # print(len(Tree))
        rdme = random.choice(list(Tree.edges()))
        rdmy = rdme[0]
        while len(list(Tree.successors(rdmy))) > 0:
            suc = list(Tree.successors(rdmy))[0]
            adjg[rdmy].append(suc)
            adjg[suc].remove(rdmy)
            rdmy = suc

    return False


def dfs(adj_list, start, marked, markedlen, visited, path, mu, k):
    visited.add(start)
    for neighbour in adj_list[start]:
        if neighbour not in visited:
            e = (start, neighbour)
            path.append(e)
            if e not in marked:
                marked.add(e)
                markedlen += 1
                if markedlen > 128 * mu * k: return path
                if random.random() <= 1.0 / (8 * mu):
                    return path
            result = dfs(adj_list, neighbour, marked, markedlen, visited, path, mu, k)
            if result is False: return False
            if result is not None:
                return result
            path.pop()
    return None


def dfs_star(adj_list, start, Tree: nx.DiGraph(), mu, k):
    if len(Tree) > 8 * mu * k:
        return Tree
    for neighbour in adj_list[start]:
        if Tree.has_node(neighbour):
            # e=(start,neighbour)
            Tree.add_edge(neighbour, start)
            result = dfs_star(adj_list, neighbour, Tree, mu, k)
            if result is not None:
                return result
    return None


def directed_flow_graph(g: nx.Graph) -> nx.DiGraph:
    """
    transform the original graph G into a directed flow graph G'.
    :param g:original graph G
    :return:a directed flow graph G'
    """
    flow_graph = nx.DiGraph()
    for node in g.nodes:
        # if node is x: continue
        flow_graph.add_edge(node, str(node) + '#', capacity=1)
    for edge in g.edges:
        # if edge[0] is x:
        #     flow_graph.add_edge(str(edge[0]), edge[1], capacity=1)
        # else:
        flow_graph.add_edge(str(edge[0]) + '#', edge[1], capacity=1)
        flow_graph.add_edge(str(edge[1]) + '#', edge[0], capacity=1)
    return flow_graph


if __name__ == '__main__':
    import networks as ns

    # G=nx.path_graph(500)
    # G.add_edge(0,499)
    nx.gomory_hu_tree()
    G = nx.karate_club_graph()
    # G=ns.facebook_combined()   # connectivity=
    # print(nx.node_connectivity(G))
    # g=directed_flow_graph(G)

    # gg=nx.to_dict_of_lists(g)
    # print(gg)
    # print(gg['2069#'])
    print(FastCheck(G, 2))

    # print(nx.node_connectivity(G))
    # print(dfs(gg,'6#',set(),0,set(),[],10,3))
    # print(LocalEC(gg,'6#',10,3))
