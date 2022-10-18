import networkx as nx

import networks
import EnumKvcc
from networkx.algorithms.flow import maximum_flow
from networkx.algorithms.flow import minimum_cut
import matplotlib.pyplot as plt
from networkx.algorithms.connectivity import cuts
import time
from networkx import k_core


def k_vcc(g: nx.Graph, k):
    EnumKvcc.strong_side_vertices(g, k)
    return EnumKvcc.kvcc_enum(g, k)
def k_vcc_random(g: nx.Graph, k):
    # EnumKvcc.strong_side_vertices(g, k)
    return EnumKvcc.kvcc_enum(g, k)

class treeNode:
    def __init__(self,V):
        self.graph=V
        self.child=[]
def build_kVC_tree(g,k):
    if len(g)<=0:
        return None
    T = treeNode(list(g.nodes()))
    if k==8: print(list(g.edges()))
    print(k-1)
    print(len(g))
    if k>1:
        # VCCs=k_vcc_random(g,k)
        VCCs=k_vcc(g,k)
        # print(VCCs)
    else:
        VCCs=[list(g.nodes())]
    for x in VCCs:
        # print(len(x))
        T.child.append(build_kVC_tree(g.subgraph(x),k+1))
    return T
def levelorder(T):
    # g=nx.Graph()
    if T:
        que=[]
        old=[T]
        sizes=[]
        count=0
        sizeDic=[]
        while len(old)>0:
            cur=old.pop(0)

            print("%d level"%count)
            print(len(cur.graph))
            sizeDic.append(len(cur.graph))
            que.extend(cur.child)
            sizes.append(len(cur.child))
            if len(old)==0:
                print(sizes)
                print(sizeDic)
                count+=1
                print("next level")
                input("input something")
                old=que[:]
                que.clear()
                sizes.clear()
                sizeDic.clear()
    return


def levegraph(T):
    G=nx.Graph()
    labels=dict()
    if T:
        que = []
        old = [(T,0)]
        # count = 0
        cnt=0
        while len(old) > 0:
            cur,idx = old.pop(0)
            # G.add_node(idx,node_label=len(cur.graph))
            labels[idx]=len(cur.graph)
            # print(len(cur.graph))
            for ch in cur.child:
                cnt+=1
                que.append((ch,cnt))
                labels[cnt] = len(ch.graph)
                # G.add_node(cnt,node_label=len(ch.graph))
                G.add_edge(idx,cnt)



            if len(old) == 0:
                # print("next level")
                # input("input something")
                old = que[:]
                que.clear()
    return G,labels

if __name__ == '__main__':
    # geng=networks.zacharyclub()
    # geng = nx.read_edgelist("dataEpoch//ca-AstroPh.txt")
    # geng = nx.read_edgelist("dataEpoch//out.moreno_blogs_blogs.txt")
    # geng=nx.read_edgelist("dataEpoch//soc-hamsterster.edges.txt")
    # geng=networks.USpowergrid_n4941()
    # geng=networks.highschool()
    # geng=ne
    # geng=networks.SFHHconf()
    # geng=networks.CaNewman()
    # geng = nx.read_edgelist("dataEpoch//fb-pages-tvshow.txt",delimiter=',')
    # geng=networks.advogado()
    # geng=networks.primarySchool()
    # geng = nx.read_edgelist("dataEpoch//ia-reality.txt")
    # geng=networks.escorts()
    # geng = nx.read_edgelist("dataEpoch//ia-crime-moreno.txt")
    geng = nx.read_edgelist("dataEpoch//M1Anonymized.txt")
    namestr="M1.pdf"
    #ca grqc 406 subgroups
    # geng=networks.lesmis()
    # geng=networks.Email_eucore()
    # geng = nx.karate_club_graph()

    # geng = networks.dolphins()
    # nx.draw(geng, with_labels=True)
    # geng=networks.infectious()
    # geng=networks.USpowergrid_n4941()
    # plt.show()
    # geng = networks.CollegeMsgNetwork()
    # geng=networks.USairport_2010()
    # geng=networks.USairport500()
    geng.remove_edges_from(nx.selfloop_edges(geng))
    largest_cc = max(nx.connected_components(geng), key=len)
    geng=geng.subgraph(largest_cc)
    g = nx.Graph()
    edges = geng.edges()
    g.add_edges_from([(str(e[0]), str(e[1])) for e in edges])
    # print(len(g))
    # print(g.edges())
    print(nx.number_of_edges(g))
    # start=time.time()
    # vcc = k_vcc(g, k)
    # for x in vcc:
        # subg = g.subgraph(x)
        # print(len(x))
        # print(x)
        # nx.draw(subg, with_labels=True)
        # plt.show()
    # print(time.time()-start)
    # print(len(k_vcc(g,2,2500)))

    # quit()


    T=build_kVC_tree(g,1)
    # newg=levelorder(T)

    newg,labels=levegraph(T.child[0])

    import pydot
    from networkx.drawing.nx_pydot import graphviz_layout
    # H = nx.convert_node_labels_to_integers(newg, label_attribute='node_label')
    # H_layout = nx.nx_pydot.pydot_layout(H, prog='dot')
    # G_layout = {H.nodes[n]['node_label']: p for n, p in H_layout.items()}

    pos = graphviz_layout(newg, prog="dot",root=0)
    nx.draw(newg, pos=pos,with_labels=False,labels=labels,node_size=0.04,width=0.02)
    # plt.show()
    for node, (x, y) in pos.items():
        plt.text(x, y, labels[node], fontsize=3, ha='right', va='center')
    plt.show()
    plt.savefig(namestr,format="pdf",figsize=(200, 80), dpi=60)