"""  
** Core library to extract hierarchical skill network. 
** Version: 0.8
** Date: April. 16. 2018
** Python version >= 3.x || >= 2.6

"""
from __future__ import print_function
import networkx as nx
import itertools, math, pickle, sys, re 


def write_edges(f, G, **kwargs):
    """
    Input: 
        1. file object: f
        2. graph object: G
    Output: 
        None
    
    Write edges (u, v) (direction: u ->v) with hierarchical measure alpha
    """
    try:
        weight_col = kwargs['weight']
    except:
        weight_col = 'weight'

    for u, v in G.edges_iter():
        f.write("%s\t%s\t%f\n" % (str(u), str(v), float(G[u][v][weight_col])) )

    f.close()

def parsimony_network(G):
    G_parsi = G.copy()

    n_list = list(G_parsi.nodes())
    direct_path_list = []

    for node in n_list:
        nn_set = set(G_parsi.successors(node))
        bfs_reachable_nodes = set()

        for nn in list(nn_set):
            bfs_dic = nx.bfs_successors(G_parsi, nn)
            bfs_reachable_nodes = bfs_reachable_nodes | set( [n for s, nn in bfs_dic for n in nn ] )

        cur_direct_path = bfs_reachable_nodes & nn_set

        for nn in list(cur_direct_path):
            direct_path_list.append((node, nn))
    
    print("Remove %d direct paths in G" % len(direct_path_list))

    G_parsi.remove_edges_from(direct_path_list)

    return G_parsi


def projection_from_list(tag_deg_list, tag_co_list, item_size, **kwargs):

    try:
        z_score_thres = float(kwargs['zscore'])
    except:
        z_score_thres = 2.0

    nofnode_pair_axis = int(item_size)
    tag_deg_dic = {}
    cp = {}
    edge_dic_filtered = {}

    for tag, deg in tag_deg_list:
        tag_deg_dic[tag] = int(deg)

    totlen = len(tag_co_list)
    thres_percent = 1
    
    for tag_node, deg in tag_deg_list:
        edge_dic_filtered[tag_node] = []
    
    for idx, (u, v, co_occur) in enumerate(tag_co_list):
        len_v_nn = float(tag_deg_dic[v])
        len_u_nn = float(tag_deg_dic[u])
        common_nn = float(co_occur)

        average = len_v_nn * len_u_nn / nofnode_pair_axis
        deviation = math.sqrt(len_v_nn * len_u_nn * (nofnode_pair_axis - len_v_nn) * (nofnode_pair_axis - len_u_nn) \
                / (nofnode_pair_axis * nofnode_pair_axis * (nofnode_pair_axis-1) )  )

        if ((common_nn - average)/deviation) >= z_score_thres:
            cp[(u, v)] = common_nn / len_v_nn
            cp[(v, u)] = common_nn / len_u_nn
            edge_dic_filtered[u].append(v)
            edge_dic_filtered[v].append(u)
        else:
            #print(u, v, (common_nn - average)/deviation) 
            pass

        if( (idx / float(totlen)*100.) > thres_percent):
            print("Progress: %d%%" % (int)(idx / float(totlen) * 100.0), end="\r")
            thres_percent += 1

    print("= Finish removing noisy-like edges & calculating conditional probability.")

    return edge_dic_filtered, cp



def projection_from_con(con, **kwargs):
    """
    Input: 
        1. list of connections [ (u, v), (w, x), ... , (k, t)]
        2. keyword arguemnt: axis= 0 or 1
           the list "con" has two dimension for its elements.
           Axis determines the dimension of node layer which 
           the bipartite network will projected onto.

           For example, 
               with the keyword of axis = 1, 
               an projection network with nodes in second dimension of 
               connection list is generated.
         
    
    Output: 
        1. dictionary of edges {v: [u, w, x, ...], w: [v, a, c, ...] }
        2. dictionary of conditional probability
            {(u,v): value, (v,u):value, (x,y):value, ... }
            Note: conditional prob. is asymmetric for commutation of two nodes.
            cp[(u,v)] != cp[(v, u)]

    """
    if kwargs["axis"] == 0:
        axis = 0
        pair_axis = 1
    elif kwargs["axis"] == 1:
        axis = 1
        pair_axis = 0
    else:
        print("axis: 0 or 1")
    
    try:
        z_score_thres = float(kwargs['zscore'])
    except:
        z_score_thres = 2.0
        
    
    print("-Make a dictionary")
    ## make a dic
    proj_dic = {}
    for elem in con:
        if elem[pair_axis] in proj_dic.keys():
            proj_dic[elem[pair_axis]].append(elem[axis])
        else:
            proj_dic[elem[pair_axis]] = [elem[axis]]
    
    print("-Ready to project the network")
    print("-Read %d nodes in the pair-axis layer" % len(proj_dic.keys()))
    print("-Projecting the biparitite network")
    #make edge list
    
    edge_list = {}
    
    thres_percent = 1
    totlen = len(proj_dic.keys())
    cnt = 0
    
    for key in proj_dic.keys():
        for node in proj_dic[key]:
            if node in edge_list.keys():
                edge_list[node] = edge_list[node].union(proj_dic[key])
            else:
                edge_list[node] = set(proj_dic[key])
        
        cnt += 1
        if( (cnt / float(totlen)*100.) > thres_percent):
            print("Progress: %d%%" % (int)(cnt / float(totlen) * 100.0), end="\r")
            thres_percent += 1
    print("=Complete to make projection network")
    
    for node in edge_list.keys():
        edge_list[node] = list(edge_list[node] - set([node]))
            
    nofnode_pair_axis = len(proj_dic.keys())
    del proj_dic
    
    print("-Make a dictionary")
    proj_dic2 = {}
    for elem in con:
        if elem[axis] in proj_dic2.keys():
            proj_dic2[elem[axis]].append(elem[pair_axis])
        else:
            proj_dic2[elem[axis]] = [elem[pair_axis]]
    
    for key in proj_dic2.keys():
        proj_dic2[key] = set(proj_dic2[key])
        
    deg_dic = {}
    cp= {}
    edge_dic_filtered = {}
    
    for u in proj_dic2.keys():
        deg_dic[u] = 0
        edge_dic_filtered[u] = []
        
    print("-Ready to remove random edges (z-score thres: %f)" % z_score_thres)
    print("-Read %d nodes in the axis layer" % len(proj_dic2.keys()))
    print("-Calculating z-score and conditional prob. for all edges in the projection network")   
    thres_percent = 1
    totlen = len(edge_list.keys())
    cnt = 0        
    
    for u in edge_list.keys():
        for v in edge_list[u]:
            len_v_nn = float(len(proj_dic2[v]))
            len_u_nn = float(len(proj_dic2[u]))
            common_nn = float(len(proj_dic2[v] & proj_dic2[u]))
        
            #print(u, v, len_v_nn, len_u_nn, nofnode_pair_axis)
            average = len_v_nn * len_u_nn / nofnode_pair_axis
            deviation = math.sqrt(len_v_nn * len_u_nn * (nofnode_pair_axis - len_v_nn) * (nofnode_pair_axis - len_u_nn) \
                              / (nofnode_pair_axis * nofnode_pair_axis * (nofnode_pair_axis-1) )  )
        
            if deviation == 0 or ((common_nn - average)/deviation) >= z_score_thres:
                cp[(u, v)] = common_nn / len_v_nn
                cp[(v, u)] = common_nn / len_u_nn
                edge_dic_filtered[u].append(v)
                edge_dic_filtered[v].append(u)
            else:
                #print(u, v, (common_nn - average)/deviation) 
                pass
        
        cnt += 1
        
        if( (cnt / float(totlen)*100.) > thres_percent):
            print("Progress: %d%%" % (int)(cnt / float(totlen) * 100.0), end="\r")
            thres_percent += 1
            
    print("=Complete to calculate conditional probability and remove random edges")

    return edge_dic_filtered, cp
    

def extract_hierarchical_subgraph_from_edges(edge_dic, cp, thres, **kwargs):
    """
    Input: 
        1. edge_dic: dictionary of edges such as {u: [u, w, x, ...] }
        2. cp: dictionary of conditional probability obtained from the function
               of projection_from_con().
        3. thres: threshold value of alpha. 
                When you input 0 to thres, it generate the whole projection 
                network.

    Output:
        Networkx DiGraph for hierarchical subgraph 
    """
    max_k = max(list(len(v) for v in edge_dic.values()))
    
    thres_percent = 1
    totlen = len(edge_dic.keys())
    cnt = 0   

    if not 'output' in kwargs.keys():
        kwargs['output'] = "graph"
    
    if kwargs['output'] == "graph":
        print("-Option: return graph object")
        print("-Calculate the hierarchy measure, alpha")
        edges_bunch = []

        for u in edge_dic.keys():
            for v in edge_dic[u]:
                min_deg = float(min([len(edge_dic[u]), len(edge_dic[v])]))
                asym = abs(cp[(u,v)] - cp[(v, u)])
                alpha = min_deg * asym / max_k
        
                if alpha > thres:
                    if cp[(u,v)] > cp[(v, u)]: # root: u
                        edges_bunch.append([u, v, {'weight': alpha}])
                    elif cp[(u,v)] < cp[(v, u)]:  #root: v
                        edges_bunch.append([v, u, {'weight': alpha}])
                    else:
                        edges_bunch.append([u, v, {'weight': alpha}])
                        edges_bunch.append([v, u, {'weight': alpha}])
    
            cnt += 1
    
            if( (cnt / float(totlen)*100.) > thres_percent):
                print("Progress: %d%%" % (int)(cnt / float(totlen) * 100.0), end="\r")
                thres_percent += 1 
        
        print("=Complete the calculation of hierarhcy measure")
        print("-Making graph object using NetworkX")
        subG = nx.DiGraph()
        subG.add_edges_from(edges_bunch)
        print("=Complete!")
        
        return subG
    
    elif kwargs['output'] == 'file':
        try:
            fn = kwargs['fn']
        except KeyError:
            print("Input file name as keywords: fn=\'alpha.txt\'")
        
        
        fs = open(fn, 'w')
        fs.write("#Direction: u --> v\n")
        fs.write("#u\tv\talpha\n")

        print("-Option: write results in given file")
        print("-Calculate the hierarchy measure, alpha")
        #for u, v in edge_pair_list:
        for u in edge_dic.keys():
            for v in edge_dic[u]:
                min_deg = float(min([len(edge_dic[u]), len(edge_dic[v])]))
                asym = abs(cp[(u,v)] - cp[(v, u)])
                alpha = min_deg * asym / max_k
                 
                if alpha > thres:
                    if cp[(u,v)] > cp[(v, u)]: # root: u
                        fs.write("%s\t%s\t%f\n" % (u, v, alpha) )
                    elif cp[(u,v)] < cp[(v, u)]:  #root: v
                        fs.write("%s\t%s\t%f\n" % (v, u, alpha) )
                    else:
                        fs.write("%s\t%s\t%f\n" % (u, v, alpha) )
                        fs.write("%s\t%s\t%f\n" % (v, u, alpha) )
    
            cnt += 1
    
            if( (cnt / float(totlen)*100.) > thres_percent):
                print("Progress: %d%%" % (int)(cnt / float(totlen) * 100.0), end="\r")
                thres_percent += 1 
        
        print("Save the result to %s" % fn)
        fs.close()


def edges_file_to_graph(f):

    hbG = nx.DiGraph()
    
    cnt = 0
    thres_cnt = 100

    while True:
        line = f.readline()
        if not line: break
        if line[0] == "#": continue
        elem = re.split("\t|\n", line)
        hbG.add_edge(elem[0], elem[1], {'weight':float(elem[2])})
        cnt+=1

        if ( cnt / 100.) > thres_cnt:
            print("Read edges: %d" % cnt, end="\r")
            thres_cnt += 100

    
    return hbG


def bipartite_projection(B, target_node_list):
    """
    Input
        B: bipartite network
        target_node_list: a bunch of nodes divided to be parallelized 
        output: instance of Queue

    Output:
        returns edge list with weight 
        [[u, v, w], [u, v, w], ... ]

    Find connections from a bunch of nodes to make projection network
    It returns edge list with weight such as [u, v, weight]

    """
    r = []
    pred=B.adj
    totlen = len(target_node_list)
    cnt=0
    thres_percent = 1
    for u in target_node_list:
        unbrs = set(B[u])
        nbrs2 = set((n for nbr in unbrs for n in B[nbr])) - set([u])
        for v in nbrs2:
            vnbrs = set(pred[v])
            common = unbrs & vnbrs
            weight = len(common)
            r.append([u, v, weight])

        cnt+=1
        if( (cnt / float(totlen) * 100.) > thres_percent):
            print("Process: %d%%,  Used memory: %d" % ((int)(cnt / float(totlen) * 100.), len(r) ))
            thres_percent += 1
    
    r2 = [ [u, v, {'weight': w}] for u, v, w in r ]
    proG = nx.Graph()
    proG.add_edges_from(r2)

    return proG


def zscore(Q, Qi, Qj, Qij):
    """
    Calculate z-score from probability of co-occurance 
    between two any skills in the given biparaite network.
    
    Inputs:
        Q  : The number of users 
        Qi : The number of users who connect to i-th skill 
        Qj : The number of users who connect to j-th skill.
        Qij: The number of users who connect to both of 
            the i-th skill and the j-th skill.
    
    Return:
        It returns z-score value.
    """

    average = float(Qj)*float(Qi)/float(Q) 
    deviation = math.sqrt(Qi*Qj*(Q-Qi)*(Q-Qj) / float(Q*Q*(Q-1)) ) 

    return (Qij - average) / deviation


def find_communities_from_hierarchy(G, **kwargs):
    newman_ziff_dic = {}
    comm_dic = {}

    #### Check the prev_hier
    if "prev_hier" in kwargs.keys() and "new_edge" in kwargs.keys():
        prev_hier = kwargs["prev_hier"]
        new_edge = kwargs["new_edge"]
        
        newman_ziff_dic = prev_hier.copy()
        
        ## new_edge: (u, v)
        ## Direction: u -> v
        ## if u node is incoming in the graph G, find its root
        if not u in prev_hier.keys() and not v in prev_hier.keys():
           newman_ziff_dic[u] = -1 
           newman_ziff_dic[v] = u 
        
        elif u in prev_hier.keys() and not v in prev_hier.keys():
           newman_ziff_dic[v] = u 
        elif not u in prev_hier.keys() and v in prev_hier.keys():
            newman_ziff_dic[u] = -1 
            
            ## update v's parent
            max_weight = 0
            for parent in G.in_edges(v):
                if G[parent][v]['weight'] > max_weight:
                    parent_node = parent
            newman_ziff_dic[v] = parent_node
        
        else:
            ## update u's parent
            max_weight = 0
            for parent in G.in_edges(u):
                if G[parent][u]['weight'] > max_weight:
                    parent_node = parent
            newman_ziff_dic[u] = parent_node
        
            ## update v's parent
            max_weight = 0
            for parent in G.in_edges(v):
                if G[parent][v]['weight'] > max_weight:
                    parent_node = parent
            newman_ziff_dic[v] = parent_node

        ### find communities
        for u in newman_ziff_dic.keys():
            comm_dic[u] = find_root_hierarchy_from_newman_ziff(newman_ziff_dic, u)                   

    else:
        ## find the most probable parent node
        for u in G.nodes_iter():
            if G.in_degree(u) == 0:
                newman_ziff_dic[u] = -1
            else:
                max_weight = 0
                for parent, cur_u in G.in_edges(u):
                    if G[parent][cur_u]['weight'] > max_weight:
                        max_weight = G[parent][cur_u]['weight'] 
                        parent_node = parent
                newman_ziff_dic[u] = parent_node
        
        ### find communities
        for u in newman_ziff_dic.keys():
            comm_dic[u] = find_root_hierarchy_from_newman_ziff(newman_ziff_dic, u)                   

    return comm_dic, newman_ziff_dic

