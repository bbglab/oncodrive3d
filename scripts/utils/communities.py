"""
Contains function to initialize a network for communities detection.
"""

import networkx as nx
import numpy as np


def get_network(nodes, mut_count_v, cmap):
    """
    Generate a network with significative pos as nodes 
    and ratio of shared mutation (Jaccard score) as edges. 
    """
    
    G = nx.Graph()
    G.add_nodes_from(nodes)

    # Iterate through each hit
    for i, ipos in enumerate(nodes):
        for j, jpos in enumerate(nodes):
            if i > j:
                
                # Add an edge if they are in contact
                ix, jx = ipos-1, jpos-1
                if cmap[ix, jx] == 1:

                    # Get the common res and their mut
                    neigh_vec_i, neigh_vec_j = cmap[ix], cmap[jx]
                    common_neigh = neigh_vec_j * neigh_vec_i
                    num_mut = np.dot(common_neigh, mut_count_v)

                    # Get the sum of the union of the mut
                    all_neigh = (neigh_vec_i + neigh_vec_j != 0).astype(int)
                    union_num_mut = np.dot(all_neigh, mut_count_v)

                    # Compute the Jaccard score or avg ratio between ij shared mut    
                    jaccard = np.round(num_mut/union_num_mut, 3)
                    G.add_edge(ipos, jpos, weight = jaccard)
    return G


def get_community_index_nx(pos_hits, communities):
    
    """
    Parse the labels returned by communities detection algorithms from NetworkX.
    """
    
    communities_mapper = {}
    for ic, c in enumerate(communities):
        for p in c:
            communities_mapper[p] = ic

    return pos_hits.map(communities_mapper).values