
import numpy as np
import networkx as nx
from string import ascii_uppercase
import itertools

SQRT_2 = np.sqrt(2)


def distance(x, x1, y, y1):
    return np.sqrt((x - x1) ** 2.0 + (y - y1) ** 2.0)


def product_gen(n):
    for r in itertools.count(1):
        for i in itertools.product(n, repeat=r):
            yield "".join(i)


def grapher(fil_arr):
    '''
    Experimental function to directly create a graph without labeling end,
    body, inter pts
    '''

    total_pix = fil_arr.sum()
    if total_pix == 0:
        raise ValueError("fil_arr contains no valid points.")
    # You probably shouldn't care about a one pixel skeleton anyways...
    elif total_pix == 1:
        G = nx.Graph()
        G.add_node(1)
        return G

    numbered_arr = fil_arr.astype(int)

    # Assign each pixel its own number
    for i, (y, x) in enumerate(zip(*np.where(fil_arr))):
        numbered_arr[y, x] = i + 1

    yshape, xshape = fil_arr.shape
    # Now we loop through the pixels and make a list of its neighbours
    G = nx.Graph()
    # Add all of the nodes to the graph
    for i in range(1, total_pix + 1):
        G.add_node(i)

    for i, (y, x) in enumerate(zip(*np.where(fil_arr))):
        slicer = (slice(max(0, y - 1), min(yshape, y + 2)),
                  slice(max(0, x - 1), min(xshape, x + 2)))
        slice_arr = numbered_arr[slicer].ravel()
        # Find neighbours that are not itself
        neighb_idx = np.where(np.logical_and(slice_arr != 0,
                                             slice_arr != i + 1))[0]

        # It is only disconnected if it is the ONLY pixel in the skeleton
        # Which means you probably shouldn't care about it, but still...
        if neighb_idx.size == 0:
            raise ValueError("Found an unconnected pixel. fil_arr must be a "
                             "set of 8-connected pixels.")

        for idx in neighb_idx:
            # Calculate the distance between the points
            y1, x1 = np.where(numbered_arr == slice_arr[idx])
            dist = distance(x, x1, y, y1)

            if dist > SQRT_2:
                raise ValueError("The distance between any two connected "
                                 "pixels cannot be larger than sqrt(2)")

            G.add_edge(i + 1, slice_arr[idx], weight=float(dist))

    # First look for corners with extra connection and remove those first
    num_conns = np.array([len(conns) for conns in G.adjacency_list()])
    gt_3_nodes = np.where(num_conns >= 3)[0] + 1
    for node in gt_3_nodes:
        G = remove_doubletriangle(node, G)

    num_conns = np.array([len(conns) for conns in G.adjacency_list()])
    gt_3_nodes = np.where(num_conns >= 3)[0] + 1
    for node in gt_3_nodes:
        G = remove_triangle(node, G)

    # Now we iterate through all of the nodes and merge those with 2 edges
    # into its neighbours.
    # Continue until none remain
    while True:
        num_conns = np.array([len(conns) for conns in G.adjacency_list()])
        two_nodes = np.where(num_conns == 2)[0] + 1
        if two_nodes.size == 0:
            break
        for node in two_nodes:
            G = merge_nodes(node, G)

    return label_graph(G)


def is_4conn(node_a, node_b, G):
    if G[node_a][node_b]['weight'] == 1.:
        return True
    return False


def is_8conn(node_a, node_b, G):
    if G[node_a][node_b]['weight'] == SQRT_2:
        return True
    return False


def remove_triangle(node, G):
    '''
    Remove cases where are corner is both 4 and 8 connected.
    Removes the 8-connection.
    '''

    # The node must have three connections
    conns = G.adjacency_list()[node - 1]
    if len(conns) != 3:
        return G

    # Must be at least one 8-conn and one 4-conn
    four_conn = []
    eight_conn = []
    for conn in conns:
        if is_4conn(node, conn, G):
            four_conn.append(conn)
        elif is_8conn(node, conn, G):
            eight_conn.append(conn)
        else:
            continue

    # Must have at least one of each
    if len(four_conn) < 1 or len(eight_conn) < 1:
        return G

    # Now check if one of the 4-conn are connected to the 8-conn
    for fconn in four_conn:
        for econn in eight_conn:
            if any(econn == G[fconn].keys()):
                G.remove_edge(node, econn)
                removal = True
                break
            else:
                removal = False
        if removal:
            break

    return G


def remove_doubletriangle(node, G):
    '''
    Remove cases where are corner is both 4 and 8 connected.
    Removes both the 8-connections.
    '''

    # The node must have three connections
    conns = G.adjacency_list()[node - 1]
    if len(conns) > 3:
        return G

    # Must be at least two 8-conn and one 4-conn
    four_conn = []
    eight_conn = []
    for conn in conns:
        if is_4conn(node, conn, G):
            four_conn.append(conn)
        elif is_8conn(node, conn, G):
            eight_conn.append(conn)
        else:
            continue

    # Must have at least one of each
    if len(four_conn) < 1 or len(eight_conn) < 2:
        return G

    # Now check if one of the 4-conn are connected to two 8-conn
    removal = []
    for fconn in four_conn:
        print(fconn)
        for econn in eight_conn:
            print(econn)
            if econn in G[fconn].keys():
                removal.append(econn)

    if len(removal) == 2:
        for econn in removal:
            G.remove_edge(node, econn)

    return G


def merge_nodes(node, G):
    '''
    Combine a node into its neighbors.
    '''

    neigb = G[node].keys()

    if len(neigb) != 2:
        return G

    new_weight = G[node][neigb[0]]['weight'] + \
        G[node][neigb[1]]['weight']

    G.remove_node(node)
    G.add_edge(neigb[0], neigb[1], weight=new_weight)

    return G


def label_graph(G):
    '''
    Letters for intersections and numbers for ends
    '''

    nodes = np.array(G.nodes())

    num_conns = np.array([len(conns) for conns in G.adjacency_list()])

    intersecs = nodes[np.where(num_conns > 1)]

    inter_labels = {inter: let for let, inter in zip(product_gen(ascii_uppercase), intersecs)}

    ends = nodes[np.where(num_conns == 1)]

    end_labels = {end: i + 1 for i, end in enumerate(ends)}

    mapping = {}
    for node in intersecs:
        mapping[node] = inter_labels[node]
    for node in ends:
        mapping[node] = end_labels[node]

    new_G = nx.relabel_nodes(G, mapping)

    return new_G
