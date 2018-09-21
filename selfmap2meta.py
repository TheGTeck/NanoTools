# coding=utf-8
#############################################################################
#                                                                           #
# Find the metareads from a read vs read mapping in paf format              #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#       Version: 0.2                                                        #
#       Author: Quentin Bonenfant                                           #
#               quentin.bonenfant@gmail.com                                 #
#############################################################################


from Read import Read
from Align import Align
import networkx as nx
import pysdam
import sys
import os

SEP = "\t"  # separator


def paf2graph(fileName):
    """ Convert a PAF file into a networkx graph and return the graph"""
    G = nx.Graph()
    print("Graphing")
    with open(fileName, "r") as file:
        for line in file.readlines():
            data = line.split("\t")

            # target and query data
            id1 = data[0]
            len1 = int(data[1])
            id2 = data[5]
            len2 = int(data[6])

            # Alignement data - Not used yet
            # alignLen = int(data[10])
            # matching = int(data[9])
            # identity = float(matching)/alignLen
            # coverage = float(alignLen) / max(len1,len2)

            # TODO: Implement CIGAR management

            G.add_node(id1, size=len1, orientation=None, relPos=None)
            G.add_node(id2, size=len2, orientation=None)
            G.add_edge(id1, id2, done=0,
                       relative=data[4],
                       mapping={id1: (data[2], data[3]),
                                id2: (data[7], data[8])},
                       direction=None)
    return(G)


def flatten(graph):
    """Return an array of array of overlaping reads,
    one array per overlap 'cluster' """

    print("Flattening")
    clusters = [[i for i in elem.nodes]
                for elem in list(nx.connected_component_subgraphs(graph))]
    return(clusters)


def getRelativePositions(graph):
    """ Figure out how reads of a subgraph are positionned,
    relatively to each other"""

    # /!\ This function needs a lot of cleaning / refactoring"

    print("Placing reads")

    # Setting first node as the + reference
    firstNode = list(graph.nodes)[0]
    graph.node[firstNode]['orientation'] = "+"
    nodes = [firstNode]

    relativPos = {n: None for n in nodes}  # relative position of nodes
    relativPos[nodes[0]] = 0     # Relative position from node 0 is 0
    graph.node[firstNode]['relPos'] = 0
    minNode = nodes[0]
    # Work until we run out of nodes
    while(nodes):
        newNodes = []  # Futur nodes to go through
        # For each node to explore
        for node in nodes:
            # Get the neighbourgs
            for ng in graph.adj[node]:

                edge = graph[ng][node]  # fetching edge data
                # Processing edges only if not already done
                if(not edge['done']):

                    edge['done'] = 1  # Marking the edge as done
                    # We will start again from this node in next iteration
                    newNodes.append(ng)
                    # Orienting the reads relative to each other
                    if edge['relative'] == "+":
                        if not graph.node[ng]['orientation']:
                            graph.node[ng]['orientation'] = "+" if graph.node[node]['orientation'] == "+" else "-"
                    elif edge['relative'] == "-":
                        if not graph.node[ng]['orientation']:
                            graph.node[ng]['orientation'] = "-" if graph.node[node]['orientation'] == "+" else "+"
                    else:  # In case there is a parsing error / incorrect PAF file or graph
                        print('Problem with the graph, aborting')
                        print(ng)
                        print(node)
                        raise ValueError("No relative position between reads")

                    # Fetching mapping regions for each nodes
                    s1, e1 = edge['mapping'][node]
                    s2, e2 = edge['mapping'][ng]

                    # getting the "real" start of the overlap for each read
                    rs1 = int(
                        s1) if graph.node[node]['orientation'] == '+' else graph.node[node]['size'] - int(e1)
                    rs2 = int(
                        s2) if graph.node[ng]['orientation'] == '+' else graph.node[ng]['size'] - int(e2)

                    # Calculating overlap offset
                    if (rs1 >= rs2):
                        edge['direction'] = (node, ng, rs1 - rs2)
                    else:
                        edge['direction'] = (ng, node, rs2 - rs1)

                    # storing relativ position (need to be optimized)
                    relPos = relativPos[node] + (rs1 - rs2)
                    relativPos[ng] = relPos
                    graph.node[ng]['relPos'] = relPos
                    minNode = minNode if graph.node[minNode]['relPos'] <= relPos else ng
        nodes = newNodes

    # Now we define each node position and mapping info relative to the left most node (minNode)
    data = []
    offset = graph.node[minNode]['relPos']
    length = 0
    for node in graph.nodes:
        start = graph.node[node]['relPos'] - offset
        nodeLen = graph.node[node]['size']
        # getting offsetted mapping of the read
        if graph.node[node]['orientation'] == "+":
            mapp = [(int(a) + start, int(b) + start)
                    for a, b in [graph[node][ng]['mapping'][node] for ng in graph.adj[node]]]
        else:
            mapp = [((nodeLen - int(b)) + start, (nodeLen - int(a)) + start)
                    for a, b in [graph[node][ng]['mapping'][node] for ng in graph.adj[node]]]
        data.append((node, start, start + nodeLen, mapp,
                     graph.node[node]['orientation']))
        length = max(length, start + nodeLen)

    return(data, length)
