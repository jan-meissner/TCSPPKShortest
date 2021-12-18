import time
from math import inf
import networkx as nx
import random
import numpy
from fibheap import *
from bisect import bisect_right, bisect_left

def custom_bisect_left(a, x, lo=0, hi=None, compare=None):
    while lo < hi:
        mid = (lo + hi) // 2
        if compare(a[mid], x):  # compare(a[mid],x): #if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo

class JITModifiedDijkstra:
    def __init__(self, graph, start):
        self.start = start
        self.graph = graph

        # init earliest arrival times
        self._t = {v: inf for v in list(nx.nodes(graph))}
        self._t[start] = 0

        # nodes not yet finalized - fibo heap
        self.S = makefheap()
        self.NodeDict = {v: Node((inf, v)) for v in list(nx.nodes(graph))}
        for node in self.NodeDict.values():
            self.S.insert(node)
        self.S.decrease_key(self.NodeDict[start], (self._t[start], start))

        # finalized nodes
        self.Fin = set()
        self.max_finalized_t = -inf

    def _cost(self, u, v):
        return self.graph.edges[u, v]["weight"]

    def _TS(self, u, v):
        return self.graph.edges[u, v]["TS"]

    def is_threshold_strictly_smaller_t_target(self, treshold, target):
        while target not in self.Fin:
            if treshold < self.max_finalized_t: return True
            self.dijkstraIteration()
        return treshold < self._t[target]

    def dijkstraIteration(self):
        curr_t, curr = self.S.extract_min().key
        self.Fin.add(curr)
        self.max_finalized_t = max(curr_t, self.max_finalized_t)
        for neighbor in self.graph.neighbors(curr):

            j = bisect_left(self._TS(curr, neighbor), self._t[curr])
            depart_time = self._TS(curr, neighbor)[j]
            path = depart_time + self._cost(curr, neighbor)

            if path < self._t[neighbor]:
                self._t[neighbor] = path
                self.S.decrease_key(self.NodeDict[neighbor], (self._t[neighbor], neighbor))

def DFSJIT_KShortest(graph, start, target, K):
    """
    O(|E| \log(|V| r) + K |E| \eta r \log r) implementation of Algorithm 3 with query based just-in-time dijkstra.
    """
    def cost(u, v):
        return graph.edges[u,v]["weight"]

    def TS(u,v):
        return graph.edges[u, v]["TS"]

    lazyDijk = JITModifiedDijkstra(graph, start)

    def R(u,v,j,pre):
        if u == start:
            KPaths.append([(u,j),] + pre)

        for pred in graph.predecessors(u):
            ipred_upper = bisect_right(TS(pred, u), TS(u, v)[j]-cost(pred, u))-1
            ipred_lower = custom_bisect_left(TS(pred, u),pred,hi=ipred_upper+1, compare=lazyDijk.is_threshold_strictly_smaller_t_target)
            for i in range(ipred_lower,ipred_upper+1):
                R(pred, u, i, [(u,j),] + pre)

                if len(KPaths) == K: return  #early stopping


    # gather candiates
    Rcanidates = [(v,j) for v in graph.predecessors(target) for j in range(len(TS(v,target)))]

    # sort canidates
    Rcanidates = sorted(Rcanidates, key=lambda x: TS(x[0],target)[x[1]] + cost(x[0],target))
    # enumerate paths
    KPaths = []
    for x in Rcanidates:
        if not lazyDijk.is_threshold_strictly_smaller_t_target(TS(x[0], target)[x[1]], x[0]):
            R(x[0], target, x[1], [])
        if len(KPaths) >= K: break

    return KPaths[:K], lazyDijk.Fin




def rndTSN(N = 30, p = 0., TS = 15):
    # Directed circle with additional random edges
    G = nx.DiGraph()
    Edges = []
    for u in range(N-1):
        Edges.append((u,u+1))
    Edges.append((N-1,0))

    n = numpy.random.binomial(N*N,p,1)[0]
    print("Graph generated with", n ,"edges")
    for _ in range(n):
        u = random.randint(0,N-1)
        v = random.randint(0,N-1)
        Edges.append((u, v))

    G.add_edges_from(Edges)
    for (u, v, w) in G.edges(data=True):
        w['weight'] = random.randint(1,5)
        w['TS'] = sorted(random.sample(range(1, TS*2), TS))
        #print(w['TS'])
    return G

def verify(graph,paths):
    print("Checking validity of paths...")
    def cost(u, v):
        return graph.edges[u,v]["weight"]

    def TS(u,v):
        return graph.edges[u, v]["TS"]

    def verifySingle(path):
        #check edges
        last = path[0][0]
        lasti = path[0][1]
        lastarrival = 0
        for (u,i) in path[1:]:
            if not graph.has_edge(last,u):
                raise Exception("Edge not in graph!")

            if lastarrival > TS(last, u)[lasti]:
                print("Time travel", last, u, TS(last, u)[lasti], cost(last, u))
                raise Exception("Time travel", last, u )
            lastarrival = TS(last, u)[lasti] + cost(last, u)
            last = u
            lasti = i

    for path in paths:
        verifySingle(path)

    print("All paths are valid!")

def printArrivalTimes(graph, paths):
    def cost(u, v):
        return graph.edges[u, v]["weight"]

    def TS(u, v):
        return graph.edges[u, v]["TS"]
    print("Arrival times:", [TS(p[-2][0], p[-1][0])[p[-2][1]]+ cost(p[-2][0], p[-1][0]) for p in paths])


if __name__ == '__main__':
    seed = int(time.time()*1000)%2**31
    print("random seed",seed)
    random.seed(seed)
    numpy.random.seed(seed)

    print("Generating graph...")
    G = rndTSN(5000,10/(5000^2),800) #
    print("Graph generated.")
    start = 1
    target = 1000
    K = 100

    startTime = time.time()
    paths, finalizedNodes = DFSJIT_KShortest(G,start,target,K)
    print("Done with DFSJIT", time.time() - startTime)
    print("Nodes finalized by DFSJIT", len(finalizedNodes))
    verify(G,paths)
    printArrivalTimes(G, [p + [(target, 0)] for p in paths])

