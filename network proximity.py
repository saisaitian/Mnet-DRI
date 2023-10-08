# -*- coding: utf-8 -*-

import os
import random
import numpy as np
import numpy.ma as ma
import networkx as nx
from matplotlib.pyplot import plot


# ========================= Distance ========================= #
def Distance(SD, modFrom, modTo=None):
    if modTo is None:
        return SD[np.ix_(modFrom, modFrom)].min(1).mean()
    else:
        return SD[np.ix_(modFrom, modTo)].min(1).mean()


def Shortest(SD, mod1, mod2):
    return SD[np.ix_(mod1, mod2)].mean()


def Closest(SD, mod1, mod2):
    sub = SD[np.ix_(mod1, mod2)]
    closest1 = sub.min(0)
    closest2 = sub.min(1)
    return (closest1.sum() + closest2.sum()) / (closest1.count() + closest2.count())


def Separation(SD, mod1, mod2):
    dAA = Distance(SD, mod1)
    dBB = Distance(SD, mod2)
    dAB = Closest(SD, mod1, mod2)
    return dAB - (dAA + dBB) / 2


def Center(SD, mod1, mod2):
    c1 = GetCenterNodes(SD, mod1)
    c2 = GetCenterNodes(SD, mod2)
    return SD[np.ix_(c1, c2)].mean()


def GetCenterNodes(SD, mod):
    sum_ = SD[np.ix_(mod, mod)].sum(1)
    return [mod[idx] for idx, isCenter in enumerate(sum_ == sum_.min()) if isCenter]  # bool(numpy.ma.core.MaskedConstant) => False


def Kernel(SD, mod1, mod2):
    sub = SD[np.ix_(mod1, mod2)]
    exp = ma.exp(-sub)
    alb = ma.log(exp.sum(1))
    bla = ma.log(exp.sum(0))
    nA = alb.count()
    nB = bla.count()
    return (nA * np.log(nB) + nB * np.log(nA) - alb.sum() - bla.sum()) / (nA + nB) + 1


# ========================= Network ========================== #
class Network(object):
    DISTANCE_MEASURE = {"DISTANCE"  : Distance,
                        "SHORTEST"  : Shortest,
                        "CLOSEST"   : Closest,
                        "SEPARATION": Separation,
                        "CENTER"    : Center,
                        "KERNEL"    : Kernel}

    def __init__(self, pathG, pathSD):
        self.G = nx.read_adjlist(pathG)
        self.SD = np.load(pathSD, allow_pickle=True)
        self.nodes = sorted(self.G.nodes())
        self.i2n = {index: node for index, node in enumerate(self.nodes)}
        self.n2i = {node: index for index, node in enumerate(self.nodes)}
        self.i2d = {}  # index: degree
        self.d2i = {}  # degree: [index1, index2, ...]
        for node, degree in self.G.degree():
            index = self.n2i[node]
            if degree in self.d2i:
                self.d2i[degree].append(index)
            else:
                self.d2i[degree] = [index]
            self.i2d[index] = degree
        self.dmin = min(self.d2i.keys())
        self.dmax = max(self.d2i.keys())
        self.d2b = {}  # degree: bin
        self.b2i = {}  # bin: [index1, index2, ...]

    def Name2Index(self, names, skipUnknown=True):
        if skipUnknown:
            return [self.n2i[n] for n in names if n in self.n2i]
        else:
            return [self.n2i[n] for n in names]

    def Index2Name(self, indexes, skipUnknown=True):
        if skipUnknown:
            return [self.i2n[i] for i in indexes if i in self.i2n]
        else:
            return [self.i2n[i] for i in indexes]

    # ------------------------------------------
    def PrepareBins(self, binSize=150):
        index = 0
        self.d2b = {}
        self.b2i[index] = []
        degrees = sorted(self.d2i.keys())
        for curr, next in zip(degrees, degrees[1:] + [degrees[-1] + 1]):
            for d in range(curr, next):
                self.d2b[d] = index
            self.b2i[index].extend(self.d2i[curr])
            if curr != degrees[-1] and len(self.b2i[index]) >= binSize:
                index += 1
                self.b2i[index] = []
        if len(self.b2i[index]) < binSize and index > 0:  # Merge last two bins if last bin < binSize
            for d in range(degrees[-1], -1, -1):
                if self.d2b[d] != index:
                    break
                self.d2b[d] = index - 1
            self.b2i[index - 1].extend(self.b2i[index])
            del self.b2i[index]

    def DegreeMimicSampling(self, indexes):
        binCount = {}
        for index in indexes:
            d = self.i2d[index]
            if d < self.dmin:
                d = self.dmin
            elif d > self.dmax:
                d = self.dmax
            b = self.d2b[d]
            if b in binCount:
                binCount[b] += 1
            else:
                binCount[b] = 1
        while True:  # TODO ValueError
            yield sum([random.sample(self.b2i[b], binCount[b]) for b in sorted(binCount.keys())], [])

    def TotalRandomSampling(self, n):
        indexes = list(range(len(self.nodes)))
        while True:
            yield random.sample(indexes, n)

    # ------------------------------------------
    def Proximity(self, mod1, mod2, method="DISTANCE"):
        return self.DISTANCE_MEASURE[method](self.SD, mod1, mod2)

    def ProximityRandom(self, mod1, mod2, method="DISTANCE", repeat=1000, seed=None):
        if seed is not None:
            random.seed(seed)
        method = self.DISTANCE_MEASURE[method]
        result = np.zeros(repeat)
        index = 0
        for mod1r, mod2r in zip(self.DegreeMimicSampling(mod1), self.DegreeMimicSampling(mod2)):
            v = method(self.SD, mod1r, mod2r)
            if not ma.is_masked(v):
                result[index] = v
                index += 1
                if index == repeat:
                    break
        return result

    def LCC(self, mod):
        return len(max(nx.connected_components(self.G.subgraph(self.Index2Name(mod))), key=len))

    def LCCRandom(self, mod, repeat=1000, seed=None):
        if seed is not None:
            random.seed(seed)
        result = np.zeros(repeat)
        index = 0
        for modr in self.DegreeMimicSampling(mod):
            result[index] = self.LCC(modr)
            index += 1
            if index == repeat:
                break
        return result


# ========================= Utility ========================== #
def Tsv2Adj(pathIn, pathOut, encoding="utf-8", removeSelfloop=True):
    G = nx.Graph()
    with open(pathIn, encoding=encoding) as fi:
        for line in fi:
            if not line.startswith("#"):
                G.add_edge(*line.strip().split("\t")[0:2])
    if removeSelfloop:
        G.remove_edges_from(nx.selfloop_edges(G))
    nx.write_adjlist(G, pathOut)
    return G


def Z_Score(real, background):
    m = background.mean()
    s = background.std(ddof=1)
    z = (real - m) / s
    p = np.mean(background < real)  # left
    return z, p


# ========================= ======== ========================= #
def ToDict(fp, Net):
    mapping = {}
    with open(fp) as fi:
        for line in fi:
            k, v = line.split()
            if v in Net.nodes:
                if k in mapping:
                    mapping[k].append(Net.n2i[v])
                else:
                    mapping[k] = [Net.n2i[v]]
    return mapping


if __name__ == "__main__":
    FILE1 = "drug.txt"
    FILE2 = "PDL1.txt"
    OUTPUT = "PDL1_result_1k.txt"
    REPEAT = 1000
    MEASURE = "CLOSEST"

    Net = Network("HumanInteractome.adj", "HumanInteractome.npy")
    Net.PrepareBins(100)
    MAPPING1 = ToDict(FILE1, Net)
    MAPPING2 = ToDict(FILE2, Net)
    with open(OUTPUT, "w") as fo:
        for k1 in sorted(MAPPING1.keys()):
            for k2 in sorted(MAPPING2.keys()):
                g1 = MAPPING1[k1]
                g2 = MAPPING2[k2]
                d = Net.Proximity(g1, g2, method=MEASURE)
                b = Net.ProximityRandom(g1, g2, method=MEASURE, repeat=REPEAT, seed=1024)
                z, p = Z_Score(d, b)
                fo.write("%s\t%s\t%.3f\t%.3f\t%.3f\n" % (k1, k2, d, z, p))
