#!/usr/bin/python3 

import traceback
import sys
import logging
import fileinput
import itertools
from pprint import pprint
from pprint import pformat

def convert_edges_to_null_notation(edges):
    for e in edges:
        e[:2] = [x-1 for x in e[:2]]

class Blossom:
    def __init__(self):
        self.stem = None
        self.charge = None


class Blossom_simple(Blossom):
    def __init__(self, vertex, charge=None, parent_blossom=None):
        self.vertex = vertex
        self.stem = vertex
        self.charge = charge or 0
        self.parent_blossom = parent_blossom

    def get_all_vertices(self):
        return [self.vertex]

    def __repr__(self):
        return "{1}[v{0}]".format(str(self.vertex), str(self.charge))


class Blossom_composite(Blossom):
    def __init__(self, blossoms, edges, charge=None, parent_blossom=None):
        # blossoms[0] is the K1
        assert(len(blossoms) % 2 == 1)
        self.stem = blossoms[0].stem
        self.blossoms = blossoms
        self.edges = edges
        self.charge = charge or 0
        self.parent_blossom = parent_blossom

    def get_all_vertices(self):
        vss = []
        for b in self.blossoms:
            vss.append(b.get_all_vertices())
        return itertools.chain(*vss)

    def __repr__(self):
        return "{}[{}]".format(self.charge, pformat(self.blossoms))

class HTNode:
    def __init__(self, blossom, children=None):
        # children = [(<edge>, <child_node>)]
        self.blossom = blossom
        self.children = children or []

    def __repr__(self):
        if len(self.children) > 0:
            return "({}->{})".format(str(self.blossom), pformat(self.children)) 
        else:
            return "({})".format(str(self.blossom))

def find_critical_bubble(node, level=0):
    res = None
    # we subtract epsilon on odd levels
    if level % 2 == 1 and isinstance(node, Blossom_composite):
        res = node.blossom.charge
    for e, child in node.children:
        subres = find_critical_bubble(child, 1 + level)
        if res is None or (subres is not None and res.charge > subres.charge):
            res = subres
    return res

def get_positive_and_negative_vertices(node, level=0):
    xpos = []
    xneg = []
    if level % 2 == 0:
        xpos.extend(node.blossom.get_all_vertices())
    else:
        xneg.extend(node.blossom.get_all_vertices())
    for e, child in node.children:
        spos, sneg = get_positive_and_negative_vertices(child, level + 1)
        xpos.extend(spos)
        xneg.extend(sneg)
    return (xpos, xneg)

def all_outer_blossoms(node):
    # pre-order
    yield node.blossom
    for e, child in node.children:
        yield from all_outer_blossoms(child)

def all_outer_blossoms_with_level(node, level=0):
    yield (node.blossom, level)
    for e, child in node.children:
        yield from all_outer_blossoms_with_level(node, level+1)

def find_min_cost_1_factor(edges, N):
    logging.warning("CORE PROCEDURE IS NOT FULLY IMPLEMENTED")
    
    infty = sum([e[2] for e in edges]) + 47

    charges = {(min(x,y), max(x,y)):0 for x,y,c in edges}
    def get_charge(x,y):
        nonlocal charges
        return charges[(min(x,y), max(x,y))]
    def add_charge(x,y, dc):
        nonlocal charges
        charges[(min(x,y), max(x,y))] += dc

    M = [] # current pairing 
    L = [] # blocking edges (edges in the hungar forest from even to odd )
    dumbbells = [] # paired blossoms
    hungar_forest = [HTNode(Blossom_simple(x)) for x in range(N)] # unpaired blossoms

    iteration_count = 0
    while True:
        iteration_count += 1
        logging.debug("=============== ITERATION: {} ================".format(str(iteration_count)))
        if iteration_count > N**3:
            logging.debug("TOO MUCH ITERATIONS: >{}**3 == {}".format(str(N), str(N**3)))
            break
        # eval post_condition
        if len(hungar_forest) == 0:
            logging.debug("NO FOREST LEFT, SUCCESS :)")
            break
        # @TODO

        logging.debug("FOREST: " + pformat(hungar_forest))
        logging.debug("DUMBBELLS: " + pformat(dumbbells))
        logging.debug("M: {}".format(M))
        logging.debug("Edge charges: {}".format(pformat(charges)))
        # we are trying to add charge to all root blossoms in the hungar forest
        # we need to evaulate maximum possible charge to add
        # find the critical spot - 4 possibilites
        # P1: composite blossom get the charge 0 - break the blossom, rebuild the tree
        # P2: an edge between a dumbbell and blossom (on the even level) is full - add to the tree
        # P3: an edge is filled between two blossoms in one tree - join the blossoms, rebuild the tree
        # P4: an edge is filled between two blossoms in the different trees

        # so two general branches: either bubble pops or an edge is filled

        # lowest charge to pop an outer bubble:
        min_blossom = None
        for root in hungar_forest:
            blossom = find_critical_bubble(root)
            logging.debug("critical blossom: " + pformat(blossom))
            if blossom is not None:
                if min_blossom is None:
                    min_blossom = blossom
                elif min_blossom.charge > blossom.charge:
                    min_blossom = blossom
        pop_eps = min_blossom.charge if min_blossom is not None else None
        logging.debug("pop_eps: {}".format(pop_eps))

        # second part: filling edges
        vpos, vneg = [], []
        for root in hungar_forest:
            xpos, xneg = get_positive_and_negative_vertices(root)
            vpos.extend(xpos)
            vneg.extend(xneg)

        vpos = set(vpos)
        vneg = set(vneg)
        set(range(N)).difference(vpos).difference(vneg)

        logging.debug("Positive vertices: {}".format(pformat(vpos)))
        logging.debug("Negative vertices: {}".format(pformat(vneg)))
        logging.debug("Neutral vertices: {}".format(pformat(set(range(N)).difference(vpos).difference(vneg))))


        min_edge = (None, None, None)
        edge_sign = {(x,y):0 for x,y,c in edges}
        for x, y, c in edges:
            logging.debug("edge: {}".format(str((x,y,c))))
            # eval current charge
            charge = get_charge(x, y)
            logging.debug("edge charge: {}".format(str(charge)))
            # eval the eps polarity of x and y
            xsgn = 0
            if x in vpos: xsgn = +1
            elif x in vneg: xsgn = -1
            ysgn = 0
            if y in vpos: ysgn = +1
            elif y in vneg: ysgn = -1 
            logging.debug("xsgn: {}, ysgn: {}".format(str(xsgn), str(ysgn)))

            s, r = sorted([xsgn, ysgn])

            # if all are positive or one is positive and one is neutral, we have problems
            # otherwise they will either cancel out or decrease

            edge_sign[(x,y)] = s + r
            #if s >= 0 and r >= 0 and s * r > 0:
            if s + r > 0:
                limit = c - charge
                if s + r == 2:
                    limit /= 2
                if min_edge[2] is None or min_edge[2] > limit:
                    min_edge = (x,y,limit)
        logging.debug("edge_eps: {}, edge: {}".format(str(min_edge[2]), str(min_edge[:2])))

        did_something = False

        if min_edge[2] is None and pop_eps is None:
            logging.warning("DAFUQ: no constraints on adding charge")
            did_something = False
        elif min_edge[2] is None or (pop_eps is not None and pop_eps <= min_edge[2]):
            logging.debug("P1: composite blossom get the charge 0")
            logging.warning("UNIMPLEMENTED")
        elif pop_eps is None or (min_edge is not None and min_edge[2] <= pop_eps):
            # filled edge
            # three possibilities: same tree, two trees or tree and a dumbbell
            x, y = min_edge[:2]
            eps = min_edge[2]

            logging.debug("filled edge: {}".format((x, y)))
            
            xti = None
            xdi = None
            xblossom = None

            yti = None
            ydi = None
            yblossom = None

            for i, root in enumerate(hungar_forest):
                for blossom in all_outer_blossoms(root):
                    if x in blossom.get_all_vertices():
                        xti = i
                        xblossom = blossom
                    if y in blossom.get_all_vertices():
                        yti = i
                        yblossom = blossom

            for i, (b1, b2) in enumerate(dumbbells):
                for b in (b1, b2):
                    if x in b.get_all_vertices():
                        xdi = i
                        xblossom = b
                    if y in b.get_all_vertices():
                        ydi = i
                        yblossom = b

            logging.debug("xtree, xdumbbell, xblossom:{}".format(pformat((xti, xdi, xblossom))))
            logging.debug("ytree, ydumbbell, yblossom:{}".format(pformat((yti, ydi, yblossom))))

            if xti is not None and yti is not None and xti == yti:
                logging.debug("P3: an edge is filled between two blossoms in one tree")
                pass
            elif xti is not None and yti is not None and xti != yti:
                logging.debug("P4: an edge is filled between two blossoms in the different trees")
                # simple case just to get going
                if len(hungar_forest[xti].children) == 0 and len(hungar_forest[yti].children) == 0:
                    logging.debug("Both trees are just nodes")
                    # adding charges to blossoms
                    for root in hungar_forest:
                        for b, level in all_outer_blossoms_with_level(root):
                            b.charge += eps * (1 if level % 2 == 0 else -1)
                    # adding charges to edges
                    for xx,yy,cc in edges:
                        add_charge(xx, yy, eps * edge_sign[(xx, yy)])
                    # creating dumbbell
                    dumbbells.append((xblossom, yblossom))
                    # removing trees
                    if xti > yti: xti, yti = yti, xti
                    hungar_forest.pop(yti)
                    hungar_forest.pop(xti)
                    # adding an edge to M
                    M.append((x,y))
                    did_something = True
                else:
                    logging.debug("Serious case")
                    pass
            else:
                logging.debug("P2: an edge between a dumbbell and blossom")
                # adding charges to blossoms
                for root in hungar_forest:
                    for b, level in all_outer_blossoms_with_level(root):
                        b.charge += eps * (1 if level % 2 == 0 else -1)
                # adding charges to edges
                for xx,yy,cc in edges:
                    add_charge(xx, yy, eps * edge_sign[(xx, yy)])

                if xti is None:
                    x, xti, xdi, xblossom, y, yti, ydi, yblossom = y, yti, ydi,yblossom, x, xti, xdi, xblossom
                # yti is none, so y is in the dumbbell
                root = hungar_forest[xti]
                db = dumbbells[ydi]
                if yblossom == db[1]:
                    db = db[::-1]
                # y is in the db[0]
                h2 = HTNode(db[1])
                h1 = HTNode(db[0], [  ((db[0].stem, db[1].stem), h2),   ])
                root.children.append(((x, y), h1))
                # remove db from dumbbells
                dumbbells.pop(ydi)
                did_something = True

        else:
            logging.warning("DAFUQ: you sholdn't be here")
            did_something = False
        
        if not did_something:
            logging.warning('DID NOTHING! STOPPING THE LOOP')
            break

    return [(x,y, get_charge(x,y)) for x,y in M]

def main():
    logging.basicConfig(level=logging.DEBUG)

    lines = [l.strip() for l in fileinput.input()]
    n, m = [int(x) for x in lines[0].split()]
    edges = [[int(x) for x in l.split()] for l in lines[1:]][:m]
    assert(len(edges) == m)
    
    convert_edges_to_null_notation(edges)

    factor = find_min_cost_1_factor(edges, n)
    logging.debug(factor)    

    print(int(sum([c for x,y,c in factor])))
    for x,y,c in factor:
        print(x+1, y+1)

if __name__ == "__main__":
    main()