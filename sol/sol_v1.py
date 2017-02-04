#!/usr/bin/python3

import traceback
import sys
import logging
import fileinput
import itertools
from pprint import pprint
from pprint import pformat

# vertex is just a number

def convert_edges_to_null_notation(edges):
    for e in edges:
        e[:2] = [x-1 for x in e[:2]]

class Graph_edge:
    def __init__(self, x, y, capacity, charge=None):
        self.x = x
        self.y = y
        self.capacity = capacity
        self.charge = charge or 0

    def __repr__(self):
        return "({}={}/{}={})".format(self.x, self.charge, self.capacity, self.y)

class Graph:
    def __init__(self, n, edges):
        self.n = n
        self.edges = []
        for x,y,c in edges:
            self.edges.append(Graph_edge(x,y,c))

    def get_edges(self):
        return self.edges

    def get_vertices(self):
        return range(self.n)

    def get_edge_by_vertices(self, x, y):
        for e in self.edges:
            if {e.x, e.y} == {x, y}:
                return e
        return None

    def __repr__(self):
        return "G:" + pformat(self.edges, compact=True)

class Matching:
    def __init__(self, graph, m=None):
        self.graph = graph
        self.m = m or []

    def get_edges(self):
        return self.m

    def add_edge(self, edge):
        logging.debug("MATCHING: adding an edge: {}".format(edge))
        assert (edge.capacity == edge.charge), "you've tried to add badly filled edge to matching: {}".format(edge)
        self.m.append(edge)

    def remove_edge(self, edge):
        logging.debug("MATCHING: remove edge: {}".format(edge))
        self.m.remove(edge)

    def __repr__(self):
        return "M:" + pformat(self.m, compact=True)


class HTNode:
    def __init__(self, blossom, parent_node=None, children_nodes=None):
        self.blossom = blossom
        self.parent_node = parent_node
        self.children_nodes = children_nodes or []

    def get_all_nodes(self):
        yield self
        for e, child_node in self.children_nodes:
            yield from child_node.get_all_nodes()

    def get_all_nodes_with_level(self, level=0):
        yield (self, level)
        for e, child_node in self.children_nodes:
            yield from child_node.get_all_nodes_with_level(1 + level)

    def get_root(self):
        node = self
        while node.parent_node is not None:
            node = node.parent_node[1]
        return node

    def get_children(self):
        return self.children_nodes

    def add_child(self, e, child):
        self.children_nodes.append((e, child))

    def add_children(self, children):
        self.children_nodes.extend(children)

    def remove_child_by_node(self, node):
        pos = None
        for i, (e, child) in enumerate(self.children_nodes):
            if child == node:
                pos == i
                break
        self.children_nodes.pop(i)

    def eval_path_to_the_root(self):
        a = [self]
        while a[-1].parent_node is not None:
            a.append(a[-1].parent_node[1])
        return a

    def __repr__(self):
        return "{" + str(self.blossom) + ("->" + str([self.children_nodes]) if len(self.children_nodes) > 0 else "") + "}"

    def verify_parent_links(self, level=0):
        for e, child in self.children_nodes:
            assert (e, self) == child.parent_node, "bad parent link: {} != {}".format((e, child), child.parent_node)
            child.verify_parent_links(level + 1)


class HTForest:
    def __init__(self, nodes=None):
        self.nodes = nodes or []

    def get_trees(self):
        return self.nodes

    def add_tree(self, node):
        self.nodes.append(node)

    def remove_tree(self, node):
        self.nodes.remove(node)

    def __repr__(self):
        return "F:" + ";\n".join([str(x) for x in self.nodes])

class Dumbbell:
    def __init__(self, b1, b2, edge=None):
        self.b1 = b1
        self.b2 = b2
        self.edge = edge

    def __repr__(self):
        return "<" + pformat(self.b1) + "|{}|".format(self.edge) + pformat(self.b2) + ">"

class Dumbbell_array:
    def __init__(self, dumbbells=None):
        self.dumbbells = dumbbells or []

    def __repr__(self):
        return "D:" + pformat(self.dumbbells, compact=True)

    def get_dumbbells(self):
        return self.dumbbells

    def add_dumbbell(self, d):
        self.dumbbells.append(d)

    def remove_dumbbell(self, d):
        self.dumbbells.remove(d)

class Blossom:
    pass

class Blossom_simple(Blossom):
    def __init__(self, v, charge=None, parent_blossom=None, node=None, dumbbell=None):
        self.v = v
        self.charge = charge or 0
        self.parent_blossom = parent_blossom
        self.node = node
        self.dumbbell = dumbbell

    def is_outer(self):
        return (self.parent_blossom is None)

    def get_node(self):
        b = self
        while b.parent_blossom is not None:
            b = b.parent_blossom
        return b.node

    def get_stem_vertex(self):
        return self.v

    def get_all_vertices(self):
        yield self.v

    def __repr__(self):
        return "{}[{}]".format(self.charge, self.v)

    def flip_path(self, edge, M):
        # M is matching
        assert (self.v == edge.x or self.v == edge.y), "Edge is not connected to the blossom"

    def verify_parent_links(self):
        pass # intentionally left blank

    def verify_charge_sign(self):
        pass # intentionally left blank

class Blossom_composite(Blossom):
    def __init__(self, blossoms, blossom_edges, stem_blossom, charge=None, parent_blossom=None, node=None, dumbbell=None):
        assert (len(blossoms) > 1 and len(blossoms) % 2 == 1), "incorrect amount blossoms in composite blossom: {}".format(len(blossoms))
        self.blossoms = blossoms
        assert (len(blossom_edges) == len(self.blossoms)), "incorrect amount of interblossom links"
        self.blossom_edges = blossom_edges
        assert (stem_blossom in self.blossoms), "stem subblossom not in this blossom"
        self.stem_blossom = stem_blossom
        self.charge = charge or 0
        self.parent_blossom = parent_blossom
        self.node = node
        self.dumbbell = dumbbell

    def is_outer(self):
        return (self.parent_blossom is None)

    def get_node(self):
        b = self
        while b.parent_blossom is not None:
            b = b.parent_blossom
        return b.node

    def get_stem_vertex(self):
        blossom = self.stem_blossom
        while not isinstance(blossom, Blossom_simple):
            blossom = blossom.stem_blossom
        return blossom.get_stem_vertex()

    def get_all_vertices(self):
        for b in self.blossoms:
            yield from b.get_all_vertices()

    def __repr__(self):
        return "{}[{}]".format(self.charge, ",".join(["({}: {})".format(b,e) if b != self.stem_blossom else "(!{}: {}!)".format(b, e)   for b, e in  zip(self.blossoms, self.blossom_edges)] ))

    def get_index_of_stem_subblossom(self):
        for i, b in enumerate(self.blossoms):
            if b == self.stem_blossom:
                return i
        return None

    def get_index_of_subblossom_coincident_with_edge(self, edge):
        for i, b in enumerate(self.blossoms):
            for v in b.get_all_vertices():
                if edge.x == v or edge.y == v:
                    return (i, v)
        return (None, None)

    def shift_stem_subblossom_left(self):
        i_stem = self.get_index_of_stem_subblossom()
        self.blossoms = self.blossoms[i_stem:] + self.blossoms[:i_stem]
        self.blossom_edges = self.blossom_edges[i_stem:] + self.blossom_edges[:i_stem]

    def shift_stem_subblossom_right(self):
        i_stem = self.get_index_of_stem_subblossom()
        self.blossoms = self.blossoms[i_stem+1:] + self.blossoms[:i_stem + 1]
        self.blossom_edges = self.blossom_edges[i_stem+1:] + self.blossom_edges[:i_stem+1]

    def flip_path(self, edge, M):
        # M is the matching
        logging.debug("flip_path: {}, edge: {}".format(self, edge))
        stem_vertex = self.get_stem_vertex()
        if stem_vertex == edge.x or stem_vertex == edge.y:
            logging.debug("easy")
            return
        else:
            i_sub, v = self.get_index_of_subblossom_coincident_with_edge(edge)
            assert(i_sub is not None), "there is no subblossom coincident with {}".format(edge)
            logging.debug("i_sub: {}, v:{}".format(i_sub, v))
            i_stem = self.get_index_of_stem_subblossom()
            assert(i_stem is not None), "blossom {} has no stem subblossom!".format(self)
            logging.debug("i_stem: {}".format(i_stem))

            if i_stem == i_sub:
                logging.debug("stem and outer edge are coincident, easy")
                self.stem_blossom.flip_path(edge, M)
            else:
                logging.debug("stem and outer edge are not coincident, hard")
                if abs(i_stem - i_sub) % 2 == 0:
                    logging.debug("an even path is in the inner cycle")
                    if i_stem < i_sub:
                        for i in range(i_stem, i_sub, 2):
                            M.add_edge(self.blossom_edges[i])
                            M.remove_edge(self.blossom_edges[i+1])
                            self.blossoms[i].flip_path(self.blossom_edges[i], M)
                            self.blossoms[i+1].flip_path(self.blossom_edges[i], M)
                        self.blossoms[i_sub].flip_path(edge, M)
                        self.stem_blossom = self.blossoms[i_sub]
                    else:
                        for i in range(i_stem, i_sub, -2):
                            M.add_edge(self.blossom_edges[i-1])
                            M.remove_edge(self.blossom_edges[i-2])
                            self.blossoms[i].flip_path(self.blossom_edges[i-1], M)
                            self.blossoms[i-1].flip_path(self.blossom_edges[i-1], M)
                        self.blossoms[i_sub].flip_path(edge, M)
                        self.stem_blossom = self.blossoms[i_sub]
                else:
                    logging.debug("an even path is in the outer cycle")
                    if i_stem < i_sub:
                        logging.debug("shifting stem blossom to the right and recursive")
                        self.blossoms = self.blossoms[i_stem + 1:] + self.blossoms[:i_stem+1]
                        self.blossom_edges = self.blossom_edges[i_stem+1:] + self.blossom_edges[:i_stem+1]
                        self.flip_path(edge, M)
                    else:
                        logging.debug("shifting stem blossom to the left and recursive")
                        self.blossoms = self.blossoms[i_stem:] + self.blossoms[:i_stem]
                        self.blossom_edges = self.blossom_edges[i_stem:] + self.blossom_edges[:i_stem]
                        self.flip_path(edge, M)
                        #raise Exception("Not tested")

    def verify_parent_links(self):
        for b in self.blossoms:
            assert(b.parent_blossom == self), "Parent blossom is set incorrectly for blossom {} ({} instead of {})".format(b, b.parent_blossom, self)
            b.verify_parent_links()

    def verify_charge_sign(self):
        assert(self.charge >= 0), "composite blossom {} has negative charge {}!".format(self, self.charge)
        for b in self.blossoms:
            b.verify_charge_sign()

class Charger:
    def __init__(self, graph, forest, dumbbell_array):
        self.graph = graph
        self.forest = forest
        self.dumbbell_array = dumbbell_array

    def find_first_blossom_to_pop(self):
        eps, blossom = None, None
        for root in self.forest.get_trees():
            for node, level in root.get_all_nodes_with_level():
                if level % 2 == 1 and isinstance(node.blossom, Blossom_composite):
                    if eps is None or eps > node.blossom.charge:
                        eps = node.blossom.charge
                        blossom = node.blossom
        return (eps, blossom)

    def find_first_edge_to_fill(self):
        v_blossom = [None for v in self.graph.get_vertices()]
        v_sign = [0 for v in self.graph.get_vertices()]
        for root in self.forest.get_trees():
            for node, level in root.get_all_nodes_with_level():
                sign = +1 if level % 2 == 0 else -1
                blossom = node.blossom
                for v in blossom.get_all_vertices():
                    v_blossom[v] = blossom
                    v_sign[v] = sign
        for dumbbell in self.dumbbell_array.get_dumbbells():
            for b in (dumbbell.b1, dumbbell.b2):
                for v in b.get_all_vertices():
                    v_blossom[v] = b
        
        min_eps, min_edge = None, None
        for edge in self.graph.get_edges():
            if v_blossom[edge.x] != v_blossom[edge.y]:
                total_sign = v_sign[edge.x] + v_sign[edge.y]
                if total_sign > 0:
                    eps = (edge.capacity - edge.charge)/total_sign
                    if min_eps is None or eps < min_eps:
                        min_eps, min_edge = eps, edge
        if min_eps is None:
            return (None, None, None)
        else:
            return (min_eps, min_edge, (v_blossom[min_edge.x], v_sign[min_edge.x], v_blossom[min_edge.y], v_sign[min_edge.y]))

    def add_charge(self, charge):
        v_blossom = [None for v in self.graph.get_vertices()]
        v_sign = [0 for v in self.graph.get_vertices()]
        for root in self.forest.get_trees():
            for node, level in root.get_all_nodes_with_level():
                sign = +1 if level % 2 == 0 else -1
                blossom = node.blossom
                blossom.charge += charge * sign # sic!
                for v in blossom.get_all_vertices():
                    v_blossom[v] = blossom
                    v_sign[v] = sign
        
        for edge in self.graph.get_edges():
            if v_blossom[edge.x] != v_blossom[edge.y]:
                total_sign = v_sign[edge.x] + v_sign[edge.y]
                edge.charge += total_sign * charge # sic!


class Solver:
    def __init__(self, graph):
        self.graph = graph
        # @TODO add customisation of state if needed

        self.M = Matching(graph)
        self.forest = HTForest()
        for v in self.graph.get_vertices():
            b = Blossom_simple(v)
            node = HTNode(b)
            b.node = node
            self.forest.add_tree(node)
        self.dumbbell_array = Dumbbell_array()

        self.charger = Charger(self.graph, self.forest, self.dumbbell_array)

    def verify_state(self):
        # edges are fine
        for edge in self.graph.get_edges():
            assert (edge.capacity >= edge.charge), "edge: {} is overcharged!".format(edge)
        logging.debug("edges are not overcharged")

        # all edges in matching are filled
        for edge in self.M.get_edges():
            assert (edge.capacity == edge.charge), "edge: {} from matching is badly charged!".format(edge)
        logging.debug("edges in matching are filled")
        # all nodes has correct parent node set
        for root in self.forest.get_trees():
            assert (root.parent_node is None), "root of the tree {} appently has a parent in {}!".format(root, root.parent_node)
            root.verify_parent_links()
        logging.debug("all blossoms in the forest has correct parent links")
        # all blossoms has correct parent blossom set
        for root in self.forest.get_trees():
            for node in root.get_all_nodes():
                assert(node.blossom.node == node), "incorrect parent node set to a blossom {} ({} instead of {})".format(node.blossom, node.blossom.node, node)
                assert(node.blossom.dumbbell is None), "incorrect parend dumbbell set to a blossom {} ({})".format(node.blossom, node.blossom.dumbbell)
        logging.debug("all blossoms in the forest has correct node link")
        # all blossoms has correct parent node set
        for root in self.forest.get_trees():
            for node in root.get_all_nodes():
                assert(node.blossom.parent_blossom is None), "outer blossom {} has a parent blossom {}!".format(node.blossom, node.blossom.parent_blossom)
                node.blossom.verify_parent_links()
                node.blossom.verify_charge_sign()
        logging.debug("All blossoms in the forest has correct parent blossom set")
        # all blossoms has correct dumbbell set
        for d in self.dumbbell_array.get_dumbbells():
            for b in (d.b1, d.b2):
                assert(b.node is None), "{} has a parent node {}!".format(b, b.node)
                assert(b.dumbbell == d), "blossom {} has incorrect dumbbell link ({} instead of {})".format(b, b.dumbbell, d)
        logging.debug("all blossoms in the dumbbells has correct links to dumbbells")
        # all blossoms in dumbbells has correct parent blossom links
        for d in self.dumbbell_array.get_dumbbells():
            for b in (d.b1, d.b2):
                assert(b.parent_blossom is None), "outer blossom {} has a parent blossom {}!".format(b, b.parent_blossom)
                b.verify_parent_links()
                b.verify_charge_sign()
        logging.debug("All blossoms in the dumbbells has correct parent blossom links")

    def decompose_tree_into_dumbbells(self, node):
        assert(len(node.get_children()) == 1), "not exaclty one child on the odd level! node: {}, children:{} ".format(node, pformat(node.get_children()))
        E, Child = node.get_children()[0]

        dumbbell = Dumbbell(node.blossom, Child.blossom, E)
        node.blossom.node = None
        node.blossom.dumbbell = dumbbell
        Child.blossom.node = None
        Child.blossom.dumbbell = dumbbell
        self.dumbbell_array.add_dumbbell(dumbbell)

        for e, child in Child.get_children():
            self.decompose_tree_into_dumbbells(child)

    def decompose_connected_trees_into_dumbbells(self, x_tree, x_blossom, y_tree, y_blossom, critical_edge):
        if x_tree == x_blossom.node and len(x_tree.children_nodes) == 0 and isinstance(x_blossom, Blossom_simple) \
                and y_tree == y_blossom.node and len(y_tree.children_nodes) == 0 and isinstance(y_blossom, Blossom_simple):
            logging.debug("Both trees are just nodes with simple blossoms, simplified routine is used")
            self.M.add_edge(critical_edge)
            dumbbell = Dumbbell(x_blossom, y_blossom, critical_edge)
            x_blossom.node = None
            x_blossom.dumbbell = dumbbell
            y_blossom.node = None
            y_blossom.dumbbell = dumbbell
            self.forest.remove_tree(x_tree)
            self.forest.remove_tree(y_tree)
            self.dumbbell_array.add_dumbbell(dumbbell)
        else:
            logging.debug("Serious case")
            xnode = x_blossom.node
            ynode = y_blossom.node
            xpath = xnode.eval_path_to_the_root()
            ypath = ynode.eval_path_to_the_root()
            logging.debug("xpath: {}".format(pformat(xpath)))
            logging.debug("ypath: {}".format(pformat(ypath)))

            # couple newly connected blossoms
            # first, remove subtrees under these blossoms
            for e, child in itertools.chain(xnode.get_children(), ynode.get_children()):
                self.decompose_tree_into_dumbbells(child)
            # we may also erase all children links and this moment, but meh
            self.M.add_edge(critical_edge)
            xnode.blossom.flip_path(critical_edge, self.M)
            ynode.blossom.flip_path(critical_edge, self.M)
            xnode.blossom.node = None
            ynode.blossom.node = None
            dumbbell = Dumbbell(xnode.blossom, ynode.blossom, critical_edge)
            xnode.blossom.dumbbell = dumbbell
            ynode.blossom.dumbbell = dumbbell
            self.dumbbell_array.add_dumbbell(dumbbell)

            for path in (xpath, ypath):
                for i in range(1, len(path), 2):
                    self.M.remove_edge(path[i-1].parent_node[0])
                    self.M.add_edge(path[i].parent_node[0])
                    edge = path[i].parent_node[0]
                    path[i].blossom.flip_path(edge, self.M)
                    path[i].blossom.flip_path(edge, self.M)
                    for e, child in path[i+1].get_children():
                        if e != edge:
                            self.decompose_tree_into_dumbbells(child)
                    path[i].blossom.flip_path(edge, self.M)
                    path[i+1].blossom.flip_path(edge, self.M)
                    path[i].blossom.node = None
                    path[i+1].blossom.node = None
                    dumbbell = Dumbbell(path[i].blossom, path[i+1].blossom, edge)
                    path[i].blossom.dumbbell = dumbbell
                    path[i+1].blossom.dumbbell = dumbbell
                    self.dumbbell_array.add_dumbbell(dumbbell)
                self.forest.remove_tree(path[-1])


    def wrap_edge_within_a_new_blossom(self, tree, x_blossom, y_blossom, critical_edge):
        xnode = x_blossom.node
        ynode = y_blossom.node        

        xpath = xnode.eval_path_to_the_root()[::-1]
        logging.debug("xpath: {}".format(xpath))
        ypath = ynode.eval_path_to_the_root()[::-1]
        logging.debug("ypath: {}".format(ypath))

        # finding the LCA
        q = 0
        while q < min(len(xpath), len(ypath)) and xpath[q] == ypath[q]:
            q += 1
        q -= 1
        logging.debug("q: {}".format(q))
        xpath = xpath[q:]
        ypath = ypath[q:]
        logging.debug("cutted xpath: {}".format(xpath))
        logging.debug("cutted ypath: {}".format(ypath))

        blossoms = []
        bedges = []
        for node in xpath:
            blossoms.append(node.blossom)
        for node in xpath[1:]:
            bedges.append(node.parent_node[0])
        bedges.append(critical_edge)
        for node in ypath[::-1][:-1]:
            blossoms.append(node.blossom)
        for node in ypath[::-1][:-1]:
            logging.debug(node)
            bedges.append(node.parent_node[0])

        logging.debug("blossoms: {}".format(pformat(blossoms)))
        logging.debug("bedges: {}".format(pformat(bedges)))
        
        # rebuilding
        newnode = xpath[0]
        blossom = Blossom_composite(blossoms=blossoms, blossom_edges=bedges, stem_blossom=blossoms[0], charge=0, node=newnode)
        for b in blossoms:
            b.parent_blossom = blossom

        nedges = []
        for node in itertools.chain(xpath[1:], ypath[1:]):
            for e, child in node.children_nodes:
                if e not in bedges:
                    nedges.append((e, child))
        for e, child in nedges:
            child.parent_node = (e, newnode)
        if len(xpath) > 1:
            newnode.remove_child_by_node(xpath[1])
        if len(ypath) > 1:
            newnode.remove_child_by_node(ypath[1])
        newnode.add_children(nedges)
        newnode.blossom = blossom

    def pop_a_bubble(self, critical_blossom):
        logging.debug("pop_a_bubble: blossom: {}".format(critical_blossom))
        upedge, upnode = critical_blossom.node.parent_node
        logging.debug("upedge: {}, upnode: {}".format(upedge, upnode))
        downedge, downnode = critical_blossom.node.children_nodes[0]
        logging.debug("downedge: {}, downnode: {}".format(downedge, downnode))

        i_stem = critical_blossom.get_index_of_stem_subblossom()
        i_sub, v = critical_blossom.get_index_of_subblossom_coincident_with_edge(upedge)
        logging.debug("i_stem: {}".format(i_stem))
        logging.debug("i_sub: {}, v: {}".format(i_sub, v))

        if i_stem == i_sub or \
            (i_stem < i_sub and (i_sub - i_stem) % 2 == 0) or \
            (i_stem > i_sub and (i_stem - i_sub) % 2 == 1):
            logging.debug("shift to the left")
            critical_blossom.shift_stem_subblossom_left()
        else:
            logging.debug("shift to the right")
            critical_blossom.shift_stem_subblossom_right()

        logging.debug("shifting the stem subblossom: {}".format(critical_blossom))

        i_stem = critical_blossom.get_index_of_stem_subblossom()
        i_sub, v = critical_blossom.get_index_of_subblossom_coincident_with_edge(upedge)
        logging.debug("i_stem: {}".format(i_stem))
        logging.debug("i_sub: {}, v: {}".format(i_sub, v))
        assert(abs(i_stem - i_sub) % 2 == 0), "the distance between blossoms is not even"

        upnode.remove_child_by_node(critical_blossom.node)

        if i_stem <= i_sub:
            logging.debug("from left to right")
            assert(i_stem == 0)

            nodes = []

            nodes.append(HTNode(critical_blossom.blossoms[i_sub]))
            critical_blossom.blossoms[i_sub].parent_blossom = None
            critical_blossom.blossoms[i_sub].node = nodes[0]
            for i in range(i_sub-1, -1, -1):
                nodes.append(HTNode(critical_blossom.blossoms[i]))
                critical_blossom.blossoms[i].parent_blossom = None
                critical_blossom.blossoms[i].node = nodes[-1]

            nodes[0].parent_node = (upedge, upnode)
            upnode.add_child(upedge, nodes[0])
            for i in range(1, len(nodes)):
                node_edge = critical_blossom.blossom_edges[i_sub - i]
                nodes[i].parent_node = (node_edge, nodes[i-1])
                nodes[i-1].add_child(node_edge, nodes[i])
            nodes[-1].add_child(downedge, downnode)
            downnode.parent_node = (downedge, nodes[-1])

            # new dumbbells
            for i in range(i_sub+1, len(critical_blossom.blossoms), 2):
                dumbbell = Dumbbell(critical_blossom.blossoms[i], critical_blossom.blossoms[i+1], critical_blossom.blossom_edges[i])
                for b in (critical_blossom.blossoms[i], critical_blossom.blossoms[i+1]):
                    b.parent_blossom = None
                    b.dumbbell = dumbbell
                    b.node = None
                self.dumbbell_array.add_dumbbell(dumbbell)
        else:
            logging.debug("from right to left")
            assert(i_stem == len(critical_blossom.blossoms) - 1)
            
            nodes = []
            nodes.append(HTNode(critical_blossom.blossoms[i_sub]))
            critical_blossom.blossoms[i_sub].parent_blossom = None
            critical_blossom.blossoms[i_sub].node = nodes[0]
            for i in range(i_sub + 1, len(critical_blossom.blossoms)):
                nodes.append(HTNode(critical_blossom.blossoms[i]))
                critical_blossom.blossoms[i].parent_blossom = None
                critical_blossom.blossoms[i].node = nodes[-1]

            nodes[0].parent_node = (upedge, upnode)
            upnode.add_child(upedge, nodes[0])
            for i in range(1, len(nodes)):
                node_edge = critical_blossom.blossom_edges[i_sub + i-1]
                nodes[i].parent_node = (node_edge, nodes[i-1])
                nodes[i-1].add_child(node_edge, nodes[i])
            nodes[-1].add_child(downedge, downnode)
            downnode.parent_node = (downedge, nodes[-1])

            # new dumbbells
            for i in range(0, i_sub, 2):
                dumbbell = Dumbbell(critical_blossom.blossoms[i], critical_blossom.blossoms[i+1], critical_blossom.blossom_edges[i])
                for b in (critical_blossom.blossoms[i], critical_blossom.blossoms[i+1]):
                    b.parent_blossom = None
                    b.node = None
                    b.dumbbell = dumbbell
                self.dumbbell_array.add_dumbbell(dumbbell)

            #raise Exception("Not implemented")

    def add_dumbbell_to_a_tree(self, tree_blossom, dumbbell_blossom, critical_edge):
        dumbbell = dumbbell_blossom.dumbbell
        logging.debug("dumbbell: {}".format(dumbbell))
        other_blossom = dumbbell.b1 if dumbbell.b1 != dumbbell_blossom else dumbbell.b2
        
        node_1 = HTNode(blossom=dumbbell_blossom, parent_node=(critical_edge, tree_blossom.node))
        dumbbell_blossom.dumbbell = None
        dumbbell_blossom.node = node_1
        
        node_2 = HTNode(blossom=other_blossom, parent_node=(dumbbell.edge, node_1))
        other_blossom.dumbbell = None
        other_blossom.node = node_2
        node_1.add_child(dumbbell.edge, node_2)

        tree_blossom.node.add_child(critical_edge, node_1)

        self.dumbbell_array.remove_dumbbell(dumbbell)


    def get_1_factor(self, max_iterations=None):
        max_iterations = max_iterations or (self.graph.n ** 2)
        iter_counter = 0
        while True:
            iter_counter += 1
            logging.debug("========================= ITERATION: {} =========================".format(iter_counter))
            logging.debug(self.graph)
            logging.debug(self.forest)
            logging.debug(self.dumbbell_array)
            logging.debug(self.M)

            self.verify_state()

            if iter_counter > max_iterations:
                logging.warning("TOO MUCH ITERATIONS, STOPPING")
                break
            if len(self.forest.get_trees()) == 0:
                logging.debug("NO TREES LEFT, SUCCESS :)")
                break

            eps_pop, critical_blossom = self.charger.find_first_blossom_to_pop()
            logging.debug("eps_pop: {}, critical blossom: {}".format(eps_pop, critical_blossom))
            eps_edge, critical_edge, edge_data = self.charger.find_first_edge_to_fill()
            logging.debug("eps_edge: {}, critical_edge: {}, edge_data: {}".format(eps_edge, critical_edge, pformat(edge_data)))

            if eps_pop is None and eps_edge is None:
                logging.warning("WTF: No constraints on adding charge!")
                break
            elif eps_pop is not None and (eps_edge is None or eps_pop <= eps_edge):
                logging.debug("P1: composite blossom get the charge 0")
                self.charger.add_charge(eps_pop)
                self.pop_a_bubble(critical_blossom)
            elif eps_edge is not None and (eps_pop is None or eps_edge < eps_pop):
                logging.debug("Edge is filled")
                self.charger.add_charge(eps_edge)
                x_blossom, x_sign, y_blossom, y_sign = edge_data
                if x_sign == +1 and y_sign == +1:
                    x_tree = x_blossom.node.get_root()
                    logging.debug("x_tree: {}".format(x_tree))
                    y_tree = y_blossom.node.get_root()
                    logging.debug("y_tree: {}".format(y_tree))
                    if x_tree == y_tree:
                        logging.debug("P3: an edge is filled between two blossoms in one tree")
                        self.wrap_edge_within_a_new_blossom(x_tree, x_blossom, y_blossom, critical_edge)
                    else:
                        logging.debug("P4: an edge is filled between two blossoms in the different trees")
                        self.decompose_connected_trees_into_dumbbells(x_tree, x_blossom, y_tree, y_blossom, critical_edge)
                else:
                    logging.debug("P2: an edge between a dumbbell and blossom")
                    if x_sign == +1:
                        self.add_dumbbell_to_a_tree(x_blossom, y_blossom, critical_edge)
                    else:
                        self.add_dumbbell_to_a_tree(y_blossom, x_blossom, critical_edge)
            else:
                logging.warning("WTF: you shouldn't be here!")
                break
        return self.M


def main():
    logging.basicConfig(level=logging.DEBUG)

    lines = [l.strip() for l in fileinput.input()]
    n, m = [int(x) for x in lines[0].split()]
    edges = [[int(x) for x in l.split()] for l in lines[1:]][:m]
    assert(len(edges) == m)
    
    convert_edges_to_null_notation(edges)

    graph = Graph(n, edges)
    solver = Solver(graph)
    matching = solver.get_1_factor()

    total_cost = sum([e.capacity for e in matching.get_edges()])
    print(total_cost)
    for e in matching.get_edges():
        print(e.x + 1, e.y + 1)


if __name__ == "__main__":
    main()