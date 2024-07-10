import sys
import os
import math
import numpy as np
import pickle
import CellModeller
import subprocess
import string
import shutil
from CellModeller.Simulator import Simulator
import networkx
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib import units
from reportlab.lib.colors import Color
from reportlab.graphics.shapes import Circle
import pandas as pd
import glob
import os


def generate_network(G, num_cells, ct_tos, n_cts, pos, celldata, cell_ids):
    for i in range(0, num_cells):  # make all nodes
        G.add_node(i, pos=pos[i], type=celldata[i], label=cell_ids[i],
                   color=(1, 0, 0) if celldata[i] == 0 else (0, 1, 0))
    for i in range(0, num_cells):  # make all edges (when all nodes are present)
        for j in range(0, n_cts[i]):  # make all edges from contact_tos
            G.add_edge(i, ct_tos[i, j], width=8, color='black')


def get_current_contacts(G, data):
    cs = data['cellStates']
    it = iter(cs)
    n = len(cs)
    cell_type = {}
    pos_dict = {}
    cell_ids = {}
    for it in cs:
        cell_ids[cs[it].idx] = cs[it].id
        cell_type[cs[it].idx] = cs[it].cellType
        pos_dict[cs[it].idx] = cs[it].pos[0:2]

    # print(cell_ids)
    modname = data['moduleName']
    moduleStr = data['moduleStr']
    sim = Simulator(modname, 0.0, moduleStr=moduleStr, saveOutput=False)
    sim.loadFromPickle(data)
    sim.phys.update_grid()  # we assume local cell_centers is current
    sim.phys.bin_cells()
    sim.phys.cell_sqs = sim.phys.cell_sqs_dev.get()  # get updated cell sqs
    sim.phys.sort_cells()
    sim.phys.sorted_ids_dev.set(sim.phys.sorted_ids)  # push changes to the device
    sim.phys.sq_inds_dev.set(sim.phys.sq_inds)
    sim.phys.find_contacts(predict=False)
    sim.phys.get_cts()
    ct_pts = sim.phys.ct_pts  # these are the points on the cell surface - they can be transformed into the global
    # coordinate system
    ct_tos = sim.phys.ct_tos  # this is the list of *some* contacts from each cell (only for lower cell_ids)
    ct_dists = sim.phys.ct_dists
    cell_cts = sim.phys.cell_n_cts  # not really all the contacts of the cell, because the tos are only defined
    # between a cell and ones with lower ids
    generate_network(G, n, ct_tos, cell_cts, pos_dict, cell_type, cell_ids)


def find_neighbors(G, time, rows):
    # List to hold the rows of the dataframe

    # Iterate over each node in the graph
    for node in G.nodes():
        # Get the label of the node
        node_label = G.nodes[node]['label']
        # Get the neighbors of the node
        neighbors = list(G.neighbors(node))
        # Append a tuple (node_label, neighbor_label) for each neighbor
        for neighbor in neighbors:
            neighbor_label = G.nodes[neighbor]['label']
            rows.append((time + 1, node_label, neighbor_label))

    return rows


def neighbor_finders(pickle_folder):

    bg_color = Color(1.0, 1.0, 1.0, alpha=1.0)
    rows = []

    for time, fname in enumerate(glob.glob(pickle_folder + '/*.pickle')):

        G = networkx.Graph()
        # fname = sys.argv[1]
        data = pickle.load(open(fname, 'rb'))
        cs = data['cellStates']
        it = iter(cs)
        n = len(cs)
        oname = fname.replace('.pickle', '_graph.pdf')

        # print(("num_cells = " + str(n)))

        cell_type = {}
        pos_dict = {}
        for it in cs:
            cell_type[cs[it].idx] = cs[it].cellType
            pos_dict[cs[it].idx] = cs[it].pos[0:2]

        get_current_contacts(G, data)

        # print(("num_contacts = " + str(networkx.number_of_edges(G))))
        # degrees = list(G.degree().values())
        # degrees = [degree for node, degree in G.degree()]
        # print(("mean degree = " + str(np.mean(degrees))))
        rows = find_neighbors(G, time, rows)
        #
        # breakpoint()

        # if list(networkx.get_edge_attributes(G, 'color').items()):
        #    edges, ecolor = list(zip(*list(networkx.get_edge_attributes(G, 'color').items())))
        # ncolor = list(networkx.get_node_attributes(G, 'color').values())

        # pdf = CellModellerPDFGenerator(oname, G, bg_color)
        # world = (120, 120)
        # page = (50, 50)
        # center = (0, 0)

        # pdf.draw_frame(oname, world, page, center)

        # from this, I want the spatial position of every contact, and the graph of the cells that each cell is touching
        # This can be all in one data structure - a graph where each vertex has a position, and each node has a cellState

    # Create a dataframe from the list of tuples
    df = pd.DataFrame(rows, columns=['Image Number', 'First Object id', 'Second Object id'])
    return df
