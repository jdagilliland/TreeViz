#!/usr/bin/env python
import csv

import ete2
import numpy as np

def print_tree(fname, **kwopts):
    """
    Print a tree read from a `fname`.

    Parameters
    ----------
    fname : str
        The file name to use to print the tree from.

    Returns
    -------
    None
    """
    tree = ete2.Tree(fname)
    return None

def color_nodes(tree, tabfile, column, dict_color=None):
    """
    Color nodes of `tree` according to data in `tabfile`.

    Parameters
    ----------
    tree : ete2.tree.TreeNode 
        The tree to color.
    tabfile : str
        The filename of the `tabfile` from which to draw data.
    column : str
        A string to indicate the column of the `tabfile` to use.
    dict_color : dict, optional
        A dictionary of colors to use based on the values in the
        selected `column`.
        If not provided, a default dictionary of appropriately spaced colors
        will be constructed for all values found in `column`.

    Returns
    -------
    tree : ete2.tree.TreeNode
        The colored `tree`.

    """
    with open(tabfile,'rb') as f:
        reader = csv.DictReader(f,delimiter='\t')
        lst_dict_entries = [row for row in reader]
    if dict_color == None:
        column_data = [entry.get(column) for entry in lst_dict_entries]
        dict_color = _get_color_dict(column_data)
    for node in tree.traverse():
        try:
            dict_entry = _get_node_entry(node.name, lst_dict_entries)
        except ValueError:
            continue
        column_data = dict_entry.get(column)
        color = dict_color.get(column_data,'none')
        node.color = color
        style = ete2.NodeStyle()
        style['fgcolor'] = color
        node.set_style(style)
    return tree

def _get_node_entry(nodename, lst_dict_entries):
    for dict_entry in lst_dict_entries:
        if dict_entry['SEQUENCE_ID'][-9:] == nodename:
            return dict_entry
        else: continue
    raise ValueError('Nodename {name} not found'.format(name=nodename))

def _get_color_dict(lst_col_data):
    set_data = set(lst_col_data)
    print(len(set_data))
#    dict_color = dict((data, ete2.random_color()) for data in set_data)
    dict_color = dict()
    for data in set_data:
        dict_color[data] = _random_color()
    print(dict_color)
    return dict_color

def _random_color():
    tpl_color = np.random.rand(3,)
    return _color_3pl2hex(tpl_color)

def _color_3pl2hex(tpl_color):
    hex_color = '#' + ''.join([hex(int(256*iI))[-2:] for iI in tpl_color])
    return hex_color

