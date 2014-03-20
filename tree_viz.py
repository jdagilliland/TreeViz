#!/usr/bin/env python
import csv

import ete2
import numpy as np

def print_tree(fname, **kwarg):
    """
    Print a tree read from a `fname`.

    Parameters
    ----------
    fname : str
        The file name to use to print the tree from.
    tabfile : str, optional
        The name of the tabfile to use for coloring.
    column : str, optional
        The heading of the title to use to color the graph.
        (default: 'COLOR_GROUP')
    outfile : str, optional
        If provided, renders the tree to an image file, rather than
        displaying it in the GUI.

    Returns
    -------
    None
    """
    treestyle = ete2.TreeStyle()
    treestyle.show_leaf_name = True
    treestyle.title.add_face(ete2.TextFace(fname), column=0)
    treestyle.force_topology = True
    treestyle.rotation = 90
    treestyle.show_branch_length = True
    
    tree = ete2.Tree(fname)
    # tabfile really may be a better positional arg
    tabfile = kwarg.pop('tabfile', None)
    if tabfile == None:
        raise Exception
    # if column specified: use, otherwise default
    column = kwarg.pop('column','COLOR_GROUP')
    
    # color nodes
    color_nodes(tree, tabfile, column)
    # if outfile specified, use, otherwise just show
    outfile = kwarg.pop('outfile',None)
    if outfile:
        tree.render(outfile, w=600, tree_style=treestyle)
    else:
        tree.show(tree_style=treestyle)
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
        style['size'] = 10
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
    dict_color = dict()
    for data in set_data:
        dict_color[data] = _random_color()
    return dict_color

def _random_color():
    tpl_color = np.random.rand(3,)
    return _color_3pl2hex(tpl_color)

def _color_3pl2hex(tpl_color):
    hex_color = '#' + ''.join([hex(int(256*iI))[-2:] for iI in tpl_color])
    return hex_color

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Visualize newick trees',
        )
    parser.add_argument('treefile',
#        nargs=1,
#        dest='treefile',
        )
    parser.add_argument('-t', '--tabfile', dest='tabfile')
    parser.add_argument('-r', '--render', dest='outfile', default=None)
    argspace = parser.parse_args()
    print_tree(argspace.treefile,
        tabfile=argspace.tabfile,
        outfile=argspace.outfile,
        )
