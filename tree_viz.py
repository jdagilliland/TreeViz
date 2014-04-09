#!/usr/bin/env python
import csv

import ete2
import numpy as np
import Bio

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
    treestyle.rotation = 0
    treestyle.show_branch_length = True
    treestyle.show_leaf_name = False
    treestyle.mode = 'c'
#    treestyle.layout_fn = _internal_layout
    
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
        tree.render(outfile, 
                w=841,
                h=1189,
                units='mm',
                tree_style=treestyle,
                )
    else:
        tree.show(tree_style=treestyle)
    return None

def nwk2nkx(fname):
    """
    Loads a Newick tree from `fname` into a networkx DiGraph.

    Parameters
    ----------
    fname : str
        The file name to use to print the tree from.

    Returns
    -------
    nkx_tree : networkx.DiGraph
        The networkx tree.
    """
    phylo_tree = Bio.Phylo.read(fname,'newick')
    nkx_tree = Bio.Phylo.to_networkx(phylo_tree)
    return nkx_tree

def cleanup_tree(nkx_tree):
    """
    Cleans up a networkx tree (`nkx_tree`) made from a PHYLIP generated
    Newick tree.
    This includes:

        - Collapsing branches of 0 distance.
        - Placing the germline node at the root.

    Parameters
    ----------
    nkx_tree : networkx.DiGraph instance.
        The networkx tree read straight from the raw PHYLIP generated tree.

    Returns
    -------
    clean_tree : networkx.DiGraph
        The cleaned networkx tree.
    """
    # no-op so far
    return clean_tree

def _collapse_null_branches(nkx_tree):
    """
    Collapses 0 distance branches of a networkx tree (`nkx_tree`).

    Parameters
    ----------
    nkx_tree : networkx.DiGraph instance.
        The networkx tree read straight from the raw PHYLIP generated tree.

    Returns
    -------
    collapsed_tree : networkx.DiGraph
        The collapsed branch networkx tree.
    """
    # no-op so far
    return collapsed_tree

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
    dict_color : dict
        The dictionary of colors used in the colorization of the tree.
        This is especially useful if it was not specified as an input
        argument, because otherwise one has no way to know how to
        correlate colors on the plot with the values they color.

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
        style['size'] = 50
        node.set_style(style)
    return dict_color

def _internal_layout(node):
    if node.is_leaf():
         # If terminal node, draws its name
         name_face = ete2.AttrFace("name")
    else:
         # If internal node, draws label with smaller font size
         name_face = ete2.AttrFace("name")
    # Adds the name face to the image at the preferred position
    ete2.faces.add_face_to_node(name_face, node, column=0,
        position="branch-right")

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
    parser.add_argument('-c', '--column', dest='column', default=None)
    argspace = parser.parse_args()
    print_tree(argspace.treefile,
        tabfile=argspace.tabfile,
        outfile=argspace.outfile,
        column=argspace.column,
        )
