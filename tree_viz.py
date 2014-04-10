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
    color_column : str, optional
        The heading of the tabfile to use to color the graph.
        (default: 'COLOR_GROUP')
    size_column : str, optional
        The heading of the tabfile to use to size the nodes.
        (default: 'COPY_NUMBER')
    outfile : str, optional
        If provided, renders the tree to an image file, rather than
        displaying it in the GUI.

    Returns
    -------
    None
    """
    ## parse kwarg
    # tabfile really may be a better positional arg
    tabfile = kwarg.pop('tabfile', None)
    if tabfile == None:
        raise Exception('No tabfile provided.')
    # if color_column specified: use, otherwise default
    color_column = kwarg.pop('color_column','COLOR_GROUP')
    # if size_column specified: use, otherwise default
    size_column = kwarg.pop('size_column','COPY_NUMBER')
    
    treestyle = ete2.TreeStyle()
    treestyle.show_leaf_name = True
    treestyle.title.add_face(ete2.TextFace(fname), column=0)
    treestyle.force_topology = True
    treestyle.rotation = 0
    treestyle.show_branch_length = True
    treestyle.show_leaf_name = False
    treestyle.mode = 'c'
    treestyle.layout_fn = _internal_layout
    treestyle.legend_position = 1
    
    tree = ete2.Tree(fname)
    
    # color, size nodes
    dict_color = format_nodes(tree, tabfile,
            color_column=color_column,
            size_column=size_column,
            )
    _add_legend(treestyle, dict_color)
    # if outfile specified, use, otherwise just show
    outfile = kwarg.pop('outfile',None)
    if outfile:
        tree.render(outfile, 
#                w=841,
#                h=1189,
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

def format_nodes(tree, tabfile,
        color_column=None,
        size_column=None,
        dict_color=None,
        ):
    """
    Color nodes of `tree` according to data in `tabfile`.

    Parameters
    ----------
    tree : ete2.tree.TreeNode 
        The tree to color.
    tabfile : str
        The filename of the `tabfile` from which to draw data.
    color_column : str, optional
        A string to indicate the column of the `tabfile` to use for
        coloration.
        (default: None)
    size_column : str, optional
        A string to indicate the column of the `tabfile` to use to
        determine node size.
        (default: None)
    dict_color : dict, optional
        A dictionary of colors to use based on the values in the
        selected `color_column`.
        If not provided, a default dictionary of appropriately spaced colors
        will be constructed for all values found in `color_column`.

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
        _data = [entry.get(color_column) for entry in lst_dict_entries]
        dict_color = _get_color_dict(_data)
    for node in tree.traverse():
        try:
            dict_entry = _get_node_entry(node.name, lst_dict_entries)
        except ValueError:
            continue
        # get color data (default: 'none')
        color_data = dict_entry.get(color_column)
        color = dict_color.get(color_data,'none')
        # get size data (default: 1, assume single copy)
        size = dict_entry.get(size_column, 1)
        if size == '':
            size = 1
        # set node style
        style = ete2.NodeStyle()
        style['fgcolor'] = color
        style['size'] = _scale_size(int(size))
        node.set_style(style)
    return dict_color

def _add_legend(treestyle, dict_legend, legend_field=None):
    """
    Add an informative legend to the treestyle based on the dictionary
    of information used to color the tree.

    Returns
    -------
    None
    """
    basesize = 20
    for data, color in dict_legend.items():
        treestyle.legend.add_face(ete2.CircleFace(basesize, color), column=1)
        treestyle.legend.add_face(ete2.TextFace(data), column=0)
    return None

def _scale_size(size,
        **kwarg):
    """
    Scale raw size data on a linear scale to a log scale (or not).

    Parameters
    ----------
    size : float
        The raw size data on a linear scale.
    logscale : bool, optional
        Boolean whether or not to log scale.
        Set `False` to leave in linear scale.
        (default: `True`)
    basesize : float, optional
        The size of a node with representing a single copy number clone.
        (default: 20)

    Notes
    -----
    This function also works on arrays of size data, which might be
    useful when economizing code.
    """
    logscale = kwarg.pop('logscale', True)
    basesize = kwarg.pop('basesize', 20)
    if logscale:
        size = int(basesize * (1 + np.log(size)))
    else:
        size = int(basesize * size)
    return size

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
    parser.add_argument('-c', '--color-column',
            dest='color',
            default='COLOR_GROUP',
            )
    parser.add_argument('-s', '--size-column',
            dest='size',
            default='COPY_NUMBER',
            )
    argspace = parser.parse_args()
    print_tree(argspace.treefile,
        tabfile=argspace.tabfile,
        outfile=argspace.outfile,
        color_column=argspace.color,
        size_column=argspace.size,
        )
