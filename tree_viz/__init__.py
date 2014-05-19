#!/usr/bin/env python
import csv
import os

import ete2
import numpy as np
import Bio

class GermTree(ete2.coretype.tree.TreeNode):
    """
    Class to encapsulate the many actions that I need to make on
    ete2.Tree instances.
    """
    def __init__(self, *args, **kwarg):
        ete2.coretype.tree.TreeNode.__init__(self, *args, **kwarg)
        # initialize value of sequence length to 1
        self.len_seq = 1

    def set_tabfile(self, tabfile):
        """
        Read in data from a tabfile to which any method may refer.
        """
        with open(tabfile,'rb') as f:
            reader = csv.DictReader(f,delimiter='\t')
            self.lst_dict_tab_entries = [row for row in reader]
        return None

    def set_phyfile(self, phyfile):
        """
        Read in data from a PHYLIP formatted *.phy file.
        This will set the sequence length properly, as well as root the
        tree at the first sequence found in the PHY file.
        """
        with open(phyfile, 'rt') as phyfileobj:
            headrow = phyfileobj.readline()
            germrow = phyfileobj.readline()
        self.len_seq = int(headrow.split()[1])
        print('Set sequence length from {file0}'.format(file0=phyfile))
        self.germname = germrow.split()[0]

    def root_tree(self):
        """
        Properly root an ete tree according the the first sequence in the *.phy
        file.
        """
        if not hasattr(self, 'germname'):
            raise ValueError(
                '''Has no attribute 'germname'. ''' +
                '''Make sure to run self.set_phyfile first.''')
        # try to root tree at germline
        try:
            self.germnode = self.search_nodes(name=self.germname)[0]
            self.set_outgroup(self.germnode)
        except Exception as exc:
            print(exc)
            print('''Could not find node named {name} in this tree, the tree
            may be inappropriately rooted.'''.format(name=self.germname))

    def find_distant_subtrees(self,
            root_node=None,
            dist_lim=4):
        """
        Splits a tree into subtrees which root nodes are more than 4
        mutations away from the overall root (generally a germline
        sequence).

        Parameters
        ----------
        root_node : ete2.tree.TreeNode
            The root node from which to asses mutation distance.
        dist_lim : float
            The minimum distance a subtree must be from the root node to be
            considered a distant subtree.

        Returns
        -------
        lst_distsub : list of ete2.tree.Treenode
            The tree will have been split up into trees, each of whose
            respective roots will be at least `dist_lim` common mutations away
            from the root nodeself.
        """
        if root_node == None:
            # If this function is called without specifying a root node,
            # assume that the user is calling it on a whole tree, and
            # intends for the root of that tree to be treated as the root
            # node.
            root_node = self
            print('No root node specified. Using provided tree.')
        distance = self.get_distance(root_node)
        print('Distance between node and root: {:0.2f}'.format(distance))
        print('dist_lim: {:0.2f}'.format(dist_lim))
        if distance >= dist_lim:
            # If this node is already far enough away from root, return it,
            # skipping all descendants.
            print('Found distant subtree, skipping descendants')
            print(self)
            return [self]
        else:
            # If this node is still too close to root, iterate through each
            # child node to check for subtrees that start far enough away
            # from the root.
            print('Attempting to find distant subtrees among children nodes')
            lst_distsub = list()
            for child in self.get_children():
                print(child)
                print(root_node)
                lst_distsub.extend(child.find_distant_subtrees(
                    root_node=root_node,
                    dist_lim=dist_lim))
        return lst_distsub

    def collapse_null_branches(self):
        """
        Collapse 0 distance branches from an ete2 tree
        (`ete2.tree.TreeNode`).
        Applies recursively.

        Returns
        -------
        None
        """
        node = self
        # Remove a node's parent iff the node's branch length is 0 and the
        # parent's name is 'NoName', that way we avoid removing named, and
        # thus possibly informative, nodes.
        if node.dist == 0 and node.up.name == 'NoName':
            parent = node.up
            # grandparent = parent.up
            # node.detach()
            for child in [child for child in parent.get_children() if
                    child != node]:
                # Remove all children except the node whose branch length is
                # 0.
                # Ensure that the non-zero distances are preserved.
                dist = child.dist
                node.add_child(child.detach(), dist=dist)
            # Make sure that the children are all detached from the parent.
            # This can be removed once tested
            # print(parent.get_children())
            # assert(len(parent.get_children()) == 0)
            # Remove the empty parent, connecting the node to the
            # grandparent with the branch length preserved.
            parent.delete(preserve_branch_length=True)
        # Recurse through all children
        for child in node.get_children():
            child.collapse_null_branches()
        return None

    def format_nodes(self,
            tabfile=None,
            color_column=None,
            size_column=None,
            dict_color=None,
            len_seq=1,
            ):
        """
        Color nodes of `tree` according to data in `tabfile`.

        Parameters
        ----------
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
        len_seq : int, optional
            The integer length of the sequences represented by this tree.
            (default: 1)

        Returns
        -------
        dict_color : dict
            The dictionary of colors used in the colorization of the tree.
            This is especially useful if it was not specified as an input
            argument, because otherwise one has no way to know how to
            correlate colors on the plot with the values they color.

        """
        if dict_color == None:
            _data = [entry.get(color_column) for entry in
                    self.lst_dict_tab_entries]
            self.dict_color = _get_color_dict(_data)
        else:
            self.dict_color = dict_color
        for node in self.traverse():
            node.dist *= len_seq
            try:
                dict_entry = _get_node_entry(node.name,
                        self.lst_dict_tab_entries)
            except ValueError:
                continue
            # get color data (default: 'none')
            color_data = dict_entry.get(color_column)
            color = self.dict_color.get(color_data, 'none')
            node.add_feature('group', color_data)
            # print(node.group)
            # get size data (default: 1, assume single copy)
            size = dict_entry.get(size_column, 1)
            if size == '':
                size = 1
            # set node style
            style = ete2.NodeStyle()
            style['fgcolor'] = color
            style['size'] = _scale_size(int(size))
            node.set_style(style)
            # scale distance to represent mutation length
        return dict_color

    @classmethod
    def print_tree(cls, fname, **kwarg):
        """
        Print a tree read from a `fname` using ete2.

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
        phyfile : str, optional
            If provided, uses the PHYLIP formatted file to determine how long
            a the sequences in the tree are, and then uses that data to convert
            edge lengths into proper mutation counts.
        outputdir : str, optional
            Destination for output tree analysis files. (default: False)

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
        # if phyfile specified: use, otherwise False
        phyfile = kwarg.pop('phyfile',False)
        # if outputdir specified: use, otherwise False
        outputdir = kwarg.pop('outputdir', False)
        # if set not to display, don't, otherwise display tree
        display = kwarg.pop('display', True)
        # print(outputdir)

        tree = cls(fname)
        tree.ete_treestyle = _get_ete_treestyle()
        tree.set_tabfile(tabfile)
        if phyfile:
            tree.set_phyfile(phyfile)
            tree.root_tree()
        else:
            len_seq = 1
            print(('''Since no *.phy file was provided, branch lengths will be''' +
                    ''' distances rather than mutation counts''').format())

        # color, size nodes
        dict_color = tree.format_nodes(
                color_column=color_column,
                size_column=size_column,
                len_seq=tree.len_seq,
                )
        tree.collapse_null_branches()
        lineages_monophyly = tree.find_pure_subtrees()
        n_lineage_monophyly = len(lineages_monophyly)
        print('Number of pure lineages identified {:d}'.format(
            n_lineage_monophyly))
        for iI, (group, lineage) in enumerate(lineages_monophyly):
            print(group)
            print(lineage.get_ascii())
            if outputdir:
                try:
                    os.listdir(outputdir)
                except:
                    os.makedirs(outputdir)
                print('Writing output files...')
                # 4-digit output file basename
                outbasename = os.path.join(outputdir, 'monophyly_' + group +
                        '_{:04d}'.format(iI))
                # print Newick tree to file
                treefname = outbasename + '.tree'
                lineage.write(outfile=treefname)
                # print ascii tree to file
                outfname = outbasename + '.asc'
                open(outfname,'wb').write(lineage.get_ascii())
        lineages_dist = tree.find_distant_subtrees(root_node=tree.germnode)
        n_lineage_dist = len(lineages_dist)
        print('Number of distant lineages identified {:d}'.format(
            n_lineage_dist))
        for iI, lineage in enumerate(lineages_dist):
            print(lineage.get_ascii())
            if outputdir:
                try:
                    os.listdir(outputdir)
                except:
                    os.makedirs(outputdir)
                print('Writing output files...')
                # 4-digit output file basename
                outbasename = os.path.join(outputdir, 'dist_' +
                        '_{:04d}'.format(iI))
                # print Newick tree to file
                treefname = outbasename + '.tree'
                lineage.write(outfile=treefname)
                # print ascii tree to file
                outfname = outbasename + '.asc'
                open(outfname,'wb').write(lineage.get_ascii())
        # show filename
        tree.ete_treestyle.title.add_face(ete2.TextFace(fname), column=0)
        _add_legend(tree.ete_treestyle, tree.dict_color)
        # if outfile specified, use, otherwise just show
        outfile = kwarg.pop('outfile',None)
        if outfile:
            tree.render(outfile,
                    # w=841,
                    # h=1189,
                    units='mm',
                    tree_style=tree.ete_treestyle,
                    )
        elif display:
            tree.show(tree_style=tree.ete_treestyle)
        else:
            return None
        return None

    def find_pure_subtrees(self):
        """
        Identify subtrees of a colored tree whose descendants are pure.

        Returns
        -------
        lst_puresub : list of ete2.tree.TreeNode
            The list of subtrees that are monophyletic by group.
        """
        set_group = {getattr(node, 'group', None) for node in
                self.traverse()}
        lst_puresub = list()
        for group in set_group:
            lineages = list(self.get_monophyletic(values=[group],
                    target_attr='group'))
            print(str(group) + '----------------------------------------------')
            for lineage in lineages:
                print(lineage.get_ascii())
            lst_puresub.extend([(group, lineage) for lineage in lineages])
        return lst_puresub

    @staticmethod
    def _internal_layout(node):
        if node.is_leaf():
            # If terminal node, draws its name
            name_face = ete2.AttrFace("name")
        elif node.name == 'NoName' or not node.name:
            return None
        else:
            # If internal node, draws label with smaller font size
            name_face = ete2.AttrFace("name")
        # Adds the name face to the image at the preferred position
        ete2.faces.add_face_to_node(name_face, node, column=0,
            position="branch-right")

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

def _get_ete_treestyle():
    """
    Returns an ete2 treestyle.

    Returns
    -------
    ete_treestyle : ete2.TreeStyle
    """
    ete_treestyle = ete2.TreeStyle()
    # do not use branch lengths to influence graphical branch lengths
    ete_treestyle.force_topology = True
    # do not rotate, this makes navigation weird
    ete_treestyle.rotation = 0
    # show the branch lengths
    ete_treestyle.show_branch_length = True
    # do not automatically show leaf names, this is handled by layout_fn
    ete_treestyle.show_leaf_name = False
    # circular layout
    ete_treestyle.mode = 'c'
    # custom layout function to be used on each node
    ete_treestyle.layout_fn = GermTree._internal_layout
    # place legend in top-left
    ete_treestyle.legend_position = 1
    return ete_treestyle

def _get_node_entry(nodename, lst_dict_entries):
    for dict_entry in lst_dict_entries:
        if dict_entry['SEQUENCE_ID'][-9:] == nodename:
            return dict_entry
        else: continue
    raise ValueError('Nodename {name} not found'.format(name=nodename))

default_color_list = [
        '#ff0000',
        '#0000ff',
        '#00ff00',
        '#7f7f00',
        '#007f7f',
        '#7f007f',
        ]

def _get_color_dict(lst_col_data):
    set_data = set(lst_col_data)
    n_data = len(set_data)
    n_rand = n_data - len(default_color_list)
    lst_color = default_color_list
    for iI in range(n_rand):
        lst_color.append(_random_color())
    dict_color = dict()
    for iI, data in enumerate(set_data):
        dict_color[data] = lst_color[iI]
    return dict_color

def _random_color():
    tpl_color = np.random.rand(3,)
    return _color_3pl2hex(tpl_color)

def _color_3pl2hex(tpl_color):
    hex_color = '#' + ''.join([hex(int(256*iI))[-2:] for iI in tpl_color])
    return hex_color

def _treeviz_main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Visualize newick trees',
        )
    parser.add_argument('treefile',
            # nargs=1,
            # dest='treefile',
        )
    parser.add_argument('-t', '--tabfile', dest='tabfile')
    parser.add_argument('-r', '--render', dest='outfile', default=None)
    parser.add_argument('-n', '--no-display', dest='display',
            action='store_false',
            help="""
            Use if you do not actually want to display the tree.
            This option is only useful if you are performing subtree
            analyses, because otherwise the whole point is to display
            the tree.
            """,
            )
    parser.add_argument('-c', '--color-column',
            dest='color',
            default='COLOR_GROUP',
            )
    parser.add_argument('-s', '--size-column',
            dest='size',
            default='COPY_NUMBER',
            )
    parser.add_argument('-p', '--phy-file',
            dest='phyfile',
            )
    parser.add_argument('-o', '--output-dir',
            dest='outputdir',
            help="""
            The directory into which to put analysis files, if not
            supplied, analysis is not performed.
            If this argument is supplied with no path, the current
            directory is used for output analysis files.
            """,
            default=None,
            const=os.getcwd(),
            nargs='?',
            )
    argspace = parser.parse_args()
    GermTree.print_tree(argspace.treefile,
        tabfile=argspace.tabfile,
        outfile=argspace.outfile,
        color_column=argspace.color,
        size_column=argspace.size,
        phyfile=argspace.phyfile,
        outputdir=argspace.outputdir,
        display=argspace.display,
        )
    return None

if __name__ == '__main__':
    _main()
