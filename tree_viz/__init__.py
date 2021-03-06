#!/usr/bin/env python
import csv
import os

import numpy as np
# import Bio
try:
    ## First try ete2.
    import ete2
except ImportError as exc:
    print('''Module 'ete2' not found, looking for 'ete_dev'...''')
    try:
        ## Since ete2 was not found, try to import ete_dev as ete2.
        import ete_dev as ete2
        print('''Using 'ete_dev'.''')
    except ImportError:
        ## If neither are found, you are in trouble.
        print('''Module 'ete_dev' not found, make sure you are using the
                right python.''')

import fasta2tab

LINEAGE_COL_NAME = 'lineage0'

def set_header_rev(header_rev):
    ## Switch some global variables based on header rev
    ## Could there be a better way to do this?
    global id_col
    global seq_col
    global clone_col
    global dmask_col
    if header_rev == 0:
        id_col = 'SEQUENCE_ID'
        seq_col = 'SEQUENCE'
        clone_col = 'CLONE'
        dmask_col = 'GERMLINE_GAP_DMASK'
    elif header_rev == 1:
        id_col = 'seqID'
        seq_col = 'sequence'
        clone_col = 'cloneID'
        dmask_col = 'germline'
    else:
        raise ValueError(
            'You have selected an invalid header_rev: {:s}'.format(header_rev))

class GermTree(ete2.coretype.tree.TreeNode):
    """
    Class to encapsulate the many actions that I need to make on
    ete2.Tree instances.
    """
    def __init__(self, *args, **kwarg):
        ete2.coretype.tree.TreeNode.__init__(self, *args, **kwarg)
        # initialize value of sequence length to 1
        self.len_seq = 1

    def set_tabfile(self, lst_tabfile):
        """
        Read in data from a tabfile to which any method may refer.
        """
        if not hasattr(self, 'lst_phy_names'):
            raise Exception('''First use 'set_phyfile'.''')
        self.lst_dict_tab_entries = list()
        for fname in lst_tabfile:
            with open(fname, 'rU') as f:
                reader = csv.DictReader(f,delimiter='\t')
                self.lst_dict_tab_entries.extend([row for row in reader
                    if row[id_col][-9:] in self.lst_phy_names])
        print("""Number of TAB entries read into memory: {:d}""".format(
            len(self.lst_dict_tab_entries)))
        return None

    def set_phyfile(self, phyfile):
        """
        Read in data from a PHYLIP formatted *.phy file.
        This will set the sequence length properly, as well as root the
        tree at the first sequence found in the PHY file.
        This also sets the PHYLIP entries used so that only the relevant
        lines from the TAB file will be read into memory.
        """
        with open(phyfile, 'rt') as phyfileobj:
            headrow = phyfileobj.readline()
            self.lst_phy_entries = [row.split() for row in phyfileobj]
        self.lst_phy_names = [tpl[0] for tpl in self.lst_phy_entries]
        self.len_seq = int(headrow.split()[1])
        print('Set sequence length from {file0}'.format(file0=phyfile))
        self.germname = self.lst_phy_names[0]

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

    def analyze_tree(self, outputdir,
            tab_out=False,
            verbose=False,
            ):
        """
        Perform analyses on the tree, write the output files in
        outputdir.

        Parameters
        ----------
        outputdir : str
            The name of the directory into which to place analysis
            files.
        tab_out : bool, optional
            Whether or not to output new TAB files based on the
            classifications arrived at by tree analysis. (default:
            False)
        verbose : bool, optional
            Whether or not to print identified subtrees to terminal
            (default: False)
        """
        self.lineages_monophyly = self.find_pure_subtrees()
        n_lineage_monophyly = len(self.lineages_monophyly)
        print('Number of pure lineages identified {:d}'.format(
            n_lineage_monophyly))
        for iI, (group, lineage) in enumerate(self.lineages_monophyly):
            if verbose:
                print(group)
                print(lineage.get_ascii())
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
        self.lineages_dist = self.find_distant_subtrees(root_node=self.germnode)
        lst_tpl_lineageid = self.tab_classification(
                self.lineages_dist,
                outputdir=outputdir,
                # outname='tabout.tab',
                )
        n_lineage_dist = len(self.lineages_dist)
        print('Number of distant lineages identified {:d}'.format(
            n_lineage_dist))
        for iI, (lineage_id, lineage) in enumerate(lst_tpl_lineageid):
            if verbose:
                print(lineage.get_ascii())
            try:
                os.listdir(outputdir)
            except:
                os.makedirs(outputdir)
            print('Writing output files...')
            # 4-digit output file basename
            outbasename = os.path.join(outputdir,
                    ('dist_' + lineage_id).format(iI))
            # print Newick tree to file
            treefname = outbasename + '.tree'
            lineage.write(outfile=treefname)
            # print ascii tree to file
            outfname = outbasename + '.asc'
            open(outfname,'wb').write(lineage.get_ascii())

    def tab_classification(self,
            lst_lineages,
            **kwarg):
        # lst_dict_tab_entries_new = list()
        outputdir = kwarg.get('outputdir', os.getcwd())
        outname = kwarg.get('outname', 'tabout.tab')
        fname_tab = os.path.join(outputdir, outname)
        classification_column = LINEAGE_COL_NAME
        lst_tpl_lineageid = list()
        for lineage in lst_lineages:
            lineage_id = get_lineage_id()
            for node in lineage:
                try:
                    node_entry = _get_node_entry(node.name,
                            self.lst_dict_tab_entries)
                    node_entry[classification_column] = lineage_id
                    pass
                except:
                    print('Skipping internal node...')
                    pass
                pass
            # Append each labelled lineage and its label list to be
            # returned
            lst_tpl_lineageid.append((lineage_id, lineage,))
            pass
        # Write tabfile
        # fasta2tab.write_tab_file(fname_tab, self.lst_dict_tab_entries)
        write_tabfile(fname_tab, self.lst_dict_tab_entries)
        return lst_tpl_lineageid

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
        # print('Distance between node and root: {:0.2f}'.format(distance))
        # print('dist_lim: {:0.2f}'.format(dist_lim))
        if distance >= dist_lim:
            # If this node is already far enough away from root, return it,
            # skipping all descendants.
            print('Found distant subtree, skipping descendants')
            # print(self)
            return [self]
        else:
            # If this node is still too close to root, iterate through each
            # child node to check for subtrees that start far enough away
            # from the root.
            print('Attempting to find distant subtrees among children nodes')
            lst_distsub = list()
            for child in self.get_children():
                # print(child)
                # print(root_node)
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
        if (node.dist == 0 and hasattr(node.up, 'name') and
            node.up.name == 'NoName'):
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
            **kwarg):
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
        logscale : bool, optional
        areascale : bool, optional

        Returns
        -------
        dict_color : dict
            The dictionary of colors used in the colorization of the tree.
            This is especially useful if it was not specified as an input
            argument, because otherwise one has no way to know how to
            correlate colors on the plot with the values they color.

        """
        logscale = kwarg.pop('logscale', True)
        areascale = kwarg.pop('areascale', False)
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
            except ValueError as exc:
                # print(""" Node info not found in tab file.""")
                # print(exc)
                node.add_feature('group', 'internal')
                continue
            # get color data (default: 'none')
            color_data = dict_entry.get(color_column)
            color = self.dict_color.get(color_data, 'none')
            # This is what is causing problems for finding monophyletic
            # lineages. The feature 'group' is overwritten after being
            # set to 'internal' if it is an inferred node.
            node.add_feature('group', color_data)
            # print(node.group)
            # get size data (default: 1, assume single copy)
            size = dict_entry.get(size_column, 1)
            if size == '':
                size = 1
            # set node style
            style = ete2.NodeStyle()
            style['fgcolor'] = color
            style['size'] = _scale_size(
                    int(size),
                    logscale=logscale,
                    areascale=areascale,
                    )
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
        correct_lengths : bool, optional
            Whether or not to use the information in the *.phy file to
            correct branch lengths expressed as distances to branch
            lengths as mutation counts. (default: False)

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
        correct_lengths = kwarg.pop('correct_lengths', False)
        # if true, display legend, else not
        legend = kwarg.pop('legend', True)
        # if true, display node names, else not
        shownames = kwarg.pop('shownames', True)
        # A few options for how to process size
        logscale = kwarg.pop('logscale', True)
        areascale = kwarg.pop('areascale', False)

        tree = cls(fname)
        tree.ete_treestyle = _get_ete_treestyle(shownames=shownames)
        tree.set_phyfile(phyfile)
        tree.set_tabfile(tabfile)
        tree.root_tree()
        if not correct_lengths:
            tree.len_seq = 1
            print("""Not correcting branch lengths...""")
            # print(('''Branch lengths will be''' +
            #         ''' distances rather than mutation counts''').format())

        # color, size nodes
        dict_color = tree.format_nodes(
                color_column=color_column,
                size_column=size_column,
                len_seq=tree.len_seq,
                logscale=logscale,
                areascale=areascale,
                )
        tree.collapse_null_branches()
        if outputdir:
            tree.analyze_tree(outputdir)
        if legend:
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
            lineages = list(self.get_monophyletic(values=[group,
            'internal'],
                    target_attr='group'))
            # print(str(group) + '----------------------------------------------')
            # for lineage in lineages:
            #     print(lineage.get_ascii())
            lst_puresub.extend([(group, lineage) for lineage in lineages])
        return lst_puresub

    @staticmethod
    def _internal_layout(node):
        # This function is currently used ONLY to draw the name on a node,
        # that way it can be specified or not depending on the user's
        # preference.
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

def write_tabfile(fname_tab, lst_dict_entries):
    """
    Write a new TAB file from a list of dict entries.
    """
    tpl_fields = list(fasta2tab.tpl_cols)
    for entry in lst_dict_entries:
        tpl_fields.extend([field for field in entry.keys()
            if field not in tpl_fields])
        pass
    with open(fname_tab, 'w') as f:
        tab_writer = csv.DictWriter(f, tpl_fields,
                extrasaction='ignore',
                delimiter='\t')
        tab_writer.writeheader()
        for entry in lst_dict_entries:
            tab_writer.writerow(entry)
            pass
        pass
    return None


    return clean_tree

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
    areascale : bool, optional
        Boolean whether or not to scale area as opposed to linear size.
        (default: `False`)
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
    areascale = kwarg.pop('areascale', False)
    if areascale:
        size = np.sqrt(size)
    if logscale:
        size = int(basesize * (1 + np.log(size)))
    else:
        size = int(basesize * size)
    return size

def _get_ete_treestyle(shownames=True):
    """
    Returns an ete2 treestyle.

    Returns
    -------
    ete_treestyle : ete2.TreeStyle
    """
    ete_treestyle = ete2.TreeStyle()
    # do not use branch lengths to influence graphical branch lengths
    ete_treestyle.force_topology = True
    # it doesn't make any sense to show a scale bar when forcing
    # topology
    ete_treestyle.show_scale = not ete_treestyle.force_topology
    # do not rotate, this makes navigation weird
    ete_treestyle.rotation = 0
    # show the branch lengths
    ete_treestyle.show_branch_length = True
    # do not automatically show leaf names, this is handled by layout_fn
    ete_treestyle.show_leaf_name = False
    # circular layout
    ete_treestyle.mode = 'c'
    # custom layout function to be used on each node
    # currently it just prints the name
    if shownames:
        ete_treestyle.layout_fn = GermTree._internal_layout
    # place legend in top-left
    ete_treestyle.legend_position = 1
    return ete_treestyle

def _get_node_entry(nodename, lst_dict_entries):
    for dict_entry in lst_dict_entries:
        if dict_entry[id_col][-9:] == nodename:
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

def _get_lineage_id(**kwarg):
    # Eventually, this will be a different function for generating lineage
    # IDs
    return None
# get_lineage_id = _get_lineage_id
get_lineage_id = fasta2tab.get_uuid

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
    parser.add_argument('-t', '--tabfile', dest='tabfile',
            nargs='*',
            help="""
            Provide a TAB file which should contain all of the relevant
            information about each sequence in the tree, especially the
            columns specified by COLOR_GROUP and COPY_NUMBER.
            Multiple TAB files can be provided.
            """,
            )
    parser.add_argument('-p', '--phy-file',
            dest='phyfile',
            help="""
            The PHY file to use to set the length of the sequences, as
            well as the root node, which is assumed to be the first
            sequence found in the PHY file.
            """,
            )
    parser.add_argument('-c', '--color-column',
            dest='color',
            default='COLOR_GROUP',
            help="""
            The column of the TAB file to use to color the tree's nodes.
            A color will be assigned to each unique value in the color
            column, so take care not to use a column that contains data
            unique to each sequence, otherwise the colors will be
            meaningless, as will be any analyses performed.
            This option should usually be provided, since the default is
            only useful if the TAB file was tailor-made to suit the
            default. (default: COLOR_GROUP)
            """,
            )
    parser.add_argument('-s', '--size-column',
            dest='size',
            default='COPY_NUMBER',
            help="""
            The column of the TAB file to use to inform the sizes of the
            nodes displayed.
            There is no good reason for a user to employ this option,
            unless they are dealing with a TAB file that has a different
            label for the COPY_NUMBER column. (default: COPY_NUMBER)
            """,
            )
    parser.add_argument('-A', '--area-scale', dest='areascale',
            action='store_true',
            help="""
            Whether to scale node sizes to area rather than
            diameter.
            """,
            )
    parser.add_argument('-i', '--linear-scale', dest='logscale',
            action='store_false',
            help="""
            Whether to scale node sizes linearly by copy number.
            """,
            )
    parser.add_argument('-r', '--render', dest='outfile', default=None,
            help="""
            If provided, render the graphical tree as a file.
            This option can be provided instead of showing the tree
            interactively.
            In future versions, this option may be able to be provided
            in addition to interactively displaying the tree.
            """,
            )
    parser.add_argument('-n', '--no-display', dest='display',
            action='store_false',
            help="""
            Use if you do not actually want to display the tree.
            This option is only useful if you are performing subtree
            analyses, because otherwise the whole point is to display
            the tree.
            """,
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
    parser.add_argument('-l', '--length-correct', dest='correct_lengths',
            action='store_true',
            help="""
            Use if branch lengths in the target tree are expressed as
            distances rather than mutation counts. (default: False)
            You must take care to use this option only when appropriate;
            misuse could result in trees with nonsensically long branch
            lengths, and worse, mistaken analyses based on mutation
            counts.
            """,
            )
    parser.add_argument('--no-legend', dest='legend',
            action='store_false',
            help="""
            Do not display the legend.
            """,
            )
    parser.add_argument('--no-names', dest='shownames',
            action='store_false',
            help="""
            Do not display the node names.
            """,
            )
    parser.add_argument('--header',
            dest='header',
            default=1,
            type=int,
            help="""
            If specified, one can change the version of headers to use.
            Old-style headers are 0, new-style headers are 1 (default).
            """,
            )
    argspace = parser.parse_args()
    ## Choose which column names to use for PHY file based on header
    ## rev.
    set_header_rev(argspace.header)
    print('''Using the following files as tab files''')
    for fname in argspace.tabfile:
        print(fname)
    GermTree.print_tree(argspace.treefile,
        tabfile=argspace.tabfile,
        outfile=argspace.outfile,
        color_column=argspace.color,
        size_column=argspace.size,
        phyfile=argspace.phyfile,
        outputdir=argspace.outputdir,
        display=argspace.display,
        correct_lengths=argspace.correct_lengths,
        legend=argspace.legend,
        shownames=argspace.shownames,
        logscale=argspace.logscale,
        areascale=argspace.areascale,
        )
    return None

if __name__ == '__main__':
    _main()
