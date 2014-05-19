"""
Test tree_viz module
"""
import unittest

import ete2

import tree_viz

class TestEteCollapse(unittest.TestCase):

    def setUp(self):
        self.tree = ete2.Tree('(A:0.0,(B:1,(C:1,D:1):0.0):0.5);')
        self.collapsed = ete2.Tree('((B:1,C:1,D:1):0.5)A:0.0;')
        tree_viz.collapse_null_branches_ete(self.tree)
        self.B = self.tree.search_nodes(name='B')[0]

    def test_collapse_length(self):
        """
        Check that all lengths check out.
        """
        self.assertEqual(self.tree.get_distance(self.B), 1.5)

    def test_collapse_topo_length(self):
        """
        Check that all lengths check out.
        """
        self.assertEqual(
                self.tree.get_distance(self.B, topology_only=True),
                1)

    def test_equal_collapsed(self):
        """
        Check that the tree equals the expected output tree.
        """
        rf_metric=self.tree.robinson_foulds(self.collapsed,
                unrooted_trees=True)[0]
        print(rf_metric)
        self.assertEqual(rf_metric, 0)

if __name__ == '__main__':
    unittest.main()

