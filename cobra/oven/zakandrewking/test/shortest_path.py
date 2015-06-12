from unittest import TestCase, TestLoader, TextTestRunner, skipIf

import sys

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model
    from cobra.oven.zakandrewking.shortest_path import (shortest_path_metabolites,
                                                        shortest_path_reactions)
    sys.path.pop(0)
# else:
#     from . import create_test_model
#     from ..shortest_path import (shortest_path_metabolites,
#                                  shortest_path_reactions)

class TestCobraShortestPath(TestCase):
    """Test the simulation functions in cobra.oven.zakandrewking.shortest_path"""

    def setUp(self):
        self.model = create_test_model(test_pickle='ecoli')

    def test_shortest_path_reactions(self):
        model = self.model
        ignore_metabolites = ['atp_c', 'adp_c', 'nadh_c', 'nad_c', 'h2o_c',
                              'h2o_p', 'h_c', 'h_p', 'nadph_c', 'nadp_c',
                              'pi_c', 'pi_p', 'nh4_c', 'co2_c']
        solution = shortest_path_reactions(model, 'LALDO2x',
                                           ignore_metabolites=ignore_metabolites)
        self.assertEqual(solution.x_dict['LALDO2x'], 0.0)
        self.assertEqual(solution.x_dict['LCARR'], 1.0)
        self.assertEqual(solution.x_dict['MGSA'], 10.0)
        self.assertEqual(solution.x_dict['12PPDRtpp'], 2.0)

    def test_shortest_path_metabolites(self):
        model = self.model
        ignore_metabolites = ['atp_c', 'adp_c', 'nadh_c', 'nad_c', 'h2o_c',
                              'h2o_p', 'h_c', 'h_p', 'nadph_c', 'nadp_c',
                              'pi_c', 'pi_p', 'nh4_c', 'co2_c']
        solution = shortest_path_metabolites(model, 'lald__D_c',
                                             ignore_metabolites=ignore_metabolites)
        self.assertEqual(solution.x_dict['lald__D_c'], 0.0)
        self.assertEqual(solution.x_dict['12ppd__R_c'], 1.0)
        self.assertEqual(solution.x_dict['mthgxl_c'], 10.0)
        self.assertEqual(solution.x_dict['12ppd__R_p'], 2.0)

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
