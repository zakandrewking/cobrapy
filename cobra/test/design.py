from unittest import TestCase, TestLoader, TextTestRunner, skipIf

import sys

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra.design import *
    from cobra.design.design_algorithms import _add_decision_variable
    from cobra.solvers import get_solver_name
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from ..design import *
    from ..design.design_algorithms import _add_decision_variable
    from ..solvers import get_solver_name

try:
    solver = get_solver_name(mip=True)
except:
    no_mip_solver = True
else:
    no_mip_solver = False


class TestDesignAlgorithms(TestCase):
    """Test functions in cobra.design"""

    def test_dual_twice(self):
        model = create_test_model("textbook")
        self.assertAlmostEqual(model.optimize("maximize").f, 0.874, places=3)
        dual = dual_problem(model)
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.874, places=3)
        dual2 = dual_problem(dual, objective_sense="minimize")
        self.assertAlmostEqual(dual2.optimize("minimize").f, -0.874, places=3)

    def test_dual_minimize(self):
        model = create_test_model("textbook")
        model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound = 0.6
        for reaction in model.reactions:
            reaction.objective_coefficient = int(reaction.id == "GAPD")
        self.assertAlmostEqual(model.optimize("minimize").f, 7.524, places=3)
        dual = dual_problem(model, objective_sense="minimize")
        self.assertAlmostEqual(dual.optimize("minimize").f, -7.524, places=3)

    def test_dual_integer_vars_as_lp(self):
        model = create_test_model("textbook")
        var = _add_decision_variable(model, "AKGDH")
        self.assertAlmostEqual(model.optimize("maximize").f, 0.874, places=3)
        # as lp: make integer continuous, set to 1
        dual = dual_problem(model, "maximize", [var.id], copy=True)
        r = dual.reactions.get_by_id(var.id)
        r.variable_kind = "continuous"
        r.lower_bound = r.upper_bound = 1
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.874, places=3)
        r.lower_bound = r.upper_bound = 0
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.858, places=3)

    @skipIf(no_mip_solver, "no MILP solver found")
    def test_dual_integer_vars_as_mip(self):
        # mip
        model = create_test_model("textbook")
        var = _add_decision_variable(model, "AKGDH")
        dual = dual_problem(model, "maximize", [var.id], copy=True)
        var_in_dual = dual.reactions.get_by_id(var.id)

        # minimization, so the optimal value state is to turn off AKGDH
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.858, places=3)

        # turn off AKGDH in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 1
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.874, places=3)

        # turn on AKGDH in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 0
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.858, places=3)

    @skipIf(no_mip_solver, "no MILP solver found")
    def test_optknock(self):
        model = create_test_model("textbook")
        model.reactions.get_by_id("EX_o2_e").lower_bound = 0
        knockable_reactions = ["ACKr", "AKGDH", "ACALD", "LDH_D"]
        optknock_problem = set_up_optknock(model, "EX_lac__D_e",
                                           knockable_reactions, n_knockouts=2,
                                           copy=False)
        solution = run_optknock(optknock_problem)
        self.assertIn("ACKr", solution.knockouts)
        self.assertIn("ACALD", solution.knockouts)
        self.assertAlmostEqual(solution.f, 17.891, places=3)

    @skipIf(no_mip_solver, "no MILP solver found")
    def test_robustknock(self):
        model = create_test_model("textbook")
        model.reactions.get_by_id("EX_o2_e").lower_bound = 0

        # add LDH_L and EX_lac__L_c
        ldh_l = Reaction("LDH_L")
        ldh_l.add_metabolites({
            model.metabolites.get_by_id("pyr_c"): -1,
            model.metabolites.get_by_id("nadh_c"): -1,
            model.metabolites.get_by_id("h_c"): -1,
            Metabolite("lac__L_c"): 1,
            model.metabolites.get_by_id("nad_c"): 1,
        })
        model.add_reaction(ldh_l)
        l_lact2 = Reaction("L_LACt2")
        l_lact2.add_metabolites({
            model.metabolites.get_by_id("lac__L_c"): -1,
            model.metabolites.get_by_id("h_c"): -1,
            Metabolite("lac__L_e"): 1,
            model.metabolites.get_by_id("h_e"): 1,
        })
        model.add_reaction(l_lact2)
        ex_lac_l = Reaction("EX_lac__L_e")
        ex_lac_l.add_metabolites({
            model.metabolites.get_by_id("lac__L_e"): -1,
        })
        model.add_reaction(ex_lac_l)

        # non-unique
        knockable_reactions = [] #"ACKr", "AKGDH", "ACALD", "GLUN", "SUCCt3"]
        robustknock_problem = set_up_robustknock(model, "EX_etoh_e", #"EX_lac__L_e",
                                                 knockable_reactions, n_knockouts=3,
                                                 copy=True)
        solution = run_robustknock(robustknock_problem, 'gurobi')
        self.assertEqual(solution.f, None)

        # unique
        knockable_reactions = ["ACKr", "AKGDH", "ACALD", "GLUN", "SUCCt3", "LDH_D"]
        robustknock_problem = set_up_robustknock(model, "EX_lac__L_e",
                                                 knockable_reactions,
                                                 n_knockouts=3, copy=True)
        solution = run_robustknock(robustknock_problem, 'gurobi')
        self.assertIn("ACKr", solution.knockouts)
        self.assertIn("ACALD", solution.knockouts)
        self.assertIn("LDH_D" in solution.knockouts)
        self.assertAlmostEqual(solution.f, 17.891, places=3)

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
