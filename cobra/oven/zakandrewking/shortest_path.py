# -*- coding: utf-8 -*-

# Author: Zachary King, 2014

from cobra import Reaction, Metabolite, Model
from cobra.manipulation.modify import convert_to_irreversible

def shortest_path_metabolites(model, origin_metabolite, max_length=10,
                              solver=None, ignore_metabolites=[]):
    """Find the shortest path between two metabolites. Returns the shortest pathway
    to each metabolite in the model.
    
    A COBRApy implementation of the algorithm for calculated minimum pathway
    distance, described by:

    1. Simeonidis, E., Rison, S. C. G., Thornton, J. M., Bogle, I. D. L. &
    Papageorgiou, L. G. Analysis of metabolic networks using a pathway distance
    metric through linear programming. Metab. Eng. 5, 211–219 (2003).

    Arguments
    ---------

    model: a cobra.Model object

    origin_metabolite: the ID or cobra.Metabolite object to start with

    max_length: the maximum pathway length

    solver: string of solver name. If None is given, the default solver will be
    used.
    
    ignore_metabolites: a list of metabolite IDs to ignore in the simulation. A
    typical list of metabolites to ignore could be: ['atp_c', 'adp_c', 'nadh_c',
    'nad_c', 'h2o_c', 'h2o_p', 'h_c', 'h_p', 'h_e', 'nadph_c', 'nadp_c', 'pi_c',
    'pi_p', 'nh4_c', 'co2_c'].

    """

    try:
        origin_id = origin_metabolite.id
    except AttributeError:
        origin_id = origin_metabolite

    # make the model
    shortest_model = Model("shortest_paths")

    # keep track of old and new reactions
    variable_dict = {}
    
    # add new reactions
    for metabolite in model.metabolites:
        if metabolite.id in ignore_metabolites:
            continue

        # Reaction is a variable
        distance = Reaction(metabolite.id)
        variable_dict[metabolite] = distance
        
        # bounds
        if metabolite.id == origin_id:
            distance.upper_bound = 0
        else:
            distance.upper_bound = max_length
        distance.lower_bound = 0
            
        # all metabolites have objective coefficient of 1
        distance.objective_coefficient = 1
        
    for reaction in model.reactions:
        # constraints on connectivity
        for D_i in [m for m, coeff in reaction.metabolites.items()
                    if m.id not in ignore_metabolites and (coeff < 0 or reaction.reversibility)]:
            for D_j in [m for m, coeff in reaction.metabolites.items()
                        if m.id not in ignore_metabolites and (coeff > 0 or reaction.reversibility)]:
                if D_i is D_j:
                    continue
                # -D_i + D_j <= 1
                # Metabolite is a constraint
                constraint = Metabolite(reaction.id + "_" + D_i.id + "_" + D_j.id)
                constraint._constraint_sense = "L"
                constraint._bound = 1
                # shortest_model.add_metabolite(constraint)
                variable_dict[D_i].add_metabolites({constraint: -1})
                variable_dict[D_j].add_metabolites({constraint: 1})
                
    # add variables
    shortest_model.add_reactions(variable_dict.values())

    solution = shortest_model.optimize(objective_sense="maximize",
                                       solver=solver)
    
    return solution
    
def shortest_path_reactions(model, origin_reaction, max_length=10, solver=None,
                            ignore_metabolites=[]):
    """Find the shortest path between two reactions. Returns the shortest pathway
    to each reaction in the model.
    
    A COBRApy implementation of the algorithm for calculated minimum pathway
    distance, described by:

    1. Simeonidis, E., Rison, S. C. G., Thornton, J. M., Bogle, I. D. L. &
    Papageorgiou, L. G. Analysis of metabolic networks using a pathway distance
    metric through linear programming. Metab. Eng. 5, 211–219 (2003).

    Arguments
    ---------

    model: a cobra.Model object

    origin_reaction: the ID or cobra.Reaction object to start with

    max_length: the maximum pathway length

    solver: string of solver name. If None is given, the default solver will be
    used.

    ignore_metabolites: a list of metabolite IDs to ignore in the simulation. A
    typical list of metabolites to ignore could be: ['atp_c', 'adp_c', 'nadh_c',
    'nad_c', 'h2o_c', 'h2o_p', 'h_c', 'h_p', 'h_e', 'nadph_c', 'nadp_c', 'pi_c',
    'pi_p', 'nh4_c', 'co2_c'].

    """

    try:
        origin_id = origin_reaction.id
    except AttributeError:
        origin_id = origin_reaction

    # make the model
    shortest_model = Model("shortest_paths")

    # keep track of old and new reactions
    variable_dict = {}
    
    # add new reactions
    for reaction in model.reactions:
        # Reaction is a variable
        distance = Reaction(reaction.id)
        variable_dict[reaction] = distance
        
        # bounds
        if reaction.id == origin_id:
            distance.upper_bound = 0
        else:
            distance.upper_bound = max_length
        distance.lower_bound = 0
            
        # all reactions have objective coefficient of 1
        distance.objective_coefficient = 1
        
    for metabolite in model.metabolites:
        if metabolite.id in ignore_metabolites:
            continue
        # constraints on connectivity
        for D_i in [r for r in metabolite.reactions if r._metabolites[metabolite] > 0 or r.reversibility]:
            for D_j in [r for r in metabolite.reactions if r._metabolites[metabolite] < 0 or r.reversibility]:
                if D_i is D_j:
                    continue
                # -D_i + D_j <= 1
                # Metabolite is a constraint
                constraint = Metabolite(metabolite.id + "_" + D_i.id + "_" + D_j.id)
                constraint._constraint_sense = "L"
                constraint._bound = 1
                # shortest_model.add_metabolite(constraint)
                variable_dict[D_i].add_metabolites({constraint: -1})
                variable_dict[D_j].add_metabolites({constraint: 1})
                
    # add variables
    shortest_model.add_reactions(variable_dict.values())

    solution = shortest_model.optimize(objective_sense="maximize", 
                                       solver=solver)
    
    return solution
