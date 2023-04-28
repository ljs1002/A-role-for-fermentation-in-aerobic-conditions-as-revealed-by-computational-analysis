from build_roots_model import *
import pandas as pd
from cobra import Reaction, Metabolite
from cobra.flux_analysis.parsimonious import pfba
from cobra.flux_analysis.variability import *
from cobra import io
import numpy as np
import pickle as pkl
# from datetime import datetime, date, time
def mainEBC2(fixphloem=1,testrange = None):
    # Load model
    rootsdir = '/Users/user/Documents/Ox/roots/'
    orootsdir = '/Users/user/Documents/Ox/roots/roots-model/'

    # modelpwy='final_roots_model.xml'
    # modelpwy='tissue_dnut_model.xml'
    # modelpwy='roots_model-1-1-1.xml'
    modelpwy='roots_model_current.xml'
    rmodel=io.read_sbml_model(rootsdir+modelpwy)
    print('loaded')
    temp=rootspFBAobj(rmodel)
    temp=rootsConstrProtons(rmodel)
    print('constrs added')
    # rmodel.objective = 'Phloem_import_reaction'
    biomass_rxns = [rxn.id for rxn in rmodel.reactions if 'Biomass_expanding_cell' in rxn.id]

    solo=rmodel.optimize()
    rmodel.reactions.Phloem_carbon_import.lower_bound=solo['Phloem_carbon_import']
    rmodel.reactions.Phloem_carbon_import.upper_bound=solo['Phloem_carbon_import']
    # fva_sol=flux_variability_analysis(rmodel)
    data={'reaction_id':[x.id for x in rmodel.reactions],
        'cell':['_'.join(x.id.split('_')[-2:]) if ('_l' in x.id or '_d' in x.id) else '' for x in rmodel.reactions],
        'subsystem':[x.notes['SUBSYSTEM']  if 'SUBSYSTEM' in str(x.notes) else '' for x in rmodel.reactions],
        'reaction':[x.reaction for x in rmodel.reactions],
        '1-1-1':[solo[x.id] for x in rmodel.reactions]}
        # 'fva_min':[fva_sol.minimum[x.id] for x in rmodel.reactions],
        # 'fva_max':[fva_sol.maximum[x.id] for x in rmodel.reactions]}
    

    # Change biomass proportions
    if not testrange:
        test_rats = [1.2,1.5]
    else:
        test_rats = np.linspace(1,testrange[0],testrange[1])
    # for ii in test_rats:
    #     for jj in test_rats:
    #         for kk in test_rats:
    #             for rxn in biomass_rxns:
    #                 new_model=rmodel.copy()
    #                 changeBiomassComposition(new_model.reactions.get_by_id(rxn),[ii,jj,kk])
    #                 sol=pfba(new_model)
    #                 data['-'.join([str(x) for x in [ii,jj,kk]])]=[sol[x.id] for x in rmodel.reactions]

    new_model=rmodel
    # with open('roots_'+label+'.pkl','wb') as file:
    #     pkl.dump(sol,file)
    original_biomass = [new_model.reactions.get_by_id(biomass_rxns[3]).get_coefficient(y.id) for y in new_model.reactions.get_by_id(biomass_rxns[3]).metabolites]
    for ii in range(3):
        for jj in test_rats:
            these_rats = [1,1,1]
            these_rats[ii]=jj
            for rxn in biomass_rxns:
                changeBiomassComposition_eBC2(new_model.reactions.get_by_id(rxn),these_rats)
            mid_biomass = [new_model.reactions.get_by_id(biomass_rxns[3]).get_coefficient(y.id) for y in new_model.reactions.get_by_id(biomass_rxns[3]).metabolites]
            mid_biomass = [mid_biomass[it]/original_biomass[it] for it in range(3)]
            sol=new_model.optimize()
            print(these_rats)
            data['-'.join([str(x) for x in these_rats])]=[sol[x.id] for x in rmodel.reactions]
            for rxn in biomass_rxns:
                resetBiomassComposition_eBC2(new_model.reactions.get_by_id(rxn),these_rats)
    new_model=add_maint_constraints(new_model)
    late_biomass=[new_model.reactions.get_by_id(biomass_rxns[3]).get_coefficient(y.id) for y in new_model.reactions.get_by_id(biomass_rxns[3]).metabolites]
    # if late_biomass != original_biomass:
    #     print('Rescale problem:')
    #     print(these_rats)
    #     print(original_biomass)
    #     print(mid_biomass)
    #     print(late_biomass)
    print('reset attempt')
    ub={}
    lb={}
    for rxn in new_model.reactions:
        ub[rxn.id]=rxn.upper_bound
        lb[rxn.id]=rxn.lower_bound
        rxn.upper_bound=solo[rxn.id]
        rxn.lower_bound=solo[rxn.id]
    sol_temp=new_model.optimize()
    for rxn in new_model.reactions:
        rxn.upper_bound=ub[rxn.id]
        rxn.lower_bound=lb[rxn.id]
    print('reset over')
        
    for ii in np.linspace(0,1.8,5):
        new_model.reactions.get_by_id("ATPase_balancer").lower_bound=ii
        sol=new_model.optimize()
        data[ii]=[sol[x.id] for x in rmodel.reactions]
    df=pd.DataFrame(data)
    # now = datetime.now().strftime('%Y_%m_%d')
    df.to_csv('/Users/user/Documents/Ox/roots/Spreadsheets/BiomassComposition2'+fixphloem+'.csv')    

def changeBiomassComposition_eBC2(rxn,ratios):
    # original_rxn = [rxn.get_coefficient(y.id) for y in rxn.metabolites]
    cellulose_met = [met for met in rxn.metabolites if 'Cellulose' in met.id][0]
    lipid_met = [met for met in rxn.metabolites if 'L_PHOSPHATIDATE' in met.id][0]
    protein_met = [met for met in rxn.metabolites if 'Protein' in met.id][0]
    mets = [cellulose_met,lipid_met,protein_met]
    # original_rxn = [rxn.get_coefficient(met) for met in mets]
    
    # print(rxn.reaction)
    # print('cBC: ',rxn.get_coefficient(cellulose_met),ratios[0],(ratios[0]-1),rxn.get_coefficient(cellulose_met)*(ratios[0]-1))
    rxn.add_metabolites({cellulose_met:rxn.get_coefficient(cellulose_met)*(ratios[0]-1),
    lipid_met:rxn.get_coefficient(lipid_met)*(ratios[1]-1),
    protein_met:rxn.get_coefficient(protein_met)*(ratios[2]-1)})
    # new_rxn = [rxn.get_coefficient(met) for met in mets]
    # print(rxn.reaction)
    # print([rxn.get_coefficient(met) for met in mets])
    # print(mets)
    # if [round(new_rxn[it]/original_rxn[it],4) for it in range(len(original_rxn))]!=ratios:
    #     print('cBC problem!',rxn.id)
    #     print([round(new_rxn[it]/original_rxn[it],4) for it in range(len(original_rxn))])
    #     print(ratios)
    #     print(original_rxn)
    #     print(new_rxn)

def resetBiomassComposition_eBC2(rxn,ratios):
    # original_rxn = [rxn.get_coefficient(y.id) for y in rxn.metabolites]
    cellulose_met = [met for met in rxn.metabolites if 'Cellulose' in met.id][0]
    lipid_met = [met for met in rxn.metabolites if 'L_PHOSPHATIDATE' in met.id][0]
    protein_met = [met for met in rxn.metabolites if 'Protein' in met.id][0]
    rxn.add_metabolites({cellulose_met:rxn.get_coefficient(cellulose_met)*(1-ratios[0])/ratios[0],
    lipid_met:rxn.get_coefficient(lipid_met)*(1-ratios[1])/ratios[1],
    protein_met:rxn.get_coefficient(protein_met)*(1-ratios[2])/ratios[2]})
    # new_rxn = [rxn.get_coefficient(cellulose_met),rxn.get_coefficient(lipid_met),rxn.get_coefficient(protein_met)]
    # if [round(original_rxn[it]/new_rxn[it],4) for it in range(len(original_rxn))]!=ratios:
    #     print('rBC problem!',rxn.id)
    #     print([round(original_rxn[it]/new_rxn[it],4) for it in range(len(original_rxn))])
    #     print(ratios)
    #     print(original_rxn)
    #     print(new_rxn)
        
def changeBiomassComposition2(model,rxn,ratios):
    cellulose_met = [met for met in rxn.metabolites if 'Cellulose' in met.id][0]
    lipid_met = [met for met in rxn.metabolites if 'L_PHOSPHATIDATE' in met.id][0]
    protein_met = [met for met in rxn.metabolites if 'Protein' in met.id][0]
    new_rxn = rxn.copy()
    new_rxn.add_metabolites({cellulose_met:rxn.get_coefficient(cellulose_met)*ratios[0],
    lipid_met:rxn.get_coefficient(lipid_met)*(ratios[1]),
    protein_met:rxn.get_coefficient(protein_met)*(ratios[2])})
    model.add_rxn(new_rxn)

def resetBiomassComposition2(model,rxn,ratios):
    cellulose_met = [met for met in rxn.metabolites if 'Cellulose' in met.id][0]
    lipid_met = [met for met in rxn.metabolites if 'L_PHOSPHATIDATE' in met.id][0]
    protein_met = [met for met in rxn.metabolites if 'Protein' in met.id][0]
    rxn.add_metabolites({cellulose_met:rxn.get_coefficient(cellulose_met)*(1-ratios[0])/ratios[0],
    lipid_met:rxn.get_coefficient(lipid_met)*(1-ratios[1])/ratios[1],
    protein_met:rxn.get_coefficient(protein_met)*(1-ratios[2])/ratios[2]})

from cleaner_roots import *
def add_maint_constraints(rmodel_in,split=0):
    if split:
        rmodel=quietcopy(rmodel_in)
    else:
        rmodel = rmodel_in

    rxns = ['ATPase_tx','NADPHoxc_tx','NADPHoxp_tx','NADPHoxm_tx']
    cells = ['_'+rxn.id.split('_')[-1] for rxn in rmodel.reactions if 'CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c' in rxn.id]
    for cell in cells:
        ATP_NAD_met = Metabolite('ATPase_NADPHoxidase_constraint_c'+cell,compartment='c'+cell)
        rmodel.add_metabolites({ATP_NAD_met})

        rmodel.reactions.get_by_id(rxns[0]+cell).add_metabolites({ATP_NAD_met:1})
        for rxn in rxns[1:]:
            rmodel.reactions.get_by_id(rxn+cell).add_metabolites({ATP_NAD_met:-3})

        # NAD_balcp_constr = rmodel.problem.Constraint(rmodel.reactions.get_by_id(rxns[1]+cell).forward_variable-rmodel.reactions.get_by_id(rxns[2]+cell).forward_variable,lb=0,ub=0,name='NAD_balcp'+cell)
        # rmodel.add_cons_vars([NAD_balcp_constr])
        NAD_balcm_constr = rmodel.problem.Constraint(rmodel.reactions.get_by_id(rxns[1]+cell).forward_variable-rmodel.reactions.get_by_id(rxns[3]+cell).forward_variable,lb=0,ub=0,name='NAD_balcm'+cell)
        rmodel.add_cons_vars([NAD_balcm_constr])
        rmodel.solver.update()       

    return rmodel