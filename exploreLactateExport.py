from build_roots_model import *
import pandas as pd
from cobra import Reaction, Metabolite
from cobra.flux_analysis.parsimonious import pfba
from cobra.flux_analysis.variability import *
# from datetime import datetime, date, time
def mainLE():
    # Load model
    rootsdir = '/Users/user/Documents/Ox/roots/'
    orootsdir = '/Users/user/Documents/Ox/roots/roots-model/'

    # modelpwy='final_roots_model.xml'
    # modelpwy='tissue_dnut_model.xml'
    modelpwy='roots_model-1-1-1.xml'
    from cobra import io
    rmodel=io.read_sbml_model(rootsdir+modelpwy)
    

    # Remove phloem metabolites not listed and add ratio constraint to remaining phloem mets
    ph_ratios = pd.read_csv(rootsdir+'Phloem_composition_z_maize.csv')
    rmodel.remove_reactions([rxn for rxn in rmodel.reactions if 'EX_X_' in rxn.id and '_ph_' in rxn.id and not any([x in rxn.id for x in list(ph_ratios['met_id'])+['PROTON']])])
    phloem_ratio_met=Metabolite('phloem_ratio_met')
    rxn=rmodel.reactions.EX_X_SUCROSE_ph_per00
    rxn.add_metabolites({phloem_ratio_met:-0.3})
    for rxn in rmodel.reactions:
        if 'EX_X_' in rxn.id and '_ph_' in rxn.id and not any([x in rxn.id for x in ['EX_X_SUCROSE_ph_per00','PROTON']]):
            # print(rxn.id)
            rxn.add_metabolites({phloem_ratio_met:0.7})

    # # Add reaction that in phloem ratios to maximise
    # # Restrict phloem import
    # # Replace all ph_per reactions with phloem import reactions
    # ph_ratios = pd.read_csv(rootsdir+'Phloem_composition_z_maize.csv')
    # ph_mets = list(ph_ratios['met_id'])
    # ph_nums = list(ph_ratios['Oshima1990'])
    # ph_imp_rxn = Reaction('Phloem_import_reaction')
    # rmodel.add_reaction(ph_imp_rxn)
    # phsuffix = '_ph_per00'
    # persuffix = '_c_per00'
    # for ii in range(1,len(ph_nums)):
    #     if not any([ph_mets[ii]+phsuffix in met.id for met in rmodel.metabolites]):
    #         newmet = Metabolite(ph_mets[ii]+phsuffix,compartment='ph_per00',name=ph_mets[ii],formula=rmodel.metabolites.get_by_id(ph_mets[ii]+persuffix).formula)
    #         metinrxn = Reaction('EX_X_'+ph_mets[ii]+phsuffix)
    #         metinrxn.add_metabolites({newmet:1})
    #         rmodel.add_reaction(metinrxn)
    #         # print(adda)
    #     # rmodel.remove_reactions([rxn for rxn in phrxns if ph_mets[ii] in rxn])
    #     ph_imp_rxn.add_metabolites({ph_mets[ii]+phsuffix:-ph_nums[ii],ph_mets[ii]+persuffix:ph_nums[ii]})
    #     # sol=rmodel.optimize()
    #     # print(sol['Phloem_import_reaction'],sol['CO2_tx_epi00'],[rxn for rxn in phrxns if ph_mets[ii] in rxn])

    # Define obj
    rootspFBAobj(rmodel)
    # rmodel.objective = 'Biomass_expanding_cell_epi00'
    # atpsink = Reaction('ATP_sink_c_per00')
    # rmodel.add_reaction(atpsink)
    # atpsink.add_metabolites({'ATP_c_per00':-0.65,'aATP_c_per00':-0.35})
    fva_sol=flux_variability_analysis(rmodel)
    data={'reaction_id':[x.id for x in rmodel.reactions],
        'cell':['_'.join(x.id.split('_')[-2:]) if ('_l' in x.id or '_d' in x.id) else '' for x in rmodel.reactions],
        'subsystem':[x.notes['SUBSYSTEM']  if 'SUBSYSTEM' in str(x.notes) else '' for x in rmodel.reactions],
        'reaction':[x.reaction for x in rmodel.reactions],
        'fva_min':[fva_sol.minimum[x.id] for x in rmodel.reactions],
        'fva_max':[fva_sol.maximum[x.id] for x in rmodel.reactions]}
    

    # Change atpsink proportions
    lbs = [0,0.01,0.1,1,10]
    for x in lbs:
        rmodel.reactions.ATPase_tx_per00.lower_bound=x
        # sol2=flux_variability_analysis(rmodel)
        # data[str(x)+'_min']=sol2.minimum
        # data[str(x)+'_max']=sol2.maximum
        sol2=rmodel.optimize()
        data[str(x)+'_opt']=[sol2[x.id] for x in rmodel.reactions]
    df=pd.DataFrame(data)
    # now = datetime.now().strftime('%Y_%m_%d')
    df.to_csv('/Users/user/Documents/Ox/roots/Spreadsheets/ATPsink.csv')    

def changeBiomassComposition(rxn,ratios):
    cellulose_met = [met for met in rxn.metabolites if 'Cellulose' in met.id][0]
    lipid_met = [met for met in rxn.metabolites if 'L_PHOSPHATIDATE' in met.id][0]
    protein_met = [met for met in rxn.metabolites if 'Protein' in met.id][0]
    rxn.add_metabolites({cellulose_met:rxn.get_coefficient(cellulose_met)*(ratios[0]-1),
    lipid_met:rxn.get_coefficient(lipid_met)*(ratios[1]-1),
    protein_met:rxn.get_coefficient(protein_met)*(ratios[2]-1)})

def resetBiomassComposition(rxn,ratios):
    cellulose_met = [met for met in rxn.metabolites if 'Cellulose' in met.id][0]
    lipid_met = [met for met in rxn.metabolites if 'L_PHOSPHATIDATE' in met.id][0]
    protein_met = [met for met in rxn.metabolites if 'Protein' in met.id][0]
    rxn.add_metabolites({cellulose_met:rxn.get_coefficient(cellulose_met)*(1-ratios[0]),
    lipid_met:rxn.get_coefficient(lipid_met)*(1-ratios[1]),
    protein_met:rxn.get_coefficient(protein_met)*(1-ratios[2])})
