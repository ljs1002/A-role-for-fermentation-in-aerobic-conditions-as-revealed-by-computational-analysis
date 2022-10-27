from build_roots_model import *
import pandas as pd
from cobra import Reaction, Metabolite
from cobra.flux_analysis.parsimonious import pfba
from cobra.flux_analysis.variability import *
from cobra import io
# from datetime import datetime, date, time
def mainEBC2(fixphloem=1):
    # Load model
    rootsdir = '/Users/user/Documents/Ox/roots/'
    orootsdir = '/Users/user/Documents/Ox/roots/roots-model/'

    # modelpwy='final_roots_model.xml'
    # modelpwy='tissue_dnut_model.xml'
    modelpwy='roots_model-1-1-1.xml'
    rmodel=io.read_sbml_model(rootsdir+modelpwy)
    
    # Remove ammonium as an osmolyte
    for rxn in rmodel.reactions:
        if 'AMMONIUM' in rxn.id and 'biomass' in rxn.id:
            # print(rxn.id, rxn.reaction)
            rmodel.remove_reactions({rxn})
        # Stop Fe transfer through apoplast
        if 'FeII' in rxn.id and 'apoplast_exchange' in rxn.id:
            rmodel.remove_reactions({rxn})
        # Prevent PRO import
        if 'PRO_ph' in rxn.id:
            rxn.upper_bound=0
    # Prevent some symp transport
    no_symp_transport = ['ATP','UTP','GTP','PROTON','UDP_GLUCOSE','CARBON_DIOXIDE','GDP_D_GLUCOSE','CPD_14553','OXYGEN_','AMP','NADH','FeII','SUPER_OXIDE','PPI']
    for rxn in rmodel.reactions:
        if 'Transfer' in rxn.id and any([x in rxn.id for x in no_symp_transport]):
            rxn.upper_bound=0
            rxn.lower_bound=0
            # print(rxn.id)
    if fixphloem:
        # Remove phloem metabolites not listed and add ratio constraint to remaining phloem mets (30% sucrose to 70% aas)
        ph_ratios = pd.read_csv(rootsdir+'Phloem_composition_z_maize.csv')
        rmodel.remove_reactions([rxn for rxn in rmodel.reactions if 'EX_X_' in rxn.id and '_ph_' in rxn.id and not any([x in rxn.id for x in list(ph_ratios['met_id'])+['PROTON']])])
        phloem_ratio_met=Metabolite('phloem_ratio_met')
        rxn=rmodel.reactions.EX_X_SUCROSE_ph_per00
        rxn.add_metabolites({phloem_ratio_met:-0.3})
        for rxn in rmodel.reactions:
            if 'EX_X_' in rxn.id and '_ph_' in rxn.id and not any([x in rxn.id for x in ['EX_X_SUCROSE_ph_per00','PROTON']]):
                # print(rxn.id)
                rxn.add_metabolites({phloem_ratio_met:0.7})
    for rxn in rmodel.reactions:
        if 'L_LACTATEDEHYDROG_RXN_c' in rxn.id:
            rxn.upper_bound=0
            rxn.lower_bound=0
    rxn=Reaction('Dummy_lact_dehydrog_cor00')
    rmodel.add_reaction(rxn)
    lrxn=rmodel.reactions.L_LACTATEDEHYDROG_RXN_c_cor00
    for met in lrxn.metabolites:
        if 'PROTON' not in met.id:
            rxn.add_metabolites({met.id:lrxn.get_coefficient(met.id)})
    rxn.upper_bound=1000
    rxn.lower_bound=-1000
    # Define obj
    rootspFBAobj(rmodel)
    # rmodel.objective = 'Phloem_import_reaction'
    biomass_rxns = [rxn.id for rxn in rmodel.reactions if 'Biomass_expanding_cell' in rxn.id]
    fva_sol=flux_variability_analysis(rmodel)
    data={'reaction_id':[x.id for x in rmodel.reactions],
        'cell':['_'.join(x.id.split('_')[-2:]) if ('_l' in x.id or '_d' in x.id) else '' for x in rmodel.reactions],
        'subsystem':[x.notes['SUBSYSTEM']  if 'SUBSYSTEM' in str(x.notes) else '' for x in rmodel.reactions],
        'reaction':[x.reaction for x in rmodel.reactions],
        'fva_min':[fva_sol.minimum[x.id] for x in rmodel.reactions],
        'fva_max':[fva_sol.maximum[x.id] for x in rmodel.reactions]}
    

    # Change biomass proportions
    test_rats = [1,0.8,0.5]
    # for ii in test_rats:
    #     for jj in test_rats:
    #         for kk in test_rats:
    #             for rxn in biomass_rxns:
    #                 new_model=rmodel.copy()
    #                 changeBiomassComposition(new_model.reactions.get_by_id(rxn),[ii,jj,kk])
    #                 sol=pfba(new_model)
    #                 data['-'.join([str(x) for x in [ii,jj,kk]])]=[sol[x.id] for x in rmodel.reactions]

    new_model=rmodel.copy()
    for ii in range(3):
        for jj in test_rats:
            these_rats = [1,1,1]
            these_rats[ii]=jj
            for rxn in biomass_rxns:
                changeBiomassComposition(new_model.reactions.get_by_id(rxn),these_rats)
            sol=pfba(new_model)
            data['-'.join([str(x) for x in these_rats])]=[sol[x.id] for x in rmodel.reactions]
            for rxn in biomass_rxns:
                resetBiomassComposition(new_model.reactions.get_by_id(rxn),these_rats)
    df=pd.DataFrame(data)
    # now = datetime.now().strftime('%Y_%m_%d')
    df.to_csv('/Users/user/Documents/Ox/roots/Spreadsheets/BiomassComposition2'+fixphloem+'.csv')    

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
