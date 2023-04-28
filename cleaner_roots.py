# from exploreBiomassComposition import *
from cobra import io,Metabolite,Reaction
from build_roots_model import *
import pandas as pd
import pickle as pkl
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm
import matplotlib.colors as mcolors
import numpy as np
from matplotlib.cm import ScalarMappable

def rebuild_roots(removefermentationbool=None):
    rootsdir = '/Users/user/Documents/Ox/roots/'
    modelpwy='roots_model-1-1-1.xml'
    rmodel=io.read_sbml_model(rootsdir+modelpwy)

    # Remove ammonium as an osmolyte
    for rxn in rmodel.reactions:
        if 'AMMONIUM' in rxn.id and any([x in rxn.id for x in ['biomass']]):
            # print(rxn.id, rxn.reaction)
            rmodel.remove_reactions({rxn})
    # Stop Fe transfer through apoplast
        if 'FeII' in rxn.id and 'apoplast_exchange' in rxn.id:
            rmodel.remove_reactions({rxn})
    # Remove Proline from phloem import
        if 'PRO_ph' in rxn.id:
            rxn.upper_bound=0
        if 'unlProtHYPO' in rxn.id:
            rxn.upper_bound=0
            rxn.lower_bound=0
    # Constrain direction of phosphate transport
        if 'Pi_c' in rxn.id and 'Transfer' in rxn.id:
            rxn.lower_bound=0
        if 'AMMONIUM_c' in rxn.id and 'Transfer' in rxn.id:
            rxn.lower_bound=0
    # Constrain minimum xylem export
    rxn = rmodel.reactions.NITRATE_c_xy_per00
    rxn.lower_bound=rmodel.reactions.Nitrate_ec_epi00.lower_bound*.19
    for met in rxn.metabolites:
        if 'PROTON' in met.id:
            rxn.add_metabolites({met.id:-rxn.get_coefficient(met.id)})
    # cells = [rxn.id.split('_')[-1] for rxn in rmodel.reactions if 'CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c' in rxn.id]
    # rmodel.reactions.H_tx_epi00.lower_bound=-1000 # Neg means drawing in protons from the soil
    # Phloem sucrose to AA ratio
    ph_ratios = pd.read_csv(rootsdir+'Phloem_composition_z_maize.csv')
    rmodel.remove_reactions([rxn for rxn in rmodel.reactions if 'EX_X_' in rxn.id and '_ph_' in rxn.id and not any([x in rxn.id for x in list(ph_ratios['met_id'])+['PROTON']])])
    phloem_ratio_met=Metabolite('phloem_ratio_met',compartment='c_per00')
    rxn=rmodel.reactions.EX_X_SUCROSE_ph_per00
    rxn.add_metabolites({phloem_ratio_met:-0.3})
    for rxn in rmodel.reactions:
        if 'EX_X_' in rxn.id and '_ph_' in rxn.id and not any([x in rxn.id for x in ['EX_X_SUCROSE_ph_per00','PROTON']]):
            # print(rxn.id)
            rxn.add_metabolites({phloem_ratio_met:0.7})
    no_symp_transport = ['ATP','UTP','GTP','PROTON','UDP_GLUCOSE','CARBON_DIOXIDE','GDP_D_GLUCOSE','CPD_14553','OXYGEN_','AMP','NADH','FeII','SUPER_OXIDE','PPI']+['Heteroglycans','CELLULOSE','Cellulose','STARCH','ACETALD','ARACHIDIC_ACID_c_','DOCOSANOATE_c_',
                'CPD_16709_c_','STEARIC_ACID_c_','OLEATE_CPD_c_','Octadecadienoate_c_','LINOLENIC_ACID_c_',
                'PALMITATE','CPD_9245_c_','CPD_17412_c_','CPD_17291_c_','CPD_674_c_','Fatty_Acids_c_','CPD_2961_']
    # no_symp_transport = ['ATP','UTP','GTP','PROTON','UDP_GLUCOSE','CARBON_DIOXIDE','GDP_D_GLUCOSE','MAL_','OXALACETIC_ACID','NADH_','NAD_','GAP','DPG','FRUCTOSE_16_DIPHOSPHATE_','DIHYDROXY_ACETONE_PHOSPHATE_']
    for rxn in rmodel.reactions:
        if 'Transfer' in rxn.id and any([x in rxn.id for x in no_symp_transport]):
            rxn.upper_bound=0
            rxn.lower_bound=0
            # print(rxn.id)

    # exudates = ['FRU_','GLC_','CIT_','SUC_','FUM_','MAL_','OXALACETIC_ACID_','L_ALPHA_ALANINE_','HIS_','GLY_','L_ASPARTATE_','LEU_','GLT_','TYR_','PHE_','GLN_','VAL_','THR_','SER_','ILE_','ASN_','LYS_','MANNOSE_','ETOH_']
    exudates = ['ETOH_']
    for metfrag in exudates:
        if metfrag+'c_epi00' in rmodel.metabolites:
            metc = rmodel.metabolites.get_by_id(metfrag+'c_epi00')
            if metfrag+'e_epi00' not in rmodel.metabolites:
                mete = Metabolite(metfrag+'e_epi00',compartment='e_epi00')
                rmodel.add_metabolites({mete})
            else:
                mete = rmodel.metabolites.get_by_id(metfrag+'e_epi00')
            hc   = rmodel.metabolites.get_by_id('PROTON_c_epi00')
            he   = rmodel.metabolites.get_by_id('PROTON_e_epi00')
            rxnec = Reaction(metfrag+'ec_epi00')
            rxntx = Reaction(metfrag+'tx_epi00')
            rmodel.add_reactions({rxnec,rxntx})
            if metfrag == 'ETOH_':
                rxnec.add_metabolites({metc:-1,mete:1})
            else:
                rxnec.add_metabolites({metc:-1,mete:1,hc:-1,he:1})
            rxntx.add_metabolites({mete:-1})
        else:
            print('False ',metfrag)

    # Recalculate nitrate/phosphate etc
    # Mean length of lateral roots = 2.2 cm = 22 000 um
    # Final length of RootSlice 300 um
    # Moreno-Ortega et al., 2017 says cells elongate from 25 um to 150 um
    # Ours elongate from 25 um to 100 um
    # So final RootSlice length should be 450 um
    # 450 / 22000 * 1666 = ~34 um / hr
    # 225 um growth -> 6.6 hrs
    # Nitrate uptake.. 14.66 pmol/cm/s York 2016
    #add epidermal cell specific uptake constraints
    sa_pergDW = sa_calculator("roots-model/model_constraints_data/cell_dimensions.xlsx")

    # At length 25um (time=0) SA = 393 um^2, at length 100um (time = 4hrs) SA = 1571um^2
    # Over length of time, we can take av. = (393+1571)/2 = 982 um^2 = 9.82e-6 cm^2
    rootslice = pd.read_excel("roots-model/model_constraints_data/cell_dimensions.xlsx", index_col=0)
    SA_per_model = 9.82e-6
    d_SA=rootslice.loc['_epi00','change_surface_area']
    f_SA=rootslice.loc['_epi00','final_surface_area']
    SA_per_model = (f_SA-d_SA/2)*1e-8 # Av. SA cm^2/model
    print('Model SA: ',SA_per_model)
    # print('New nitrate uptake: ',14.66*1e-6*60*60*6.6) # umol/cm^2/6.6hrs
    # print('Old nitrate uptake: ',14.66*1e-6*60*60*3.05) # umol/cm^2/6.6hrs
    new_nitrate_uptake = 14.66*1e-6*60*60*6.6*1e-11*sa_pergDW
    print('New nitrate uptake: ',new_nitrate_uptake) # mmol/gDW/6.6hrs
    print('Old nitrate uptake: ',14.66*1e-6*60*60*3.05*1e-11*sa_pergDW) # mmol/gDW/3.05hrs


    # Phosphate uptake 1e-6 umol/cm/s
    # print(1e-6*24*60*60) #umol/cm/day
    # Change to Stefan's 0.055 μmol/cm2/day
    new_Pi_uptake = 0.055/24*6.6*1e-11*sa_pergDW
    print('New phosphate uptake: ',new_Pi_uptake) # mmol/gDW/6.6hrs
    rmodel.reactions.Nitrate_ec_epi00.lower_bound = -1000
    rmodel.reactions.Nitrate_ec_epi00.upper_bound = float(new_nitrate_uptake)
    rmodel.reactions.Nitrate_ec_epi00.lower_bound = float(new_nitrate_uptake)

    rmodel.reactions.Pi_ec_epi00.lower_bound = 0
    rmodel.reactions.Pi_ec_epi00.upper_bound = 1000
    # rmodel.reactions.Pi_ec_epi00.lower_bound = float(new_Pi_uptake)
    # Constrain minimum xylem export
    rxn = rmodel.reactions.NITRATE_c_xy_per00
    rxn.lower_bound=rmodel.reactions.Nitrate_ec_epi00.lower_bound*.19

    # Add biomass ratio
    biomass_rxn = rmodel.reactions.soluble_biomass_tissue
    nit_accum_rxns = [rxn for rxn in rmodel.reactions if 'NITRATE_b_biomass' in rxn.id]
    met = Metabolite('Nitrate_b_biomass_ratio',compartment='b_epi00')
    biomass_rxn.add_metabolites({met:-0.2})
    for rxn in nit_accum_rxns:
        rxn.add_metabolites({met:1})
    if removefermentationbool:
        rmodel=remove_fermentation(rmodel)
    return rmodel

def remove_fermentation(rmodel,newmodel=None):
    if newmodel:
        this_model=quietcopy(rmodel)
    else:
        this_model=rmodel
    ferm_rxns=['L_LACTATEDEHYDROG_RXN_c','RXN_6161_c','ALANINE_AMINOTRANSFERASE_RXN_c','ALCOHOL_DEHYDROG_RXN']
    for rxn in this_model.reactions:
        if any([x in rxn.id for x in ferm_rxns]):
            rxn.upper_bound=0
            rxn.lower_bound=0
    return this_model


def just_solve(model):
    this_model=quietcopy(model)
    sol1=rootspFBAobj(this_model)
    sol=this_model.optimize()
    return sol
def solve_and_pickle(rmodel,label='nitrate_accum_ratio_'):
    this_model = quietcopy(rmodel)
    sol1=rootspFBAobj(this_model)
    sol=rmodel.optimize()
    fva_sol=flux_variability_analysis(rmodel)
    date = datetime.now().strftime('%Y_%m_%d')
    label += date
    with open('roots_'+label+'.pkl','wb') as file:
        pkl.dump(sol,file)
    with open('roots_fva_'+label+'.pkl','wb') as file:
        pkl.dump(fva_sol,file)
    return this_model,sol,fva_sol, label

def unpickle(file_label):
    with open('roots_'+file_label+'.pkl','rb') as file:
        sol=pkl.load(file)
    with open('roots_fva_'+file_label+'.pkl','rb') as file:
        fva_sol=pkl.load(file)
    return sol, fva_sol

def plot_gradients(rmodel,sol,fva_sol,savebool,model_time = 6.6,log_on=1,minlog=1e-3,cb_on=None):
    
    cells = [rxn.id.split('_')[-1] for rxn in rmodel.reactions if 'CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c' in rxn.id]
    cell_labels = ['EPI1','COR1','COR2','COR3','COR4','COR5','COR6','COR7','COR8','END1','PER1']
    
    glyc_rxns = ['SUCROSE_SYNTHASE_RXN_c_cor03','3_PERIOD_2_PERIOD_1_PERIOD_48_RXN_v_cor02','CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c_epi00','FRUCTOKINASE_RXN_c_epi00','2_PERIOD_7_PERIOD_1_PERIOD_90_RXN_c_cor06','F16ALDOLASE_RXN_c_per00','3PGAREARR_RXN_p_cor00','2PGADEHYDRAT_RXN_p_cor06','PEPDEPHOS_RXN_c_cor00','PYRUVDEH_RXN_m_epi00','L_LACTATEDEHYDROG_RXN_c_end00','GLYC3PDEHYDROGBIOSYN_RXN_c_epi00','ACETYL_COA_CARBOXYLTRANSFER_RXN_p_epi00','ALCOHOL_DEHYDROG_RXN_c_epi00','GLUC1PURIDYLTRANS_RXN_c_epi00','2_PERIOD_3_PERIOD_1_PERIOD_180_RXN_c_epi00']
    glyc_rxns =['_'.join(x.split('_')[:-2]) for x in glyc_rxns]
    glyc_rxn_dict = {'SUCROSE_SYNTHASE_RXN':'Sucrose Synthase','3_PERIOD_2_PERIOD_1_PERIOD_48_RXN':'Invertase',
    'CELLULOSE_SYNTHASE_UDP_FORMING_RXN':'Cellulose Synthase','FRUCTOKINASE_RXN':'Fructose Kinase','2_PERIOD_7_PERIOD_1_PERIOD_90_RXN':'F6P Phosphotransferase','F16ALDOLASE_RXN':'FBP Aldolase','3PGAREARR_RXN':'Phosphoglycerate mutase','2PGADEHYDRAT_RXN':'Enolase','PEPDEPHOS_RXN':'Pyruvate Kinase','PYRUVDEH_RXN':'Pyruvate Dehydrogenase','L_LACTATEDEHYDROG_RXN':'Lactate Dehydrogenase','GLYC3PDEHYDROGBIOSYN_RXN':'Glycerol-3P Dehydrogenase','ACETYL_COA_CARBOXYLTRANSFER_RXN':'Acetyl-CoA Carboxyltransfer','2_PERIOD_3_PERIOD_1_PERIOD_180_RXN':'Lipid Biosynthesis','CITSYN_RXN':'Citrate Synthase','ALCOHOL_DEHYDROG_RXN':'Alcohol Dehydrogenase','FRUCTOKINASE_RXN':'Fructose Kinase','GLUC1PURIDYLTRANS_RXN':'UGPase'}

    TCA_rxns = ['MALATE_DEH_RXN','CITSYN_RXN','ACONITATEDEHYDR_RXN','ACONITATEHYDR_RXN','ISOCITDEH_RXN','ISOCITRATE_DEHYDROGENASE_NAD_RXN','2OXOGLUTARATEDEH_RXN','SUCCINYL_COA_HYDROLASE_RXN','SUCCCOASYN_RXN','FUMHYDR_RXN','SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN','Mitochondrial_ATP_Synthase','GLUTAMATE_DEHYDROGENASE_RXN']+['GABATRANSAM_RXN','SUCCINATE_SEMIALDEHYDE_DEHYDROGENASE_RXN']
    tca_rxn_dict = {'MALATE_DEH_RXN':'Malate Dehydrogenase','CITSYN_RXN':'Citrate Synthase','ACONITATEDEHYDR_RXN':'Aconitate Dehydrogenase','ACONITATEHYDR_RXN':'Aconitase','ISOCITDEH_RXN':'Isocitrate Dehydrogenase','ISOCITRATE_DEHYDROGENASE_NAD_RXN':'Isocitrate Dehydrogenase NAD','2OXOGLUTARATEDEH_RXN':'2-Ketoglutarate dehydrogenase','SUCCINYL_COA_HYDROLASE_RXN':'Succinyl-CoA Hydrolase','SUCCCOASYN_RXN':'Succinyl-CoA Synthetase','FUMHYDR_RXN':'Fumarase','SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN':'Succinate Dehydrogenase','Mitochondrial_ATP_Synthase':'Mitochondrial ATP Synthase','GLUTAMATE_DEHYDROGENASE_RXN':'Glutamate Dehydrogenase','GABATRANSAM_RXN':'GABA Transaminase','SUCCINATE_SEMIALDEHYDE_DEHYDROGENASE_RXN':'Succinate-Semialdehyde Dehydrogenase'}

    these_rxns = glyc_rxns+TCA_rxns
    this_rxn_dict = glyc_rxn_dict | tca_rxn_dict

    glyc_dict={}
    for cell in cells:
        glyc_dict[cell]={}
        for rxn in these_rxns:
            glyc_dict[cell][rxn]=sum([sol[x.id] for x in rmodel.reactions if cell in x.id and rxn in x.id])/model_time
    # temp=pd.DataFrame(glyc_dict)
    # temp1=[x for y in temp.values for x in y]

    glyc_dict={}
    for cell in cells:
        glyc_dict[cell]={}
        for rxn in these_rxns:
            glyc_dict[cell][rxn]=round(sum([fva_sol.maximum[x.id] for x in rmodel.reactions if cell in x.id and rxn in x.id])/model_time,4)
    tempmax=pd.DataFrame(glyc_dict)
    temp1max=[x for y in tempmax.values for x in y]
    # tmax=max(temp1max)
    # tmin=min(temp1max)
    glyc_dict={}
    for cell in cells:
        glyc_dict[cell]={}
        for rxn in these_rxns:
            glyc_dict[cell][rxn]=round(sum([fva_sol.minimum[x.id] for x in rmodel.reactions if cell in x.id and rxn in x.id])/model_time,4)
    tempmin=pd.DataFrame(glyc_dict)
    # temp1min=[x for y in tempmin.values for x in y]


    colors1 = plt.cm.afmhot_r(np.linspace(0., 1, 256))
    colors2 = plt.cm.Blues_r(np.linspace(0, 1, 256))
    # combine them and build a new colormap
    colors = np.vstack((colors2, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    #new
    tmax=np.ceil(max(temp1max))
    cmo = mymap#plt.get_cmap('afmhot').reversed()

    if log_on:
        sm=ScalarMappable(cmap=cmo, norm=LogNorm(minlog,tmax))
    else:
        sm = ScalarMappable(cmap=cmo, norm=plt.Normalize(0,tmax))
        sm.set_array([])
    if not cb_on:
        plt.figure(figsize=(10,1))
        plt.bar(range(1),[0 for x in range(1)],width=1,edgecolor='none')
        plt.axis('off')
        # if log_on:
            # cbar=plt.colorbar(ScalarMappable(cmap=cmo, norm=LogNorm(minlog,tmax)),orientation='horizontal',pad=0.3,aspect=35)
            # cbar=plt.colorbar(ScalarMappable(cmap=cmo, norm=SymLogNorm(minlog,vmin=-tmax,vmax=tmax)),orientation='horizontal',pad=0.3,aspect=35)
        # else:
        #     cbar=plt.colorbar(sm,orientation='horizontal',pad=0.3,aspect=35)
        if savebool:
            plt.savefig('Figures/barcode_colorbar.eps',bbox_inches='tight')
    for ii in tempmin.index:
        # colorsmax = [[1-x*y/255/tmax for y in [208, 66, 255]] for x in tempmax.loc[ii]]
        # colorsmin = [[1-x*y/255/tmax for y in [208, 66, 255]] for x in temp1min.loc[ii]]

        cm = sm.get_cmap()
        if log_on:
            colorsmax = [cm(0.5+np.log10(max(abs(x),minlog)/minlog)/np.log10(tmax/minlog)/2) if x>0 else cm(0.5) if abs(x)<minlog else cm(0.5-np.log10(max(abs(x),minlog)/minlog)/np.log10(tmax/minlog)/2) for x in tempmax.loc[ii]]
            colorsmin = [cm(0.5+np.log10(max(abs(x),minlog)/minlog)/np.log10(tmax/minlog)/2) if x>0 else cm(0.5) if abs(x)<minlog else cm(0.5-np.log10(max(abs(x),minlog)/minlog)/np.log10(tmax/minlog)/2) for x in tempmin.loc[ii]]
        else:
            colorsmax = [cm(x/tmax) for x in tempmax.loc[ii]]
            colorsmin = [cm(x/tmax) for x in tempmin.loc[ii]]

        # print(this_rxn_dict[ii],[x for x in tempmax.loc[ii]])
        # print([x for x in temp1min.loc[ii]])
        # plt.figure(figsize=(10,1))

        plt.figure(figsize=(1,1))
        plt.bar(range(1),[0.5 for x in range(1)],color=colorsmin[0],width=1,edgecolor='none')
        plt.bar(range(1),[0.5 for x in range(1)],bottom = [0.5 for x in range(1)],color=colorsmax[0],width=1,edgecolor='none')
        plt.bar(range(1),[0.02 for x in range(1)],color='black',width=1)
        plt.bar(range(1),[0.02 for x in range(1)],bottom = [0.49 for x in range(1)],color='black',width=1)
        plt.bar(range(1),[0.02 for x in range(1)],bottom = [0.99 for x in range(1)],color='black',width=1)
        plt.bar(np.arange(1)-0.5,[1 for x in range(1)],color='none',width=0.001,edgecolor='black')
        plt.bar(np.arange(1)+0.5,[1 for x in range(1)],color='none',width=0.001,edgecolor='black')
        plt.axis('off')
        for i in range(1):
            plt.text(i-0.4,-0.3,cell_labels[i],fontsize=16)
        plt.title(this_rxn_dict[ii],fontsize=24)
        if savebool:
            plt.savefig('Figures/barcode_epi_'+savebool+this_rxn_dict[ii]+'.png',bbox_inches='tight')
        if cb_on:
            plt.figure(figsize=(10,1.8))
        else:
            plt.figure(figsize=(10,1))
        plt.bar(range(10),[0.5 for x in range(10)],color=colorsmin[1:],width=1,edgecolor='none')
        plt.bar(range(10),[0.5 for x in range(10)],bottom = [0.5 for x in range(10)],color=colorsmax[1:],width=1,edgecolor='none')
        plt.bar(range(10),[0.02 for x in range(10)],color='black',width=1)
        plt.bar(range(10),[0.02 for x in range(10)],bottom = [0.49 for x in range(10)],color='black',width=1)
        plt.bar(range(10),[0.02 for x in range(10)],bottom = [0.99 for x in range(10)],color='black',width=1)
        plt.bar(np.arange(10)-0.5,[1 for x in range(10)],color='none',width=0.0001,edgecolor='black')
        plt.bar(np.arange(10)+0.5,[1 for x in range(10)],color='none',width=0.0001,edgecolor='black')
        plt.axis('off')
        plt.title(this_rxn_dict[ii],fontsize=24)
        if cb_on:
            for i in range(1,11):
                plt.text((i-1)-0.25,-0.35,cell_labels[i])
            if log_on:
                plt.colorbar(ScalarMappable(cmap=cmo, norm=LogNorm(minlog,tmax)))
            else:
                sm = ScalarMappable(cmap=cm, norm=plt.Normalize(0,tmax))
                sm.set_array([])
                cbar=plt.colorbar(sm,orientation='horizontal',pad=0.3,aspect=35)
        else:
            for i in range(1,11):
                plt.text((i-1)-0.4,-0.3,cell_labels[i],fontsize=16)
        if savebool:
            plt.savefig('Figures/barcode_'+savebool+this_rxn_dict[ii]+'.png',bbox_inches='tight')
        # plt.axhline(0,color='black')
    return tmax

def generate_glyc_table(rmodel,sol,solf,label='latest'):
    glyc_rxns = ['SUCROSE_SYNTHASE_RXN_c_cor03','3_PERIOD_2_PERIOD_1_PERIOD_48_RXN_v_cor02','CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c_epi00','FRUCTOKINASE_RXN_c_epi00','2_PERIOD_7_PERIOD_1_PERIOD_90_RXN_c_cor06','F16ALDOLASE_RXN_c_per00','3PGAREARR_RXN_p_cor00','2PGADEHYDRAT_RXN_p_cor06','PEPDEPHOS_RXN_c_cor00','PYRUVDEH_RXN_m_epi00','L_LACTATEDEHYDROG_RXN_c_end00','GLYC3PDEHYDROGBIOSYN_RXN_c_epi00','ACETYL_COA_CARBOXYLTRANSFER_RXN_p_epi00','2_PERIOD_3_PERIOD_1_PERIOD_180_RXN_p_epi00','CITSYN_RXN_m_epi00']
    glyc_rxns =['_'.join(x.split('_')[:-2]) for x in glyc_rxns]
    TCA_rxns = ['MALATE_DEH_RXN','CITSYN_RXN','ACONITATEDEHYDR_RXN','ACONITATEHYDR_RXN','ISOCITDEH_RXN','ISOCITRATE_DEHYDROGENASE_NAD_RXN','2OXOGLUTARATEDEH_RXN','SUCCINYL_COA_HYDROLASE_RXN','SUCCCOASYN_RXN','FUMHYDR_RXN','SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN','Mitochondrial_ATP_Synthase','GLUTAMATE_DEHYDROGENASE_RXN']+['GABATRANSAM_RXN','SUCCINATE_SEMIALDEHYDE_DEHYDROGENASE_RXN']
    cells = [rxn.id.split('_')[-1] for rxn in rmodel.reactions if 'CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c' in rxn.id]

    mets = ['SUCROSE','FRU_','GLC_','UDP_GLUCOSE','FRUCTOSE_6P','FRUCTOSE_16_DIPHOSPHATE','GAP_','G3P_','3PGA_','DIHYDROXY_ACETONE_PHOSPHATE','PHOSPHO_ENOL_PYRUVATE','PYRUVATE_','ACETYL_COA','CIT','2_KETOGL','MAL_','OXALOACE','FUM_','SUC_']
    model_time=6.6
    table_dict={}
    for rxnroot in glyc_rxns+['CELLULOSE_SYNTHASE_GDP_FORMING_RXN','TRIOSEPISOMERIZATION_RXN','GLUC1PURIDYLTRANS']+TCA_rxns+list(set(['_'.join(x.id.split('_')[:-1]) for x in rmodel.reactions if '_biomass_' in x.id and any([y in x.id for y in mets])])):#TCA_rxns:#
        rxns = [x for x in rmodel.reactions if rxnroot in x.id and any([y in x.id for y in cells[1:]])]
        epi_rxns = [x for x in rmodel.reactions if rxnroot in x.id and 'epi00' in x.id]
        table_dict[rxnroot]=[epi_rxns[0].reaction,sum([sol[x.id] for x in epi_rxns])/model_time,sum([sol[x.id] for x in rxns])/10/model_time,sum([solf[x.id] for x in epi_rxns])/model_time,sum([solf[x.id] for x in rxns])/10/model_time]
    transfers = ['_'.join(x.id.split('_')[-5:]) for x in rmodel.reactions if 'MAL_c_' in x.id and 'Transfer' in x.id]
    transfers.pop(7)
    transfers=['epi00_Transfer_epi00_to_cor00']+transfers
    transfer_rxns = list(set(['_'.join(x.id.split('_')[:-5])+'_' for x in rmodel.reactions if '_Transfer_' in x.id and '_r_' not in x.id]))
    for rxnroot in transfer_rxns:
        rsum=sum([sol[rxnroot+transfers[ii]]-sol[rxnroot+transfers[ii+1]] for ii in range(len(transfers)-1)])/10
        rsumf=sum([solf[rxnroot+transfers[ii]]-solf[rxnroot+transfers[ii+1]] for ii in range(len(transfers)-1)])/10
        if any([-sol[rxnroot+transfers[0]],rsum,-solf[rxnroot+transfers[0]],rsumf]):
            table_dict[rxnroot] = ['',-sol[rxnroot+transfers[0]]/model_time,rsum/model_time,-solf[rxnroot+transfers[0]]/model_time,rsumf/model_time]
    for key in list(table_dict.keys()):
        table_dict[key]+=[30*x for x in table_dict[key][1:]]
    flux_table = pd.DataFrame.from_dict(table_dict,orient='index')
    flux_table.columns = ['Stoichiometry','EPI','Other','No Ferm EPI','No Ferm Other','30x EPI','30x Other','30x No Ferm EPI','30x No Ferm Other']
    flux_table.head(2)
    flux_table.to_excel('Spreadsheets/Fluxes_for_glyc_diag_'+label+'.xlsx')

def plot_symp_transport(rmodel,thissol,savebool='test',partitions = 6,partition_size=6,ferm=None):
    col = [ '#1f78b4', '#b2df8a', '#8ae1ed', '#bebada', '#fdb462', '#d68aed', '#f5737a', '#fee08b', '#e69595', '#66c2a5', '#d690b2', '#3288bd', '#ffff99', '#a6cee3', '#ffffb3', '#6a3d9a', '#ff7f00', '#d9d9d9', '#9e0142', '#b3de69', '#33a02c', '#f46d43', '#abdda4', '#d53e4f', '#fccde5', '#fb8072', '#ccebc5', '#ffffbf', '#e6f598', '#bc80bd', '#5e4fa2', '#e89f82', '#80b1d3', '#ffed6f']
    met_dict={'PHOSPHO_ENOL_PYRUVATE':'PEP', 'UDP_GLUCURONATE':'UDP-Glucuronate', 'GLC_1_P':'Glucose-1P', 'FRU':'Fructose', 'G3P':'3-PGA', 'CPD_510':'Glucuronate-1P', 'L_LACTATE':'Lactate', 'FRUCTOSE_6P':'Fructose-6P', 'NITRATE':r'NO$_3$', 'FRUCTOSE_16_DIPHOSPHATE':'Fructose-1,6P', 'PYRUVATE':'Pyruvate', 'AMMONIUM':r'NH$_4$', 'SUCROSE':'Sucrose', 'GAP':'GAP', '2_PG':'2-PG', 'ETOH':'Ethanol', 'MGII':r'Mg$^{2+}$', 'DPG':'DPGA', 'Pi':'Pi','SER':'Serine','GLC_6_P':'Glucose-6P','DIHYDROXY_ACETONE_PHOSPHATE':'DHAP','UDP':'UDP','MAL':'Malate','ARG':'Arginine','OXALACETIC_ACID':'OAA','GLT':'GLT','L_ASPARTATE':'Aspartate','MET':'Methionine','TYR':'Tyrosine','PHOSPHATIDYL_CHOLINE':'Phosphatidyl Choline','GLYCEROL_3P_c':'Glycerol-3P'}
    transfers = set(['_'.join(x.id.split('_')[-5:]) for x in rmodel.reactions if 'Transfer' in x.id and not 'WATER' in x.id])
    large_symp = sorted([rxn for rxn in rmodel.reactions if any([x in rxn.id for x in list(transfers)])],key=lambda x:abs(thissol[x.id]),reverse=True)
    subset_symp = set(['_'.join(x.id.split('_')[:-5]) for x in large_symp[:150]])

    transfers_list=sorted(list(transfers),key=lambda x:x)
    transfers_list =[transfers_list[-1]]+transfers_list[:-1]
    # transfers_list 

    large_symp = sorted([rxn for rxn in rmodel.reactions if any([x in rxn.id for x in list(transfers)]) and not rxn.id[0]=='a' and not 'WATER' in rxn.id],key=lambda x:abs(thissol[x.id]),reverse=True)
    subset_symp = set(['_'.join(x.id.split('_')[:-5]) for x in large_symp[:]])
    subs_order = {}
    for met in subset_symp:
        subs_order[met]=max([abs(thissol[x.id]) for x in large_symp if met in x.id])
    
    sym_partitions={}
    temp_symp=sorted(list(subset_symp),key = lambda x: subs_order[x],reverse=True)
    for ii in range(partitions):
        if not ferm:
            if 'FRU_c' in temp_symp:
                ref_met = 'FRU_c'
            elif 'NITRATE_c' in temp_symp:
                ref_met = 'NITRATE_c'
            elif 'PYRUVATE_c' in temp_symp:
                ref_met='PYRUVATE_c'
            elif 'L_LACTATE_c' in temp_symp:
                ref_met = 'L_LACTATE_c'
            elif 'MAL_c' in temp_symp:
                ref_met = 'MAL_c'
            elif 'CPD_510_c' in temp_symp:
                ref_met = 'CPD_510_c'
            # elif 'PYRUVATE_c' in temp_symp:
            #     ref_met = 'PYRUVATE_c'
            else:
                ref_met = sorted(temp_symp,key=lambda x: formula_bkdown(rmodel.metabolites.get_by_id(x+'_cor01').formula)[0])[-1]
            if ref_met=='FRU_c':
                sym_partitions[ii] = [x for x in temp_symp if x[:-2] in list(met_dict.keys()) and 'ose' in met_dict[x[:-2]]]
                sym_partitions[ii]=sym_partitions[ii][:min(len(sym_partitions[ii]),partition_size)]
            # elif ref_met=='NITRATE_c':
            #     sym_partitions[ii] = [x for x in ['NITRATE_c','AMMONIUM_c','MGII_c','CAII_c','HCO3_c'] if x in temp_symp]+[x+'_c' for x in aas if x+'_c' in temp_symp]
            # elif ref_met=='DIHYDROXY_ACETONE_PHOSPHATE_c':
            #     sym_partitions[ii] = [x for x in ['GAP_c','G3P_c','DIHYDROXY_ACETONE_PHOSPHATE_c','MAL_c','OXALACETIC_ACID_c','DPG_c'] if x in temp_symp]
            # elif ref_met=='PYRUVATE_c':
            #     sym_partitions[ii] = ['PYRUVATE_c','L_LACTATE_c','ETOH_c','PHOSPHO_ENOL_PYRUVATE_c']
            # For nitrate accum ratio:
            elif ref_met=='NITRATE_c':
                sym_partitions[ii] = [x for x in ['NITRATE_c','AMMONIUM_c','HCO3_c','GLT_c'] if x in temp_symp]#+[x+'_c' for x in aas if x+'_c' in temp_symp]
            elif ref_met=='MAL_c':
                sym_partitions[ii] = [x for x in ['MAL_c','OXALACETIC_ACID_c'] if x in subset_symp]
            elif ref_met=='PYRUVATE_c':
                sym_partitions[ii] = [x for x in ['G3P_c','DPG_c','2_PG_c','PHOSPHO_ENOL_PYRUVATE_c','PYRUVATE_c'] if x in subset_symp]
            # elif ref_met=='PYRUVATE_c':
            #     sym_partitions[ii] = ['PYRUVATE_c','PHOSPHO_ENOL_PYRUVATE_c','2_PG_c','G3P_c','DIHYDROXY_ACETONE_PHOSPHATE_c','DPG_c']
            elif ref_met=='CPD_510_c':
                sym_partitions[ii] = ['CPD_510_c','UDP_GLUCURONATE_c']
            elif ref_met == 'L_LACTATE_c':
                sym_partitions[ii] = ['L_LACTATE_c','ETOH_c']
            elif ref_met== 'GLYCEROL_3P_c':
                sym_partitions[ii] = ['GLYCEROL_3P_c','DIHYDROXY_ACETONE_PHOSPHATE_c','GAP_c','FRUCTOSE_16_DIPHOSPHATE_c']
            # elif len(temp_symp)>partition_size:
            #     sym_partitions[ii] = sorted(temp_symp,key=lambda x: dls(rmodel.metabolites.get_by_id(ref_met+'_cor01'),rmodel.metabolites.get_by_id(x+'_cor01')))[:partition_size]
            else:
                sym_partitions[ii]=temp_symp[:partition_size]
        else:
            if 'FRU_c' in temp_symp:
                ref_met = 'FRU_c'
            elif 'NITRATE_c' in temp_symp:
                ref_met = 'NITRATE_c'
            elif 'PYRUVATE_c' in temp_symp:
                ref_met='PYRUVATE_c'            
            elif 'MAL_c' in temp_symp:
                ref_met = 'MAL_c'
            # elif 'GLYCEROL_3P_c' in temp_symp:
            #     ref_met='GLYCEROL_3P_c'
            elif 'CIT_c' in temp_symp:
                ref_met = 'CIT_c'
            elif 'CPD_510_c' in temp_symp:
                ref_met = 'CPD_510_c'
            # elif 'PYRUVATE_c' in temp_symp:
            #     ref_met = 'PYRUVATE_c'
            else:
                ref_met = sorted(temp_symp,key=lambda x: formula_bkdown(rmodel.metabolites.get_by_id(x+'_cor01').formula)[0])[-1]
            if ref_met=='FRU_c':
                sym_partitions[ii] = [x for x in temp_symp if x[:-2] in list(met_dict.keys()) and 'ose' in met_dict[x[:-2]]]
                sym_partitions[ii]=sym_partitions[ii][:min(len(sym_partitions[ii]),partition_size)]
            elif ref_met=='NITRATE_c':
                sym_partitions[ii] = [x for x in ['NITRATE_c','AMMONIUM_c','GLT_c','SER_c'] if x in temp_symp]#+[x+'_c' for x in aas if x+'_c' in temp_symp]
            elif ref_met=='MAL_c':
                sym_partitions[ii] = [x for x in ['MAL_c','OXALACETIC_ACID_c','FUM_c'] if x in subset_symp]
            elif ref_met=='PYRUVATE_c':
                sym_partitions[ii] = [x for x in ['G3P_c','DPG_c','2_PG_c','PHOSPHO_ENOL_PYRUVATE_c','PYRUVATE_c'] if x in subset_symp]
            elif ref_met=='CPD_510_c':
                sym_partitions[ii] = ['CPD_510_c','UDP_GLUCURONATE_c']
            elif ref_met == 'L_LACTATE_c':
                sym_partitions[ii] = ['L_LACTATE_c','ETOH_c']
            elif ref_met== 'GLYCEROL_3P_c':
                sym_partitions[ii] = ['GLYCEROL_3P_c','DIHYDROXY_ACETONE_PHOSPHATE_c','GAP_c','FRUCTOSE_16_DIPHOSPHATE_c']
            elif ref_met=='CIT_c':
                sym_partitions[ii] = ['CIT_c','CIS_ACONITATE_c','THREO_DS_ISO_CITRATE_c']
            # elif len(temp_symp)>partition_size:
            #     sym_partitions[ii] = sorted(temp_symp,key=lambda x: dls(rmodel.metabolites.get_by_id(ref_met+'_cor01'),rmodel.metabolites.get_by_id(x+'_cor01')))[:partition_size]
            else:
                sym_partitions[ii]=temp_symp[:partition_size]
        
        sym_partitions[ii] = sorted(list(set(sym_partitions[ii])),key = lambda x: subs_order[x],reverse=True)

        if ii==partitions-1:
            plot_sym_transport_roots(thissol,sym_partitions[ii],rmodel,col,'Figures/SymplasticTransport'+str(ii)+savebool+'.pdf')
        else:
            plot_sym_transport_roots(thissol,sym_partitions[ii],rmodel,col,'Figures/SymplasticTransport'+str(ii)+savebool+'.pdf',xlab=0,xtick=0)
        # print(sorted(temp_symp,key=lambda x: dls(rmodel.metabolites.get_by_id(ref_met+'_cor01'),rmodel.metabolites.get_by_id(x+'_cor01'))))
        temp_symp = sorted(list(set(temp_symp).difference(set(sym_partitions[ii]))),key = lambda x: subs_order[x],reverse=True)


import re
import numpy as np
def formula_bkdown(chem_formula):
    elems=['C','H','O','P','N']
    formula_dict = [0,0,0,0,0]
    if chem_formula:
        temp1 = re.findall(r'[A-Z][a-z]?',chem_formula)
        for it,elt in enumerate(elems):
            if elt in temp1:
                it1=temp1.index(elt)
                if chem_formula.index(elt)+1==len(chem_formula):
                    formula_dict[it]=1
                elif it1+1==len(temp1):
                    formula_dict[it]=float(chem_formula[(chem_formula.index(elt)+1):])
                elif chem_formula.index(elt)==(chem_formula.index(temp1[it1+1])-1):
                    formula_dict[it]=1
                else:
                    formula_dict[it]=float(chem_formula[(chem_formula.index(elt)+1):(chem_formula.index(temp1[it1+1]))])
    return formula_dict

import matplotlib.pyplot as plt
import numpy as np
def plot_bdry(rmodel,sol,it=1,tot_num=1,rxns=[],ylim=None):
    exp_mets = {}
    exbit='EX_X_'
    for rxn in [x for x in rmodel.reactions if (x.boundary and abs(sol[x.id])>1e-4 and not '_b_' in x.id and not 'TotalSoluteConstraint' in x.id and not 'Phloem_carbon' in x.id) or x.id in rxns]:
            rxns+=[rxn.id]
            rxnroot = '_'.join(rxn.id.split('_')[:-1])
            if exbit in rxnroot:
                rxnroot = '_'.join(rxnroot.split('_')[2:])
            if rxnroot not in list(exp_mets.keys()):
                exp_mets[rxnroot]=0
            exp_mets[rxnroot]+=sol[rxn.id]
            # print(round(sol[rxn.id],4),rxn.id ,rxn.reaction)
    width = 0.9/(tot_num)
    plt.bar(np.arange(len(exp_mets))+it*width,exp_mets.values(),width=width)
    t1=plt.xticks(range(len(exp_mets)),exp_mets.keys(),rotation=90)
    if ylim:
        plt.ylim(0,ylim)
    return rxns

import matplotlib.pyplot as plt

def plot_sym_transport_roots(sol,subset_symp,rmodel,col,saveas=None,ylab=1,xlab=1,xtick=1,model_time=6.6):
    if subset_symp:
        transfers = set(['_'.join(x.id.split('_')[-5:]) for x in rmodel.reactions if 'Transfer' in x.id])
        transfers_list=list(sorted(list(transfers),key=lambda x:x))
        transfers_list =[transfers_list[-1]]+transfers_list[:-1]
        x_ticks = {'cor03_Transfer_03_to_04':'COR 4', 'cor06_Transfer_06_to_07':'COR 7', 'epi00_Transfer_epi00_to_cor00': 'EPI 1', 'cor01_Transfer_01_to_02': 'COR 2', 'cor05_Transfer_05_to_06':'COR 6', 'cor00_Transfer_00_to_01':'COR 1', 'end00_Transfer_end00_to_per00':'END 1', 'cor07_Transfer_cor07_to_end00':'COR 8', 'cor02_Transfer_02_to_03':'COR 3', 'cor04_Transfer_04_to_05':'COR 5'}
        x_ticks2 = {'cor03_Transfer_03_to_04':'COR 5', 'cor06_Transfer_06_to_07':'COR 8', 'epi00_Transfer_epi00_to_cor00': 'COR 1', 'cor01_Transfer_01_to_02': 'COR 3', 'cor05_Transfer_05_to_06':'COR 7', 'cor00_Transfer_00_to_01':'COR 2', 'end00_Transfer_end00_to_per00':'PER 1', 'cor07_Transfer_cor07_to_end00':'END 1', 'cor02_Transfer_02_to_03':'COR 4', 'cor04_Transfer_04_to_05':'COR 6'}
        met_dict={'PHOSPHO_ENOL_PYRUVATE':'PEP', 'UDP_GLUCURONATE':'UDP-Glucuronate', 'GLC_1_P':'Glucose-1P', 'FRU':'Fructose', 'G3P':'3-PGA', 'CPD_510':'Glucuronate-1P', 'L_LACTATE':'Lactate', 'FRUCTOSE_6P':'Fructose-6P', 'NITRATE':r'NO$_3$', 'FRUCTOSE_16_DIPHOSPHATE':'Fructose-1,6P', 'PYRUVATE':'Pyruvate', 'AMMONIUM':r'NH$_4$', 'SUCROSE':'Sucrose', 'GAP':'GAP', '2_PG':'2-PG', 'ETOH':'Ethanol', 'MGII':r'Mg$^{2+}$', 'DPG':'DPGA', 'Pi':'Pi','SER':'Serine','GLC_6_P':'Glucose-6P','DIHYDROXY_ACETONE_PHOSPHATE':'DHAP','UDP':'UDP','MAL':'Malate','ARG':'Arginine','OXALACETIC_ACID':'OAA','GLT':'GLT','L_ASPARTATE':'Aspartate','MET':'Methionine','TYR':'Tyrosine','PHOSPHATIDYL_CHOLINE':'Phosphatidyl Choline'}
        symdict = {}
        plts=[]
        plted=[]
        met_labels = [met_dict[x[:-2]] if x[:-2] in (met_dict.keys()) else x for x in subset_symp]
        # plt.figure(figsize=(6,8))
        fig,(ax1,ax2) = plt.subplots(1,2, sharex=False,figsize=(12,5),sharey=True)
        for it,met in enumerate(list(subset_symp)):
            symdict[met]=[sol[met+'_'+x] for x in transfers_list]
            mets = met[:-2]
            if mets=='NITRATE':
                mets='Nitrate'
            txrxn = mets+'_tx_epi00'
            # xyrxn='EX_X_'+mets+'_xy_per00'
            xyrxn=met+'_xy_per00'
            phrxn = mets+'_ph_c_per00'
            # print(xyrxn)
            symdict[met] = [sol[txrxn]*rmodel.reactions.get_by_id(txrxn).get_coefficient(met[:-2]+'_e_epi00') if txrxn in rmodel.reactions else 0]+symdict[met]+[sum([sol[xyrxn]*rmodel.reactions.get_by_id(xyrxn).get_coefficient(met[:-2]+'_xy_per00') if xyrxn in rmodel.reactions else 0]+[sol[phrxn]*rmodel.reactions.get_by_id(phrxn).get_coefficient(mets+'_ph_per00') if phrxn in rmodel.reactions else 0])]
            if 1:#sum(symdict[met])>=0:
                ax1.stackplot(range(len(symdict[met])),[x/model_time if x>=0 else 0 for x in symdict[met]],labels = transfers_list, colors = [col[it]], baseline='zero',alpha=0.5)
                plts+=ax2.stackplot(range(len(symdict[met])),[-x/model_time if x<0 else 0 for x in symdict[met][-1::-1]],labels = transfers_list, colors = [col[it]], baseline='zero',alpha=0.5) 
                plted+=[it]
        ax2.legend(plts,[list(met_labels)[x] for x in plted],bbox_to_anchor=(1,0.5),loc='center left',fontsize=18)
        if xtick:
            ax1.set_xticks(range(len(symdict[met])),labels=['Soil']+[x_ticks[x] for x in transfers_list]+['PER 1'],rotation=90,fontsize=16)
            ax2.set_xticks(range(len(symdict[met])),labels=['Xy/Ph']+[x_ticks2[x] for x in transfers_list[-1::-1]]+['EPI 1'],rotation=90,fontsize=16)
        else:
            ax1.set_xticks(range(len(symdict[met])),labels=['' for x in symdict[met]],rotation=90,fontsize=16)
            ax2.set_xticks(range(len(symdict[met])),labels=['' for x in symdict[met]],rotation=90,fontsize=16)
        ax1.tick_params(axis='y',labelsize=16)
        ax2.tick_params(axis='y',labelsize=16)
        # ax1.set_ylim(0,sol['SUCROSE_ph_c_per00']*1.1)
        # ax1.set_ylim(0,sol['Nitrate_ec_epi00']*1.1/model_time)
        if xlab:
            ax1.set_xlabel('Transferred from',fontsize=20)
            ax2.set_xlabel('Transferred from',fontsize=20)
        if ylab:
            ax1.set_ylabel('Flux to next cell (mmol/h)',fontsize=20)
        if saveas:
            plt.savefig(saveas, bbox_inches='tight')
        plt.show()

def dls(met1,met2):
    f1 = formula_bkdown(met1.formula)
    f2 = formula_bkdown(met2.formula)
    if f1[1]!=0:
        return np.sqrt(sum([pow(f1[x]/f1[1]-f2[x]/max(f2[1],1),2) for x in range(len(f1))]))
    else:
        print('Odd dls for: ',met1,met2)
        return 1

def print_all_met_reactions(rmodel, sol,metroot = 'GLYCEROL_3P_c',sign='=',model_time=6.6):
    rxns=[]
    rxn_sols=[]
    for met in [x.id for x in rmodel.metabolites if metroot in x.id]:# and 'epi00' in x.id]:
        for rxn in rmodel.metabolites.get_by_id(met).reactions:
            if sign=='=':
                if sol[rxn.id]*rxn.get_coefficient(met)!=0 and not 'Transfer' in rxn.id:
                    rxns +=[rxn]
                    rxn_sols +=[sol[rxn.id]*rxn.get_coefficient(met)]
            elif sign=='<':
                if sol[rxn.id]*rxn.get_coefficient(met)<0 and not 'Transfer' in rxn.id:
                    rxns +=[rxn]
                    rxn_sols +=[sol[rxn.id]*rxn.get_coefficient(met)]
            elif sign=='>':
                if sol[rxn.id]*rxn.get_coefficient(met)>0 and not 'Transfer' in rxn.id:
                    rxns +=[rxn]
                    rxn_sols +=[sol[rxn.id]*rxn.get_coefficient(met)]
            else:
                rxns +=[rxn]
                rxn_sols +=[sol[rxn.id]*rxn.get_coefficient(met)]

    for rxn in sorted(rxns,key=lambda x:abs(rxn_sols[rxns.index(x)]),reverse=True):
        print(round(rxn_sols[rxns.index(rxn)],4)/model_time,rxn.id,rxn.reaction)

import io as IO
from contextlib import redirect_stdout
def quietcopy(obj):
    text_trap=IO.StringIO()
    with redirect_stdout(text_trap):
        newcopy=obj.copy()
    return newcopy

################ BUDGETS ###############

import matplotlib.pyplot as plt

def generate_process_chart_all(model, solution1,split_epi=None,savebool='',thresh=0,NADbool=1,verbose=None):
    reaction_atp_process = {
    'Mitochondrial_ATP_Synthase_m_': 'Mitochondrial ATP Synthase',
    'PHOSGLYPHOS_RXN_c_': 'Cytosolic glycolysis' ,
    'PEPDEPHOS_RXN_p_': 'Plastidic glycolysis',
    'PHOSPHORIBULOSEKINASE_RXN_p_': 'Calvin cycle',
    'PEPDEPHOS_RXN_c_': 'Cytosolic glycolysis',
    'PHOSGLYPHOS_RXN_p_': 'Plastidic glycolysis',
    'ACETYL_COA_CARBOXYLTRANSFER_RXN_p_': 'Lipid synthesis',
    'FRUCTOKINASE_RXN_c_': 'Cytosolic glycolysis',
    'PROTON_ATPase_c_': 'Plasma membrane proton ATPase',
    'PROTONATP_rev_vc':'Tonoplast proton ATPase',
    'PEPCARBOXYKIN_RXN_c_': 'PEP carboxykinase',
    'GLUC1PURIDYLTRANS_RXN_c_': 'Cellulose Synthesis',
    'SUCROSE_SYNTHASE_RXN_c_': 'Cellulose Synthesis',
    '2_PERIOD_7_PERIOD_7_PERIOD_44_RXN_c_':'UDP-glucuronate-glucuronate-1P shuttle',
    '2_PERIOD_7_PERIOD_7_PERIOD_34_RXN_c_': 'Cellulose Synthesis',
    'Protein_Polymerisation_c_': 'Protein polymerisation',
    '2_PERIOD_7_PERIOD_7_PERIOD_13_RXN_c_': 'Cytosolic glycolysis',
    'Protein_Translocation_c': 'Protein polymerisation',
    'Protein_Processing_c': 'Protein polymerisation',
    'ADENYL_KIN': 'Adenylate kinase',
    'TRNA_LIGASE_RXN':'Protein polymerisation',
    'GLURS_RXN':'Protein polymerisation',
    'SUCCCOASYN_RXN': 'TCA cycle',
    '6PFRUCTPHOS_RXN_c': 'Cytosolic glycolysis',
    '6PFRUCTPHOS_RXN_p': 'Plastidic glycolysis',
    'CARBAMATE_KINASE_RXN_p':'Amino acid degradation',
    'GLUTKIN_RXN_c':'Amino acid degradation',
    'PRPPSYN_RXN_p': 'Amino acid degradation',
    'SHIKIMATE_KINASE_RXN_p':'Amino acid degradation',
    'GLUC1PADENYLTRANS_RXN_p':'Starch synthesis',
    'ATPase_tx': 'Cell maintenance'
    # 'GLUTKIN_RXN':
    
    }
    col = ['#1f78b4','#b2df8a','#8ae1ed','#bebada','#fdb462','#d68aed','#f5737a','#fee08b','#e69595','#66c2a5','#d690b2','#3288bd','#ffff99','#a6cee3','#ffffb3','#6a3d9a','#ff7f00','#d9d9d9','#9e0142','#b3de69','#33a02c','#f46d43','#abdda4','#d53e4f','#fccde5','#fb8072','#ccebc5','#ffffbf','#e6f598','#bc80bd','#5e4fa2','#e89f82','#80b1d3','#ffed6f']

    atp_coloursdict = {
    'Mitochondrial ATP Synthase': col[0],
    'Cytosolic glycolysis': col[1] ,
    'Plastidic glycolysis': col[10],
    'Calvin cycle': col[3],
    'Lipid synthesis': col[4],
    'Cellulose Synthesis': col[5],
    'Amino acid degradation': col[24],
    'Protein polymerisation': col[6],
    'Starch synthesis': col[25],
    'Plasma membrane proton ATPase': col[18],
    'PEP carboxykinase': col[8],
    'Tonoplast proton ATPase': col[9],
    'TCA cycle': col[2],
    'UDP-glucuronate-glucuronate-1P shuttle': col[31],
    'Other': '#000000',
    'Adenylate kinase':col[23],
    'Cell maintenance':'#1A85FF'}

    nad_coloursdict = {
    'NADH dehydrogenase': col[0],
    'Glycolysis': col[1] ,
    'TCA cycle': col[2],
    'Pentose phosphate pathway': col[3],
    'Lipid synthesis': col[4],
    # 'Cellulose Synthesis': '#26C6DA',
    'Amino acid degradation': col[24],
    'Fe reduction/oxidation': col[6],
    'Fermentation': col[7],
    'Tonoplast proton ATPase': col[8],
    'Malate dehydrogenase (x)':col[13],
    'Malate dehydrogenase (c)':col[14],
    'Malate dehydrogenase (p)':col[9],
    'Malate dehydrogenase (m)':col[15],
    'NADP- malic reaction (c)':col[30],
    'NADP- malic reaction (p)':col[28],
    'NAD- malic reaction (m)':col[29],
    'Valine biosynthesis':col[10],
    'Isocitrate dehydrogenase':col[11],
    'Cell maintenance':'#1A85FF',
    'Oxidative phosphorylation':col[21],
    'GABA bypass':col[26],
    'Other': '#000000',}

    model_time = 6.6
    outdict={}
    cell_names = ["Epidermis", "Cortical 1", "Cortical 2", "Cortical 3", "Cortical 4",
                    "Cortical 5", "Cortical 6", "Cortical 7", "Cortical 8", "Endodermal", "Pericycle"]
    celltags = ['_epi00','_cor00','_cor01','_cor02','_cor03','_cor04','_cor05','_cor06','_cor07','_end00','_per00']
    # atp_rxns_fluxes={}
    if NADbool:
        fig, ((ax1_epi,ax1),(temp1,temp2),(ax2_epi,ax2),(temp5,temp6),(ax1n_epi,ax1n_oth),(temp3,temp4),(ax2n_epi,ax2n_oth)) = plt.subplots(nrows=7, ncols=2,height_ratios=[6,1,6,1,6,1,6],width_ratios=[1.1,10],figsize=(3,8),dpi=800)#(12.8,9.6))
        temp1.set_visible(False)
        temp2.set_visible(False)
        temp5.set_visible(False)
        temp6.set_visible(False)
    else:
        fig, ((ax1_epi,temp3,ax1),(temp1,temp4,temp2),(ax2_epi,temp5,ax2)) = plt.subplots(nrows=3, ncols=3,height_ratios=[6,1,6],width_ratios=[1.1,1,10],figsize=(3,4),dpi=800)#(12.8,9.6))
        temp1.set_visible(False)
        temp2.set_visible(False)
        temp3.set_visible(False)
        temp4.set_visible(False)
        temp5.set_visible(False)
    poskeykeep=[]
    negkeykeep=[]
    p1={}
    p2={}
    if split_epi:
        start=1
    else:
        start=0
    for it,celltag in enumerate(celltags[start:]):
        outdict[celltag]={}
        temp, atp_rxns_fluxes= generateATPbudgetforonecell_ed(model,solution1.fluxes,celltag, percentage=False,thresh=0)
        process_rxn_list = list(reaction_atp_process.keys())
        posdict = {}
        negdict={}
        for rxn in atp_rxns_fluxes.keys():
            matching_proc_keys=[x for x in process_rxn_list if x in rxn]
            if matching_proc_keys:
                if len(matching_proc_keys)>1 and verbose:
                    print(matching_proc_keys)
                rxn_process = reaction_atp_process[matching_proc_keys[0]]
            else:
                rxn_process = 'Other'
                if atp_rxns_fluxes[rxn] and verbose:
                    print('Oao: ',round(atp_rxns_fluxes[rxn],4),rxn)
            if atp_rxns_fluxes[rxn]>0:
                if rxn_process not in list(posdict.keys()):
                    posdict[rxn_process]=0
                posdict[rxn_process]+=atp_rxns_fluxes[rxn]
            elif atp_rxns_fluxes[rxn]<0:
                if rxn_process not in list(negdict.keys()):
                    negdict[rxn_process]=0
                negdict[rxn_process]+=-atp_rxns_fluxes[rxn]
        pos_keys = list(posdict.keys())
        neg_keys = list(negdict.keys())
        pos_values = list(posdict.values())
        neg_values = list(negdict.values())
        plt.style.reload_library()
        my_colors = ["r","g","b","k","y","m","c"]  #red, green, blue, black, etc.
        possum=0
        for itp,ii in enumerate(pos_keys):
            if ii not in poskeykeep and pos_values[itp]>1e-5:
                poskeykeep+=[ii]
        for i in range(len(poskeykeep)):
            if poskeykeep[i] in pos_keys:
                posval=pos_values[pos_keys.index(poskeykeep[i])]/model_time
                if posval>0:
                    p1[i]=ax1.bar(cell_names[it+start], posval, bottom=possum,color=atp_coloursdict[poskeykeep[i]])
                    if celltag == '_cor01a' and verbose:
                        print(poskeykeep[i],posval)
                    possum+=posval
        for itn,ii in enumerate(neg_keys):
            if ii not in negkeykeep and neg_values[itn]>1e-5:
                negkeykeep+=[ii]
        negsum=0
        for i in range(0,len(negkeykeep)):
            if negkeykeep[i] in neg_keys:
                negval=neg_values[neg_keys.index(negkeykeep[i])]/model_time
                if negval>0:
                    p2[i]=ax2.bar(cell_names[it+start], negval, bottom=negsum,color=atp_coloursdict[negkeykeep[i]])
                    if celltag == '_cor01a'and verbose:
                        print(negkeykeep[i],negval)
                    negsum+=negval
        # ax2.tick_params("x", labelrotation=50)
        # ax1.tick_params("x", labelrotation=50)
        ax1.set_title("ATP generating processes",x=0.3)
        ax2.set_title("ATP consuming processes",x=0.3)
        ax1.set_xticks(range(len(cell_names[it+start])),labels=['' for x in cell_names[it+start]])
        if NADbool:
            ax2.set_xticks(range(len(cell_names[it+start])),labels=['' for x in cell_names[it+start]])
    for i in range(0,len(poskeykeep)):
        p1[i]=ax1.bar(cell_names[it-1+start], 0, bottom=possum,color=atp_coloursdict[poskeykeep[i]])
    for i in range(0,len(negkeykeep)):
        p2[i]=ax2.bar(cell_names[it-1+start], 0, bottom=negsum,color=atp_coloursdict[negkeykeep[i]])

    if not split_epi and savebool:
        plt.savefig('Figures/ATP_budget2_'+savebool+'.pdf',bbox_inches='tight')

    if split_epi:
        poskeykeep_epi=[]
        negkeykeep_epi=[]
        p1_epi={}
        p2_epi={}
        width = 0.7
        for it,celltag in enumerate([celltags[0]]):
            outdict[celltag]={}
            temp, atp_rxns_fluxes= generateATPbudgetforonecell_ed(model,solution1.fluxes,celltag, percentage=False,thresh=0)
            process_rxn_list = list(reaction_atp_process.keys())
            posdict = {}
            negdict={}
            for rxn in atp_rxns_fluxes.keys():
                matching_proc_keys=[x for x in process_rxn_list if x in rxn]
                if matching_proc_keys:
                    if len(matching_proc_keys)>1 and verbose:
                        print(matching_proc_keys)
                    rxn_process = reaction_atp_process[matching_proc_keys[0]]
                else:
                    rxn_process = 'Other'
                    if atp_rxns_fluxes[rxn] and verbose:
                        print('Oae: ',round(atp_rxns_fluxes[rxn],4),rxn)
                if atp_rxns_fluxes[rxn]>0:
                    if rxn_process not in list(posdict.keys()):
                        posdict[rxn_process]=0
                    posdict[rxn_process]+=atp_rxns_fluxes[rxn]
                elif atp_rxns_fluxes[rxn]<0:
                    if rxn_process not in list(negdict.keys()):
                        negdict[rxn_process]=0
                    negdict[rxn_process]+=-atp_rxns_fluxes[rxn]
            pos_keys = list(posdict.keys())
            neg_keys = list(negdict.keys())
            pos_values = list(posdict.values())
            neg_values = list(negdict.values())
            
            possum_epi=0
            for itp,ii in enumerate(pos_keys):
                if ii not in poskeykeep_epi and pos_values[itp]>1e-5:
                    poskeykeep_epi+=[ii]
            for i in range(len(poskeykeep_epi)):
                if poskeykeep_epi[i] in pos_keys:
                    posval=pos_values[pos_keys.index(poskeykeep_epi[i])]/model_time
                    if posval>0:
                        p1_epi[i]=ax1_epi.bar(cell_names[it], posval, bottom=possum_epi,color=atp_coloursdict[poskeykeep_epi[i]],width=width)
                        possum_epi+=posval
                        ax1_epi.set_xlim(-0.5,0.5)
            for itn,ii in enumerate(neg_keys):
                if ii not in negkeykeep_epi and neg_values[itn]>1e-5:
                    negkeykeep_epi+=[ii]
            negsum_epi=0
            for i in range(0,len(negkeykeep_epi)):
                if negkeykeep_epi[i] in neg_keys:
                    negval=neg_values[neg_keys.index(negkeykeep_epi[i])]/model_time
                    if negval>0:
                        p2_epi[i]=ax2_epi.bar(cell_names[it], negval, bottom=negsum_epi,color=atp_coloursdict[negkeykeep_epi[i]],width=width)
                        ax2_epi.set_xlim(-0.5,0.5)
                        negsum_epi+=negval
            # ax2_epi.tick_params("x", labelrotation=50)
            # ax1_epi.tick_params("x", labelrotation=50)
            ax1_epi.set_ylabel("ATP biosynthesis" '\n' "(mmol/h)")
            ax2_epi.set_ylabel("ATP usage" '\n' "(mmol/h)")
            ax2_epi.set_ylim(0,negsum_epi*1.05)
            ax1_epi.set_xticks([])
            if NADbool:
                ax2_epi.set_xticks([])
        
        for i in range(len(poskeykeep),len(poskeykeep)+len(poskeykeep_epi)):
            if poskeykeep_epi[i-len(poskeykeep)] not in poskeykeep:
                p1[i]=ax1.bar(cell_names[-2], 0, bottom=possum,color=atp_coloursdict[poskeykeep_epi[i-len(poskeykeep)]],width=width)
        # ax1.legend(p1.values(),[poskeykeep[x] for x in range(len(poskeykeep))]+[poskeykeep_epi[x] for x in range(len(poskeykeep_epi)) if poskeykeep_epi[x] not in (poskeykeep)],bbox_to_anchor=(1,1.32))
    
        for i in range(len(negkeykeep),len(negkeykeep)+len(negkeykeep_epi)):
            if negkeykeep_epi[i-len(negkeykeep)] not in negkeykeep:
                p2[i]=ax2.bar(cell_names[-2], 0, bottom=negsum,color=atp_coloursdict[negkeykeep_epi[i-len(negkeykeep)]],width=width)
        # ax2.legend(p2.values(),[negkeykeep[x] for x in range(len(negkeykeep))]+[negkeykeep_epi[x] for x in range(len(negkeykeep_epi)) if negkeykeep_epi[x] not in (negkeykeep)],bbox_to_anchor=(1,1.55))
        leg1=[poskeykeep[x] for x in range(len(poskeykeep))]+[poskeykeep_epi[x] for x in range(len(poskeykeep_epi)) if poskeykeep_epi[x] not in (poskeykeep)]
        leg2 = [negkeykeep[x] for x in range(len(negkeykeep))]+[negkeykeep_epi[x] for x in range(len(negkeykeep_epi)) if negkeykeep_epi[x] not in (negkeykeep)]
        leg_its = [x for x in range(len(leg2)) if leg2[x] not in leg1]
        if NADbool:
            ax2.legend(list(p1.values())+[list(p2.values())[x] for x in leg_its],
                        leg1+[leg2[x] for x in leg_its],bbox_to_anchor=(1,3))
        else:
            ax2.legend(list(p1.values())+[list(p2.values())[x] for x in leg_its],
                        leg1+[leg2[x] for x in leg_its],bbox_to_anchor=(1,2.5))
    model_time = 6.6
    if NADbool:
        reaction_nad_process = {'RXN_9532_p': 'Lipid synthesis', 
        'NADH_DEHYDROG_A_RXN_mi': 'Oxidative phosphorylation', 
        'RXN0_5330_NAD_mi': 'Oxidative phosphorylation', 
        'RXN_9663_p': 'Lipid synthesis', 
        'RXN_9661_p': 'Lipid synthesis',
        'RXN_9536_p': 'Lipid synthesis',
        'PYRROLINECARBDEHYDROG_RXN_NADP_m': 'Amino acid degradation', 
        'ISOCITDEH_RXN_m': 'TCA cycle', 
        'RXN_9658_p': 'Lipid synthesis',
        'RXN_9514_p': 'Lipid synthesis',
        'RXN_9660_p': 'Lipid synthesis',
        'PYRUVDEH_RXN_p': 'Glycolysis', 
        'PYRUVDEH_RXN_m': 'Glycolysis', 
        'MALATE_DEH_RXN_x': 'Malate dehydrogenase (x)',  
        'MALATE_DEH_RXN_c': 'Malate dehydrogenase (c)',  
        'MALATE_DEH_RXN_m': 'Malate dehydrogenase (m)',
        'MALATE_DEH_RXN_p': 'Malate dehydrogenase (p)',
        'RXN_9518_p': 'Lipid synthesis',
        'RXN_9662_p': 'Lipid synthesis',
        '2OXOGLUTARATEDEH_RXN_m': 'TCA cycle', 
        'RXN_9540_p': 'Lipid synthesis',
        'RXN_9524_p': 'Lipid synthesis',
        'RXN_9657_p': 'Lipid synthesis',
        'GLUTAMATE_DEHYDROGENASE_RXN_m': 'Amino acid degradation', 
        'RXN_9528_p': 'Lipid synthesis',
        '6PGLUCONDEHYDROG_RXN_p': 'Pentose phosphate pathway', 
        'RXN_9659_p': 'Lipid synthesis',
        'GLU6PDEHYDROG_RXN_p':'Pentose phosphate pathway', 
        'GAPOXNPHOSPHN_RXN_': 'Glycolysis', 
        'FERRIC_CHELATE_REDUCTASE_RXN_c': 'Fe reduction/oxidation',
        'ALCOHOL_DEHYDROG_RXN_c': 'Fermentation', 
        'L_LACTATEDEHYDROG_RXN_c':'Fermentation',
        'GLYC3PDEHYDROGBIOSYN_RXN_c': 'Lipid synthesis',
        'MALIC_NADP_RXN_c':'NADP- malic reaction (c)',
        'MALIC_NADP_RXN_p':'NADP- malic reaction (p)',
        '1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m':'NAD- malic reaction (m)',
        'ACETOLACTREDUCTOISOM_RXN_p':'Other',#'Valine biosynthesis',
        'NADPHox': 'Cell maintenance',
        'ISOCITRATE_DEHYDROGENASE_NAD_RXN':'Isocitrate dehydrogenase',
        'ISOCITDEH_RXN_c': 'TCA cycle',
        'PYRROLINECARBREDUCT_RXN_NAD_c':'Amino acid degradation',
        'GLUTSEMIALDEHYDROG_RXN_c':'Amino acid degradation',
        'SHIKIMATE_5_DEHYDROGENASE_RXN_p': 'Amino acid degradation',
        'SUCCINATE_SEMIALDEHYDE_DEHYDROGENASE_RXN_m':'GABA bypass',
        'GCVMULTI_RXN_m':'Amino acid degradation',
        # 'NITRATE_REDUCTASE_NADH_RXN':'',
        # '1_PERIOD_1_PERIOD_1_PERIOD_255_RXN':'',
        # 'VANILLICACID_BIOSYNTHESIS':'',
        # 'NADH_KINASE_RXN':'',
        # 'PYRROLINECARBREDUCT_RXN_NAD':'',
        # 'ALLYSINE_DEHYDROG_RXN':'',
        # '1_PERIOD_5_PERIOD_1_PERIOD_9_RXN':'',
        'Beta_Oxidation':'Lipid synthesis'
        }
        
        outdict={}
        cell_names = ["Epidermis", "Cortical 1", "Cortical 2", "Cortical 3", "Cortical 4",
                        "Cortical 5", "Cortical 6", "Cortical 7", "Cortical 8", "Endodermal", "Pericycle"]
        celltags = ['_epi00','_cor00','_cor01','_cor02','_cor03','_cor04','_cor05','_cor06','_cor07','_end00','_per00']
        # atp_rxns_fluxes={}

        temp3.set_visible(False)
        temp4.set_visible(False)
        poskeykeep=[]
        negkeykeep=[]
        p1={}
        p2={}
        maxpossum=0
        maxnegsum=0
        if split_epi:
            start=1
        else:
            start=0
        for it,celltag in enumerate(celltags[start:]):
            outdict[celltag]={}
            total_atp_dict, atp_rxns_fluxes= generateNADbudgetforonecell_ed(model,solution1.fluxes,celltag, percentage=False,thresh=thresh)
            process_rxn_list = list(reaction_nad_process.keys())
            posdict = {}
            negdict={}
            for rxn in atp_rxns_fluxes.keys():
                matching_proc_keys=[x for x in process_rxn_list if x in rxn]
                if matching_proc_keys:
                    if len(matching_proc_keys)>1 and verbose:
                        print(matching_proc_keys)
                    rxn_process = reaction_nad_process[matching_proc_keys[0]]
                else:
                    rxn_process = 'Other'
                    if round(atp_rxns_fluxes[rxn],4) and verbose:
                        print('Ono: ',rxn,atp_rxns_fluxes[rxn])
                if atp_rxns_fluxes[rxn]>0:
                    if rxn_process not in list(posdict.keys()):
                        posdict[rxn_process]=0
                    posdict[rxn_process]+=atp_rxns_fluxes[rxn]
                    if celltag == '_cor01' and verbose:
                            print(rxn_process, rxn, atp_rxns_fluxes[rxn])
                elif atp_rxns_fluxes[rxn]<0:
                    if rxn_process not in list(negdict.keys()):
                        negdict[rxn_process]=0
                    negdict[rxn_process]+=-atp_rxns_fluxes[rxn]
                    if celltag == '_cor01' and verbose:
                            print(rxn_process, rxn, atp_rxns_fluxes[rxn])
            pos_keys = list(posdict.keys())
            neg_keys = list(negdict.keys())
            pos_values = list(posdict.values())
            neg_values = list(negdict.values())
            plt.style.reload_library()
            possum=0
            for itp,ii in enumerate(pos_keys):
                if ii not in poskeykeep and pos_values[itp]>1e-5:
                    poskeykeep+=[ii]
            for i in range(len(poskeykeep)):
                if poskeykeep[i] in pos_keys:
                    posval = pos_values[pos_keys.index(poskeykeep[i])]
                    if posval>0:
                        p1[i]=ax1n_oth.bar(cell_names[it+start], posval/model_time, bottom=possum,color=nad_coloursdict[poskeykeep[i]])
                        # if celltag=='_cor04':
                        #     print(pos_keys[i],pos_val,nad_coloursdict[pos_keys[i]])
                        possum+=posval/model_time
            for itn,ii in enumerate(neg_keys):
                if ii not in negkeykeep and neg_values[itn]>1e-5:
                    negkeykeep+=[ii]
            negsum=0
            # print('Neg lengths: ',len(neg_keys),len(neg_processes))
            for i in range(0,len(negkeykeep)):
                if negkeykeep[i] in neg_keys:
                    negval=neg_values[neg_keys.index(negkeykeep[i])]
                    if negval>0:
                        p2[i]=ax2n_oth.bar(cell_names[it+start], negval/model_time, bottom=negsum,color=nad_coloursdict[negkeykeep[i]])
                        # if celltag=='_cor04':
                        #     print(negkeykeep[i],negval,nad_coloursdict[negkeykeep[i]])
                        negsum+=negval/model_time
            plt.setp(ax2n_oth.xaxis.get_majorticklabels(),rotation=50,ha='right',rotation_mode='anchor')
            # ax2n_oth.tick_params("x", labelrotation=50,ha='left',rotation_mode='anchor')
            # ax1n_oth.tick_params("x", labelrotation=50)
            ax1n_oth.set_title("NAD(P)H generating processes",x=0.3)
            # ax1n_oth.set_ylabel("NAD(P)H biosynthesis (mmol/h)")
            ax2n_oth.set_title("NAD(P)H consuming processes",x=0.3)
            # ax2n_oth.set_ylabel("NAD(P)H usage (mmol/h)")
            ax1n_oth.set_xticks(range(len(cell_names[it+start])),labels=['' for x in cell_names[it+start]])
            if possum>maxpossum:
                maxpossum=possum
            if negsum>maxnegsum:
                maxnegsum=negsum
        
        for i in range(len(poskeykeep)):
            p1[i]=ax1n_oth.bar(cell_names[it+start], 0, bottom=possum,color=nad_coloursdict[poskeykeep[i]],width=width)
        # ax1n_oth.legend(p1.values(),[poskeykeep[x] for x in list(p1.keys())])
        for i in range(0,len(negkeykeep)):
            p2[i]=ax2n_oth.bar(cell_names[it+start], 0, bottom=negsum,color=nad_coloursdict[negkeykeep[i]],width=width)
        # ax2n_oth.legend(p2.values(),[negkeykeep[x] for x in list(p2.keys())])
        ax1n_oth.set_ylim(0,maxpossum*1.05)
        ax2n_oth.set_ylim(0,maxnegsum*1.05)
        if split_epi:
            poskeykeep_epi=[]
            negkeykeep_epi=[]
            p1_epi={}
            p2_epi={}
            for it,celltag in enumerate([celltags[0]]):
                outdict[celltag]={}
                total_atp_dict, atp_rxns_fluxes= generateNADbudgetforonecell_ed(model,solution1.fluxes,celltag, percentage=False,thresh=thresh)
                # [print(t1,atp_rxns_fluxes[t1]) for t1 in list(atp_rxns_fluxes.keys()) if atp_rxns_fluxes[t1]]
                process_rxn_list = list(reaction_nad_process.keys())
                posdict = {}
                negdict={}
                for rxn in atp_rxns_fluxes.keys():
                    matching_proc_keys=[x for x in process_rxn_list if x in rxn]
                    if matching_proc_keys:
                        if len(matching_proc_keys)>1 and verbose:
                            print(matching_proc_keys)
                        rxn_process = reaction_nad_process[matching_proc_keys[0]]
                    else:
                        rxn_process = 'Other'
                        if round(atp_rxns_fluxes[rxn],4) and verbose:
                            print('One: ',rxn,atp_rxns_fluxes[rxn])
                    if atp_rxns_fluxes[rxn]>0:
                        if rxn_process not in list(posdict.keys()):
                            posdict[rxn_process]=0
                        posdict[rxn_process]+=atp_rxns_fluxes[rxn]
                        if celltag == '_cor01' and verbose:
                                print(rxn_process, rxn, atp_rxns_fluxes[rxn])
                    elif atp_rxns_fluxes[rxn]<0:
                        if rxn_process not in list(negdict.keys()):
                            negdict[rxn_process]=0
                        negdict[rxn_process]+=-atp_rxns_fluxes[rxn]
                        if celltag == '_cor01' and verbose:
                                print(rxn_process, rxn, atp_rxns_fluxes[rxn])
                pos_keys = list(posdict.keys())
                neg_keys = list(negdict.keys())
                pos_values = list(posdict.values())
                neg_values = list(negdict.values())
                plt.style.reload_library()
                possum=0
                for itp,ii in enumerate(pos_keys):
                    if ii not in poskeykeep_epi and pos_values[itp]>1e-5:
                        poskeykeep_epi+=[ii]
                for i in range(len(poskeykeep_epi)):
                    if poskeykeep_epi[i] in pos_keys:
                        posval = pos_values[pos_keys.index(poskeykeep_epi[i])]
                        if posval>0:
                            p1_epi[i]=ax1n_epi.bar(cell_names[it], posval/model_time, bottom=possum,color=nad_coloursdict[poskeykeep_epi[i]],width=width)
                            ax1n_epi.set_xlim(-0.5,0.5)
                            possum+=posval/model_time
                # ax1n_epi.tick_params("x", labelrotation=50)
                ax1n_epi.set_xticks([])
                # ax1n_epi.set_ylim(0, 4)
                # ax1n_epi.set_title("NAD(P)H generating processes")
                ax1n_epi.set_ylabel('NAD(P)H biosynthesis' '\n' '(mmol/gDW/h)')
                ax1n_epi.set_ylim(0,1.05*possum)
                for itn,ii in enumerate(neg_keys):
                    if ii not in negkeykeep_epi and neg_values[itn]>1e-5:
                        negkeykeep_epi+=[ii]
                negsum=0
                # print('Neg lengths: ',len(neg_keys),len(neg_processes))
                for i in range(0,len(negkeykeep_epi)):
                    if negkeykeep_epi[i] in neg_keys:
                        negval=neg_values[neg_keys.index(negkeykeep_epi[i])]
                        if negval>0:
                            p2_epi[i]=ax2n_epi.bar(cell_names[it], negval/model_time, bottom=negsum,color=nad_coloursdict[negkeykeep_epi[i]],width=width)
                            ax2n_epi.set_xlim(-0.5,0.5)
                            negsum+=negval/model_time
                ax2n_epi.set_ylabel('NAD(P)H usage' '\n' '(mmol/gDW/h)')
                ax2n_epi.set_ylim(0,1.05*negsum)
                # ax2n_epi.tick_params("x", labelrotation=50,ha='left',rotation_mode='anchor')
                plt.setp(ax2n_epi.xaxis.get_majorticklabels(),rotation=50,ha='right',rotation_mode='anchor')
            for i in range(len(poskeykeep),len(poskeykeep)+len(poskeykeep_epi)):
                if poskeykeep_epi[i-len(poskeykeep)] not in poskeykeep:
                    p1[i]=ax1n_oth.bar(cell_names[-2], 0, bottom=possum,color=nad_coloursdict[poskeykeep_epi[i-len(poskeykeep)]],width=width)
            # ax1n_oth.legend(p1.values(),[poskeykeep[x] for x in range(len(poskeykeep))]+[poskeykeep_epi[x] for x in range(len(poskeykeep_epi)) if poskeykeep_epi[x] not in (poskeykeep)],bbox_to_anchor=(1,1.))
            for i in range(len(negkeykeep),len(negkeykeep)+len(negkeykeep_epi)):
                if negkeykeep_epi[i-len(negkeykeep)] not in negkeykeep:
                    p2[i]=ax2n_oth.bar(cell_names[-2], 0, bottom=negsum,color=nad_coloursdict[negkeykeep_epi[i-len(negkeykeep)]],width=width)
            # ax2n_oth.legend(p2.values(),[negkeykeep[x] for x in range(len(negkeykeep))]+[negkeykeep_epi[x] for x in range(len(negkeykeep_epi)) if negkeykeep_epi[x] not in (negkeykeep)],bbox_to_anchor=(1,1.42))
            leg1=[poskeykeep[x] for x in range(len(poskeykeep))]+[poskeykeep_epi[x] for x in range(len(poskeykeep_epi)) if poskeykeep_epi[x] not in (poskeykeep)]
            leg2 = [negkeykeep[x] for x in range(len(negkeykeep))]+[negkeykeep_epi[x] for x in range(len(negkeykeep_epi)) if negkeykeep_epi[x] not in (negkeykeep)]
            leg_its = [x for x in range(len(leg2)) if leg2[x] not in leg1]
            ax2n_oth.legend(list(p1.values())+[list(p2.values())[x] for x in leg_its],
                            leg1+[leg2[x] for x in leg_its],bbox_to_anchor=(1,3))
            plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.4)
            dx = 3/72.; dy = 0/72. 
            offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)

            for label in ax2n_oth.yaxis.get_majorticklabels():
                label.set_transform(label.get_transform() + offset)
            for label in ax1n_oth.yaxis.get_majorticklabels():
                label.set_transform(label.get_transform() + offset)
            plt.savefig('Figures/Budget_all'+savebool+'.jpeg',bbox_inches='tight')
    else:
        plt.setp(ax2.xaxis.get_majorticklabels(),rotation=50,ha='right',rotation_mode='anchor')
        plt.savefig('Figures/Budget_all'+savebool+'_no_nad.jpeg',bbox_inches='tight')
    return outdict

import matplotlib

def generateATPbudgetforonecell_ed(model,solution,tag = "_epi00", outfile="",show_plot=False,percentage=False,save_plot_to="",thresh=0.005):

    #find colour hex codes here: https://www.colorhexa.com/color-names
    atp_coloursdict = {
        'Mitochondrial_ATP_Synthase_m_': '#5C6BC0',
        'PHOSGLYPHOS_RXN_c_':'#8D6E63' ,
        'PEPDEPHOS_RXN_p_': '#AB47BC',
        'PHOSPHORIBULOSEKINASE_RXN_p_': '#2E7D32',
        'PEPDEPHOS_RXN_c_': '#EF5350',
        'PHOSGLYPHOS_RXN_c_': '#FF7043',
        'GLUC1PADENYLTRANS_RXN_c_': '#26C6DA',
        'ACETYL_COA_CARBOXYLTRANSFER_RXN_p_': '#9CCC65',
        'UDPKIN_RXN_c_': '#BDBDBD',
        'GDPKIN_RXN_c_': '#fbceb1',
        'FRUCTOKINASE_RXN_c_': '#cd9575',
        'PROTON_ATPase_c_': '#000000',
        'Others-pos': '#000000',}
    colourDict = atp_coloursdict

    total_atp_dict = {}
    model_tag = tag
    if outfile!="":
        fout = open(outfile,"w")
    ATPdict = dict()
    total = 0
    excl = ['UDPKIN_RXN_c_','GDPKIN_RXN_c_','_mc_','_xc_','_pc_','Transfer']
    energymets=['ATP_','GTP_','UTP_']
    metcomps={'ATP_':("c","p","m","x"),'GTP_':("c","p","m"),'UTP_':("c","p")}
    
    
    for thismet in energymets:
        for p in metcomps[thismet]:
            met=model.metabolites.get_by_id(thismet+p+model_tag)
            met1=model.metabolites.get_by_id("a"+thismet+p+model_tag)
            for rxn in met.reactions:
                if any([rxn.id.__contains__(x) for x in excl]):
                    continue
                sto=rxn.metabolites.get(met)
                sto1=rxn.metabolites.get(met1)
                try:
                    if outfile!="":
                        fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
                    ATPdict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
                except:
                    continue
                try:
                    if solution.get(rxn.id)*(sto+sto1) > 0:
                        total = total + (solution.get(rxn.id)*(sto+sto1))
                except:
                    continue
    if outfile!="":
        fout.close()

    tempDict = dict()
    for rxn in ATPdict.keys():
        tempDict[rxn]=abs(ATPdict[rxn])

    #sort ATPdict by values
    import operator
    sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)

    ATPdict2 = dict()
    ATPdict2["Others-pos"]=0
    ATPdict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    for TEMP in sorted_by_value:
            rxn = TEMP[0]
            if ATPdict[rxn]>0:
                if ATPdict[rxn] < total*thresh:
                    if percentage:
                        ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+float(ATPdict[rxn]*100)/total
                    else:
                        ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+ATPdict[rxn]
                    continue
                base = pos_base
                if percentage:
                    ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
                    pos_base = pos_base + float(ATPdict[rxn]*100)/total
                else:
                    pos_base = pos_base + ATPdict[rxn]
                    ATPdict2[rxn]=ATPdict[rxn]
            else:
                if abs(ATPdict[rxn]) < total*thresh:
                    if percentage:
                        ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+float(ATPdict[rxn]*100)/total
                    else:
                        ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+ATPdict[rxn]
                    continue
                base = neg_base
                if percentage:
                    ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
                    neg_base = neg_base + float(ATPdict[rxn]*100)/total
                else:
                    neg_base = neg_base + ATPdict[rxn]
                    ATPdict2[rxn]=ATPdict[rxn]
            i=i+1
            baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    #obtaining total ATP expenditure
    total = 0
    for value in ATPdict2.values():
        if value > 0:
            total += value
    # print(model_tag, " total ATP consumed: ", total)
    total_atp_dict[model_tag] = total


    #making plots
    if show_plot:
        import matplotlib.pyplot as plt
        #plt.suptitle(title)
        plt.rcParams.update({'font.size': 10}) #sets a global fontsize
        plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['axes.linewidth']=2 # makes axes line thicker
        plt.figure(figsize=(3,4))
        for rxn in ATPdict2.keys():
            if colourDict.keys().__contains__(rxn[:-5]):
                plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn[:-5]])
            else:
                plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
        plt.xlim(0.8,1.2)
        if percentage:
            plt.ylabel("ATP produced/consumed (%)")
        else:
            plt.ylabel("ATP produced/consumed (in moles)")
        handles, labels = plt.gca().get_legend_handles_labels()
        labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
        handles2=[handles[labels.index(i)] for i in labels2]
        lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
        plt.axhline(0,linestyle="--",color="black")
        plt.tight_layout
        plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')

    return total_atp_dict, ATPdict2



import matplotlib.pyplot as plt

def generateNADbudgetforonecell_ed(model,solution,tag = "_epi00", outfile="",percentage=False,save_plot_to="",thresh=0):

    total_atp_dict = {}
    model_tag = tag
    if outfile!="":
        fout = open(outfile,"w")
    NADdict = dict()
    total = 0
    excl = ['UDPKIN_RXN_c_','GDPKIN_RXN_c_','_mc_','_xc_','_pc_','Transfer']
    energymets=['NADH_','NADPH_']
    metcomps={'NADH_':("c","p","m","x"),'NADPH_':("c","x","p","m")}
    
    
    for thismet in energymets:
        for p in metcomps[thismet]:
            met=model.metabolites.get_by_id(thismet+p+model_tag)
            for rxn in met.reactions:
                if any([rxn.id.__contains__(x) for x in excl]):
                    continue
                sto=rxn.metabolites.get(met)
                try:
                    if outfile!="":
                        fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto))+"\t"+met.compartment+"\n")
                    NADdict[rxn.id]=solution.get(rxn.id)*(sto)
                except:
                    continue
                try:
                    if solution.get(rxn.id)*(sto) > 0:
                        total = total + (solution.get(rxn.id)*(sto))
                except:
                    continue
    if outfile!="":
        fout.close()

    tempDict = dict()
    for rxn in NADdict.keys():
        tempDict[rxn]=abs(NADdict[rxn])

    #sort NADdict by values
    import operator
    sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)

    NADdict2 = dict()
    NADdict2["Others-pos"]=0
    NADdict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    for TEMP in sorted_by_value:
            rxn = TEMP[0]
            if NADdict[rxn]>0:
                if NADdict[rxn] < total*thresh:
                    if percentage:
                        NADdict2["Others-pos"]=NADdict2["Others-pos"]+float(NADdict[rxn]*100)/total
                    else:
                        NADdict2["Others-pos"]=NADdict2["Others-pos"]+NADdict[rxn]
                    continue
                base = pos_base
                if percentage:
                    NADdict2[rxn]=float(NADdict[rxn]*100)/total
                    pos_base = pos_base + float(NADdict[rxn]*100)/total
                else:
                    pos_base = pos_base + NADdict[rxn]
                    NADdict2[rxn]=NADdict[rxn]
            else:
                if abs(NADdict[rxn]) < total*thresh:
                    if percentage:
                        NADdict2["Others-neg"]=NADdict2["Others-neg"]+float(NADdict[rxn]*100)/total
                    else:
                        NADdict2["Others-neg"]=NADdict2["Others-neg"]+NADdict[rxn]
                    continue
                base = neg_base
                if percentage:
                    NADdict2[rxn]=float(NADdict[rxn]*100)/total
                    neg_base = neg_base + float(NADdict[rxn]*100)/total
                else:
                    neg_base = neg_base + NADdict[rxn]
                    NADdict2[rxn]=NADdict[rxn]
            i=i+1
            baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    #obtaining total ATP expenditure
    total = 0
    for value in NADdict2.values():
        if value > 0:
            total += value
    total_atp_dict[model_tag] = total

    return total_atp_dict, NADdict2

def vol_vs_change(rmodel, cell_dimensions, save_loc=None):
    col=['#1f78b4', '#b2df8a', '#8ae1ed', '#bebada', '#fdb462', '#d68aed', '#f5737a', '#fee08b', '#e69595', '#66c2a5', '#d690b2', '#3288bd', '#ffff99', '#a6cee3', '#ffffb3', '#6a3d9a', '#ff7f00', '#d9d9d9', '#9e0142', '#b3de69', '#33a02c', '#f46d43', '#abdda4', '#d53e4f', '#fccde5', '#fb8072', '#ccebc5', '#ffffbf', '#e6f598', '#bc80bd', '#5e4fa2', '#e89f82', '#80b1d3', '#ffed6f']
    model_time = 6.6
    #import rootslice data
    cells = [rxn.id.split('_')[-1] for rxn in rmodel.reactions if 'CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c' in rxn.id]
    rootslice = pd.read_excel(cell_dimensions, index_col=0)
    
    #calculating gDW per model -- this is the single value I should use for all scaling!!!
    whole_model_vol = rootslice.loc["whole_model", "avg_tot_dnut_vol"] #this value is the average volume of whole model across cell lengths
    gDW_per_tissue = whole_model_vol/1.06383E+13 # 2.541e+13 is um3/gDW according to https://doi.org/10.1371/journal.pone.0239075
    cell_vol={}
    axis_label_size=20
    legend_size = 18
    tick_size = 16
    scaling = 100
    for model_tag in ['_'+x for x in cells]:
        #cell wall volume
        cell_wall_vol = rootslice.loc[model_tag, "change_cellwall_dnut_vol"]
        ##cellulose content
        p_cellulose = 1.54*(10**-12)              # unit = g/um3 (density of cellulose in cell wall)
        g_cellulose = cell_wall_vol * p_cellulose           # unit = g
        mmol_cellulose = (g_cellulose/180.16) * 1000 / gDW_per_tissue                 # unit = mmol/gDW *
        #cell membrane volume
        cell_mem_vol = rootslice.loc[model_tag, "change_pm_dnut_vol"]
        ##phospholipid content
        p_phospholipid = 0.975*(10**-12)              	# unit = g/um3
        g_phospholipid = cell_mem_vol * p_phospholipid                      		# unit = g/cell
        mmol_phospholipid = (g_phospholipid/689.978637415773) * 1000 / gDW_per_tissue	# unit = mmol/gDW
        #cytoplasm volume
        cell_cyto_vol = rootslice.loc[model_tag, "change_cytoplasm_dnut_vol"]

        ##total protein content. Data borrowed from Sanu -- NEED TO FIND MAIZE ROOT SPECIFIC DATA!
        protein_conc = 21458.1747597 * (10**-18)     # unit = mmol/um3
        mmol_protein = protein_conc * cell_cyto_vol / gDW_per_tissue     #mmol/gDW

        Vv = rootslice.loc[model_tag, "change_vac_dnut_vol"]
        Vc = rootslice.loc[model_tag, "change_cytoplasm_dnut_vol"]
        Vcell = rootslice.loc[model_tag, "change_tot_dnut_vol"]
        av_vol= rootslice.loc[model_tag, "avg_tot_dnut_vol"]

        C_cell = 250                      #Bloom et al. 2012; units = mol/m3
        C_cell = 250 * (10**-15)          #units = mmol/um3
        mmol_osmolytes=C_cell*Vcell/gDW_per_tissue

        cell_vol[model_tag[1:]]=[av_vol,mmol_cellulose,mmol_phospholipid*scaling,mmol_protein*scaling,mmol_osmolytes,Vcell,Vc,Vv]

    data=pd.DataFrame(cell_vol)
    # print(data)
    data=data.applymap(lambda x: x/model_time)
    # print(data)
    # rootslice = rootslice[:11]
    rootslice.rename(columns={'avg_tot_dnut_vol': "Total torus volume"}, inplace=True)
    rootslice.rename(columns={'cell_id': "Cell type"}, inplace=True)
    Total_torus_volume = rootslice[['Total torus volume']]

    #create chart
    #v useful guide to dataframe.plot() https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
    #plot of all volumes
        #rootslice[['avg_tot_dnut_vol','avg_vac_dnut_vol','avg_pm_dnut_vol', "avg_cellwall_dnut_vol", "avg_cytoplasm_dnut_vol"]].plot(kind='bar', width = 1, figsize=(12,6), colormap="viridis")
    #rootslice[['avg_pm_dnut_vol']].plot(kind='bar', width = 1, figsize=(12,6), colormap="viridis")

    cell_names = ["Epidermis", "Cortical 1", "Cortical 2", "Cortical 3", "Cortical 4",
                    "Cortical 5", "Cortical 6", "Cortical 7", "Cortical 8", "Endodermal", "Pericycle"]

    
    # p1=data.loc[0].plot(kind='line', linewidth=5, figsize=(10,10),color=col[0],label="Initial cell volume")
    p1=data.loc[5].plot(kind='line', linewidth=5, figsize=(6,6),color=col[6],label=r"Cell (μm$^3$/h)")
    p1=data.loc[6].plot(kind='line', linewidth=5, color=col[7],label=r"Cytoplasm (μm$^3$/h)")
    p1=data.loc[7].plot(kind='line', linewidth=5, color=col[8],label=r"Vacuole (μm$^3$/h)")
    plt.ylabel(r"Increase in volume (μm$^3$/h)",fontsize=axis_label_size)
    p2=data.loc[1].plot(secondary_y=True, color=col[1], linewidth=5, label="Cellulose biosynthesis (mmol/gDW/h)")
    p3=data.loc[2].plot(secondary_y=True,color=col[2], linewidth=5, label="Lipid biosynthesis (x"+str(scaling)+") (mmol/gDW/h)")
    p4=data.loc[3].plot(secondary_y=True,color=col[3], linewidth=5, label="Protein biosynthesis (x"+str(scaling)+") (mmol/gDW/h)")
    p5=data.loc[4].plot(secondary_y=True,color=col[4], linewidth=5, label="Osmolyte accumulation (x"+str(scaling)+") (mmol/gDW/h)")
    # data.loc[4].plot(secondary_y=True, linewidth=5,color=col[5],label="Increase in cell volume (as a fraction of initial cell volume)")
    plt.ylabel("Calculated rate of... (mmol/gDW/h)",fontsize=axis_label_size)
    
    
    plt.ylim(bottom=0)
    # plt.legend(loc="center right")
    # handles=[]
    # labels=[]
    # for ax in [p1,p2]:#fig.axes:
    #     ax.set_ylim(bottom=0)
    #     ax.set_xticks(range(len(cell_names)), cell_names, rotation = 90,fontsize=14)
    #     ax.tick_params(axis='y',labelsize=tick_size)
    #     ax.tick_params(axis='x',labelsize=tick_size)
    #     for h,l in zip(*ax.get_legend_handles_labels()):
    #         handles.append(h)
    #         labels.append(l)
    # plt.legend(handles,labels,fontsize=legend_size,bbox_to_anchor=(0.5,1.2),loc='center')

    for ax in [p1,p2]:#fig.axes:
        ax.set_ylim(bottom=0)
        ax.set_xticks(range(len(cell_names)), cell_names, rotation = 90,fontsize=14)
        ax.tick_params(axis='y',labelsize=tick_size)
        ax.tick_params(axis='x',labelsize=tick_size)
    lines, labels = p1.get_legend_handles_labels()
    lines2, labels2 = p2.get_legend_handles_labels()
    labels2 = [x[:-8] for x in labels2]
    plt.legend(lines + lines2, labels + labels2,fontsize=legend_size,bbox_to_anchor=(0.5,1.4),loc='center')
    # print(data)
    #plt.title("Plot of total doughnut volume vs NAD and ATP demand")
    if save_loc:
        plt.savefig(save_loc, bbox_inches='tight')
    return

def generate_partitions(rmodel,thissol):
    met_dict={'PHOSPHO_ENOL_PYRUVATE':'PEP', 'UDP_GLUCURONATE':'UDP-Glucuronate', 'GLC_1_P':'Glucose-1P', 'FRU':'Fructose', 'G3P':'3-PGA', 'CPD_510':'Glucuronate-1P', 'L_LACTATE':'Lactate', 'FRUCTOSE_6P':'Fructose-6P', 'NITRATE':r'NO$_3$', 'FRUCTOSE_16_DIPHOSPHATE':'Fructose-1,6P', 'PYRUVATE':'Pyruvate', 'AMMONIUM':r'NH$_4$', 'SUCROSE':'Sucrose', 'GAP':'GAP', '2_PG':'2-PG', 'ETOH':'Ethanol', 'MGII':r'Mg$^{2+}$', 'DPG':'DPGA', 'Pi':'Pi','SER':'Serine','GLC_6_P':'Glucose-6P','DIHYDROXY_ACETONE_PHOSPHATE':'DHAP','UDP':'UDP','MAL':'Malate','ARG':'Arginine','OXALACETIC_ACID':'OAA','GLT':'GLT','L_ASPARTATE':'Aspartate','MET':'Methionine','TYR':'Tyrosine','PHOSPHATIDYL_CHOLINE':'Phosphatidyl Choline'}
    transfers = set(['_'.join(x.id.split('_')[-5:]) for x in rmodel.reactions if 'Transfer' in x.id and not 'WATER' in x.id])
    large_symp = sorted([rxn for rxn in rmodel.reactions if any([x in rxn.id for x in list(transfers)])],key=lambda x:abs(thissol[x.id]),reverse=True)
    subset_symp = set(['_'.join(x.id.split('_')[:-5]) for x in large_symp[:150]])
    large_symp = sorted([rxn for rxn in rmodel.reactions if any([x in rxn.id for x in list(transfers)]) and not rxn.id[0]=='a' and not 'WATER' in rxn.id],key=lambda x:abs(thissol[x.id]),reverse=True)
    subset_symp = set(['_'.join(x.id.split('_')[:-5]) for x in large_symp[:]])
    subs_order = {}
    for met in subset_symp:
        subs_order[met]=max([abs(thissol[x.id]) for x in large_symp if met in x.id])
    partitions = 6
    partition_size=6
    sym_partitions={}
    # temp_symp=sorted(list(subset_symp)[:partition_size*partitions],key = lambda x: subs_order[x],reverse=True)
    temp_symp=sorted(list(subset_symp),key = lambda x: subs_order[x],reverse=True)
    for ii in range(partitions):
        if 'FRU_c' in temp_symp:
            ref_met = 'FRU_c'
        elif 'NITRATE_c' in temp_symp:
            ref_met = 'NITRATE_c'
        elif 'PYRUVATE_c' in temp_symp:
            ref_met='PYRUVATE_c'
        elif 'L_LACTATE_c' in temp_symp:
            ref_met = 'L_LACTATE_c'
        elif 'MAL_c' in temp_symp:
            ref_met = 'MAL_c'
        elif 'CPD_510_c' in temp_symp:
            ref_met = 'CPD_510_c'
        # elif 'PYRUVATE_c' in temp_symp:
        #     ref_met = 'PYRUVATE_c'
        else:
            ref_met = sorted(temp_symp,key=lambda x: formula_bkdown(rmodel.metabolites.get_by_id(x+'_cor01').formula)[0])[-1]
        if ref_met=='FRU_c':
            sym_partitions[ii] = [x for x in temp_symp if x[:-2] in list(met_dict.keys()) and 'ose' in met_dict[x[:-2]]]
            sym_partitions[ii]=sym_partitions[ii][:min(len(sym_partitions[ii]),partition_size)]
        # For nitrate accum ratio:
        elif ref_met=='NITRATE_c':
            sym_partitions[ii] = [x for x in ['NITRATE_c','AMMONIUM_c','GLT_c'] if x in temp_symp]#+[x+'_c' for x in aas if x+'_c' in temp_symp]
        elif ref_met=='MAL_c':
            sym_partitions[ii] = [x for x in ['MAL_c','OXALACETIC_ACID_c'] if x in subset_symp]
        elif ref_met=='PYRUVATE_c':
            sym_partitions[ii] = [x for x in ['G3P_c','DPG_c','2_PG_c','PHOSPHO_ENOL_PYRUVATE_c','PYRUVATE_c'] if x in subset_symp]
        # elif ref_met=='PYRUVATE_c':
        #     sym_partitions[ii] = ['PYRUVATE_c','PHOSPHO_ENOL_PYRUVATE_c','2_PG_c','G3P_c','DIHYDROXY_ACETONE_PHOSPHATE_c','DPG_c']
        elif ref_met=='CPD_510_c':
            sym_partitions[ii] = ['CPD_510_c','UDP_GLUCURONATE_c']
        elif ref_met == 'L_LACTATE_c':
            sym_partitions[ii] = ['L_LACTATE_c','ETOH_c']
        # elif len(temp_symp)>partition_size:
        #     sym_partitions[ii] = sorted(temp_symp,key=lambda x: dls(rmodel.metabolites.get_by_id(ref_met+'_cor01'),rmodel.metabolites.get_by_id(x+'_cor01')))[:partition_size]
        else:
            sym_partitions[ii]=temp_symp[:partition_size]
        
        sym_partitions[ii] = sorted(list(set(sym_partitions[ii])),key = lambda x: subs_order[x],reverse=True)
        temp_symp = sorted(list(set(temp_symp).difference(set(sym_partitions[ii]))),key = lambda x: subs_order[x],reverse=True)
    return sym_partitions

# import matplotlib.pyplot as plt
def plot_all_sym_transport(sol,rmodel,col,saveas=None,xtick=1):
    met_dict={'PHOSPHO_ENOL_PYRUVATE':'PEP', 'UDP_GLUCURONATE':'UDP-Glucuronate', 'GLC_1_P':'Glucose-1P', 'FRU':'Fructose', 'G3P':'3-PGA', 'CPD_510':'Glucuronate-1P', 'L_LACTATE':'Lactate', 'FRUCTOSE_6P':'Fructose-6P', 'NITRATE':r'NO$_3$', 'FRUCTOSE_16_DIPHOSPHATE':'Fructose-1,6P', 'PYRUVATE':'Pyruvate', 'AMMONIUM':r'NH$_4$', 'SUCROSE':'Sucrose', 'GAP':'GAP', '2_PG':'2-PG', 'ETOH':'Ethanol', 'MGII':r'Mg$^{2+}$', 'DPG':'DPGA', 'Pi':'Pi','SER':'Serine','GLC_6_P':'Glucose-6P','DIHYDROXY_ACETONE_PHOSPHATE':'DHAP','UDP':'UDP','MAL':'Malate','ARG':'Arginine','OXALACETIC_ACID':'OAA','GLT':'GLT','L_ASPARTATE':'Aspartate','MET':'Methionine','TYR':'Tyrosine','PHOSPHATIDYL_CHOLINE':'Phosphatidyl Choline'}
    model_time=6.6
    transfers = set(['_'.join(x.id.split('_')[-5:]) for x in rmodel.reactions if 'Transfer' in x.id])
    transfers_list=list(sorted(list(transfers),key=lambda x:x))
    transfers_list =[transfers_list[-1]]+transfers_list[:-1]
    x_ticks = {'cor03_Transfer_03_to_04':'COR 4', 'cor06_Transfer_06_to_07':'COR 7', 'epi00_Transfer_epi00_to_cor00': 'EPI 1', 'cor01_Transfer_01_to_02': 'COR 2', 'cor05_Transfer_05_to_06':'COR 6', 'cor00_Transfer_00_to_01':'COR 1', 'end00_Transfer_end00_to_per00':'END 1', 'cor07_Transfer_cor07_to_end00':'COR 8', 'cor02_Transfer_02_to_03':'COR 3', 'cor04_Transfer_04_to_05':'COR 5'}
    x_ticks2 = {'cor03_Transfer_03_to_04':'COR 5', 'cor06_Transfer_06_to_07':'COR 8', 'epi00_Transfer_epi00_to_cor00': 'COR 1', 'cor01_Transfer_01_to_02': 'COR 3', 'cor05_Transfer_05_to_06':'COR 7', 'cor00_Transfer_00_to_01':'COR 2', 'end00_Transfer_end00_to_per00':'PER 1', 'cor07_Transfer_cor07_to_end00':'END 1', 'cor02_Transfer_02_to_03':'COR 4', 'cor04_Transfer_04_to_05':'COR 6'}
    sym_partitions=generate_partitions(rmodel,sol)
    nitrate_it = [it for it in range(len(sym_partitions)) if 'NITRATE_c' in sym_partitions[it]][0]
    temp = sym_partitions[nitrate_it]
    sym_partitions[nitrate_it]=sym_partitions[0]
    sym_partitions[0]=temp

    # plt.figure(figsize=(6,8))
    fig,axes = plt.subplots(len(sym_partitions),2, figsize=(10,27),sharey=False,sharex=False,dpi=1000)
    for frow in range(len(sym_partitions)):
        subset_symp=sym_partitions[frow]
        symdict = {}
        plts=[]
        plted=[]
        met_labels = [met_dict[x[:-2]] if x[:-2] in (met_dict.keys()) else x for x in subset_symp]
        max_met =0
        for it,met in enumerate(list(subset_symp)):
            symdict[met]=[sol[met+'_'+x] for x in transfers_list]
            mets = met[:-2]
            if mets=='NITRATE':
                mets='Nitrate'
            txrxn = mets+'_tx_epi00'
            # xyrxn='EX_X_'+mets+'_xy_per00'
            xyrxn=met+'_xy_per00'
            phrxn = mets+'_ph_c_per00'
            # print(xyrxn)
            symdict[met] = [sol[txrxn]*rmodel.reactions.get_by_id(txrxn).get_coefficient(met[:-2]+'_e_epi00') if txrxn in rmodel.reactions else 0]+symdict[met]+[sum([sol[xyrxn]*rmodel.reactions.get_by_id(xyrxn).get_coefficient(met[:-2]+'_xy_per00') if xyrxn in rmodel.reactions else 0]+[sol[phrxn]*rmodel.reactions.get_by_id(phrxn).get_coefficient(mets+'_ph_per00') if phrxn in rmodel.reactions else 0])]
            if 1:#sum(symdict[met])>=0:
                axes[frow][0].stackplot(range(len(symdict[met])),[x/model_time if x>=0 else 0 for x in symdict[met]],labels = transfers_list, colors = [col[it]], baseline='zero',alpha=0.5)
                plts+=axes[frow][1].stackplot(range(len(symdict[met])),[-x/model_time if x<0 else 0 for x in symdict[met][-1::-1]],labels = transfers_list, colors = [col[it]], baseline='zero',alpha=0.5) 
                plted+=[it]
            max_met = max([max_met]+[abs(x) for x in symdict[met]])
        axes[frow][1].legend(plts,[list(met_labels)[x] for x in plted],bbox_to_anchor=(1,0.5),loc='center left',fontsize=18)
        # if xtick:
        #     axes[frow][0].set_xticks(range(len(symdict[met])),labels=['Soil']+[x_ticks[x] for x in transfers_list]+['PER 1'],rotation=90,fontsize=16)
        #     axes[frow][1].set_xticks(range(len(symdict[met])),labels=['Xy/Ph']+[x_ticks2[x] for x in transfers_list[-1::-1]]+['EPI 1'],rotation=90,fontsize=16)
        # else:
        #     axes[frow][0].set_xticks(range(len(symdict[met])),labels=['' for x in symdict[met]],rotation=90,fontsize=16)
        #     axes[frow][1].set_xticks(range(len(symdict[met])),labels=['' for x in symdict[met]],rotation=90,fontsize=16)
        axes[frow][0].tick_params(axis='y',labelsize=16)
        axes[frow][1].tick_params(axis='y',labelsize=16)
        axes[frow][0].set_ylim(0,max_met/model_time*1.1)
        axes[frow][1].set_ylim(0,max_met/model_time*1.1)
        if frow==(len(sym_partitions)-1):
            axes[frow][0].set_xlabel('Transferred from',fontsize=20)
            axes[frow][1].set_xlabel('Transferred from',fontsize=20)
            axes[frow][0].set_xticks(range(len(symdict[met])),labels=['Soil']+[x_ticks[x] for x in transfers_list]+['PER 1'],rotation=90,fontsize=16)
            axes[frow][1].set_xticks(range(len(symdict[met])),labels=['Xy/Ph']+[x_ticks2[x] for x in transfers_list[-1::-1]]+['EPI 1'],rotation=90,fontsize=16)
            these_ticks=axes[frow][1].get_yticks()
            axes[frow][0].set_yticklabels([round(x,2) for x in these_ticks])
            axes[frow][1].set_yticklabels([round(x,2) for x in these_ticks])
        else:
            axes[frow][0].set_xticks(range(len(symdict[met])),labels=['' for x in symdict[met]],rotation=90,fontsize=16)
            axes[frow][1].set_xticks(range(len(symdict[met])),labels=['' for x in symdict[met]],rotation=90,fontsize=16)
        if frow == 3:
            axes[frow][0].set_ylabel('Flux to next cell (mmol/h)',fontsize=20)
            axes[frow][0].yaxis.set_label_coords(-0.3,1.15)
    if saveas:
        plt.savefig(saveas, bbox_inches='tight')
    plt.show()

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
def resetBiomassComposition_eBC2(rxn,ratios):
    # original_rxn = [rxn.get_coefficient(y.id) for y in rxn.metabolites]
    cellulose_met = [met for met in rxn.metabolites if 'Cellulose' in met.id][0]
    lipid_met = [met for met in rxn.metabolites if 'L_PHOSPHATIDATE' in met.id][0]
    protein_met = [met for met in rxn.metabolites if 'Protein' in met.id][0]
    rxn.add_metabolites({cellulose_met:rxn.get_coefficient(cellulose_met)*(1-ratios[0])/ratios[0],
    lipid_met:rxn.get_coefficient(lipid_met)*(1-ratios[1])/ratios[1],
    protein_met:rxn.get_coefficient(protein_met)*(1-ratios[2])/ratios[2]})
def vary_biomass_composition(rmodel2,solo=None,cheung_maint=None):
    if not cheung_maint:
        cheung_maint = rmodel2.reactions.get_by_id("ATPase_balancer").lower_bound
    # carb 1.52, lipid 4.9, protein 32
    biomass_limits = [1.52,4.9,32]
    rmodel2.reactions.get_by_id("ATPase_balancer").lower_bound=cheung_maint
    if not solo:
        solo=rmodel2.optimize()
    sols={}
    rmodel2.reactions.Phloem_carbon_import.lower_bound=solo['Phloem_carbon_import']
    rmodel2.reactions.Phloem_carbon_import.upper_bound=solo['Phloem_carbon_import']
    biomass_rxns = [rxn.id for rxn in rmodel2.reactions if 'Biomass_expanding_cell' in rxn.id]
    for ii in range(3):
        sols[0]={}
        for jj in np.linspace(1,biomass_limits[ii],5):
            these_rats = [1,1,1]
            these_rats[ii]=jj
            for rxn in biomass_rxns:
                changeBiomassComposition_eBC2(rmodel2.reactions.get_by_id(rxn),these_rats)
            print(jj)
            sols[0][jj]=rmodel2.optimize()
            print(sols[0][jj]['L_LACTATE_ec_epi00'],sols[0][jj]['ETOH_ec_epi00'])
            for rxn in biomass_rxns:
                resetBiomassComposition_eBC2(rmodel2.reactions.get_by_id(rxn),these_rats)
    rmodel2.reactions.Phloem_carbon_import.lower_bound=0
    rmodel2.reactions.Phloem_carbon_import.upper_bound=1000

    rmodel2.reactions.Phloem_carbon_import.lower_bound=solo['Phloem_carbon_import']
    rmodel2.reactions.Phloem_carbon_import.upper_bound=solo['Phloem_carbon_import']
    maint_vals = list(np.linspace(0,cheung_maint,15))+list(np.linspace(cheung_maint,2.5,20))
    sols[3]={}
    for maint in maint_vals:
        rmodel2.reactions.get_by_id("ATPase_balancer").lower_bound=maint
        sols[3][maint]=rmodel2.optimize()
        print(sols[3][maint]['L_LACTATE_ec_epi00'],sols[3][maint]['ETOH_ec_epi00'])
    rmodel2.reactions.get_by_id("ATPase_balancer").lower_bound=cheung_maint
    rmodel2.reactions.Phloem_carbon_import.lower_bound=0
    rmodel2.reactions.Phloem_carbon_import.upper_bound=1000
    return sols # c,l,p,m


def fermentation_export(rmodel,sols=None,savebool=None):
    default_cols=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
    '#7f7f7f', '#bcbd22', '#17becf']
    model_time=6.6
    if not sols:
        sols=vary_biomass_composition(rmodel)
    itc= list(sols[0].keys())
    itl=list(sols[1].keys())
    itp=list(sols[2].keys())
    # itc=itp
    # itl=itp
    itm=list(sols[3].keys())
    atpase_bal_sum=-sum([rmodel.reactions.ATP_balancer.get_coefficient(x.id) for x in rmodel.reactions.ATP_balancer.metabolites])
    xc = [100*(x-1) for x in itc]#[-(x-1)*cellulose_flux for x in itc]
    xl = [100*(x-1) for x in itl]#[-(x-1)*lipid_flux for x in itl]
    xp = [100*(x-1) for x in itp]#[-(x-1)*protein_flux for x in itp]
    xm = [x*atpase_bal_sum for x in itm]
    rxn1 = 'L_LACTATE_tx_epi00'
    rxn2 = 'ETOH_tx_epi00'
    fig,ax=plt.subplots(dpi=1000)
    ax2=ax.twiny()
    # plt.plot(lbits,[(sol_dict[str(x)+'-1-1'][rxn1]+sol_dict[str(x)+'-1-1'][rxn2])*1000/model_time for x in bits],linewidth=3,color=default_cols[0])
    ax.plot(xc,[(sols[0][x][rxn1]+sols[0][x][rxn2])*1000/model_time for x in itc],linewidth=3,color=default_cols[0])
    ax.plot(xl,[(sols[1][x][rxn1]+sols[1][x][rxn2])*1000/model_time for x in itl],linewidth=3,color=default_cols[1])
    ax.plot(xp,[(sols[2][x][rxn1]+sols[2][x][rxn2])*1000/model_time for x in itp],linewidth=3,color=default_cols[2])
    ax.plot(xm,[np.nan]*len(xm),linewidth=3,color=default_cols[3])
    ax2.plot(xm,[(sols[3][x][rxn1]+sols[3][x][rxn2])*1000/model_time for x in itm],linewidth=3,color=default_cols[3])
    # ax2.plot(sols.keys(),[(sols[key][rxn1]+sols[key][rxn2])*1000/model_time*total_maint_atp for key in list(sols.keys())],linewidth=3,color=default_cols[3])
    # ax2.plot([round(float(x)*total_maint_atp,2) for x in list(sol_dict.keys())[7:]],[(sol_dict[key][rxn1]+sol_dict[key][rxn2])*1000/model_time for key in list(sol_dict.keys())[7:]],linewidth=3,color=default_cols[3])
    ax.legend(['Cellulose','Lipids','Protein','Maintenance'],fontsize=16,bbox_to_anchor=(1.55,0.75))
    ax.set_xlim(0,max(xc))
    ylims=ax.get_ylim()
    ax.set_ylim(0,ylims[1])
    ax2.set_ylim(0,ylims[1])
    # temp=ax.set_xlabel(r'Amount component was ' '\n' 'increased in biomass''\n''(mmol/gDW/h)',fontsize=20)
    temp=ax.set_xlabel(r'Percentage component was ' '\n' 'increased in biomass',fontsize=20)
    temp=ax2.set_xlabel('Total tissue maintenance' '\n' '(mmol/gDW/h)',fontsize=20)
    ax.set_ylabel('Lactate+ethanol export' '\n' '(umol/hr)',fontsize=20)
    ax.tick_params(labelsize=16)
    ax2.tick_params(labelsize=16)
    plt.yticks(fontsize=16)
# plt.show()
    if savebool:
        plt.savefig('Figures/Biomass_comp_fermentation_fixed_ph_sum.png',bbox_inches='tight')
    return

def make_prot_dict(rmodel,sol,cell='_cor',excl=[],thresh=1e-5,metroot = 'PROTON_c'):
    label = 'roots_prot_rxn_pwy_dict'
    with open('roots_'+label+'.pkl','rb') as file:
        prot_rxn_pwy_dict=pkl.load(file)
    mets = [met for met in rmodel.metabolites if metroot in met.id and cell in met.id and not any([x in met.id for x in ['_e_','_ph_','_xy_']+excl])]
    # print(mets)
    sol_prot_dict={}
    rxns=[]
    for met in mets:
        for rxn in met.reactions:
            if sol[rxn.id]:
                if rxn.id not in sol_prot_dict.keys():
                    sol_prot_dict[rxn.id]=0
                sol_prot_dict[rxn.id]+=sol[rxn.id]*rxn.get_coefficient(met.id)
                rxns+=[rxn]    
    sol_prot_dict_nosol={}
    for met in mets:
        for rxn in met.reactions:
            if sol[rxn.id]:
                if rxn.id not in sol_prot_dict_nosol.keys():
                    sol_prot_dict_nosol[rxn.id]=0
                sol_prot_dict_nosol[rxn.id]+=rxn.get_coefficient(met.id)
            
    sol_prot_dict_root={}
    rxns=list(set(rxns))
    for rxn in rxns:
        rxnroot = return_rxn_root(rxn.id)
        if rxnroot not in sol_prot_dict_root.keys():
            sol_prot_dict_root[rxnroot]=0
        sol_prot_dict_root[rxnroot]+=sol_prot_dict[rxn.id]
    
    prot_budget_dict={}
    for key in list(sol_prot_dict_root.keys()):
        pwy = prot_rxn_pwy_dict[key]
        if pwy not in prot_budget_dict.keys():
            prot_budget_dict[pwy]=0
        prot_budget_dict[pwy]+=sol_prot_dict_root[key]
    # print(prot_budget_dict)
    for key in list(prot_budget_dict.keys()):
        if abs(prot_budget_dict[key])<thresh:
            prot_budget_dict.pop(key)
    return prot_budget_dict

def return_rxn_root(rxnid):
    rxnroot='_'.join(rxnid.split('_')[:-1])
    return rxnroot

import matplotlib as mpl
def plot_protons(prot_budget_dict,model_time=6.6,breakthresh = 0.1,savebool='',split=1,sumbool=0,ylim=None):
    sorted_dict=sorted(prot_budget_dict.items(),key=lambda x:(x[1]),reverse=True)
    keys = [x[0] for x in sorted_dict]
    values = [x[1] for x in sorted_dict]
    if split:
        fig, (axt,axm,axb) = plt.subplots(3,1,sharex=True,figsize=(len(prot_budget_dict)*6.4/20,4.8))
        mpl.rcParams['savefig.dpi']=1000
        mpl.rcParams['font.size']=16
        mpl.rcParams['legend.fontsize']=16
        mpl.rcParams['lines.linewidth']=3
        axt.bar(range(len(prot_budget_dict)),[x/model_time for x in values])
        axm.bar(range(len(prot_budget_dict)),[x/model_time for x in values])
        axb.bar(range(len(prot_budget_dict)),[x/model_time for x in values])
        if ylim:
            maxlim = ylim
        else:
            maxlim = max(max(values),abs(min(values)))/model_time*1.05
        axm.set_ylim(-breakthresh,breakthresh)
        axt.set_ylim(breakthresh,maxlim)
        axb.set_ylim(-maxlim,-breakthresh)
        temp=axb.set_xticks(range(len(prot_budget_dict)), keys,rotation=315,ha='left',rotation_mode='anchor',fontsize=14)
        axm.spines['top'].set_visible(False)
        axm.spines['bottom'].set_visible(False)
        axt.spines['bottom'].set_visible(False)
        axb.spines['top'].set_visible(False)
        axm.tick_params(bottom=False)
        axm.tick_params(top=False)
        # axt.xaxis.tick_top()
        axt.tick_params(bottom=False)
        axt.tick_params(labeltop=False)
        axb.xaxis.tick_bottom()
        # axt.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        # axm.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        # axb.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        axm.set_ylabel('Protons released' '\n' '(µmol/gDW/hr)',fontsize=20)
        if sumbool:
            # axt2=axt.twinx()
            axm2=axm.twinx()
            # axb2=axb.twinx()
            axm2.plot(range(len(prot_budget_dict)),[sum(values[:(x+1)]) for x in range(len(values))],color = '#8ae1ed')
            axm2ymax=max([abs(sum(values[:(x+1)])) for x in range(len(values))])
            axm2.set_ylim(-axm2ymax,axm2ymax)
            axm2.spines['top'].set_visible(False)
            axm2.spines['bottom'].set_visible(False)
            axm2.tick_params(bottom=False)
            axm2.set_ylabel('Cumulative protons released (µmol/gDW/hr)',fontsize=20)
    else:
        fig = plt.figure(figsize=(len(prot_budget_dict)*6.4/20,4.8))
        plt.bar(range(len(prot_budget_dict)),[x/model_time for x in values])
        maxlim = max(max(values),abs(min(values)))/model_time*1.05
        plt.ylim(-maxlim,maxlim)
        temp=plt.xticks(range(len(prot_budget_dict)), keys,rotation='vertical')
        plt.ylabel('Protons released (µmol/gDW/hr)')
    if savebool:
        plt.savefig('Figures/Protons_'+savebool+'.png',bbox_inches='tight')