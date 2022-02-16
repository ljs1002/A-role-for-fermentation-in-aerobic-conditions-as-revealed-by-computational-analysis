### function to set defaults for pyplots
def set_pyplot_defaults():
    import matplotlib.pyplot as plt
    from matplotlib import cycler
    colors = cycler('color',
                    ['#EE6666', '#3388BB', '#9988DD',
                     '#EECC55', '#88BB44', '#FFBBBB'])
    #plt.rc('axes', facecolor='#E6E6E6', edgecolor='white',
          # axisbelow=True, grid=True, prop_cycle=colors)
    plt.rc('grid', color='w', linestyle='solid')
    plt.rc('xtick', direction='out', color='gray')
    plt.rc('ytick', direction='out', color='gray')
    plt.rc('patch', edgecolor='#E6E6E6')
    plt.rc('lines', linewidth=2)

    #change font sizes
    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 30

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



### function to convert names in an escher json to clean names for publication
### takes map, cleans names, and saves map as json with same name i.e. overwrites old map
#args: json map
def map_renamer(map_raw, new_rxn_names, new_met_names):
    import json
    def getjsonfrommap(map):
        with open(map, 'r+') as f:
            data = json.load(f)
        return data

    def writejson(data, map):
        with open(map, 'w') as f:
            f.seek(0)
            json.dump(data, f, indent=4)
            f.truncate()     # remove remaining part

    map = getjsonfrommap(map_raw)

    for k in map[1]["reactions"].keys():
        for r in new_rxn_names:
            if r in map[1]["reactions"][k]["bigg_id"]:
                map[1]["reactions"][k]["name"] = new_rxn_names[r]
            else:
                continue
        if ":" in map[1]["reactions"][k]["name"]:
            clean_name = map[1]["reactions"][k]["name"].split(":")[1]
            map[1]["reactions"][k]["name"] = clean_name
        else:
            continue

    for n in map[1]["nodes"].keys():
        if "bigg_id" in map[1]["nodes"][n]:
            for m in new_met_names.keys():
                if m in map[1]["nodes"][n]["bigg_id"]:
                    map[1]["nodes"][n]["name"] = new_met_names[m]
                else:
                    continue
        else:
            continue

    writejson(map, map_raw)


### function to compare fluxes between two models with different constraints
#args: tissue_model1 (no constraint), tissue_model2(with constraints), rxn id as string
def constraint_comparer(model1, model2, rxn_id):
    for rxn in model1.reactions:
        if rxn_id in rxn.id:
            flux2 = model2.reactions.get_by_id(rxn.id).x
            print(rxn.id, "\n", "unconstrained: ",rxn.x,"\n", "constrained: ", flux2)


### function to print reactions with flux greater than and less than a certain amount
#args: model, key terms in rxn.id e.g. "c_xy", flux cutoff e.g. 0.01
def rxns_w_flux(model, key_term, flux_cutoff):
    for rxn in model.reactions:
        if key_term in rxn.id:
            if rxn.x < -flux_cutoff or rxn.x > flux_cutoff:
                print(rxn.compartments, rxn, "\n", rxn.flux)


###plot dnut volumes vs atp and nad use
def vol_vs_energy(model, cell_dimensions, total_nad_dict, total_atp_dict, save_loc):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    #convert nad and atp dictionaries to dataframes
    nad_df = pd.DataFrame.from_dict(total_nad_dict, orient='index')
    atp_df = pd.DataFrame.from_dict(total_atp_dict, orient='index')

    #import rootslice data
    rootslice = pd.read_excel(cell_dimensions, index_col=0)
    rootslice = rootslice[:11]
    rootslice.rename(columns={'avg_tot_dnut_vol': "Total torus volume"}, inplace=True)
    rootslice.rename(columns={'cell_id': "Cell type"}, inplace=True)
    Total_torus_volume = rootslice[['Total torus volume']]

    #create chart
    #v useful guide to dataframe.plot() https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
    #plot of all volumes
        #rootslice[['avg_tot_dnut_vol','avg_vac_dnut_vol','avg_pm_dnut_vol', "avg_cellwall_dnut_vol", "avg_cytoplasm_dnut_vol"]].plot(kind='bar', width = 1, figsize=(12,6), colormap="viridis")
    #rootslice[['avg_pm_dnut_vol']].plot(kind='bar', width = 1, figsize=(12,6), colormap="viridis")

    cell_names = ["Epidermal", "Cortical 1", "Cortical 2", "Cortical 3", "Cortical 4",
                    "Cortical 5", "Cortical 6", "Cortical 7", "Cortical 8", "Endodermal", "Pericycle"]

    Total_torus_volume.plot(kind='line', linewidth=5, figsize=(16,6))
    plt.ylabel("volume (μm^3)")
    nad_df[0].plot(secondary_y=True, color='r', linewidth=5, label="NAD(P)H demand")
    atp_df[0].plot(secondary_y=True,color='g', linewidth=5, label="ATP demand")
    plt.ylabel("mmol/gDW")
    plt.xlabel("")
    plt.xticks(range(len(cell_names)), cell_names, rotation = 10)
    plt.ylim(bottom=0)
    plt.legend(loc="center right")
    #plt.title("Plot of total doughnut volume vs NAD and ATP demand")

    plt.savefig(save_loc, bbox_inches='tight')

### plot of carbon usage in model
def plot_carbon_uptake(stem_model,psol):
    import matplotlib.pyplot as plt
    import numpy as np

    carbon_dict={}
    for rxn in stem_model.reactions:
        if any([('CARBON_DIOXIDE' in x.id or 'HCO3' in x.id) for x in rxn.metabolites])and abs(psol[rxn.id])>0:
            carbon_dict[rxn.id]=0
            for met in rxn.metabolites:
                if ('CARBON_DIOXIDE' in met.id or 'HCO3' in met.id):
                    carbon_dict[rxn.id]+=psol[rxn.id]*rxn.get_coefficient(met.id)
    thresh=0.1
    labels=[x for x in carbon_dict.keys()if abs(carbon_dict[x])>thresh]
    data1=[carbon_dict[x] if '_tx' in x else -carbon_dict[x] for x in labels]
    plt.figure(figsize=(7,15))
    width=0.6
    plt.yticks(range(len(data1)), labels,rotation=0)
    plt.barh(np.arange(len(data1)), data1)
    # plt.bar(np.arange(len(data2))+ width, data2, width=width)
    plt.legend(['carbon -> model'])
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()



###plot all fluxes of a particular subsytem
def plot_subsys_btwn_cells(cons_model, csol,ss=['glycolysis'],thresh='',plttitle=''):
    import matplotlib.pyplot as plt
    import numpy as np

    if not plttitle:
        plttitle=ss[0]
    if not thresh:
        thresh=0

    all_subsys_rxns=[x.id for x in cons_model.reactions if any([y in str(x.notes)for y in ss])]
    _epi00_rxns=[x for x in all_subsys_rxns if '_epi00' in x]
    _cor00_rxns=[x.replace('_epi00','_cor00') if x.replace('_epi00','_cor00') in all_subsys_rxns else '' for x in _epi00_rxns]
    _cor01_rxns=[x.replace('_epi00','_cor01') if x.replace('_epi00','_cor01') in all_subsys_rxns else '' for x in _epi00_rxns]
    _cor02_rxns=[x.replace('_epi00','_cor02') if x.replace('_epi00','_cor02') in all_subsys_rxns else '' for x in _epi00_rxns]
    _cor03_rxns=[x.replace('_epi00','_cor03') if x.replace('_epi00','_cor03') in all_subsys_rxns else '' for x in _epi00_rxns]
    _cor04_rxns=[x.replace('_epi00','_cor04') if x.replace('_epi00','_cor04') in all_subsys_rxns else '' for x in _epi00_rxns]
    _cor05_rxns=[x.replace('_epi00','_cor05') if x.replace('_epi00','_cor05') in all_subsys_rxns else '' for x in _epi00_rxns]
    _cor06_rxns=[x.replace('_epi00','_cor06') if x.replace('_epi00','_cor06') in all_subsys_rxns else '' for x in _epi00_rxns]
    _cor07_rxns=[x.replace('_epi00','_cor07') if x.replace('_epi00','_cor07') in all_subsys_rxns else '' for x in _epi00_rxns]
    _end00_rxns=[x.replace('_epi00','_end00') if x.replace('_epi00','_end00') in all_subsys_rxns else '' for x in _epi00_rxns]
    _per00_rxns=[x.replace('_epi00','_per00') if x.replace('_epi00','_per00') in all_subsys_rxns else '' for x in _epi00_rxns]

    data_epi00=[csol[x] if x else 0 for x in _epi00_rxns]
    data_cor00=[csol[x] if x else 0 for x in _cor00_rxns]
    data_cor01=[csol[x] if x else 0 for x in _cor01_rxns]
    data_cor02=[csol[x] if x else 0 for x in _cor02_rxns]
    data_cor03=[csol[x] if x else 0 for x in _cor03_rxns]
    data_cor04=[csol[x] if x else 0 for x in _cor04_rxns]
    data_cor05=[csol[x] if x else 0 for x in _cor05_rxns]
    data_cor06=[csol[x] if x else 0 for x in _cor06_rxns]
    data_cor07=[csol[x] if x else 0 for x in _cor07_rxns]
    data_end00=[csol[x] if x else 0 for x in _end00_rxns]
    data_per00=[csol[x] if x else 0 for x in _per00_rxns]

    no_flux=[]
    for ii in range(len(data_epi00)):
        if abs(data_epi00[ii])<=thresh and abs(data_cor00[ii])<=thresh and abs(data_cor01[ii])<=thresh and abs(data_cor02[ii])<=thresh and abs(data_cor03[ii])<=thresh and abs(data_cor04[ii])<=thresh and abs(data_cor05[ii])<=thresh and abs(data_cor06[ii])<=thresh and abs(data_cor07[ii])<=thresh and abs(data_end00[ii])<=thresh and abs(data_per00[ii])<=thresh:
            no_flux.append(ii)

    _epi00_rxns=[_epi00_rxns[x] for x in range(len(_epi00_rxns)) if x not in no_flux]
    _cor00_rxns=[_cor00_rxns[x] for x in range(len(_cor00_rxns)) if x not in no_flux]
    _cor01_rxns=[_cor01_rxns[x] for x in range(len(_cor01_rxns)) if x not in no_flux]
    _cor02_rxns=[_cor02_rxns[x] for x in range(len(_cor02_rxns)) if x not in no_flux]
    _cor03_rxns=[_cor03_rxns[x] for x in range(len(_cor03_rxns)) if x not in no_flux]
    _cor04_rxns=[_cor04_rxns[x] for x in range(len(_cor04_rxns)) if x not in no_flux]
    _cor05_rxns=[_cor05_rxns[x] for x in range(len(_cor05_rxns)) if x not in no_flux]
    _cor06_rxns=[_cor06_rxns[x] for x in range(len(_cor06_rxns)) if x not in no_flux]
    _cor07_rxns=[_cor07_rxns[x] for x in range(len(_cor07_rxns)) if x not in no_flux]
    _end00_rxns=[_end00_rxns[x] for x in range(len(_end00_rxns)) if x not in no_flux]
    _per00_rxns=[_per00_rxns[x] for x in range(len(_per00_rxns)) if x not in no_flux]

    data_epi00=[data_epi00[x] for x in range(len(data_epi00)) if x not in no_flux]
    data_cor00=[data_cor00[x] for x in range(len(data_cor00)) if x not in no_flux]
    data_cor01=[data_cor01[x] for x in range(len(data_cor01)) if x not in no_flux]
    data_cor02=[data_cor02[x] for x in range(len(data_cor02)) if x not in no_flux]
    data_cor03=[data_cor03[x] for x in range(len(data_cor03)) if x not in no_flux]
    data_cor04=[data_cor04[x] for x in range(len(data_cor04)) if x not in no_flux]
    data_cor05=[data_cor05[x] for x in range(len(data_cor05)) if x not in no_flux]
    data_cor06=[data_cor06[x] for x in range(len(data_cor06)) if x not in no_flux]
    data_cor07=[data_cor07[x] for x in range(len(data_cor07)) if x not in no_flux]
    data_end00=[data_end00[x] for x in range(len(data_end00)) if x not in no_flux]
    data_per00=[data_per00[x] for x in range(len(data_per00)) if x not in no_flux]

    #plotting figure
    labels=_epi00_rxns.copy()
    for x in range(len(labels)):
        #temp=labels[x].split('_')
        #labels[x]='_'.join(temp[:-2]+[temp[-1]])
        temp = labels[x]
        labels[x] = temp[:-6]

    width=0.08
    plt.figure(figsize=(9,20))
    plt.yticks([i+width*5 for i in range(len(data_epi00))], labels, rotation=0)

    plt.barh(np.arange(len(data_epi00)), data_epi00, width, label="epi00")
    plt.barh(np.arange(len(data_cor00))+ width, data_cor00, width, label="cor00")
    plt.barh(np.arange(len(data_cor01))+ 2*width, data_cor01, width, label="cor01")
    plt.barh(np.arange(len(data_cor02))+ 3*width, data_cor02, width, label="cor02")
    plt.barh(np.arange(len(data_cor03))+ 4*width, data_cor03, width, label="cor03")
    plt.barh(np.arange(len(data_cor04))+ 5*width, data_cor04, width, label="cor04")
    plt.barh(np.arange(len(data_cor05))+ 6*width, data_cor05, width, label="cor05")
    plt.barh(np.arange(len(data_cor06))+ 7*width, data_cor06, width, label="cor06")
    plt.barh(np.arange(len(data_cor07))+ 8*width, data_cor07, width, label="cor07")
    plt.barh(np.arange(len(data_end00))+ 9*width, data_end00, width, label="end00")
    plt.barh(np.arange(len(data_per00))+ 10*width, data_per00, width, label="per00")

    #adjusting figure
    plt.xlabel("flux in mmol/gDW")
    #reversing order of legend
    current_handles, current_labels = plt.gca().get_legend_handles_labels()
    reversed_handles = list(reversed(current_handles))
    reversed_labels = list(reversed(current_labels))
    plt.legend(reversed_handles,reversed_labels)
    #final bits
    plt.title(plttitle)
    plt.tick_params(axis='y', direction='inout', length=40)
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.savefig(f"fluxes_across_celltypes/{plttitle}_fluxesAcrossCellTypes", transparent=True)
    plt.show()



### function to plot all transfer reactions of a certain metabolites
#args: model, term to search e.g. "CIT"
def linker_plotter(model, met):
    import collections
    transf_rxns = collections.OrderedDict()
    for rxn in model.reactions:
        if "Transfer" in rxn.id:
            if met in rxn.id:
                transf_rxns[rxn.id] = rxn.x

    copy = transf_rxns.copy()
    for key in copy.keys():
        if "_epi00_" in key:
            transf_rxns.move_to_end(key, last=False)
        else:
            continue

    import matplotlib.pyplot as plt
    plt.barh(range(len(transf_rxns)), list(transf_rxns.values()), align='center')
    plt.yticks(range(len(transf_rxns)), list(transf_rxns.keys()), rotation = 0)
    plt.xlabel("flux in mmol/gDW")
    #plt.title(met + " linker reactions")


### function to scan reaction notes for a subsystem
#args: model, "term" to search in subsystem, model_tag e.g. "epi00"
#e.g. subsystem_search(tissue_model, "Calvin", "epi00")
# if no model_tag given it'll search all cells
def subsystem_search(model, subsystem, model_tag=""):
    for rxn in model.reactions:
        if model_tag == "":
            if subsystem in str(rxn.notes):
                print(rxn, "\n", rxn.notes, "\n")
        else:
            if rxn.id[-5:] == model_tag:
                if subsystem in str(rxn.notes):
                    print(rxn, "\n", rxn.notes, "\n")


### function to compare fluxes between cell types and return reactions with major differences
#args: base_model, tissue_model, flux difference cutoff e.g. 1
# use like so: flux_diffs = flux_comparer1(model, tissue_model, 1)
# then can do flux_diffs.to_excel("flux_differences.xlsx")
def flux_comparer(base_model, tissue_model, flux_cutoff):
    import pandas as pd
    from pandas import DataFrame
    list_tags = ["epi00","cor00","cor01","cor02","cor03","cor04","cor05","cor06","cor07","end00", "per00"]
    list_reactions = []
    base_reactions = base_model.reactions
    for rxn in base_reactions:
        for i in range(0,10):
            try:
                a = tissue_model.reactions.get_by_id(str(rxn.id)+"_"+list_tags[i])
                b = tissue_model.reactions.get_by_id(str(rxn.id)+"_"+list_tags[i+1])
                c = a.flux - b.flux
            except:
                continue
            if c > flux_cutoff or c < -flux_cutoff:
                dict1 = {}
                data1 = {"Reaction_id": a.id[:-5]}
                data2 = {"Reaction": a.reaction}
                for tag in list_tags:
                    z = tissue_model.reactions.get_by_id(str(rxn.id)+"_"+tag)
                    data_tagged = {"Flux_in_"+tag: z.flux}
                    dict1.update(data_tagged)
                dict1.update(data1)
                dict1.update(data2)
                list_reactions.append(dict1)

    flux_differences = pd.DataFrame(list_reactions, columns=["Reaction_id", "Reaction", "Flux_in_epi00", "Flux_in_cor00", "Flux_in_cor01", "Flux_in_cor02", "Flux_in_cor03",
                                               "Flux_in_cor04", "Flux_in_cor05", "Flux_in_cor06", "Flux_in_cor07", "Flux_in_end00", "Flux_in_per00"])
    flux_differences = flux_differences.drop_duplicates(subset=['Reaction_id'])

    return(flux_differences)



###Function to constraint sum of fluxes when performing FBA
#args: 1) a cobra model, 2) a python list of reactions to leave out from constrai-
#-nt, 3) the float value that sum of fluxes must be constrained to & 4) value obj-
#-ective function needs to be constraint to (provide "" to avoid constraining obj-
#ective function)
#output: a cobra model with sum of fluxes constrained to
def constrainSumOfFluxes(cobra_model, rxn2avoid,SFvalue,objvalue):
        temp=cobra_model.copy()
        SFMet = Metabolite("SFMet",name="Sum of fluxes pseudometabolite",compartment="c2")
        for rxn in cobra_model.reactions:
                if not rxn2avoid.__contains__(rxn.id):
                        if rxn.id.__contains__("reverse"):
                                temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:-1})
                else:
                        temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:1})
                        SFRxn = Reaction("SFRxn",name="Sum of fluxes pseudoreaction")
                        SFRxn.add_metabolites({SFMet:-1})
                        SFRxn.lower_bound=SFvalue
                        SFRxn.upper_bound=SFvalue
                        temp.add_reaction(SFRxn)
        if (not objvalue=="") and (len(temp.objective) == 1):
                for rxn in temp.objective.keys():
                        rxn.lower_bound=objvalue
                        rxn.upper_bound=objvalue
        return temp


##generate ATP budget
#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose whether to use percentage or absolute values in the plot. 6)#
# Provide a day or night indicator tag to specify day or night ATP    #
# summary 7) a destination file to save plot to 8) a dictionary to    #
# specify colour for fluxes in plot                                                                 #
#####################################################################
def generateATPbudget(model,solution,outfile="",show_plot=True,percentage=False,save_plot_to="temp.png",colourDict=None):
    if not colourDict:
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
    for model_tag in ["_epi00","_cor00","_cor01","_cor02","_cor03","_cor04","_cor05","_cor06","_cor07","_end00", "_per00"]:
        if outfile!="":
            fout = open(outfile,"w")
        ATPdict = dict()
        total = 0
        for p in ("c","p","m","x"):
            met=model.metabolites.get_by_id("ATP_"+p+model_tag)
            met1=model.metabolites.get_by_id("aATP_"+p+model_tag)
            for rxn in met.reactions:
                if rxn.id.__contains__("ATP_AMP_mc") or rxn.id.__contains__("ATP_ADP_mc") or rxn.id.__contains__("ATP_pc") or rxn.id.__contains__("AMP_ATP_xc") or rxn.id.__contains__("ATP_ADP_Pi_pc"):
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
                if ATPdict[rxn] < total*0.005:
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
                if abs(ATPdict[rxn]) < total*0.005:
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
        print(model_tag, " total ATP consumed: ", total)
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

def generateATPbudgetforonecell(model,solution,tag = "_epi00", outfile="",show_plot=False,percentage=False,save_plot_to=""):

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
    for p in ("c","p","m","x"):
        met=model.metabolites.get_by_id("ATP_"+p+model_tag)
        met1=model.metabolites.get_by_id("aATP_"+p+model_tag)
        for rxn in met.reactions:
            if rxn.id.__contains__("ATP_AMP_mc") or rxn.id.__contains__("ATP_ADP_mc") or rxn.id.__contains__("ATP_pc") or rxn.id.__contains__("AMP_ATP_xc") or rxn.id.__contains__("ATP_ADP_Pi_pc"):
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
                if ATPdict[rxn] < total*0.005:
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
                if abs(ATPdict[rxn]) < total*0.005:
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
    print(model_tag, " total ATP consumed: ", total)
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



### generate budget organised by subsystem rather than by reaction
def generateATPbudgetsubsystem(model,solution,outfile="",show_plot=True,percentage=False,save_plot_to="temp.png",colourDict={}):

    #find colour hex codes here: https://www.colorhexa.com/color-names
    nad_coloursdict = {
        'MALATE_DEH_RXN_m_': '#5C6BC0',
        '2OXOGLUTARATEDEH_RXN_m_':'#8D6E63' ,
        'ISOCITRATE_DEHYDROGENASE_NAD_RXN_m_': '#AB47BC',
        'PYRUVDEH_RXN_m_': '#2E7D32',
        'GAPOXNPHOSPHN_RXN_c_': '#EF5350',
        'MALATE_DEH_RXN_c_': '#FF7043',
        'NADH_DEHYDROG_A_RXN_mi_': '#26C6DA',
        'GLC': '#9CCC65',
        'Other': '#BDBDBD',
        'ISOCITDEH_RXN_m_': '#e9d66b',
        'MALATE_DEH_RXN_p_': '#fdee00',
        '6PGLUCONDEHYDROG_RXN_p_': '#f4c2c2',}

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

    total_atp_dict = {}
    for model_tag in ["_epi00","_cor00","_cor01","_cor02","_cor03","_cor04","_cor05","_cor06","_cor07","_end00", "_per00"]:
        if outfile!="":
            fout = open(outfile,"w")
        ATPdict = dict()
        total = 0
        for p in ("c","p","m","x"):
            met=model.metabolites.get_by_id("ATP_"+p+model_tag)
            met1=model.metabolites.get_by_id("aATP_"+p+model_tag)
            for rxn in met.reactions:
                if rxn.id.__contains__("ATP_AMP_mc") or rxn.id.__contains__("ATP_ADP_mc") or rxn.id.__contains__("ATP_pc") or rxn.id.__contains__("AMP_ATP_xc") or rxn.id.__contains__("ATP_ADP_Pi_pc"):
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
                if ATPdict[rxn] < total*0.005:
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
                if abs(ATPdict[rxn]) < total*0.05:
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
        print(model_tag, " total ATP consumed: ", total)
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

    return total_atp_dict


##generate NADP(H) budget
def generateNADHNADPHbudget(model,solution,outfile="",show_plot=True,percentage=False,save_plot_to="temp",colourDict=None):
    if not colourDict:
        #find colour hex codes here: https://www.colorhexa.com/color-names
        nad_coloursdict = {
            'MALATE_DEH_RXN_m_': '#5C6BC0',
            '2OXOGLUTARATEDEH_RXN_m_':'#8D6E63' ,
            'ISOCITRATE_DEHYDROGENASE_NAD_RXN_m_': '#AB47BC',
            'PYRUVDEH_RXN_m_': '#2E7D32',
            'GAPOXNPHOSPHN_RXN_c_': '#EF5350',
            'MALATE_DEH_RXN_c_': '#FF7043',
            'NADH_DEHYDROG_A_RXN_mi_': '#26C6DA',
            'GLC': '#9CCC65',
            'Other': '#BDBDBD',
            'ISOCITDEH_RXN_m_': '#e9d66b',
            'MALATE_DEH_RXN_p_': '#fdee00',
            '6PGLUCONDEHYDROG_RXN_p_': '#f4c2c2',}
        colourDict = nad_coloursdict

    total_nad_dict = {}
    for model_tag in ["_epi00","_cor00","_cor01","_cor02","_cor03","_cor04","_cor05","_cor06","_cor07","_end00", "_per00"]:
        if outfile!="":
                fout = open(outfile,"w")
        Reddict = dict()
        total = 0
        for red in ["NADPH","NADH"]:
                for p in ("c","p","m","x"):
                        if len(model.metabolites.query(red+"_"+p+model_tag))==0:
                                continue
                        met=model.metabolites.get_by_id(red+"_"+p+model_tag)
                        for rxn in met.reactions:
                                sto=rxn.metabolites.get(met)
                                sto1=0#rxn.metabolites.get(met1)
                                if outfile!="":
                                        fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
                                Reddict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
                                if solution.get(rxn.id)*(sto+sto1) > 0:
                                        total = total + (solution.get(rxn.id)*(sto+sto1))
        if outfile!="":
                fout.close()

        tempDict = dict()
        for rxn in Reddict.keys():
            tempDict[rxn]=abs(Reddict[rxn])

        #sort by values
        import operator
        sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)


        Reddict2 = dict()
        Reddict2["Others-pos"]=0
        Reddict2["Others-neg"]=0
        baseline = dict()
        pos_base=0
        neg_base=0
        i=0
        for TEMP in sorted_by_value:
                rxn = TEMP[0]
                if Reddict[rxn]>0:
                        if Reddict[rxn] < total*0.005:
                                if percentage:
                                        Reddict2["Others-pos"]=Reddict2["Others-pos"]+float(Reddict[rxn]*100)/total
                                else:
                                        Reddict2["Others-pos"]=Reddict2["Others-pos"]+Reddict[rxn]
                                continue
                        base = pos_base
                        if percentage:
                                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                                pos_base = pos_base + float(Reddict[rxn]*100)/total
                        else:
                                pos_base = pos_base + Reddict[rxn]
                                Reddict2[rxn]=Reddict[rxn]
                else:
                        if abs(Reddict[rxn]) < total*0.005:
                                if percentage:
                                        Reddict2["Others-neg"]=Reddict2["Others-neg"]+float(Reddict[rxn]*100)/total
                                else:
                                        Reddict2["Others-neg"]=Reddict2["Others-neg"]+Reddict[rxn]
                                continue
                        base = neg_base
                        if percentage:
                                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                                neg_base = neg_base + float(Reddict[rxn]*100)/total
                        else:
                                neg_base = neg_base + Reddict[rxn]
                                Reddict2[rxn]=Reddict[rxn]
                i=i+1
                baseline[rxn]=base
        baseline["Others-pos"]=pos_base
        baseline["Others-neg"]=neg_base

        #obtaining total nad expenditure
        total = 0
        for value in Reddict2.values():
            if value > 0:
                total += value
        print(model_tag, " total NAD(P)H consumed: ", total)
        total_nad_dict[model_tag] = total

        #making plots
        if show_plot:
                import matplotlib.pyplot as plt
                plt.rcParams.update({'font.size': 10}) #sets a global fontsize
                plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
                plt.rcParams['xtick.major.width'] = 1
                plt.rcParams['ytick.major.size'] = 5
                plt.rcParams['ytick.major.width'] = 1
                plt.rcParams['axes.linewidth']=2 # makes axes line thicker
                plt.figure(figsize=(3,4))
                for rxn in Reddict2.keys():
                    if colourDict.keys().__contains__(rxn[:-5]):
                        plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn[:-5]])
                    else:
                        plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
                plt.xlim(0.8,1.2)
                if percentage:
                        plt.ylabel("NAD(P)H produced/consumed (%)")
                else:
                        plt.ylabel("NAD(P)H produced/consumed (in moles)")
                handles, labels = plt.gca().get_legend_handles_labels()
                labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
                handles2=[handles[labels.index(i)] for i in labels2]
                lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
                plt.axhline(0,linestyle="--",color="black")
                plt.tight_layout
                plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')

    return total_nad_dict, Reddict2


    total_nad_dict = {}
    model_tag = tag
    Reddict = dict()
    total = 0
    for red in ["NADPH","NADH"]:
            for p in ("c","p","m","x"):
                    if len(model.metabolites.query(red+"_"+p+model_tag))==0:
                            continue
                    met=model.metabolites.get_by_id(red+"_"+p+model_tag)
                    for rxn in met.reactions:
                            sto=rxn.metabolites.get(met)
                            sto1=0#rxn.metabolites.get(met1)
                            Reddict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
                            if solution.get(rxn.id)*(sto+sto1) > 0:
                                    total = total + (solution.get(rxn.id)*(sto+sto1))

    tempDict = dict()
    for rxn in Reddict.keys():
        tempDict[rxn]=abs(Reddict[rxn])

    #sort by values
    import operator
    sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)


    Reddict2 = dict()
    Reddict2["Others-pos"]=0
    Reddict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    for TEMP in sorted_by_value:
            rxn = TEMP[0]
            if Reddict[rxn]>0:
                    if Reddict[rxn] < total*0.005:
                            if percentage:
                                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+float(Reddict[rxn]*100)/total
                            else:
                                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+Reddict[rxn]
                            continue
                    base = pos_base
                    if percentage:
                            Reddict2[rxn]=float(Reddict[rxn]*100)/total
                            pos_base = pos_base + float(Reddict[rxn]*100)/total
                    else:
                            pos_base = pos_base + Reddict[rxn]
                            Reddict2[rxn]=Reddict[rxn]
            else:
                    if abs(Reddict[rxn]) < total*0.005:
                            if percentage:
                                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+float(Reddict[rxn]*100)/total
                            else:
                                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+Reddict[rxn]
                            continue
                    base = neg_base
                    if percentage:
                            Reddict2[rxn]=float(Reddict[rxn]*100)/total
                            neg_base = neg_base + float(Reddict[rxn]*100)/total
                    else:
                            neg_base = neg_base + Reddict[rxn]
                            Reddict2[rxn]=Reddict[rxn]
            i=i+1
            baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    #obtaining total nad expenditure
    total = 0
    for value in Reddict2.values():
        if value > 0:
            total += value
    print(model_tag, " total NAD(P)H consumed: ", total)
    total_nad_dict[model_tag] = total

    return total_nad_dict, Reddict2


def generateNADHNADPHbudgetforonecell(model,solution,tag = "_epi00", outfile="",show_plot=False,percentage=False,save_plot_to=""):

    #find colour hex codes here: https://www.colorhexa.com/color-names
    nad_coloursdict = {
        'MALATE_DEH_RXN_m_': '#5C6BC0',
        '2OXOGLUTARATEDEH_RXN_m_':'#8D6E63' ,
        'ISOCITRATE_DEHYDROGENASE_NAD_RXN_m_': '#AB47BC',
        'PYRUVDEH_RXN_m_': '#2E7D32',
        'GAPOXNPHOSPHN_RXN_c_': '#EF5350',
        'MALATE_DEH_RXN_c_': '#FF7043',
        'NADH_DEHYDROG_A_RXN_mi_': '#26C6DA',
        'GLC': '#9CCC65',
        'Other': '#BDBDBD',
        'ISOCITDEH_RXN_m_': '#e9d66b',
        'MALATE_DEH_RXN_p_': '#fdee00',
        '6PGLUCONDEHYDROG_RXN_p_': '#f4c2c2',}
    colourDict = nad_coloursdict

    total_nad_dict = {}
    model_tag = tag
    if outfile!="":
            fout = open(outfile,"w")
    Reddict = dict()
    total = 0
    for red in ["NADPH","NADH"]:
            for p in ("c","p","m","x"):
                    if len(model.metabolites.query(red+"_"+p+model_tag))==0:
                            continue
                    met=model.metabolites.get_by_id(red+"_"+p+model_tag)
                    for rxn in met.reactions:
                            sto=rxn.metabolites.get(met)
                            sto1=0#rxn.metabolites.get(met1)
                            if outfile!="":
                                    fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
                            Reddict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
                            if solution.get(rxn.id)*(sto+sto1) > 0:
                                    total = total + (solution.get(rxn.id)*(sto+sto1))
    if outfile!="":
            fout.close()

    tempDict = dict()
    for rxn in Reddict.keys():
        tempDict[rxn]=abs(Reddict[rxn])

    #sort by values
    import operator
    sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)


    Reddict2 = dict()
    Reddict2["Others-pos"]=0
    Reddict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    for TEMP in sorted_by_value:
            rxn = TEMP[0]
            if Reddict[rxn]>0:
                    if Reddict[rxn] < total*0.005:
                            if percentage:
                                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+float(Reddict[rxn]*100)/total
                            else:
                                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+Reddict[rxn]
                            continue
                    base = pos_base
                    if percentage:
                            Reddict2[rxn]=float(Reddict[rxn]*100)/total
                            pos_base = pos_base + float(Reddict[rxn]*100)/total
                    else:
                            pos_base = pos_base + Reddict[rxn]
                            Reddict2[rxn]=Reddict[rxn]
            else:
                    if abs(Reddict[rxn]) < total*0.005:
                            if percentage:
                                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+float(Reddict[rxn]*100)/total
                            else:
                                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+Reddict[rxn]
                            continue
                    base = neg_base
                    if percentage:
                            Reddict2[rxn]=float(Reddict[rxn]*100)/total
                            neg_base = neg_base + float(Reddict[rxn]*100)/total
                    else:
                            neg_base = neg_base + Reddict[rxn]
                            Reddict2[rxn]=Reddict[rxn]
            i=i+1
            baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    #obtaining total nad expenditure
    total = 0
    for value in Reddict2.values():
        if value > 0:
            total += value
    print(model_tag, " total NAD(P)H consumed: ", total)
    total_nad_dict[model_tag] = total

    #making plots
    if show_plot:
            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 10}) #sets a global fontsize
            plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
            plt.rcParams['xtick.major.width'] = 1
            plt.rcParams['ytick.major.size'] = 5
            plt.rcParams['ytick.major.width'] = 1
            plt.rcParams['axes.linewidth']=2 # makes axes line thicker
            plt.figure(figsize=(3,4))
            for rxn in Reddict2.keys():
                if colourDict.keys().__contains__(rxn[:-5]):
                    plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn[:-5]])
                else:
                    plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
            plt.xlim(0.8,1.2)
            if percentage:
                    plt.ylabel("NAD(P)H produced/consumed (%)")
            else:
                    plt.ylabel("NAD(P)H produced/consumed (in moles)")
            handles, labels = plt.gca().get_legend_handles_labels()
            labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
            handles2=[handles[labels.index(i)] for i in labels2]
            lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
            plt.axhline(0,linestyle="--",color="black")
            plt.tight_layout
            plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')

    return total_nad_dict, Reddict2



def generate_atp_process_chart1(model, solution1, cell_tag):
    reaction_atp_process = {
    'Mitochondrial_ATP_Synthase_m_': 'Mitochondrial ATP Synthase',
    'PHOSGLYPHOS_RXN_c_': 'Cytosolic glycolysis' ,
    'PEPDEPHOS_RXN_p_': 'Plastidic glycolysis',
    'PHOSPHORIBULOSEKINASE_RXN_p_': 'Calvin cycle',
    'PEPDEPHOS_RXN_c_': 'Cytosolic glycolysis',
    'PHOSGLYPHOS_RXN_p_': 'Plastidic glycolysis',
    'ACETYL_COA_CARBOXYLTRANSFER_RXN_p_': 'Lipid synthesis',
    'UDPKIN_RXN_c_': 'Cellulose Synthesis',
    'GDPKIN_RXN_c_': 'GDP kinase',
    'FRUCTOKINASE_RXN_c_': 'Cytosolic glycolysis',
    'PROTON_ATPase_c_': 'Plasma membrane proton ATPase',
    'PEPCARBOXYKIN_RXN_c_': 'PEP carboxykinase',}

    atp_coloursdict = {
    'Mitochondrial ATP Synthase': '#5C6BC0',
    'Cytosolic glycolysis':'#8D6E63' ,
    'Plastidic glycolysis': '#AB47BC',
    'Calvin cycle': '#2E7D32',
    'Lipid synthesis': '#FF7043',
    'Cellulose Synthesis': '#26C6DA',
    'GDP kinase': '#9CCC65',
    'Plasma membrane proton ATPase': '#fbceb1',
    'PEP carboxykinase': '#cd9575',
    'Others-pos': '#000000',}


    total_atp_dict, atp_rxns_fluxes = generateATPbudgetforonecell(model,solution1.fluxes,cell_tag, percentage=False)

    process_fluxes = {}
    for v in reaction_atp_process.values():
        process_fluxes[v] = 0

    for k in atp_rxns_fluxes.keys():
        for rxn_id in reaction_atp_process.keys():
            if rxn_id in k:
                process = reaction_atp_process[rxn_id]
                rxn_flux = atp_rxns_fluxes[k]
                process_fluxes[process] = process_fluxes[process] + rxn_flux

    pos_processes = process_fluxes.copy()
    to_del_pos = []
    for k in pos_processes.keys():
        if pos_processes[k] <= 0:
            to_del_pos.append(k)

    for k in to_del_pos:
        del pos_processes[k]

    neg_processes = process_fluxes.copy()
    to_del_neg = []
    for k in neg_processes.keys():
        if neg_processes[k] >= 0:
            to_del_neg.append(k)
    for k in to_del_neg:
        del neg_processes[k]

    pos_keys = [*pos_processes.keys()]
    pos_values = [*pos_processes.values()]
    neg_keys = [*neg_processes.keys()]
    neg_values_old = [*neg_processes.values()]
    neg_values = []
    for value in neg_values_old:
        neg_values.append(value*-1)

    import matplotlib.pyplot as plt
    plt.style.reload_library()
    with plt.style.context('style1'):
        fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2)

        fig.set_size_inches(14,4)
        my_colors = ["r","g","b","k","y","m","c"]  #red, green, blue, black, etc.

        for i in range(0,len(pos_keys)):
            ax1.bar(pos_keys[i], pos_values[i], color=atp_coloursdict[pos_keys[i]])
        ax1.tick_params("x", labelrotation=80)
        ax1.set_ylim(0, 1.8)
        ax1.set_title("ATP generating processes in "+cell_tag)
        ax1.set_ylabel("ATP in mmol")

        for i in range(0,len(neg_keys)):
            ax2.bar(neg_keys[i], neg_values[i], color=atp_coloursdict[neg_keys[i]])
        ax2.tick_params("x", labelrotation=80)
        ax2.set_ylim(0, 1.8)
        ax2.set_title("ATP consuming processes in "+cell_tag)
        ax2.set_ylabel("ATP in mmol")
        #fig.legend()

        plt.savefig("ATPorNAD production and usage/ATP/"+"ATP_"+cell_tag+"production+usage", bbox_inches='tight')



def generate_nad_process_chart(model, solution1, cell_tag):
    reaction_nad_process = {
        'MALATE_DEH_RXN_m_': 'TCA cycle',
        '2OXOGLUTARATEDEH_RXN_m_':'TCA cycle' ,
        'ISOCITRATE_DEHYDROGENASE_NAD_RXN_m_': 'TCA cycle',
        'PYRUVDEH_RXN_m_': 'TCA cycle',
        'GAPOXNPHOSPHN_RXN_c_': 'Cytosolic glycolysis',
        'GAPOXNPHOSPHN_RXN_p_': 'Plastidic glycolysis',
        'MALATE_DEH_RXN_c_': 'Cytosolic malate dehydrogenase',
        'MALATE_DEH_RXN_p_': 'Plastidic malate dehydrogenase',
        'NADH_DEHYDROG_A_RXN_mi_': 'Electron transport chain',
        'ISOCITDEH_RXN_m_': 'TCA cycle',
        'GLU6PDEHYDROG_RXN_p_': 'Plastidic OPPP',
        '6PGLUCONDEHYDROG_RXN_p_': 'Plastidic OPPP',
        'RXN_9524_p_': 'Lipid synthesis',
        'RXN_9532_p_': 'Lipid synthesis',
        'RXN_9536_p_': 'Lipid synthesis',
        'RXN_9540_p_': 'Lipid synthesis',
        'RXN_9528_p_': 'Lipid synthesis',
        'RXN_9514_p_': 'Lipid synthesis',
        'RXN_9518_p_': 'Lipid synthesis',
        'RXN_9660_p_': 'Lipid synthesis',
        'RXN_9657_p_': 'Lipid synthesis',
        'RXN_9661_p_': 'Lipid synthesis',
        'RXN_9659_p_': 'Lipid synthesis',
        'RXN_9658_p_': 'Lipid synthesis',
        'RXN_9662_p_': 'Lipid synthesis',
        'RXN_9663_p_': 'Lipid synthesis',
        'PYRUVDEH_RXN_p_': 'Plastidic pyruvate dehydrogenase',
        'L_LACTATEDEHYDROG_RXN_c_': 'Fermentation',
        'GLUTAMATE_DEHYDROGENASE_RXN_m_': 'Glutamate dehydrogenase',
        'PYRROLINECARBDEHYDROG_RXN_NADP_m_': '1-Pyrroline-5-carboxylate dehydrogenase',
        'MALIC_NADP_RXN_p_': 'Malate dehydrogenase (oxaloacetate decarboxylating)',}

    nad_coloursdict = {
        'TCA cycle': '#5C6BC0',
        'Cytosolic glycolysis':'#8D6E63' ,
        'Plastidic glycolysis': '#AB47BC',
        'Cytosolic malate dehydrogenase': '#2E7D32',
        'Plastidic malate dehydrogenase': '#EF5350',
        'Lipid synthesis': '#FF7043',
        'Plastidic pyruvate dehydrogenase': '#26C6DA',
        'Fermentation': '#9CCC65',
        'Glutamate dehydrogenase': '#BDBDBD',
        '1-Pyrroline-5-carboxylate dehydrogenase': '#fbceb1',
        'Malate dehydrogenase (oxaloacetate decarboxylating)': '#cd9575',
        'Plastidic OPPP': '#808080',
        'Electron transport chain': '#000000',}


    total_nad_dict, nad_rxns_fluxes = generateNADHNADPHbudgetforonecell(model,solution1.fluxes,tag=cell_tag,
    percentage=False)

    process_fluxes = {}
    for v in reaction_nad_process.values():
        process_fluxes[v] = 0

    for k in nad_rxns_fluxes.keys():
        for rxn_id in reaction_nad_process.keys():
            if rxn_id in k:
                process = reaction_nad_process[rxn_id]
                rxn_flux = nad_rxns_fluxes[k]
                process_fluxes[process] = process_fluxes[process] + rxn_flux

    pos_processes = process_fluxes.copy()
    to_del_pos = []
    for k in pos_processes.keys():
        if pos_processes[k] <= 0:
            to_del_pos.append(k)

    for k in to_del_pos:
        del pos_processes[k]

    neg_processes = process_fluxes.copy()
    to_del_neg = []
    for k in neg_processes.keys():
        if neg_processes[k] >= 0:
            to_del_neg.append(k)
    for k in to_del_neg:
        del neg_processes[k]

    pos_keys = [*pos_processes.keys()]
    pos_values = [*pos_processes.values()]
    neg_keys = [*neg_processes.keys()]
    neg_values_old = neg_processes.values()
    neg_values = []
    for value in neg_values_old:
        neg_values.append(value*-1)

    import matplotlib.pyplot as plt
    with plt.style.context('style1'):
        fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2)

        fig.set_size_inches(14,6)
        my_colors = ["r","g","b","k","y","m","c"]  #red, green, blue, black, etc.

        for i in range(0,len(pos_keys)):
            ax1.bar(pos_keys[i], pos_values[i], color=nad_coloursdict[pos_keys[i]])
        ax1.tick_params("x", labelrotation=80)
        ax1.set_ylim(0, 0.5)
        ax1.set_title("NAD(P)H generating processes in "+cell_tag)
        ax1.set_ylabel("NAD(P)H in mmol")

        for i in range(0,len(neg_keys)):
            ax2.bar(neg_keys[i], neg_values[i], color=nad_coloursdict[neg_keys[i]])
        ax2.tick_params("x", labelrotation=80)
        ax2.set_ylim(0, 0.5)
        ax2.set_title("NAD(P)H consuming processes in "+cell_tag)
        ax2.set_ylabel("NAD(P)H in mmol")

        plt.savefig("ATPorNAD production and usage/NAD(P)H/"+"NAD(P)H_"+cell_tag+"production+usage", bbox_inches='tight')


###
def writeSolutionFluxesToFile(sol,outfile,model):
        import pandas as pd
        rxnList = list()
        eqnList = list()
        fluxList = list()
        EClist = list()
        for rxn in model.reactions:
                rxnList.append(rxn.id)
                eqnList.append(rxn.reaction)
                if not "PROTEIN CLASS" in rxn.notes.keys():
                        EClist.append("")
                else:
                        EClist.append(rxn.notes.get("PROTEIN CLASS")[0])
                fluxList.append(sol.fluxes.get(rxn.id))
        df = pd.DataFrame(data={"ID":rxnList,"EC number":EClist,"reaction":eqnList,"flux":fluxList})
        df = df[['ID','EC number', 'reaction', 'flux']]
        df.to_csv(outfile)
        return


### modified version of pFBA to allow linker reactions to be excluded b/ they don't have an enzymatic cost
import logging
from warnings import warn
from itertools import chain

from optlang.symbolics import Zero

from cobra.util import solver as sutil
from cobra.core.solution import get_solution

def pfba_Weighted(model, weightings, fraction_of_optimum=1.0, objective=None, reactions=None):
    """Perform basic pFBA (parsimonious Enzyme Usage Flux Balance Analysis)
    to minimize total flux.
    pFBA [1] adds the minimization of all fluxes the the objective of the
    model. This approach is motivated by the idea that high fluxes have a
    higher enzyme turn-over and that since producing enzymes is costly,
    the cell will try to minimize overall flux while still maximizing the
    original objective function, e.g. the growth rate.
    Parameters
    ----------
    model : cobra.Model
        The model
    fraction_of_optimum : float, optional
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    objective : dict or model.problem.Objective
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives.
    reactions : iterable
        List of reactions or reaction identifiers. Implies `return_frame` to
        be true. Only return fluxes for the given reactions. Faster than
        fetching all fluxes if only a few are needed.
    Returns
    -------
    cobra.Solution
        The solution object to the optimized model with pFBA constraints added.
    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
       Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
       from evolved E. coli are consistent with computed optimal growth from
       genome-scale models. Molecular Systems Biology, 6,
       390. doi:10.1038/msb.2010.47
    """
    reactions = model.reactions if reactions is None \
        else model.reactions.get_by_any(reactions)
    with model as m:
        add_pfba_Weighted(m, weightings, objective=objective,
                 fraction_of_optimum=fraction_of_optimum)
        m.slim_optimize(error_value=None)
        solution = get_solution(m, reactions=reactions)
    return solution


#################################################################################
# This function is a modified version of cobrapy add_pfba function			#
#										#
#################################################################################

def add_pfba_Weighted(model, weightings, objective=None, fraction_of_optimum=1.0):
    """Add pFBA objective
    Add objective to minimize the summed flux of all reactions to the
    current objective.
    See Also
    -------
    pfba
    Parameters
    ----------
    model : cobra.Model
        The model to add the objective to
    objective :
        An objective to set in combination with the pFBA objective.
    fraction_of_optimum : float
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    """
    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction='min', sloppy=True, name="_pfba_objective")
    #print([v for v in variables])
    tempDict = dict()
    for v in variables:
        w = str(v).split("=")[1].replace(" ","").replace("<","")
        found=False
        for rxn in weightings.keys():
            if w.__contains__(rxn):
                #print(v)
                #print(rxn)
                tempDict[v]=weightings[rxn]
                found=True
                break
        if not found:
            #print("Weightings for reaction "+w+" not found, so assuming weighting = 1")
            tempDict[v] = 1
    model.objective.set_linear_coefficients(tempDict)




###attempt to repurpose subsystem plotter function to plot linker fluxes. Doesn't currently work
def plot_linkers_btwn_cells(cons_model, csol, thresh):
   # [x.id for x in cons_model.reactions if any([y in str(x.notes)for y in ss])]

    #all_subsys_rxns=[x.id for x in cons_model.reactions if "Transfer" in x.id and ss in x.id]

    all_subsys_rxns = []
    for rxn in cons_model.reactions:
        if "Transfer" in rxn.id:
            if rxn.x < -thresh or rxn.x > thresh:
                #print(rxn)
                all_subsys_rxns.append(rxn.id)


    #all_subsys_rxns=[x.id for x in cons_model.reactions if any([y in str(x.id)for y in ss and "Transfer" in x.id])]
    _epi00_rxns=[x for x in all_subsys_rxns if 'c_epi00' in x]
    _cor00_rxns=[x for x in all_subsys_rxns if 'c_cor00' in x]
    _cor01_rxns=[x for x in all_subsys_rxns if 'c_cor01' in x]
    _cor02_rxns=[x for x in all_subsys_rxns if 'c_cor02' in x]
    _cor03_rxns=[x for x in all_subsys_rxns if 'c_cor03' in x]
    _cor04_rxns=[x for x in all_subsys_rxns if 'c_cor04' in x]
    _cor05_rxns=[x for x in all_subsys_rxns if 'c_cor05' in x]
    _cor06_rxns=[x for x in all_subsys_rxns if 'c_cor06' in x]
    _cor07_rxns=[x for x in all_subsys_rxns if 'c_cor07' in x]
    _end00_rxns=[x for x in all_subsys_rxns if 'c_end00' in x]
    print(_cor00_rxns)
    print(_cor07_rxns)

    data_epi00=[csol[x] if x else 0 for x in _epi00_rxns]
    data_cor00=[csol[x] if x else 0 for x in _cor00_rxns]
    data_cor01=[csol[x] if x else 0 for x in _cor01_rxns]
    data_cor02=[csol[x] if x else 0 for x in _cor02_rxns]
    data_cor03=[csol[x] if x else 0 for x in _cor03_rxns]
    data_cor04=[csol[x] if x else 0 for x in _cor04_rxns]
    data_cor05=[csol[x] if x else 0 for x in _cor05_rxns]
    data_cor06=[csol[x] if x else 0 for x in _cor06_rxns]
    data_cor07=[csol[x] if x else 0 for x in _cor07_rxns]
    data_end00=[csol[x] if x else 0 for x in _end00_rxns]
    print(data_cor00)
    print(data_cor07)

    """no_flux=[]
    for ii in range(len(data_epi00)):
        if abs(data_epi00[ii])<=thresh and abs(data_cor00[ii])<=thresh and abs(data_cor01[ii])<=thresh and abs(data_cor02[ii])<=thresh and abs(data_cor03[ii])<=thresh and abs(data_cor04[ii])<=thresh and abs(data_cor05[ii])<=thresh and abs(data_cor06[ii])<=thresh and abs(data_cor07[ii])<=thresh and abs(data_end00[ii])<=thresh:
            no_flux.append(ii)"""

    """_epi00_rxns=[_epi00_rxns[x] for x in range(len(_epi00_rxns))]
    _cor00_rxns=[_cor00_rxns[x] for x in range(len(_cor00_rxns))]
    _cor01_rxns=[_cor01_rxns[x] for x in range(len(_cor01_rxns))]
    _cor02_rxns=[_cor02_rxns[x] for x in range(len(_cor02_rxns))]
    _cor03_rxns=[_cor03_rxns[x] for x in range(len(_cor03_rxns))]
    _cor04_rxns=[_cor04_rxns[x] for x in range(len(_cor04_rxns))]
    _cor05_rxns=[_cor05_rxns[x] for x in range(len(_cor05_rxns))]
    _cor06_rxns=[_cor06_rxns[x] for x in range(len(_cor06_rxns))]
    _cor07_rxns=[_cor07_rxns[x] for x in range(len(_cor07_rxns))]
    _end00_rxns=[_end00_rxns[x] for x in range(len(_end00_rxns))]

    data_epi00=[data_epi00[x] for x in range(len(data_epi00))]
    data_cor00=[data_cor00[x] for x in range(len(data_cor00))]
    data_cor01=[data_cor01[x] for x in range(len(data_cor01))]
    data_cor02=[data_cor02[x] for x in range(len(data_cor02))]
    data_cor03=[data_cor03[x] for x in range(len(data_cor03))]
    data_cor04=[data_cor04[x] for x in range(len(data_cor04))]
    data_cor05=[data_cor05[x] for x in range(len(data_cor05))]
    data_cor06=[data_cor06[x] for x in range(len(data_cor06))]
    data_cor07=[data_cor07[x] for x in range(len(data_cor07))]
    data_end00=[data_end00[x] for x in range(len(data_end00))]"""

    #plotting figure
    labels=all_subsys_rxns.copy()
    for x in range(len(labels)):
        #temp=labels[x].split('_')
        #labels[x]='_'.join(temp[:-2]+[temp[-1]])
        temp = labels[x]
        labels[x] = temp[:-1]

    width=0.08
    plt.figure(figsize=(9,20))
    plt.yticks([i+width*5 for i in range(len(data_epi00))], labels, rotation=0)

    plt.barh(np.arange(len(data_epi00)), data_epi00, width, label="epi00")
    plt.barh(np.arange(len(data_cor00))+ width, data_cor00, width, label="cor00")
    plt.barh(np.arange(len(data_cor01))+ 2*width, data_cor01, width, label="cor01")
    plt.barh(np.arange(len(data_cor02))+ 3*width, data_cor02, width, label="cor02")
    plt.barh(np.arange(len(data_cor03))+ 4*width, data_cor03, width, label="cor03")
    plt.barh(np.arange(len(data_cor04))+ 5*width, data_cor04, width, label="cor04")
    plt.barh(np.arange(len(data_cor05))+ 6*width, data_cor05, width, label="cor05")
    plt.barh(np.arange(len(data_cor06))+ 7*width, data_cor06, width, label="cor06")
    plt.barh(np.arange(len(data_cor07))+ 8*width, data_cor07, width, label="cor07")
    plt.barh(np.arange(len(data_end00))+ 9*width, data_end00, width, label="end00")


    #adjusting figure
    plt.xlabel("flux in mmol/gDW")
    #reversing order of legend
    current_handles, current_labels = plt.gca().get_legend_handles_labels()
    reversed_handles = list(reversed(current_handles))
    reversed_labels = list(reversed(current_labels))
    plt.legend(reversed_handles,reversed_labels)
    #final bits
    #plt.title(plttitle)
    plt.tick_params(axis='y', direction='inout', length=40)
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    #plt.savefig(f"fluxes_across_celltypes/{plttitle}_fluxesAcrossCellTypes", transparent=True)
    plt.show()
