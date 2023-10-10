import pandas as pd
import cobra
from cobra import Reaction, Metabolite

###shift header of dataframe down one row
def header_shifter(excelfile):
    rootslice = pd.read_excel("Data from Lynch lab/cell_dimensions.xlsx")
    header = rootslice.iloc[0]
    rootslice = rootslice[1:]
    rootslice.columns = header
    rootslice

###remove boundary reactions e.g. CO2_tx
#args: model
def remove_boundary_rxns(model):
    reactions_to_remove = []
    for rxn in model.reactions:
        if rxn.boundary == True:
            reactions_to_remove.append(rxn)
    model.remove_reactions(reactions_to_remove)


###create multi-cell model
#args: model, number of cells as integer, cell type e.g. "cortical"
#outputs: dictionary of cells
def set_up_multicell_model(model, number_cells, cell_type):
    # making setting up the separate time phases
    cell_number_dict= {}
    cell_number_list = []
    models = 0
    for i in range(0, number_cells):
        comp_dict = {}
        temp_model = model.copy()
        #generate new compartments for every cell
        for item in temp_model.compartments:
            comp_dict[item+"_"+(cell_type[:3] + str(i).zfill(2))] = temp_model.compartments[item]+"_"+(cell_type[:3] + str(i).zfill(2))
        temp_model.compartments = comp_dict
        temp_model.compartments.update()
        #generate new metabolites for every cell
        for met in temp_model.metabolites:
            met.id = str(met.id) + "_"+(cell_type[:3] + str(i).zfill(2))
            met.compartment = str(met.compartment) + "_"+(cell_type[:3] + str(i).zfill(2))
        # the part of updating the name of the compartment of the metabolite is essential to actually update
        # the name compartment dictionary in the model, and the names of the compartments set above have to be the same
        # as the names set here
        for reaction in temp_model.reactions:
            reaction.id = str(reaction.id) + "_"+(cell_type[:3] + str(i).zfill(2))

        exec(f'cell_number_model_{str(i).zfill(2)} = temp_model')
        cell_number_dict[i] = eval("cell_number_model_%s" % (str(i).zfill(2)))
        print("Model created for cell: "+ cell_type[:3] + str(i))
        models +=1

    # combining the models into to a single model
    multicell_model = cell_number_dict[0]
    """for model_number in cell_number_dict:
        if model_number != 0:
            multicell_model += cell_number_dict[model_number]"""
    # using this method, metabolites that are not part of any reactions will not exist in all time phases except for the first

    return multicell_model, cell_number_dict


###generating linker reactions between cells
#args: multicell model, total_cells as integer e.g. 9, cor_tag e.g. "cor"
def create_linker_fluxes(multicell_model, total_cells, core_tag):
    no_symp_transport = ['ATP','UTP','GTP','PROTON','UDP_GLUCOSE','CARBON_DIOXIDE','GDP_D_GLUCOSE','CPD_14553','AMP','NADH','FeII','SUPER_OXIDE','PPI']+['Heteroglycans','CELLULOSE','Cellulose','STARCH','ACETALD','ARACHIDIC_ACID_c_','DOCOSANOATE_c_','CPD_16709_c_','STEARIC_ACID_c_','OLEATE_CPD_c_','Octadecadienoate_c_','LINOLENIC_ACID_c_','PALMITATE','CPD_9245_c_','CPD_17412_c_','CPD_17291_c_','CPD_674_c_','Fatty_Acids_c_','CPD_2961_']
    for i in range(0, total_cells):
        for met in multicell_model.metabolites:
            #Exchange cytosolic metabolites
            met_num = "c_" + core_tag + str(i).zfill(2)
            if met_num in met.compartment:
                if "CELLULOSE" not in met.name: #exclude cellulose
                    if met.formula_weight != None:   #excludes osmolytes
                        if met.formula_weight < 900 and met.formula_weight > 0:     #vary depending on the plasmodesmata pore size
                            rxn = Reaction(str(met) + "_Transfer_" + str(i).zfill(2) + "_" + "to" + "_" + str(i+1).zfill(2))
                            rxn.name = str(met) + "_Transfer_" + str(i).zfill(2) + "_" + "to" + "_" + str(i+1).zfill(2)
                            met2 = str(met)[:-2]+ str(i+1).zfill(2)
                            if met2 not in multicell_model.metabolites:
                                new_met = met.copy()
                                new_met.id = met2
                                new_met.compartment = met_num
                                multicell_model.add_metabolites({new_met})
                            rxn.add_metabolites({multicell_model.metabolites.get_by_id(str(met)):-1, multicell_model.metabolites.get_by_id(met2):1})
                            rxn.lower_bound = -1000
                            rxn.upper_bound = 1000
                            multicell_model.add_reactions([rxn])
                            if any([y in met.id for y in no_symp_transport]):
                                rxn.lower_bound=0
                                rxn.upper_bound=0
            #exchange ER metabolites
            met_num = "r_" + str(i).zfill(2) #exchange cytosolic metabolites
            if met_num in met.compartment:
                if met.formula_weight != None:   #excludes osmolytes
                    if met.formula_weight < 900 and met.formula_weight > 0:     #vary depending on the plasmodesmata pore size
                        rxn = Reaction(str(met) + "_Transfer_" + str(i).zfill(2) + "_" + "to" + "_" + str(i+1).zfill(2))
                        rxn.name = str(met) + "_Transfer_" + str(i).zfill(2) + "_" + "to" + "_" + str(i+1).zfill(2)
                        rxn.add_metabolites({multicell_model.metabolites.get_by_id(str(met)):-1, multicell_model.metabolites.get_by_id(str(met)[:-2]+ str(i+1).zfill(2)):1})
                        rxn.lower_bound = -1000
                        rxn.upper_bound = 1000
                        multicell_model.add_reactions([rxn])


###generates one multicell model from multiple models and connects them with linker reactions
#args: model 1, model 2, label of model1 e.g. "epi00", label of model2 e.g. "cor00"
def link_two_models(model1, model2, model1_tag, model2_tag):
    no_symp_transport = ['ATP','UTP','GTP','PROTON','UDP_GLUCOSE','CARBON_DIOXIDE','GDP_D_GLUCOSE','CPD_14553','AMP','NADH','FeII','SUPER_OXIDE','PPI']+['Heteroglycans','CELLULOSE','Cellulose','STARCH','ACETALD','ARACHIDIC_ACID_c_','DOCOSANOATE_c_','CPD_16709_c_','STEARIC_ACID_c_','OLEATE_CPD_c_','Octadecadienoate_c_','LINOLENIC_ACID_c_','PALMITATE','CPD_9245_c_','CPD_17412_c_','CPD_17291_c_','CPD_674_c_','Fatty_Acids_c_','CPD_2961_']
    tissue_model = model1.merge(model2)
    for met in tissue_model.metabolites:
        #exchange cytosolic metabolites
        met_num = "c_" + model1_tag
        if met_num in met.compartment:
            if "CELLULOSE" not in met.name: #exclude cellulose
                if met.formula_weight != None:   #excludes osmolytes
                    if met.formula_weight != None:
                        if met.formula_weight < 900 and met.formula_weight > 0:     #vary depending on the plasmodesmata pore size
                            rxn = Reaction(str(met) + "_Transfer_" + model1_tag + "_" + "to" + "_" + model2_tag)
                            rxn.name = str(met) + "_Transfer_" + model1_tag + "_" + "to" + "_" + model2_tag
                            rxn.add_metabolites({tissue_model.metabolites.get_by_id(str(met)):-1, tissue_model.metabolites.get_by_id(str(met)[:-5]+ model2_tag):1})
                            rxn.lower_bound = -1000
                            rxn.upper_bound = 1000
                            tissue_model.add_reactions([rxn])
                            if any([y in met.id for y in no_symp_transport]):
                                rxn.lower_bound=0
                                rxn.upper_bound=0
                            #print(rxn)
        #exchange ER metabolites
        met_num = "r_" + model1_tag
        if met_num in met.compartment:
            if met.formula_weight != None:
                if met.formula_weight < 900 and met.formula_weight > 0:     #vary depending on whether there's an ER restriction
                    rxn = Reaction(str(met) + "_Transfer_" + model1_tag + "_" + "to" + "_" + model2_tag)
                    rxn.name = str(met) + "_Transfer_" + model1_tag + "_" + "to" + "_" + model2_tag
                    rxn.add_metabolites({tissue_model.metabolites.get_by_id(str(met)):-1, tissue_model.metabolites.get_by_id(str(met)[:-5]+ model2_tag):1})
                    rxn.lower_bound = -1000
                    rxn.upper_bound = 1000
                    tissue_model.add_reactions([rxn])
                    #print(rxn)
    return tissue_model


###delete all reactions in a compartment function
#args: model and compartment in form e.g. "c_"
def remove_compartment(model, compartment, model_tag):
    rxns_to_remove = []
    for rxn in model.reactions:
        if rxn.id[:-6] not in ["ATPase_tx","NADPHoxm_tx","NADPHoxc_tx","NADPHoxp_tx","CO2_tx","O2_tx"]:
            if compartment in str(rxn.compartments):
                rxns_to_remove.append(rxn)
    model.remove_reactions(rxns_to_remove)
    print(compartment + "compartment no longer contains any reactions")


###import constraints from excel
#args: model and excel file in form e.g. "Model_constraints.xlsx"
def import_constraints(model, excel_file):
    import pandas as pd
    constraints = pd.read_excel(excel_file)
    for label, row in constraints.iterrows():
        if str(row["rxn_id"]) != "nan":
            model.reactions.get_by_id(row["rxn_id"]).lower_bound = row["lb"]
            lb = model.reactions.get_by_id(row["rxn_id"]).lower_bound
            model.reactions.get_by_id(row["rxn_id"]).upper_bound = row["ub"]
            ub = model.reactions.get_by_id(row["rxn_id"]).upper_bound
            print(model.reactions.get_by_id(row["rxn_id"])," ub = ", ub, " lb = ", lb)


###set maintenance costs by iterating through ATPase_tx values
#args: model, respiration rate as float e.g. -7.28 (-ve value b/ signifies CO2 leaving the cell)
#make sure units of respiration_rate are correct!
#mean respiration rate: 8.76E-14 umol CO2 per sec per um3 of cytoplasm
def set_maintenance(model, respiration_rate, starting_value, model_tag):
    from cobra import io,flux_analysis
    CO2_tx_flux = 0
    i = starting_value
    while CO2_tx_flux > respiration_rate:
        i=i+1
        tempModel = model.copy()
        tempModel.reactions.get_by_id("ATPase_tx" + model_tag).lower_bound = i
        tempModel.reactions.get_by_id("ATPase_tx" + model_tag).upper_bound = i
        solution = flux_analysis.parsimonious.pfba(tempModel)
        CO2_tx_flux = solution.fluxes.get("CO2_tx" + model_tag)
        print("NGAM ATPase flux ="+str(i))
        print("CO2_tx_flux ="+str(CO2_tx_flux))
        print("---------------------------")

    print("Assumed maintenance cost, ATPase_tx =")
    print(tempModel.reactions.get_by_id("ATPase_tx" + model_tag).flux)
    print(tempModel.reactions.get_by_id("CO2_tx" + model_tag).flux)

    #set ATPase_tx value in model to value determined above
    model.reactions.get_by_id('ATPase_tx' + model_tag).upper_bound = tempModel.reactions.get_by_id("CO2_tx" + model_tag).flux
    model.reactions.get_by_id('ATPase_tx' + model_tag).lower_bound = tempModel.reactions.get_by_id("CO2_tx" + model_tag).flux

###calculate surface area
def sa_calculator(cell_dimensions, model_tag = "_epi00"):
    import pandas as pd
    #import rootslice data
    rootslice = pd.read_excel(cell_dimensions, index_col=0)

    dnut_sa_avg = rootslice.loc[model_tag, "change_surface_area"]
    dnut_vol_tot = rootslice.loc[model_tag, "change_tot_dnut_vol"]

    dnut_pergDW = 1.06383E+13/dnut_vol_tot # 2.541e+13 is um3/gDW according to https://doi.org/10.1371/journal.pone.0239075
    sa_pergDW = dnut_sa_avg*dnut_pergDW # um2/gDW # do i need to times by 4 hours? no because it's the change in surface area over the period of the model that we're using
    return sa_pergDW


###generating expanding cells
#Calculating demand for cellulose, phospholipid and protein based on volumes of cell wall, cell membrane, cytosol from rootslice
#args: model, cell dimension data from Lynch lab in form "filename.xlsl",
#excel file of soluble metabolites in form "file.xlsx",
#excel file of osmotic constraints in form "file.xlsx", model_tag in form e.g. "_cor00"
def make_expanding_cell(model, cell_dimensions, sol_mets, osmotic_constraints, model_tag):    

    #import rootslice data
    rootslice = pd.read_excel(cell_dimensions, index_col=0)

    #calculating gDW per model -- this is the single value I should use for all scaling!!!
    whole_model_vol = rootslice.loc["whole_model", "avg_tot_dnut_vol"] #this value is the average volume of whole model across cell lengths
    gDW_per_tissue = whole_model_vol/1.06383E+13 # 2.541e+13 is um3/gDW according to https://doi.org/10.1371/journal.pone.0239075


    #total cell volume (old approach to calculating scaling factor)
    #cell_vol_tot = rootslice.loc[model_tag, "change_tot_dnut_vol"]
    #cell_pergDW = 2.541e+13/cell_vol_tot # 2.541e+13 is um3/gDW according to https://doi.org/10.1371/journal.pone.0239075


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

    ##Adjusting the Biomass_tx reaction to produce a generic protein-producing reactions that turns all the AAs into Protein_b:
    rxn = model.reactions.get_by_id("Biomass_tx" + model_tag)

    #remove Ca, K and Mg from the Biomass reaction so it only uses AAs:
    for ion in ["Ca_b", "K_b", "Mg_b"]:
        met = model.metabolites.get_by_id(ion + model_tag)
        coeff = rxn.metabolites.get(met)
        rxn.add_metabolites({met:-1*coeff})

    #add protein_b metabolite to represent all protein in cell
    try:
        met = Metabolite("Protein_b" + model_tag,name="Protein_biomass", compartment = "b"+model_tag)
        formula_dict = rxn.check_mass_balance()
        met.formula = "".join([atom+str(formula_dict[atom]*-1) for atom in formula_dict.keys() if atom != "charge"])
        met.charge = int(formula_dict["charge"]*-1)
        rxn.add_metabolites({met:1})
    except:
        print("already done")

    ##Add new expanding biomass drain reaction of cellulose + lipid + protein in exactly the
    #right ratios as determined by earlier calculations based on rootslice data
    #(ratio represented by e.g.  cellulose_acc)
    #This is fixed at a given flux because we know how much it should produce based on outputs from rootslice.
    rxn = Reaction("Biomass_expanding_cell" + model_tag)
    rxn.add_metabolites({model.metabolites.get_by_id("Cellulose_b" + model_tag):-1*mmol_cellulose})
    rxn.add_metabolites({model.metabolites.get_by_id("L_PHOSPHATIDATE_p" + model_tag):-1*mmol_phospholipid})
    rxn.add_metabolites({model.metabolites.get_by_id("Protein_b" + model_tag):-1*mmol_protein})
    model.add_reactions([rxn])
    #fix flux at certain rate:
    rxn.lower_bound = 1
    rxn.upper_bound = 1
    #check biomass reactions
    #print(rxn.reaction)


    #setting up pseudo-metabolites representing osmolytes to enter osmotic constraints into model
    VO = Metabolite("VO" + model_tag,name="vacuolar osmolytes", compartment = "v"+model_tag, charge = 0)
    VC = Metabolite("VC" + model_tag,name="vacuolar charge", compartment = "v"+model_tag, charge = 0)
    CO = Metabolite("CO" + model_tag,name="cytosolic osmolytes", compartment = "c"+model_tag, charge = 0)
    CC = Metabolite("CC" + model_tag,name="cytosolic charge", compartment = "c"+model_tag, charge = 0)

    #importing soluble metabolite ratios and osmotically active metabolites
    osm_constraints = pd.read_excel(osmotic_constraints)
    #sol_mets_df = pd.read_excel(sol_mets)

    #adding soluble_biomass reaction
    rxn = Reaction("soluble_biomass" + model_tag)
    model.add_reactions([rxn])

    #organic pseudometabolites in vacuole
    for label, row in osm_constraints.iterrows():
        METB = Metabolite(row["met_id"]+"_b" + model_tag, compartment = "b"+model_tag, charge = model.metabolites.get_by_id(row["met_id"]+"_c" + model_tag).charge)
        if row["charge_form"] != "none":
            METV = model.metabolites.get_by_id(row["met_id"]+"_v" + model_tag)
            aMETV = model.metabolites.get_by_id("a"+row["met_id"]+"_v" + model_tag)
            charge = (METV.charge*row["met_stoyk"]*-1)+(aMETV.charge*row["amet_stoyk"]*-1)
            rxn = Reaction(row["met_id"]+"_v_biomass" + model_tag)
            rxn.add_metabolites({METV:row["met_stoyk"],aMETV:row["amet_stoyk"],VC:charge,VO:1,METB:1})
            model.add_reactions([rxn])

        elif row["met_id"] == "SUC" or row["met_id"] == "STARCH":
            continue

        else:
            METV = model.metabolites.get_by_id(row["met_id"]+"_v" + model_tag)
            rxn = Reaction(row["met_id"]+"_v_biomass" + model_tag)
            rxn.add_metabolites({METV:-1,VC:METV.charge,VO:1,METB:1})
            model.add_reactions([rxn])

        #adding soluble metabolites in correct ratios
        #rxn = model.reactions.get_by_id("soluble_biomass" + model_tag)
        #sol_mets_df.loc[sol_mets_df["met"] == row["met_id"], ["amount"]]
        #amount = sol_mets_df.loc[sol_mets_df["met"] == row["met_id"], ["amosunt"]]
        #if amount.empty == False:
        #    amount_float = float(amount.to_numpy()[0])
        #    rxn.add_metabolites({METB:-1*amount_float})

    #organic pseudometabolites in cytosol (pH7 => don't worry about charge state)
    for label, row in osm_constraints.iterrows():
        METC = model.metabolites.get_by_id(row["met_id"]+"_c" + model_tag)
        METB = model.metabolites.get_by_id(row["met_id"]+"_b" + model_tag)
        rxn = Reaction(row["met_id"]+"_c_biomass" + model_tag)
        rxn.add_metabolites({METC:-1,CC:METC.charge,CO:1,METB:1})
        model.add_reactions([rxn])
    ##giving ions a biomass drain in each cell##
        if row["biomass"] != "none":
            rxn = Reaction(row["met_id"]+"_b_biomass"+model_tag)
            rxn.add_metabolites({METB:-1})
            model.add_reactions([rxn])

    #Set Vv and Vc equal to vacuolar and cytoplasmic volumes at the end of the time period (like Sanu did)
    Vv = rootslice.loc[model_tag, "final_vac_dnut_vol"]
    Vc = rootslice.loc[model_tag, "final_cytoplasm_dnut_vol"]
    Vcell = rootslice.loc[model_tag, "final_tot_dnut_vol"]

    #sets ratio based on eqn 2 in paper
    rxn = Reaction("VacCytRatio" + model_tag)
    WO = Metabolite("TotalSolute" + model_tag, charge = 0, compartment = "c"+model_tag)    #totalsolute placed in cytosol for now -- seems best place for it -- removed from compartment!
    rxn.add_metabolites({CO:-1,VO:-1*(Vv/Vc),WO:1+(Vv/Vc)})
    rxn.lower_bound = 0
    rxn.upper_bound = 1000
    model.add_reactions([rxn])

    C_cell = 250                      #Bloom et al. 2012; units = mol/m3
    C_cell = 250 * (10**-15)          #units = mmol/um3

    #based on eqn 1 in paper. constrains reactions based on overall change in volume that you'll get from rootslice
    rxn = Reaction("TotalSoluteConstraint" + model_tag)
    rxn.add_metabolites({WO:-1})
    rxn.lower_bound = float(C_cell*Vcell/gDW_per_tissue) #mmol/gDW
    rxn.upper_bound = float(C_cell*Vcell/gDW_per_tissue) #mmol/gDW
    model.add_reactions([rxn])

    #print reactions to check outputs
    #print(model.reactions.get_by_id("TotalSoluteConstraint" + model_tag))
    #print(model.reactions.get_by_id("soluble_biomass" + model_tag))


###add apoplast metabolites
def apoplast_metabolite_adder(model, comp):
    from cobra import Reaction, Metabolite

    list_ions = ["FeIII", "MGII", "CAII", "Pi", "KI", "FeII", "NITRATE", "WATER", "CARBON_DIOXIDE", "PROTON"]
    for i in list_ions:
        try:
            ion_details = model.metabolites.get_by_id(i + "_c_epi00")
            ion_met = Metabolite(
            f"{i}_{comp}",
            formula = ion_details.formula,
            name = f"{i}[{comp}]",
            compartment = comp)
            model.add_metabolites(ion_met)
            #print(ion_met)
        except:
            f"{i} has already been added to {comp}"

###create apoplast uptake reactions
#args: model e.g. tissue_model, compartment e.g. "ap1", model_tag e.g. "epi00"
def apoplast_transport_reactions(model, comp, model_tags):
    from cobra import Reaction, Metabolite
    #add apoplast metabolites
    list_ions = ["FeIII", "MGII", "CAII", "Pi", "KI", "FeII", "NITRATE", "WATER", "CARBON_DIOXIDE", "PROTON"]
    for i in list_ions:
        try:
            ion_details = model.metabolites.get_by_id(i + "_c_epi00")
            ion_met = Metabolite(
            f"{i}_{comp}",
            formula = ion_details.formula,
            name = f"{i}[{comp}]",
            compartment = comp)
            model.add_metabolites(ion_met)
            #print(ion_met)
        except:
            f"{i} has already been added to {comp}"

    # add proton ATPase and transport reaction from each cell type to apoplast
    for tag in model_tags:
        #add to all cells a proton ATPase to pump protons into apoplast
        rxn = Reaction("PROTON_ATPase_ap_"+tag)
        rxn.name = "PROTON_ATPase_ap_"+tag
        rxn.add_metabolites({
            model.metabolites.get_by_id("ATP_c_"+tag): -1.0,
            model.metabolites.get_by_id("PROTON_c_"+tag): -1.0,
            model.metabolites.get_by_id("WATER_c_"+tag): 1.0,
            model.metabolites.get_by_id("aATP_c_"+tag): 1.0,
            model.metabolites.get_by_id("ADP_c_"+tag): 1.0,
            model.metabolites.get_by_id("PROTON_"+comp): 1.0,
            model.metabolites.get_by_id("Pi_c_"+tag): 1.0,
            model.metabolites.get_by_id("aADP_c_"+tag): 1.0,
            model.metabolites.get_by_id("aPi_c_"+tag): 1.0,
            })
        model.add_reactions([rxn])

        #if no uptake cost
        for i in ["FeIII", "MGII", "CAII", "FeII", "WATER", "CARBON_DIOXIDE"]:
            ion_ap = model.metabolites.get_by_id(f"{i}_{comp}")
            ion = f"{i}_c_{tag}"
            ion_met = model.metabolites.get_by_id(ion)
            reaction = Reaction(ion + "_apoplast_exchange_" + comp)
            reaction.name = ion + " exchange with apoplast"
            reaction.add_metabolites({
            ion_met: -1.0,
            ion_ap: 1.0,
            })
            reaction.lower_bound = -1000
            reaction.upper_bound = 1000
            model.add_reactions([reaction])
            #print(reaction)

        #if uptake cost
        for i in ["Pi", "KI", "NITRATE"]:
            #laoding into apoplast
            ion_ap = model.metabolites.get_by_id(f"{i}_{comp}")
            ion = f"{i}_c_{tag}"
            ion_met = model.metabolites.get_by_id(ion)
            reaction = Reaction(ion + "_apoplast_loading_" + comp)
            reaction.name = ion + " loading into apoplast"
            reaction.add_metabolites({
            ion_met: -1.0,
            ion_ap: 1.0,
            })
            reaction.lower_bound = 0
            reaction.upper_bound = 1000
            model.add_reactions([reaction])

            #unloading from apoplast
            ion_ap = model.metabolites.get_by_id(f"{i}_{comp}")
            ion = f"{i}_c_{tag}"
            ion_met = model.metabolites.get_by_id(ion)
            rxn = Reaction(ion + "_apoplast_unloading_" + comp)
            reaction.name = ion + " unloading from apoplast"
            reaction.add_metabolites({
            ion_ap: -1.0,
            model.metabolites.get_by_id("PROTON_"+comp): -1.0,
            ion_met: 1.0,
            model.metabolites.get_by_id("PROTON_c_"+tag): 1.0
            })
            reaction.lower_bound = 0
            reaction.upper_bound = 1000
            model.add_reactions([rxn])





###adds phloem mets and reactions
#args: model, excel_file of phloem metabolites in form "file", model tag e.g. per00
"""def add_phloem_rxns(model, excel_file, model_tag):
    import pandas as pd
    from cobra import Metabolite, Reaction
    phloem_met_df = pd.read_excel(excel_file)
    phloem_met = []
    for label, row in phloem_met_df.iterrows():
        met_name = str(row["phloem_met"]) + model_tag
        print(met_name)
        met_tagged = model.metabolites.get_by_id(met_name)
        phloem_met.append(met_tagged)

    phloem_EX_X_rxns = []
    proton_ph_tagged = "PROTON_ph" + model_tag
    proton_ph = Metabolite(str(model.metabolites.get_by_id(proton_ph_tagged),
                          formula = model.metabolites.get_by_id("PROTON_c"+model_tag).formula,
                          charge = model.metabolites.get_by_id("PROTON_c"+model_tag).charge, compartment = "ph")
    model.add_metabolites(proton_ph)

    rxn = Reaction(str('EX_X_') + str(proton_ph))
    rxn.add_metabolites({proton_ph: 1.0})
    model.add_reactions([rxn])

    #phloem_rxns.append(rxn)
    #takes metabolites in phloem (phloem_met), converts them from _c to _ph, generates new import reaction with phloem components as reactants, and cytosolic version as products with proton symport
    all_phloem_rxns = []
    for i in phloem_met:
        metabolite_name = Metabolite(str(i).rsplit('_',1)[0] + str('_ph'), compartment = 'ph',
                                    formula = model.metabolites.get_by_id(str(i).rsplit('_',1)[0] + str('_c')).formula,
                                     charge = model.metabolites.get_by_id(str(i).rsplit('_',1)[0] + str('_c')).charge)
        model.add_metabolites(metabolite_name)
        rxn1 = Reaction(str(metabolite_name) + str('_c'))
        rxn1.upper_bound = 1000
        rxn1.lower_bound = 0
        rxn1.add_metabolites({metabolite_name: -1.0, model.metabolites.get_by_id("PROTON_c"): 1.0, i : 1.0, model.metabolites.get_by_id("PROTON_ph"): -1.0})
        model.add_reaction(rxn1)
        all_phloem_rxns.append(rxn1)

    #all_phloem_rxns are exchange rxns where phloem metabolites enter the model
    for i in phloem_met:
        metabolite_name = str(model.metabolites.get_by_id(str(i).rsplit('_',1)[0] + str('_ph')))
        rxn2 = Reaction(str('EX_X_') + str(metabolite_name))
        model.add_reaction(rxn2)
        rxn2.add_metabolites({metabolite_name: 1.0})

        phloem_EX_X_rxns.append(rxn2)
        all_phloem_rxns.append(rxn2)


    #add proton ATPase to generate proton gradient for uptake from phloem
    reaction_name = Reaction('PROTON_ATPase_c_ph')

    reaction_name.add_metabolites({model.metabolites.get_by_id("PROTON_c"): -0.45, model.metabolites.get_by_id("WATER_c"): -1,
                                  model.metabolites.get_by_id("aATP_c"): -0.35, model.metabolites.get_by_id("ATP_c"): -0.65,
                                 model.metabolites.get_by_id("PROTON_ph"): 1, model.metabolites.get_by_id("ADP_c"): 0.5,
                                  model.metabolites.get_by_id("Pi_c"): 0.7, model.metabolites.get_by_id("aADP_c"): 0.5,
                                  model.metabolites.get_by_id("aPi_c"): 0.3})
    model.add_reaction(reaction_name)
    print(reaction_name.reaction)
"""
