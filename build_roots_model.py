# Build roots model
#import necessary packages
import cobra
from cobra import io, Reaction, Metabolite,flux_analysis
import pandas as pd
#import all my functions
from root_study_functions import *
from cobra.flux_analysis.parsimonious import pfba
from scipy.optimize import curve_fit

def main(saveas='',clp=[1,1,1],orootsdir = '/Users/user/Documents/Ox/roots/roots-model/'):
    model = cobra.io.read_sbml_model(orootsdir+"base_cell_model.xml") #import base model
    epi00_model, epi_dict = set_up_multicell_model(model, 1, "epidermal")
    #remove inappropriate compartments
    remove_compartment(epi00_model, "ph_", "_epi00")
    remove_compartment(epi00_model, "xy_", "_epi00")

    #create lactate cytosolic metabolite
    lactate = Metabolite(
        'L_LACTATE_c',
        formula='C3H5O3',
        name='lactate',
        charge=-1,
        compartment='c')
    model.add_metabolites(lactate)

    #create lactate soil metabolite
    lactate = Metabolite(
        'L_LACTATE_e',
        formula='C3H6O3',
        name='lactate',
        charge=-1,
        compartment='e')
    model.add_metabolites(lactate)

    #add pyruvate to lactate rxn
    rxn = Reaction("L_LACTATEDEHYDROG_RXN_c")
    rxn.name = "Fermentation to lactate"
    rxn.subsystem = "fermentation"
    rxn.add_metabolites({
        model.metabolites.PYRUVATE_c: -1.0,
        model.metabolites.NADH_c: -1.0,
        model.metabolites.PROTON_c: -1.0,
        model.metabolites.L_LACTATE_c: 1.0,
        model.metabolites.NAD_c: 1.0,
        })
    model.add_reaction(rxn)

    #add epidermal cell specific uptake constraints
    sa_pergDW = sa_calculator(orootsdir+"model_constraints_data/cell_dimensions.xlsx")

    #from OpenSimRoot (see Jagdeep email)
    #nitrogen uptake: Imax (uMol/cm2/day) = 1.27 
    #phosphorus uptake: Imax (uMol/cm2/day) = 0.0555 

    n_uptake = 1.27 #uMol/cm2/day
    n_uptake = n_uptake/24 #uMol/cm2/hour
    n_uptake = n_uptake/1000 #mmol/cm2/hour
    n_uptake = n_uptake/1e+8 #mmol/um2/hour
    n_uptake = n_uptake * sa_pergDW #mmol/gDW/hour
    print('N bound: ',n_uptake)
    epi00_model.reactions.Nitrate_ec_epi00.upper_bound = float(n_uptake)
    epi00_model.reactions.Nitrate_ec_epi00.lower_bound = float(n_uptake)

    p_uptake = 0.0555 #uMol/cm2/day
    p_uptake = p_uptake/24 #uMol/cm2/hour
    p_uptake = p_uptake/1000 #mmol/cm2/hour
    p_uptake = p_uptake/1e+8 #mmol/um2/hour
    p_uptake = p_uptake * sa_pergDW #mmol/gDW/hour
    p_uptake
    #adding phosphate constraint makes the model infeasible!

    #epi00_model.reactions.Pi_ec_epi00.upper_bound = float(p_uptake)
    #epi00_model.reactions.Pi_ec_epi00.lower_bound = float(p_uptake)

    #make cell into expanding cell using its specific cell dimensions 
    cell_dimensions_pwy=orootsdir+"model_constraints_data/cell_dimensions.xlsx"
    osmotic_constraints_pwy=orootsdir+"model_constraints_data/osmotic_constraints.xlsx"
    make_expanding_cell(epi00_model, cell_dimensions_pwy, osmotic_constraints_pwy, "_epi00",clp)

    placeholder, cor_models_dict = set_up_multicell_model(model, 8, "cortical")

    #remove inappropriate compartments from all 8 cortical models
    for i in range(0,8):
        remove_compartment(cor_models_dict[i], "ph_", f"_cor0{str(i)}")
        remove_compartment(cor_models_dict[i], "xy_", f"_cor0{str(i)}")
        remove_compartment(cor_models_dict[i], "e_", f"_cor0{str(i)}")
        
    #remove boundary reactions
    #for i in range(0,8):
    #    remove_boundary_rxns(cor_models_dict[i])

    #make cell into expanding cell using its specific cell dimensions 
    for i in range(0,8):
        make_expanding_cell(cor_models_dict[i], cell_dimensions_pwy, osmotic_constraints_pwy, f"_cor0{str(i)}",clp)

    #create multi_cor_model by combining the models in the dictionary adjusted models
    a = cor_models_dict[0].copy()
    multi_cor_model = a
    for i in range(1,8):
        multi_cor_model += cor_models_dict[i]

    end00_model, end_dict = set_up_multicell_model(model, 1, "endodermal")

    #remove inappropriate compartments
    remove_compartment(end00_model, "ph_", "_end00")
    remove_compartment(end00_model, "xy_", "_end00")
    remove_compartment(end00_model, "e_", "_end00") #CAREFUL if just set to "e" IT REMOVES EVERYTHING WITH "e" IN IT SO ALL COMPARTMENTS B/ THEY HAVE "E"ND00

    #make cell into expanding cell using its specific cell dimensions 
    make_expanding_cell(end00_model, cell_dimensions_pwy, osmotic_constraints_pwy, "_end00",clp)

    per00_model, peri_dict = set_up_multicell_model(model, 1, "pericycle")
    #remove inappropriate compartments
    remove_compartment(per00_model, "e_", "_per00")

    #make cell into expanding cell using its specific cell dimensions 
    make_expanding_cell(per00_model, cell_dimensions_pwy, osmotic_constraints_pwy, "_per00",clp)

    #linking cortical cells to each other
    create_linker_fluxes(multi_cor_model, 7, "cor")

    #linking epidermal epi00 to cortical cor00
    tissue_model1 = link_two_models(epi00_model, multi_cor_model, "epi00", "cor00")

    #linking cortical cor07 to endodermal end00
    tissue_model2 = link_two_models(tissue_model1, end00_model, "cor07", "end00")

    #linking endodermal end00 to pericycle per00
    #note the end00_model has already been fully incoporated into tissue_model2
    tissue_model = link_two_models(tissue_model2, per00_model, "end00", "per00")

    tissue_model.compartments = {"ap1": "apoplast_epi_to_end", "ap2": "apoplast_end_to_stele"}

    #adding apoplast from epidermis up to casparian strip
    apoplast_transport_reactions(tissue_model, "ap1", ["epi00","cor00","cor01","cor02","cor03","cor04","cor05","cor06","cor07","end00"])

    #adding apoplast post-casparian strip from endodermis to stele
    apoplast_transport_reactions(tissue_model, "ap2", ["end00","per00"])

    #create whole tissue soluble biomass reaction
    rxn = Reaction("soluble_biomass_tissue")
    tissue_model.add_reaction(rxn)

    #load soluble metabolites from excel
    sol_mets_df = pd.read_excel(orootsdir+"model_constraints_data/soluble_metabolite_composition.xlsx")
    cellList = ["_epi00","_cor00","_cor01","_cor02","_cor03","_cor04","_cor05","_cor06","_cor07","_end00", "_per00"]
    #add metabolites to tissue wide sol biomass reaction
    for index, row in sol_mets_df.iterrows():
        met = row["met"]
        met_b_tissue = Metabolite(met+"_b_tissue", compartment = "b_tissue")    
        amount_float = float(row["amount"])
        rxn = tissue_model.reactions.get_by_id("soluble_biomass_tissue")
        rxn.add_metabolites({met_b_tissue:-1*amount_float})
        
        for model_tag in cellList:
            if met == "STARCH":
                met_b_celltype = tissue_model.metabolites.get_by_id("Starch_b"+model_tag)
            else:
                met_b_celltype = tissue_model.metabolites.get_by_id(met+"_b"+model_tag)
            met_b_tissue = tissue_model.metabolites.get_by_id(met+"_b_tissue")
            rxn_to_tissue = Reaction(met+model_tag+"_tissue_biomass")
            rxn_to_tissue.add_metabolites({met_b_celltype:-1,met_b_tissue:1})
            tissue_model.add_reaction(rxn_to_tissue)
            
    ### MAINTENANCE        
    #generic ATPase and NADPH oxidase
    rootslice = pd.read_excel(orootsdir+"model_constraints_data/cell_dimensions.xlsx", index_col=0)
    Maintenance_constraints = {}
    ratio_cyt_dnut_vol={}
    reaction = Reaction('ATPase_balancer')
    reaction.name = 'ATPase_balancer'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000 
    for cell in cellList:
        Maintenance_constraints[cell] = Metabolite("ATPase_NADPHoxidase_constraint"+cell,name =  "ATPase_NADPHoxidase_constraint"+cell, compartment = "c"+cell)
        tissue_model.reactions.get_by_id("ATPase_tx"+cell).add_metabolites({Maintenance_constraints[cell]:1})
        ratio_cyt_dnut_vol[cell] = rootslice.loc[cell, "ratio_cytoplasm_dnut_vol"]
        reaction.add_metabolites({Maintenance_constraints[cell]: -1/ratio_cyt_dnut_vol[cell]})

    tissue_model.add_reaction(reaction)


    # ### STORAGE This doesn't make sense??
    # #determine ratio of cytoplasmic volumes of each cell relative to cell with smallest cyt vol (per00)
    # rootslice = pd.read_excel("model_constraints_data/cell_dimensions.xlsx", index_col=0)
    
    # a1 = rootslice.loc["_per00", "ratio_cytoplasm_dnut_vol"]
    # b1 = rootslice.loc["_end00", "ratio_cytoplasm_dnut_vol"]
    # c1 = rootslice.loc["_cor07", "ratio_cytoplasm_dnut_vol"]
    # d1 = rootslice.loc["_cor06", "ratio_cytoplasm_dnut_vol"]
    # e1 = rootslice.loc["_cor05", "ratio_cytoplasm_dnut_vol"]
    # f1 = rootslice.loc["_cor04", "ratio_cytoplasm_dnut_vol"]
    # g1 = rootslice.loc["_cor03", "ratio_cytoplasm_dnut_vol"]
    # h1 = rootslice.loc["_cor02", "ratio_cytoplasm_dnut_vol"]
    # i1 = rootslice.loc["_cor01", "ratio_cytoplasm_dnut_vol"]
    # j1 = rootslice.loc["_cor00", "ratio_cytoplasm_dnut_vol"]
    # k1 = rootslice.loc["_epi00", "ratio_cytoplasm_dnut_vol"]

    # #constraining ATPases of each cell in ratio given their cytoplasmic volumes
    # a = tissue_model.reactions.get_by_id("ATPase_tx_per00")
    # b = tissue_model.reactions.get_by_id("ATPase_tx_end00")
    # c = tissue_model.reactions.get_by_id("ATPase_tx_cor07")
    # d = tissue_model.reactions.get_by_id("ATPase_tx_cor06")
    # e = tissue_model.reactions.get_by_id("ATPase_tx_cor05")
    # f = tissue_model.reactions.get_by_id("ATPase_tx_cor04")
    # g = tissue_model.reactions.get_by_id("ATPase_tx_cor03")
    # h = tissue_model.reactions.get_by_id("ATPase_tx_cor02")
    # i = tissue_model.reactions.get_by_id("ATPase_tx_cor01")
    # j = tissue_model.reactions.get_by_id("ATPase_tx_cor00")
    # k = tissue_model.reactions.get_by_id("ATPase_tx_epi00")

    # ATPase_cells_flux  = tissue_model.problem.Constraint(
    #     (a.flux_expression) -  b1*(b.flux_expression) -  c1*(c.flux_expression) -  d1*(d.flux_expression) -  e1*(e.flux_expression) -  f1*(f.flux_expression)
    #     -  g1*(g.flux_expression) -  h1*(h.flux_expression) -  i1*(i.flux_expression) -  j1*(j.flux_expression) -  k1*(k.flux_expression), 
    #     name = "ATPase_cell_flux",
    #     lb = 0,
    #     ub= 0)
    # tissue_model.add_cons_vars(ATPase_cells_flux)

    for rxn in tissue_model.reactions:
        if 'L_LACTATE_ec_' in rxn.id:
            metlist=[x.id for x in rxn.metabolites]
            for met in metlist:
                if 'PROTON' in met:
                    rxn.add_metabolites({met:-rxn.get_coefficient(met)})

    resp_rate = 8.76E-14 #umol CO2 per sec per um3 of cytoplasm
    resp_rate = resp_rate*60*60 #umol/hour/um3 cytoplasm
    resp_rate = resp_rate*(1e-3) #mmol/hour/um3 cytoplasm
    #resp_rate = resp_rate*4 #mmol/4 hours/um3 cytoplasm

    # set_maintenance_hilary(tissue_model, resp_rate, orootsdir+"model_constraints_data/cell_dimensions.xlsx")
    # tissue_model.reactions.ATPase_tx_per00=0.15
    tissue_model.reactions.ATPase_tx_per00=0

    rootspFBAobj(tissue_model,exclude_from_pFBA=[])
    # Add phloem carbon import restriction
    ph_rxns = [rxn for rxn in tissue_model.reactions if 'EX_X' in rxn.id and '_ph_' in rxn.id]
    carboncountmet = Metabolite('Phloem_carbon_in',compartment='ph_per00')
    for rxn in ph_rxns:
        carbon_translocated = sum([count_carbon(met) for met in rxn.products])
        rxn.add_metabolites({carboncountmet:carbon_translocated})
    ph_carbon_restr = Reaction('Phloem_carbon_import')
    tissue_model.add_reaction(ph_carbon_restr)
    # This number is from Lohaus 2000. Flux through ph_caron_restr tells us how much carbon is being imported relative to this estimate.
    ph_carbon_restr.add_metabolites({carboncountmet:-3.79682498*4.2})
    ph_carbon_restr.upper_bound=1000
    ph_carbon_restr.lower_bound=0

    # Solve and save
    solution = pfba(tissue_model)
    io.write_sbml_model(tissue_model,'roots_model'+saveas+'-'+str(clp[0])+'-'+str(clp[1])+'-'+str(clp[2])+'.xml')
    return tissue_model,solution

def count_carbon(met):
    carbons=0
    metformula=met.formula
    if metformula and 'C' in metformula:
        formula_nondigits=[x for x in range(len(metformula)) if not metformula[x].isdigit()]
        cind = metformula.index('C')
        cnext = [x for x in formula_nondigits if x>cind+1]
        carbons = int(metformula[(cind+1):(cnext[0])])
    return carbons

def set_maintenance(model, respiration_rate, cell_dimensions):
    whole_model_vol = 6216571.895 #this value is the average volume of whole model across cell lengths
    gDW_per_tissue = whole_model_vol/2.541e+13 # 2.541e+13 is um3/gDW according to https://doi.org/10.1371/journal.pone.0239075

    rootslice = pd.read_excel(cell_dimensions, index_col=0)
    
    whole_model_cytoplasm_dnut_vol = rootslice.loc["whole_model", "avg_cytoplasm_dnut_vol"]
    # whole_model_tot_dnut_vol = rootslice.loc["whole_model", "avg_tot_dnut_vol"]
       
    
    #scale respiration rate to cytoplasmic volume of whole model
    vol_adj_resp_rate = respiration_rate * -1 * whole_model_cytoplasm_dnut_vol #mmol/model; value is -ve because CO2 leaving the cell
    vol_adj_resp_rate = vol_adj_resp_rate / gDW_per_tissue #mmol/gDW
    
    print(f"the respiration rate adjusted for per gDW was {vol_adj_resp_rate}")
    
    from cobra import io,flux_analysis
    tempModel = model.copy()
    CO2_tx_flux = 0
    i = 0
    while CO2_tx_flux > vol_adj_resp_rate:
        i=i+0.01
        tempModel.reactions.get_by_id("ATPase_tx_per00").lower_bound = i
        tempModel.reactions.get_by_id("ATPase_tx_per00").upper_bound = i
        solution = flux_analysis.parsimonious.pfba(tempModel)
        CO2_tx_flux = solution.fluxes.get("CO2_tx_epi00")
        print("NGAM ATPase flux = "+str(i))
        print("CO2_tx_flux = "+str(CO2_tx_flux))
        print("---------------------------")

    print("Assumed maintenance cost, ATPase_tx =")
    print(tempModel.reactions.get_by_id("ATPase_tx_per00").flux)
    print("CO2 flux: ", tempModel.reactions.get_by_id("CO2_tx_epi00").flux)

    #set ATPase_tx value in model to value determined above
    model.reactions.get_by_id("ATPase_tx_per00").upper_bound = tempModel.reactions.get_by_id("ATPase_tx_per00").flux
    model.reactions.get_by_id("ATPase_tx_per00").lower_bound = tempModel.reactions.get_by_id("ATPase_tx_per00").flux
    for rxn in model.reactions.query("ATPase_tx_"):
        print(rxn.x, rxn)        

def set_maintenance_hilary(model, respiration_rate, cell_dimensions,verbose=True):
    whole_model_vol = 6216571.895 #this value is the average volume of whole model across cell lengths
    gDW_per_tissue = whole_model_vol/2.541e+13 # 2.541e+13 is um3/gDW according to https://doi.org/10.1371/journal.pone.0239075

    rootslice = pd.read_excel(cell_dimensions, index_col=0)
    
    whole_model_cytoplasm_dnut_vol = rootslice.loc["whole_model", "avg_cytoplasm_dnut_vol"]
       
    
    #scale respiration rate to cytoplasmic volume of whole model
    vol_adj_resp_rate = respiration_rate * -1 * whole_model_cytoplasm_dnut_vol #mmol/model; value is -ve because CO2 leaving the cell
    vol_adj_resp_rate = vol_adj_resp_rate / gDW_per_tissue #mmol/gDW
    if verbose:
        print(f"the respiration rate adjusted for per gDW was {vol_adj_resp_rate}")
    
    tempModel = model.copy()
    CO2_tx_flux = 0
    i = 0
    co2_vals={}
    for it in [0.1,0.3,0.3]:
        i=i+it
        tempModel.reactions.get_by_id("ATPase_tx_per00").lower_bound = i
        # tempModel.reactions.get_by_id("ATPase_tx_per00").upper_bound = i
        solution = tempModel.optimize()
        CO2_tx_flux = solution.fluxes.get("CO2_tx_epi00")
        co2_vals[i]=(vol_adj_resp_rate - CO2_tx_flux)/vol_adj_resp_rate
        if verbose:
            print("NGAM ATPase flux = "+str(i))
            print("CO2_tx_flux = "+str(CO2_tx_flux))
            print("---------------------------")


    popt, pcov = curve_fit(Hyp_func, list(co2_vals.keys()), list(co2_vals.values()),maxfev=10000)
    i = (-popt[1])/popt[0]
    if i<0:
        i=0
        print('Warning: Predicted ATPase_tx flux is ', i)
    model_feasible=0
    try:
        tempModel.reactions.get_by_id("ATPase_tx_per00").lower_bound = i
        # model.reactions.get_by_id("ATPase_tx_per00").upper_bound = i
        sol = flux_analysis.parsimonious.pfba(tempModel)
        model_feasible=1
        count+=1
    except:
        print('Warning: Predicted ATPase_tx value infeasible!')
        while not model_feasible:
            try:
                i = (-popt[1])/popt[0]+ ((vol_adj_resp_rate -sol.fluxes['ATPase_tx_per00']))/4
                tempModel.reactions.get_by_id("ATPase_tx_per00").lower_bound = i
                sol = flux_analysis.parsimonious.pfba(tempModel)
                model_feasible=1
                count+=1
            except:
                print('Warning: Predicted ATPase_tx value infeasible again!')
    if verbose:
        print("Assumed maintenance cost, ATPase_tx =")
        print(sol["ATPase_tx_per00"])
        print("CO2 flux: ", sol["CO2_tx_epi00"])

    #set ATPase_tx value in model to value determined above
    # model.reactions.get_by_id("ATPase_tx_per00").upper_bound = tempModel.reactions.get_by_id("ATPase_tx_per00").flux
    if 'ATPase_tx_per00' not in model.reactions:
        print([rxn.id for rxn in model.reactions if 'ATPase_tx' in rxn.id])
    model.reactions.get_by_id("ATPase_tx_per00").lower_bound = i
    for rxn in model.reactions.query("ATPase_tx_"):
        print(rxn.x, rxn)        
def Hyp_func(x, a,b):
    return a*x+b

###generating expanding cells
#Calculating demand for cellulose, phospholipid and protein based on volumes of cell wall, cell membrane, cytosol from rootslice
#args: model, cell dimension data from Lynch lab in form "filename.xlsl",
#excel file of soluble metabolites in form "file.xlsx",
#excel file of osmotic constraints in form "file.xlsx", model_tag in form e.g. "_cor00"
def make_expanding_cell(model, cell_dimensions, osmotic_constraints, model_tag, bio_components_clp=[1,1,1]):

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
    rxn.add_metabolites({model.metabolites.get_by_id("Cellulose_b" + model_tag):-bio_components_clp[0]*mmol_cellulose})
    rxn.add_metabolites({model.metabolites.get_by_id("L_PHOSPHATIDATE_p" + model_tag):-bio_components_clp[1]*mmol_phospholipid})
    rxn.add_metabolites({model.metabolites.get_by_id("Protein_b" + model_tag):-bio_components_clp[2]*mmol_protein})
    model.add_reaction(rxn)
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
    model.add_reaction(rxn)

    #organic pseudometabolites in vacuole
    for label, row in osm_constraints.iterrows():
        METB = Metabolite(row["met_id"]+"_b" + model_tag, compartment = "b"+model_tag, charge = model.metabolites.get_by_id(row["met_id"]+"_c" + model_tag).charge)
        if row["charge_form"] != "none":
            METV = model.metabolites.get_by_id(row["met_id"]+"_v" + model_tag)
            aMETV = model.metabolites.get_by_id("a"+row["met_id"]+"_v" + model_tag)
            charge = (METV.charge*row["met_stoyk"]*-1)+(aMETV.charge*row["amet_stoyk"]*-1)
            rxn = Reaction(row["met_id"]+"_v_biomass" + model_tag)
            rxn.add_metabolites({METV:row["met_stoyk"],aMETV:row["amet_stoyk"],VC:charge,VO:1,METB:1})
            model.add_reaction(rxn)

        elif row["met_id"] == "SUC" or row["met_id"] == "STARCH":
            continue

        else:
            METV = model.metabolites.get_by_id(row["met_id"]+"_v" + model_tag)
            rxn = Reaction(row["met_id"]+"_v_biomass" + model_tag)
            rxn.add_metabolites({METV:-1,VC:METV.charge,VO:1,METB:1})
            model.add_reaction(rxn)

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
        model.add_reaction(rxn)
    ##giving ions a biomass drain in each cell##
        if row["biomass"] != "none":
            rxn = Reaction(row["met_id"]+"_b_biomass"+model_tag)
            rxn.add_metabolites({METB:-1})
            model.add_reaction(rxn)

    #Set Vv and Vc equal to vacuolar and cytoplasmic volumes at the end of the time period (like Sanu did)
    # This wrong. Changed to change in vol
    Vv = rootslice.loc[model_tag, "change_vac_dnut_vol"]
    Vc = rootslice.loc[model_tag, "change_cytoplasm_dnut_vol"]
    Vcell = rootslice.loc[model_tag, "change_tot_dnut_vol"]

    #sets ratio based on eqn 2 in paper
    rxn = Reaction("VacCytRatio" + model_tag)
    WO = Metabolite("TotalSolute" + model_tag, charge = 0, compartment = "c"+model_tag)    #totalsolute placed in cytosol for now -- seems best place for it -- removed from compartment!
    rxn.add_metabolites({CO:-1,VO:-1*(Vv/Vc),WO:1+(Vv/Vc)})
    rxn.lower_bound = 0
    rxn.upper_bound = 1000
    model.add_reaction(rxn)

    C_cell = 250                      #Bloom et al. 2012; units = mol/m3
    C_cell = 250 * (10**-15)          #units = mmol/um3

    #based on eqn 1 in paper. constrains reactions based on overall change in volume that you'll get from rootslice
    rxn = Reaction("TotalSoluteConstraint" + model_tag)
    rxn.add_metabolites({WO:-1})
    rxn.lower_bound = float(C_cell*Vcell/gDW_per_tissue) #mmol/gDW
    rxn.upper_bound = float(C_cell*Vcell/gDW_per_tissue) #mmol/gDW
    model.add_reaction(rxn)

    #print reactions to check outputs
    #print(model.reactions.get_by_id("TotalSoluteConstraint" + model_tag))
    #print(model.reactions.get_by_id("soluble_biomass" + model_tag))

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
        model.add_reaction(rxn)

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
            model.add_reaction(reaction)
            #print(reaction)

        #if uptake cost
        for i in ["Pi", "KI", "NITRATE"]:
            #laoding into apoplast
            ion_ap = model.metabolites.get_by_id(f"{i}_{comp}")
            ion = f"{i}_c_{tag}"
            ion_met = model.metabolites.get_by_id(ion)
            reaction = Reaction(ion + "_apoplast_loading_" + comp)
            model.add_reaction(reaction)
            reaction.name = ion + " loading into apoplast"
            reaction.add_metabolites({
            ion_met: -1.0,
            ion_ap: 1.0,
            })
            reaction.lower_bound = 0
            reaction.upper_bound = 1000

            #unloading from apoplast
            ion_ap = model.metabolites.get_by_id(f"{i}_{comp}")
            ion = f"{i}_c_{tag}"
            ion_met = model.metabolites.get_by_id(ion)
            reaction = Reaction(ion + "_apoplast_unloading_" + comp)
            model.add_reaction(reaction)
            reaction.name = ion + " unloading from apoplast"
            reaction.add_metabolites({
            ion_ap: -1.0,
            model.metabolites.get_by_id("PROTON_"+comp): -1.0,
            ion_met: 1.0,
            model.metabolites.get_by_id("PROTON_c_"+tag): 1.0
            })
            reaction.lower_bound = 0
            reaction.upper_bound = 1000

# Add pFBA obj to model that doesn't include...
def rootspFBAobj(tissue_model,exclude_from_pFBA=[]):
    #exclude sympastic transfer reactions and apoplast exchange rxns from pFBA minimisation
    for rxn in tissue_model.reactions:
        if "Transfer" in rxn.id:
            exclude_from_pFBA+=[rxn.id] #therefore these fluxes won't be minimised
        elif "apoplast_exchange" in rxn.id:
            exclude_from_pFBA+=[rxn.id] #therefore these fluxes won't be minimised
    #exclude_from_pFBA

    Dcoeffs = {}
    sum_var = tissue_model.problem.Variable('sumvar')
    tissue_model.add_cons_vars(sum_var)
    for ii in [rxn.id for rxn in tissue_model.reactions if rxn.id not in exclude_from_pFBA]:
        Dcoeffs[tissue_model.reactions.get_by_id(ii).forward_variable]=1
        Dcoeffs[tissue_model.reactions.get_by_id(ii).reverse_variable]=1
    Dcoeffs[sum_var]=-1
    Dcons = tissue_model.problem.Constraint(0,
        lb=0,
        ub=0)
    tissue_model.add_cons_vars([Dcons])
    tissue_model.solver.update()
    Dcons.set_linear_coefficients(coefficients=Dcoeffs)
    pfbaobj=tissue_model.problem.Objective(sum_var,
            direction='min')
    tissue_model.objective = pfbaobj