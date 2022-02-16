from cobra import io
from cobra.flux_analysis.variability import *
from cobra.flux_analysis.parsimonious import pfba
import pandas as pd

def main(modelfile='final_roots_model.xml',path='./',savepath='./'):
    # If no preexisting files, generate from stem model
    # stem_modelm,vfva_sol=loadModel(modelpath)
    # Dm, Dstdm = loadTranscriptome(stem_modelm,transcriptome_data,gene_mapping,minstd)
    
    stem_modelm=io.read_sbml_model(path+modelfile)
    vfva_sol=flux_variability_analysis(stem_modelm)
    solmax = pfba(stem_modelm)
    data={'reaction_id':[x.id for x in stem_modelm.reactions],
        'cell':['_'.join(x.id.split('_')[-2:]) if ('_l' in x.id or '_d' in x.id) else '' for x in stem_modelm.reactions],
        'subsystem':[x.notes['SUBSYSTEM']  if 'SUBSYSTEM' in str(x.notes) else '' for x in stem_modelm.reactions],
        'reaction':[x.reaction for x in stem_modelm.reactions],
        'dielmax':[solmax[x.id] for x in stem_modelm.reactions],
        'pfva_min':[vfva_sol.minimum[x.id] for x in stem_modelm.reactions],
        'pfva_max':[vfva_sol.maximum[x.id] for x in stem_modelm.reactions]}
    df=pd.DataFrame(data)
    from datetime import datetime, date, time
    now = datetime.now().strftime('%Y_%m_%d')
    df.to_csv(savepath+'FVA'+modelfile.split('.')[0]+now+'.csv')