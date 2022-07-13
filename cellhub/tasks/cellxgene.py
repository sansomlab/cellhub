import pandas as pd
import numpy as np

def get_range(x):
    
    return np.max(x) - np.min(x)

def facet_layout(adata, 
                 layout="X_umap",
                 name=None,
                 x_factor=None, 
                 y_factor=None, 
                 x_levels=None, 
                 y_levels=None):
    '''
        Facet a UMAP by given factors. The factors are assumed to be in the correct order.
    '''
    
    if name is None:
        name = layout + "_faceted"
    
    x = pd.DataFrame(adata.obsm[layout]).copy()
    x.index = adata.obs.index
    x.columns = ["x","y"]
    
    x_range = get_range(x["x"]) * 1.05 # add 5% space between facets
    y_range = get_range(x["y"]) * 1.05 
    

    if not x_factor is None:
        
        if x_levels is None:    
            x_levels = set(adata.obs[x_factor])
        else:
            x_levels = [x.strip() for x in x_levels.split(",")]
            
        # for each level shift
        m = 0
        for xl in x_levels:
            
            
            cells = adata.obs.index[adata.obs[x_factor]==xl]
            
            
            x.loc[cells, "x"] = x.loc[cells, "x"] + x_range * m 
             
            m += 1
    
    if not y_factor is None:
        
        if y_levels is None:    
            y_levels = set(adata.obs[y_factor])
        else:
            y_levels = [x.strip() for x in y_levels.split(",")]
            
        # for each level shift
        m = 0
        for yl in y_levels:
            
            cells = adata.obs.index[adata.obs[y_factor]==yl]
            
            x.loc[cells, "y"] = x.loc[cells, "y"] + y_range * m 
            
            m += 1
            
    adata.obsm[name] = np.array(x)
    
    return adata
    
def clip(x, lq = 0.01, uq = 0.9):
    #print(len(x))
    #print(x[1:10])
    #l = np.quantile(x,lq)
    u = np.quantile(x,uq)
    #x[x<l] = l
    x[x>u] = u
    return x