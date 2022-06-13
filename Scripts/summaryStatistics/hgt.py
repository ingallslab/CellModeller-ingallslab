import numpy as np

def count_hgt_events(cells):
    """
    Sums the cumulative HGT events
    
    @param  cells           cellStates dict
    @return sum_hgt_events  total HGT events within the cellStates dict
    """
    sum_hgt_events = 0
    for cell in cells.values():
        sum_hgt_events += cell.hgt_events
        
    return sum_hgt_events      
   
def calc_hgt_efficiency(populations, recip_types=[1], trans_types=[2]):
    """
    Calculates conjugation efficiency as transconjugants per total recipient pool
    
    @param  populations total population counts
    @param  recips      recipient cellTypes
    @param  trans       transconjugant cellTypes
    @return 
    """
    recips = np.sum(populations[recip_types])
    trans = np.sum(populations[trans_types])
    hgt_efficiency = trans/(recips + trans)
    
    return hgt_efficiency
