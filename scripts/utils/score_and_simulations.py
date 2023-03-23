"""
Contains function to assign clustering anomaly score and perform simulations
"""

import numpy as np
import pandas as pd
from scipy import stats


def get_anomaly_score(vec_mut_in_vol, gene_mut, vec_vol_miss_mut_prob):          
    """
    Compute a metric that scores the anomaly of observing a certain 
    number of mutations in the volume of a residue.
    It takes into account the volume and the mutation rate of the codon 
    of each residue within that volume.
    
    Score: loglik equal or larger mut_count / loglik(N)
    """
    
    den = stats.binom.logpmf(k=gene_mut, n=gene_mut, p=vec_vol_miss_mut_prob)

    return stats.binom.logsf(k=vec_mut_in_vol-1, n=gene_mut, p=vec_vol_miss_mut_prob) / den


def simulate_mutations(n_mutations, p, size):
    """
    Simulate the mutations given the mutation rate of a cohort.
    """

    rng = np.random.default_rng()
    samples = rng.multinomial(n_mutations, p, size=size)
    
    return samples


def get_sim_anomaly_score(mut_count, 
                            cmap,
                            gene_miss_prob,
                            vol_missense_mut_prob,
                            num_iteration=1000):
    """
    Simulated mutations following the mutation profile of the cohort.
    Compute the log-likelihood of observing k or more mutation in the 
    volume and compare it with the corresponding simualted rank.
    """
    
    # Generate x sets of random mutation distributed following the mut 
    # profile of the cohort, each with the same size of the observed mut   
    mut_sim = simulate_mutations(mut_count, gene_miss_prob, num_iteration)
    
    # Get the density of mutations of each position in each iteration
    density_sim = np.einsum('ij,jk->ki', cmap, mut_sim.T.astype(float), optimize=True)
    
    # Compute the ranked score of the densities obtained at each iteration
    # sign is used to sort in descending order
    loglik_plus = -np.sort(-get_anomaly_score(density_sim, mut_count, vol_missense_mut_prob))
    
    return pd.DataFrame(loglik_plus).T