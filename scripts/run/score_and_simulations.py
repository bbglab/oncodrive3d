"""
Contains function to assign clustering anomaly score and perform simulations
"""

import daiquiri
import numpy as np
import pandas as pd
from scipy import stats
from decimal import Decimal, getcontext
from functools import reduce
import operator

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".run.score_and_simulations")


def dcm_factorial(n):
    """
    Compute factorial.
    """
    return reduce(operator.mul, [Decimal(i) for i in range(1, int(n)+1)], Decimal(1))


def dcm_binom_coeff(n, k):
    """
    Compute binomial coefficient.
    """
    
    return dcm_factorial(n) / (dcm_factorial(k) * dcm_factorial(n - k))


def dcm_binom_cdf(k, n, p):
    """
    Compute binomial cumulative distribution function (CDF).
    """
    
    p = Decimal(p)
    q = Decimal(1) - p
    cdf = Decimal(0)
    for i in range(int(k) + 1):
        cdf += dcm_binom_coeff(n, i) * (p ** i) * (q ** (n - i))
        
    return cdf


def dcm_binom_sf(k, n, p):
    """
    Compute binomial survival function (SF).
    """
    
    return Decimal('1') - dcm_binom_cdf(k, n, p)


def dcm_binom_logsf(k, n, p):
    """
    Compute log binomial survival function.
    """
    
    sf = dcm_binom_sf(k, n, p)
    if sf <= 0:
        return np.inf
    return sf.ln()


def get_dcm_anomaly_score(k, n, p, decimal=600):          
    """
    Use the decimal package to compute the anomaly score 
    with high precision to avoid approximation of the 
    numerator.
    
    Score: loglik equal or larger mut_count / loglik(N)
    """

    getcontext().prec = decimal
    num = dcm_binom_logsf(k-1, n, p)
    den = stats.binom.logpmf(k=n, n=n, p=p)

    return float(num / Decimal(den))


def recompute_inf_score(result_pos_df, gene_mut, vol_missense_mut_prob):
    """
    Use high precision calculation to recompute the score that 
    were approximated to inf by scipy.
    
    The issue happens in extreme cases when the numerator of the score
    is so small that is approximated to 0.
    """
    
    inf_ix = np.isinf(result_pos_df.Score)
    if sum(inf_ix) > 0:
        for ix, k, n, p in zip(np.where(inf_ix)[0], 
                                result_pos_df.Mut_in_vol[inf_ix], 
                                np.repeat(gene_mut, sum(inf_ix)), 
                                vol_missense_mut_prob[inf_ix]):
            
            if np.isinf(result_pos_df.iloc[ix].Score):
                result_pos_df.loc[ix, "Score"] = get_dcm_anomaly_score(k, n, p)
            else:
                logger.warning("Trying to overwrite a non-inf score: Skipping..")
                
    return result_pos_df


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


def simulate_mutations(n_mutations, p, size, seed=None):
    """
    Simulate the mutations given the mutation rate of a cohort.
    """

    rng = np.random.default_rng(seed=seed)
    samples = rng.multinomial(n_mutations, p, size=size)
    
    return samples


def get_sim_anomaly_score(mut_count, 
                          cmap,
                          gene_miss_prob,
                          vol_missense_mut_prob,
                          num_iteration=1000,
                          seed=None):
    """
    Simulated mutations following the mutation profile of the cohort.
    Compute the log-likelihood of observing k or more mutation in the 
    volume and compare it with the corresponding simualted rank.
    """
    
    # Generate x sets of random mutation distributed following the mut 
    # profile of the cohort, each with the same size of the observed mut   
    mut_sim = simulate_mutations(mut_count, gene_miss_prob, num_iteration, seed)
    
    # Get the density of mutations of each position in each iteration
    density_sim = np.einsum('ij,jk->ki', cmap, mut_sim.T.astype(float), optimize=True)
    
    # Compute the ranked score of the densities obtained at each iteration
    # sign is used to sort in descending order
    loglik_plus = -np.sort(-get_anomaly_score(density_sim, mut_count, vol_missense_mut_prob))
    
    return pd.DataFrame(loglik_plus).T