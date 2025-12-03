## Complexity with Placebo Returns 
I simulate ($i.i.d.$) Returns from a Normal($\hat\mu(R_t),\hat\sigma(R_t)$), where $\hat\mu(R_t)$ and $\hat\sigma(R_t)$ are the sample mean and standard deviation of the empirical standardized excess returns $R_t$. 

I then use the simulated $R_t$ series with the empirical GW predictor data and KMZ's complexity code to generate VoC curves based on the $i.i.d.$ simulated data to show that the same VoC pattern arises. 

### RNG Seed
The code sets the seed to 1000+1, seeds from 1 to 1000 are used in the RFF weight generation, so it's better to not generate the Normal draws using the same seed. 