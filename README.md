## Complexity with Placebo Returns 
I simulate ($i.i.d.$) Returns from a Normal $\hat\mu(R_t) \hat\sigma(R_t)$, where $\hat\mu(R_t)$ and $\hat\sigma(R_t)$ are the sample mean and standard deviation of the empirical standardized excess returns $R_t$. 

I then use the simulated $R_t$ series with the empirical GW predictor data and KMZ's complexity code to generate VoC curves based on the $i.i.d.$ simulated data to show that the same VoC pattern arises. 

## RNG Seed
The code sets the seed to 1000+1, seeds from 1 to 1000 are used in the RFF weight generation, so it's better to not generate the Normal draws using the same seed. 

Run the code in order:
1) o1_ --> simulates the placebo returns and estimates RFF and GW regressions
2) o2_ --> combines the RFF output and creates portfolio measures
3) o3_ --> plots the results

- the file `Y_placebo_seed_1001.mat` contains the simulated/fake Y returns data.

- the main function files are in the directory: `local.functions`. 

### Function files have been modified for the placebo data to be used are:

- `KMZ_GW_benchmark_function_placebo_Ret.m`
- `KMZ_tryrff_v2_function_for_each_sim_placebo_Ret.m`


