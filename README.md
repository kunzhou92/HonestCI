# Honest confidence sets for high-dimensional regression by projection and shrinkage

This repository contains the code for the paper [Honest confidence sets for high-dimensional regression by projection and shrinkage](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1938581). 

**Abstract:** The issue of honesty in constructing confidence sets arises in nonparametric regression. While optimal rate in nonparametric estimation can be achieved and utilized to construct sharp confidence sets, severe degradation of confidence level often happens after estimating the degree of smoothness. Similarly, for high-dimensional regression, oracle inequalities for sparse estimators could be utilized to construct sharp confidence sets. Yet the degree of sparsity itself is unknown and needs to be estimated, causing the honesty problem. To resolve this issue, we develop a novel method to construct honest confidence sets for sparse high-dimensional linear regression. The key idea in our construction is to separate signals into a strong and a weak group, and then construct confidence sets for each group separately. This is achieved by a projection and shrinkage approach, the latter implemented via Stein estimation and the associated Stein unbiased risk estimate. Our confidence set is honest over the full parameter space without any sparsity constraints, while its diameter adapts to the optimal rate of $n^{-1/4}$ when the true parameter is indeed sparse. Through extensive numerical comparisons, we demonstrate that our method outperforms other competitors with big margins for finite samples, including oracle methods built upon the true sparsity of the underlying model.


Here are 8 simulation scenarios corresponding 8 folders under mainline. Please check #5-8 for the simulation presented in the paper.
* `1_n100p400` contains the initial implementation for the simulation study.
*  `2_n100p400_test_differenet_sigma` has the same setup as `1_n100p400` except the variance of noise is set to a smaller number.
*  `3_n100p400_test2_differenet_sigma` has the same setup as `1_n100p400` except the variance of noise is set to a larger number.
*  `4_n200p800` is the setup higher-dimensional scenario.
*  `5_n200p800_new_c_alpha` updates the computation logic of $c_(\alpha)$.
*  `6_n200p800_unknown_sigma` is the setup where `\sigma` is unknown and is estimated.
*  `7_n200p800_against_normality_homogeneity` is the setup where the variance of noise ($\sigma$) is correlated to coefficients ($\beta$) 
*  `8_real_data_simulation` is the setup for real-data simulation.

For each folder, please check `implement.R` for how to create/process data and implement our method. Check `summary_figures.ipynb` and `summary_paper.ipynb` for visualization of experiment results.

**If you have questions**, please reach out to Kun Zhou (k.zhou@ucla.edu) and Qing Zhou (zhou@stat.ucla.edu).

