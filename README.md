
# The conjugate beta estimator statistical algorithm

This repo contains a reference implementation for a statistical algorithm called
Conjugate Beta Estimator (CBE) for computing
CIs for population means using (weighted) sample means and
potentially noisy labels.

The basic formula for CBE is

    alpha = mu*n + alpha_prior
    beta = (1-mu)*n + beta_prior
    ci = [ppf(0.05, alpha, beta), ppf(0.95, alpha, beta)

where mu is the mean of the (weighted) sample labels and n is the sample size in bits

Both mu and n can be adjusted to account for label noise.

mu should be adjusted using the [Rogan Gladen](https://en.wikipedia.org/wiki/Beth_Gladen) (RG) estimator for the sample mean:

    rg(mu, sensitivity, specificity) = (mu + specificity - 1) / (sensitivity + specificity - 1)

n should be adjusted using the following formula:

    num_bits_per_label = (1 - entropy((sensitivity + specificity) / 2))
    n_modified = num_bits_per_label * n

The reason for the num_bits_per_label formula is that the rg formula is increasingly unstable when
the mean of sensitivity and specificity approaches 0.5 (the max entropy value) and the
rg denominator approaches 0. The rg formula
is maximally stable when sensitivity = specificity = 1, which is the case when labels are perfectly
accurate. Therefore, we want the CI derived from the Beta distribution to grow wider as
(sensitivity + specificity)/2 approeaches 0.5 from 0 or 1.


See the [CONTRIBUTING](CONTRIBUTING.md) file for how to help out.

## License
Conjugate Estimators is MIT licensed, as found in the LICENSE file.
