
# The conjugate Beta estimator statistical algorithm


## Introduction
This repo contains a reference implementation for a statistical algorithm called
Conjugate Beta Estimator (CBE) for computing
confidence intervals (CIs) for population means using (weighted) sample means and
potentially noisy labels.

## Basic algorithm
The basic formula for CBE is

    alpha = mu*n + alpha_prior
    beta = (1-mu)*n + beta_prior
    ci = [ppf(0.05, alpha, beta), ppf(0.95, alpha, beta)]

where mu is the mean of the (weighted) sample labels, n is the sample size in bits, and ppf is the
[inverse CDF](https://en.wikipedia.org/wiki/Quantile_function) function.

CBE is an unbiased estimator of the population mean because CBE's mean is alpha/(alpha + beta).
As n->inf, alpha_prior and beta_prior wash out, and the expression simplifies to mu.

The idea behind CBE -- a distribution whose mean converges to the population mean and whose
variance shrinks with the number of labels -- can be generalized to non-Beta distributions
if you need a different support than [0,1]. For example, the Gamme distribution can be parameterized as

    Gamma(alpha=1+(n*mu), beta=1+n)

The Gaussian distribution can be parameterized as

    Gaussian(mean=(prior + n*mu)/(n+1), var=1/(n+1))


Consider using the [Kumaraswamy](https://en.wikipedia.org/wiki/Kumaraswamy_distribution) distribution
instead of Beta for a more efficient analytical inverse CDF function.

CBE can be applied in streaming applications, where mu can represent a moving average.

## Accounting for label noise
Both mu and n can be adjusted to account for label noise.

mu should be adjusted using the [Rogan Gladen](https://en.wikipedia.org/wiki/Beth_Gladen) (RG) estimator for the sample mean:

    mu_modified = (mu + specificity - 1) / (sensitivity + specificity - 1)

n should be adjusted using the following formula:

    num_bits_per_label = (1 - entropy((sensitivity + specificity) / 2))
    n_modified = num_bits_per_label * n

the accuracy adjusted CBE formula is therefore

    alpha_modified = mu_modified*n_modified + alpha_prior
    beta_modified = (1-mu_modified)*n_modified + beta_prior
    ci = [ppf(0.05, alpha_modified, beta_modified), ppf(0.95, alpha_modified, beta_modified)

The reason for the num_bits_per_label formula is that the rg formula is increasingly unstable when
the mean of sensitivity and specificity approaches 0.5 (the max entropy value) and the
rg denominator approaches 0. The rg formula
is maximally stable when sensitivity = specificity = 1, which is the case when labels are perfectly
accurate. Therefore, we want the CI derived from the Beta distribution to grow wider as
(sensitivity + specificity)/2 approaches 0.5 from 0 or 1.

The relationship between average accuracy and # of bits per label is visualized in the graph below.

<img width="695" alt="Screenshot 2024-11-05 at 12 44 57 PM" src="https://github.com/user-attachments/assets/975f7141-6ed6-4327-9035-052b419fbc51">

## Mathematical derivation

The following derivation shows that the 1-H(X) formula, where H is entropy and X is (sensitivity + specificity)/2, can be expressed as the Bayesian information gain, or KL divergence
between the posterior P and uniform prior Q:

$$D_{\text{KL}}(P \parallel Q) = p * \log_2\left(\frac{p}{q}\right) + (1-p) *\log_2\left(\frac{1-p}{1-q}\right) =$$
$$p *\log_2\left(\frac{p}{0.5}\right) + (1-p) * \log_2\left(\frac{1-p}{0.5}\right) =$$
$$p *\left(\log_2(p) - \log_2(0.5)\right) + (1-p) * \left(\log_2(1-p) - \log_2(0.5)\right) =$$
$$p *\log_2(p) + p + (1-p) *\log_2(1-p) + (1-p) =$$
$$p *\log_2(p) + (1-p) * \log_2(1-p) + 1 = $$
$$ 1 - H(X) $$

## Future Directions

If you want to introduce uncertainty about sensitivity and specificity, it might be possible to mulitply
num_bits_per_label by another 1-H(X) term, where X is the entropy of the sensitivity or specificity distribution.
The intuition is that an infinite number of labels whose accuracy is unknown (uniformly distributed on [0,1])
are as informative as zero labels.

## License
Conjugate Estimators is MIT licensed, as found in the LICENSE file.
