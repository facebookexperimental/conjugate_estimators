# Copyright (c) Meta Platforms, Inc. and affiliates.

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly import express as px
from scipy.stats import beta, gamma


##### Initialize constants

# fix the random seed for reproducibility
np.random.seed(0)

# How many rows are in the population table
population_size = 1000

# What % of population rows will get sampled
sample_rate = 0.01

# The number of population rows that get sampled
sample_size = int(population_size * sample_rate)
sample_size = 100

# The beta parameters that are used for sampling ground truth labels
# for each row.
beta_prior_negative = 50
beta_prior_positive = 10


###### Create the population

# Sample the p(positive) score for each population row
pop_positive_probs = np.random.beta(
    beta_prior_positive, beta_prior_negative, size=(population_size)
)

# Sample a ground truth label for each population row
pop_labels = np.random.rand((pop_positive_probs.size)) < pop_positive_probs
pop_mean = pop_labels.mean()
print("Population mean label:", pop_mean)

num_trials = 100

# Add some noise to the probabilities to simulate importance sampling with an imperfect classifier
noise_mean = 0
noise_std_dev = 0.1
classifier_preds = (pop_positive_probs) + np.random.normal(
    noise_mean, noise_std_dev, pop_positive_probs.shape
)
classifier_preds = np.clip(classifier_preds, 0, 1) + 1e-9


# take the sqrt because that's what we do today
# TODO verify necessity
sample_weights = np.sqrt(classifier_preds)

pop_sampling_probs = sample_weights / sample_weights.sum()


sampled_indices = np.random.choice(
    np.arange(population_size),
    size=(num_trials, sample_size),
    replace=True,
    p=pop_sampling_probs,
)

# the probability each sample is picked
sample_probs = pop_sampling_probs[sampled_indices]
sample_labels = pop_labels[sampled_indices]
sample_rep_weights = (sample_size / sample_probs) / (sample_size / sample_probs).sum(
    -1, keepdims=True
)


# This function implements the conjugate beta estimator algorithm.
# Its forumala is
#   alpha = mu*n + alpha_prior
#   beta = (1-mu)*n + beta_prior
#   ci = [ppf(0.05, alpha, beta), ppf(0.95, alpha, beta]
# where mu is the mean of the (weighted) sample labels and n is the sample size in bits
# (see below).
def cbe(sample_mean, sample_size):
    alpha_param = sample_mean * sample_size + 1e-9  # optional: add an alpha prior
    beta_param = (1 - sample_mean) * sample_size + 1e-9  # optional: add a beta prior
    lower_bound = beta.ppf(0.025, alpha_param, beta_param)
    median = beta.ppf(0.5, alpha_param, beta_param)
    upper_bound = beta.ppf(0.975, alpha_param, beta_param)
    return lower_bound.mean(), median.mean(), upper_bound.mean()


### Perfect labels section
sample_mean = (sample_labels * sample_rep_weights).sum(-1)
res = cbe(sample_mean, sample_size)
print("CI with perfect labels:", res)


### Noisy labels section

# Generate random noise
rns = np.random.uniform(0, 1, size=sample_labels.shape)
noisy_labels = sample_labels.copy().astype(int)

sensitivity = 0.85
specificity = 0.93

# Add random noise to negative samples according to sensitivity
noisy_labels[sample_labels == 0] += rns[sample_labels == 0] > sensitivity

# Add random noise to positive samples according to specificity
noisy_labels[sample_labels == 1] = (
    noisy_labels[sample_labels == 1]
    + (rns[sample_labels == 1] > specificity).astype(int)
) % 2


# Rogan Gladen estimator for the sample mean. See: https://en.wikipedia.org/wiki/Beth_Gladen
# To derive it, we start with the following equation, where o is the observed mean, p is the true mean,
# and tpr and fpr are the true positive and false positive rates.
#   o = tpr*p + fpr*(1-p)
#   o = tpr*p + fpr - fpr*p
#   p*(tpr - fpr) = o - fpr
#   p = (o - fpr)/(tpr - fpr)
def rg(mean):
    return (mean + specificity - 1) / (sensitivity + specificity - 1)


noisy_mean = (noisy_labels * sample_rep_weights).sum(-1)
print("Noisy mean:", noisy_mean.mean())
res = cbe(noisy_mean, sample_size)
print("CI with noisy labels:", res)

rg_mean = rg((noisy_mean))
print("RG sample mean:", rg_mean.mean())

res = cbe(rg_mean, sample_size)
print("CI with RG mean and standard sample size", res)


def entropy(val):
    return -val * np.log2(val) - (1 - val) * np.log2(1 - val)


# When the mean of sensitivity and specificity approachs 0.5, the denominator
# in the Rogan Gladen formula, (mean + specificity - 1)/(sensitivity + specificity - 1)
# approaches 0. This causes the RG mean estimate to be highly unstable to small
# fluctuations in sensitivity and/or specificity. To counter this instability, we want
# to widen the CI. The entropy in this case is 1 (because entropy(0.5) = 1).
# In the formula below, an entropy value of 1 leads to 0 bits per sample, which
# amounts to having no samples when those samples are basically random noise.
# However, as the mean of sensitivity and specificity approaches 1, the entropy
# approaches 0, leading to about 1 bit per sample, which is the behavior
# of the original CBE algorithm.
num_bits = (1 - entropy((sensitivity + specificity) / 2)) * sample_size
res = cbe(rg_mean, num_bits)
print("CI with RG mean and entropy weighted sample size:", res)
