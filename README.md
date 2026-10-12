# When Does the Central Limit Theorem Kick In?

**Math 167R, Project 1: a Monte Carlo investigation of the CLT** · Dhruv Shah

> **Is the number 30?** There is no single sample size that guarantees an approximately normal sampling distribution. The required $n$ depends on the population’s skewness, discreteness, dependence, and variance. Therefore, $n=30$ is only a rough rule of thumb, not a universal answer.

📄 **Full report:** [`report/clt_project_1.md`](report/clt_project_1.md) (readable on GitHub, with every figure) · [`report/clt_project_1.html`](report/clt_project_1.html) (self-contained HTML; download and open in a browser)

---

## What I did

For each population I generated **10,000 independent samples** at each of **11 sample sizes** (n = 2, 5, 10, 20, 30, 50, 100, 250, 500, 1000, 2000), computed every sample mean, and asked whether the resulting sampling distribution was approximately normal.

**Prediction first.** For i.i.d. data the skewness of $\bar X$ is $\gamma_1/\sqrt{n}$ and its excess kurtosis is $\gamma_2/n$. That predicts how large n must be *before* simulating anything, and the simulation then tests the prediction.

**My definition of "approximately normal".** It was fixed in advance and applied identically to every population. All four checks must pass:

| Check | Rule |
|---|---|
| Symmetry | \|skewness of the means\| ≤ 0.20 |
| Tails | \|excess kurtosis of the means\| ≤ 0.30 |
| Spread | simulated SD / theoretical SE between 0.97 and 1.03 |
| Usefulness | 4%–6% of means fall outside μ ± 1.96 SE (a normal gives 5%) |

Each population also gets histograms with the CLT normal curve, normal Q–Q plots, a theory-vs-simulation table, and a Shapiro–Wilk test (reported, but not used for the verdict; see the report for why).

## Final table

<!-- SUMMARY_TABLE_START -->
|Population                        |CLT assumptions                |Predicted n |Observed n |Skew @ n=30 |SD/SE @ n=30 |Tail rate @ n=30 |Normal @ n=30? |
|:---------------------------------|:------------------------------|:-----------|:----------|:-----------|:------------|:----------------|:--------------|
|Standard Normal                   |None (i.i.d., finite variance) |1           |2          |-0.01       |1.01         |5.2%             |Yes            |
|Six-Sided Die                     |None (i.i.d., finite variance) |5           |10         |0.00        |1.00         |4.8%             |Yes            |
|Exponential (rate = 1)            |None (i.i.d., finite variance) |100         |100        |0.40        |0.99         |4.8%             |No             |
|Binomial (5, 0.10)                |None (i.i.d., finite variance) |36          |50         |0.21        |1.00         |4.4%             |No             |
|Cauchy (0, 1)                     |Breaks: finite variance        |n/a         |> 2000     |undefined   |undefined    |undefined        |No             |
|Dependent Machine Failure         |Breaks: independent            |n/a         |1000       |0.75        |1.01         |4.4%             |No             |
|Non-identical Bernoulli (Fatigue) |Breaks: identical              |15          |100        |-0.14       |0.99         |5.9%             |Yes            |
|Age at Death (Mixture)            |None (i.i.d., finite variance) |19          |50         |-0.20       |1.01         |5.1%             |No             |
<!-- SUMMARY_TABLE_END -->

*Predicted n*: from theory, before simulating. *Observed n*: the smallest simulated n from which all four checks pass and keep passing. *SD/SE*: simulated spread divided by the correct theoretical standard error. *Tail rate*: share of sample means beyond ±1.96 SE.

## Conclusions

**Does the CLT get the center and spread right at $n=30$?**
Yes. For every finite-variance population, the simulated standard deviation was very close to the theoretical standard error. The main differences at $n=30$ were in the shape, especially skewness and tail behavior.

1. **What does the skewness table show at $n=30$? Did it match the prediction?**
   Yes, the results mostly matched the prediction that skewness decreases as $\gamma_1/\sqrt n$. For example, the exponential distribution had predicted skewness about 0.37 and simulated skewness about 0.40. The binomial distribution had predicted skewness about 0.22 and simulated skewness about 0.21.
2. **Does skewness determine how quickly the sampling distribution becomes normal?**
   Generally, yes. The symmetric normal and die distributions became approximately normal quickly. The more skewed exponential and binomial distributions needed larger sample sizes. However, skewness is not the only factor because discreteness and dependence also affect convergence.
3. **What happens with discrete populations?**
   Discrete populations can have sampling distributions that are close to normal in spread but still show irregular tail probabilities. This happened with the binomial and fatigue models because their sample means can only take certain discrete values.
4. **What happens with the Cauchy population?**
   The Cauchy sample mean does not become approximately normal, even when $n$ increases to 2,000. The Cauchy distribution has no finite variance, so the usual Central Limit Theorem does not apply.
5. **What happens when observations are dependent?**
   Dependence makes the usual standard-error formula incorrect. For the machine-failure process, negative autocorrelation reduced the true standard error, so the formula $\sigma/\sqrt n$ overestimated the spread. The distribution also remained skewed for much longer.
6. **Can non-identical but independent observations still follow a CLT?**
   Yes. The fatigue model became approximately normal even though its success probability changed between observations. However, its standard error had to be calculated using the individual probabilities rather than simply using the average probability.
7. **Was $n=30$ enough for every population?**
   No. It was enough for the normal, die, and fatigue populations, but not for the exponential, binomial, Cauchy, dependent machine-failure, or age-at-death populations.
8. **What is the overall conclusion?**
   There is no single sample size that guarantees an approximately normal sampling distribution. The required $n$ depends on the population’s skewness, discreteness, dependence, and variance. Therefore, $n=30$ is only a rough rule of thumb, not a universal answer.

## Populations studied

| # | Population | Why it is interesting | Required / extra credit |
|---|---|---|---|
| 1 | Standard Normal | Control: the mean is exactly normal at every n | Required |
| 2 | Six-Sided Die | Symmetric, discrete, flat | Required |
| 3 | Exponential (λ = 1) | Strong right skew | Required |
| 4 | Binomial (5, 0.10) | Highly discrete and skewed | Required |
| 5 | Cauchy (0, 1) | No finite mean or variance | Required |
| 6 | Dependent machine failure | Breaks **independence** | Extra credit |
| 7 | Non-identical Bernoulli (fatigue) | Breaks **identical distribution** | Extra credit |
| 8 | Age at death (mixture) | Bimodal, real-world-shaped | Extra credit |

## Reproduce it

Requirements: R (≥ 4.1, base R only, no extra packages for the analysis), `knitr` and `rmarkdown`, and [Quarto](https://quarto.org) (≥ 1.4).

```bash
git clone https://github.com/Dhruv-cs50/clt-monte-carlo-project1.git
cd clt-monte-carlo-project1
./render.sh          # about 6 minutes; regenerates every number, table and figure
```

The seed is the last four digits of my student ID (**6877**), set once in [`R/config.R`](R/config.R). Each population uses `SEED + k`, so any section reproduces on its own. Every setting (replications, sample sizes, normality thresholds) is in that same file.

## Repository layout

```
├── README.md              ← this executive summary
├── AI_PROMPTS.md          ← AI prompts used, and how the output was checked and corrected
├── render.sh              ← one command to reproduce everything
├── R/
│   ├── config.R           ← seed, replications, n grid, normality thresholds
│   ├── populations.R      ← the 8 populations + their theoretical moments
│   ├── simulate.R         ← Monte Carlo engine and normality diagnostics
│   └── plots.R            ← reusable plotting functions
├── report/
│   ├── clt_project_1.qmd  ← the full analysis (Quarto source)
│   ├── clt_project_1.md   ← rendered for GitHub, with figures
│   └── clt_project_1.html ← rendered, self-contained HTML
└── output/
    ├── summary_table.csv  ← the final table above
    └── all_diagnostics.csv← every statistic for every population × n
```

## AI use

This project was built with help from an AI assistant (Claude). The prompts, the design decisions, how every result was checked against theory, and the errors that were caught and corrected are documented in [`AI_PROMPTS.md`](AI_PROMPTS.md).
