# AI usage log

The course asks for a record of the AI prompts used, the follow-up prompts, and how the resulting code was checked or corrected. This project was built with Claude (Anthropic), working as a coding assistant in an R + Quarto environment.

## 1. Prompts

### Initial prompt

> Using all this information build the project for me as asked. Also make a GitHub repo for the same and make good commits for the same too.

Attached with the prompt were all of the course materials: the Project 1 assignment text, `Project_1_Student_Guidance.qmd` (with the reduced 9/22 population list), the Class 05 CLT lesson, the Class 06 statistical distributions lesson, and the three extra-credit handouts (dependent failure process, basketball fatigue, custom age-at-death distribution).

### Follow-up prompts

> My student ID number is 015216877. I made a new repo: https://github.com/Dhruv-cs50/clt-monte-carlo-project1

This set the seed to **6877** (the last four digits, as the guidance requires) and gave the push destination.

## 2. Design decisions made with the AI

| Decision | Choice | Reason |
|---|---|---|
| Populations | All 5 required + all 3 extra credit | Guidance list after the 9/22 update |
| Replications | 10,000 per sample size | Same as the class lessons; enough to resolve skewness to about ±0.025 |
| Sample sizes | 2, 5, 10, 20, 30, 50, 100, 250, 500, 1000, 2000 | Brackets 30, and reaches far enough to see the slow cases converge |
| "Approximately normal" | All four checks must pass: \|skew\| ≤ 0.20, \|excess kurtosis\| ≤ 0.30, SD/SE in [0.97, 1.03], tail rate in [4%, 6%] | Each threshold is several Monte Carlo SEs wide, so noise cannot decide the verdict |
| Prediction | skew(X̄) = γ₁/√n, exkurt(X̄) = γ₂/n | Gives a testable prediction before simulating |
| Dependent process | 50-step burn-in; SE from autocorrelations | Measures the process's long-run behaviour, and uses a correct SE instead of σ/√n |
| Fatigue model | A sample of size n is n shots with p spread from 0.90 to 0.50 | Keeps p̄ = 0.70 at every n, and the exact Poisson-binomial SE is available |
| Shapiro–Wilk | Reported, not used for the verdict | It rejects tiny, irrelevant deviations at 5,000 values |

## 3. How the output was checked

Every piece of generated code was run, and the results were compared with known theory before being trusted:

1. **Control case.** The standard normal must pass every check at every n, because the mean of normals is exactly normal. It did, with skewness and kurtosis within noise of 0 and SD/SE within about 1%. This shows the criteria do not produce false failures.
2. **Theoretical moments.** Population means, SDs, skewness and kurtosis were checked against textbook formulas: die σ = √(35/12); Binomial(5, 0.1) skewness = 0.8/√0.45; exponential skewness 2 and excess kurtosis 6. The age-at-death mean (49.98) and SD (27.63) match the values in the course handout.
3. **Shape vs theory.** Simulated skewness of the sample mean tracks γ₁/√n for every i.i.d. population. For the exponential it was 0.18 at n = 100, against 0.20 predicted.
4. **Cauchy.** Theory says the mean of n Cauchy draws is again Cauchy(0, 1). The simulated IQR of the means stays at about 2 for every n, and the share of means beyond ±10 matches 2·P(X < −10) ≈ 6.3%.
5. **Seed robustness.** The whole study was repeated with five seeds to make sure the conclusions do not depend on one lucky draw.

## 4. Corrections made along the way

These are real issues found while checking the AI's first drafts. Each one was fixed before the final render.

| Problem found | How it was caught | Fix |
|---|---|---|
| The dependent-process SE first used the *asymptotic* long-run variance σ²(1 + 2Σρₖ)/n. At n = 2 the SD/SE ratio came out at 0.94, failing the spread check for the wrong reason. | The SD/SE column was not ≈ 1 at small n | Switched to the exact finite-n formula (σ²/n)[1 + 2Σ(1 − k/n)ρₖ] |
| The n grid originally stopped at 500, where the exponential and dependent cases were still converging | The convergence table had no stable pass at the top of the grid | Extended to n = 1000 and 2000, and processed large n in blocks to limit memory |
| `plot_population()` crashed with `invalid 'xlim' value` for continuous populations | First Quarto render failed | Compute a default x-range when none is given |
| The robustness table's "Seed 6877" column did not reproduce the main results | Age at death showed 30 in one table and 50 in the other | The robustness loop now applies the same per-population seed offset as the main run |
| Several sentences in the first draft hardcoded numbers (e.g. "age at death needs about 20") that did not hold for seed 6877 | Read the rendered text against the result tables | Replaced with inline R values, so the prose always matches the current seed |
| The summary table showed a Cauchy "skewness" of about 99, which is meaningless because the moment does not exist | Reviewed the summary table | Displayed as "undefined", with quantile-based evidence (IQR) used instead |
| The GitHub-markdown figures disappeared after rendering | The `.md` referenced images that no longer existed | HTML's self-contained mode deletes the figure folder, so `render.sh` now renders HTML first and GitHub markdown second |

## 5. What I would verify myself before relying on it

- Re-run `./render.sh` and confirm the numbers match the committed outputs. They should be identical, because the seed is fixed.
- Change `STUDENT_ID_LAST4` in `R/config.R` to a different value and confirm the *conclusions* (not the exact numbers) stay the same.
- Tighten the skewness threshold to 0.10 in `R/config.R` and confirm every required n grows by roughly a factor of 4, as γ₁/√n predicts.
