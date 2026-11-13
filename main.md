**GitHub
repository***https://github.com/lwo-owl/DSA4212-Project-/tree/main*

# Introduction

Black-box optimisation seeks a good solution when an objective can be
evaluated but its mathematical form or derivatives are unavailable to
the algorithm. Evaluations may involve expensive computations or
experiments, making efficient use of a limited evaluation budget
important. Bayesian optimisation (BO) addresses this problem by
constructing a probabilistic surrogate and selecting evaluations through
an acquisition function [@frazier2018bo]. The covariance matrix
adaptation evolution strategy (CMA-ES) samples populations of candidate
solutions and adapts its search distribution using their rankings
[@hansen2016cma]. These different mechanisms suggest that the methods
may respond differently to objective geometry, observation noise and the
interaction between them.

Observation noise may make geometric challenges harder to resolve.
Differences in objective value between candidates can be small relative
to measurement noise, particularly when candidates differ mainly along a
direction of weak curvature. The extent of the disruption by noise may
depend on both the valley geometry and the optimisation method.
Comparing noisy performance with each method's noiseless baseline can
reveal whether particular geometries make one method more susceptible to
noise than the other.

Noise also changes how an algorithm can spend its evaluation budget.
Evaluating new locations may uncover better solutions, while repeating
evaluations can help distinguish promising candidates already found.
Repetition uses calls that could otherwise support further search. Its
value therefore depends on whether the main difficulty is finding a good
solution or recognising one that has already been evaluated.

This report investigates how valley geometry affects noiseless search,
how noise changes search and final selection and whether reevaluation or
a change to the recommendation rule improves the returned solution
within a fixed budget.

# Literature review

#### Valley geometry.

A Gaussian process (GP) with separate length scales for each coordinate
can represent different rates of variation along the coordinate axes.
Rotating an anisotropic valley generally produces a distance metric with
cross-coordinate terms that this kernel cannot express. This gives a
reason to expect BO to be sensitive to orientation. An isotropic kernel
provides a useful control because its distances are unchanged by
rotation, although a single length scale may be a poor fit for a narrow
valley [@rasmussen2006gpml]. Full-covariance CMA-ES can adapt its
sampling distribution to the valley's principal directions. When the
objective, feasible region and initial search distribution are rotated
together, its behaviour is equivalent in distribution. Learning the
valley's shape still takes evaluations, so performance can depend on
conditioning even when it is insensitive to orientation
[@hansen2016cma].

#### Noise during search and selection.

Noisy observations can change where BO searches by altering the
surrogate, while noise can change how CMA-ES updates its sampling
distribution by reversing candidate rankings
[@frazier2018bo; @hansen2009uncertainty]. Noise can also cause either
method to return an inferior point whose observed value was unusually
favourable [@smith2006curse]. This motivates measuring both the true
objective value of the best point visited and that of the final
recommendation [@astetemorales2016regret]. Their difference, when the
recommendation is a visited point, measures how much was lost in
selection. BO's posterior mean also provides a way to test whether a
different recommendation rule improves the final choice using the same
observations [@frazier2018bo].

#### Continued search versus verification.

@boesel2003cleanup described a two-stage approach to selecting among
candidates visited during a noisy search. The first stage uses existing
observations to screen candidates, and the second uses fresh evaluations
to distinguish those that remain. In a fixed-budget setting, those fresh
evaluations use calls that could otherwise be spent on new points.
Therefore, verification was compared with continued search under the
same total budget. To interpret the comparison, the loss from excluding
better candidates at shortlisting was separated from the loss from
selecting a worse point among the finalists.

# Geometry under noiseless evaluations {#sec:stage-a}

#### Experimental design.

This section investigates how conditioning, rotation and valley shape
affect search before observation noise is introduced. A synthetic valley
objective is minimised on the two-dimensional unit disk $\mathcal D$.
For an optimum location $c$ in the unrotated coordinate system, it is
defined $$(u,v)^\top = R_\theta^\top x-c,\qquad
  f_{\kappa,b,\theta}(x)
  = \frac{u^2+\kappa(v-bu^2)^2}{Z},$$ where $Z$ is the expected value of
the numerator at a uniformly sampled point in $\mathcal D$. The minimum
is zero, and one uniform draw has an expected gap of one. Three
condition numbers $\kappa\in\{1,10,100\}$, two bend values
$b\in\{0,0.75\}$ and two orientations $\theta\in\{0^\circ,45^\circ\}$
were compared. The Hessian at the optimum has condition number $\kappa$
for either bend value. Bending changes the geometry away from the
optimum by making the valley direction vary with location. Rotation
preserves the feasible disk and the objective geometry while changing
their representation in the search coordinates.

#### Optimisers and evaluation budget.

Uniform random search, full-covariance CMA-ES and Gaussian-process BO
with a squared-exponential kernel were compared against one another. BO
used either a separate length scale per coordinate (ARD) or a shared
length scale (isotropic) in every setting. Both variants minimised the
lower confidence bound $\mu(x)-2s(x)$ over 1,024 fresh uniform
candidates per call. The refined versions also ran constrained local
optimisation from up to two promising candidates. This refinement used
only the surrogate and incurred no extra objective calls. Kernel
hyperparameters were refitted by marginal likelihood every six calls.
The main experiment used 24 independent seed blocks. Within each block
the unrotated optimum was held fixed and orientations were paired by
rotating the corresponding random draws. Each BO run used the first six
CMA-ES points as its initial design. These points counted towards the
budget of 96 calls, which allowed 16 complete CMA-ES generations of six
candidates. CMA-ES started with an isotropic search distribution and
resampled infeasible points. The disk and the paired sampling procedure
made it possible to check rotation invariance without introducing a
change in the boundary or in initial information.

#### Performance measures.

The median best visited gap was plotted across seeds. At 96 calls, the
rotation effect was measured within each seed as
$D=\log_{10}(g_{45}/g_0)$, where $g_\theta$ is the terminal gap. A
positive $D$ suggests that rotation worsened performance. To assess the
effect of bending, we compared $D$ between curved and straight valleys
within the same seed. Means and pointwise 95% percentile intervals were
reported from 10,000 bootstrap resamples of complete seed blocks.

<figure id="fig:stage-a-convergence">
<span class="image placeholder"
data-original-image-src="straight-vs-curved-convergence.png"
data-original-image-title="" width="80%" height="30%"></span>
<figcaption>Median best normalised objective gap across 24 seed blocks.
Rows show straight and curved valleys, and columns show condition
numbers. Solid and dotted lines denote <span
class="math inline">0<sup>∘</sup></span> and <span
class="math inline">45<sup>∘</sup></span>, respectively. The vertical
axis is logarithmic.</figcaption>
</figure>

#### Results.

At 96 calls, both refined BO variants achieved lower median gaps than
CMA-ES in every tested setting, and CMA-ES improved on random search
(Figure [1](#fig:stage-a-convergence){reference-type="ref"
reference="fig:stage-a-convergence"}).

::: {#tab:stage-a-rotation}
    $\kappa$           Straight                      Curved               Curved minus straight
  ---------- ---------------------------- ---------------------------- ----------------------------
           1  $+0.231\,[-0.261,\,0.735]$   $-0.023\,[-0.727,\,0.584]$   $-0.254\,[-1.160,\,0.630]$
          10  $+0.499\,[0.042,\,0.958]$    $-0.127\,[-0.510,\,0.300]$   $-0.626\,[-1.288,\,0.069]$
         100  $-0.071\,[-0.586,\,0.446]$   $+0.586\,[0.271,\,0.904]$    $+0.657\,[0.046,\,1.280]$

  : Mean paired rotation effect for refined ARD BO at 96 calls, with
  pointwise 95% bootstrap intervals. The last column is the within-seed
  change in the rotation effect from straight to curved.
:::

CMA-ES matched its paired rotated runs to numerical precision, and
refined isotropic BO showed small rotation effects across all six
geometries. Refined ARD BO was more sensitive to orientation, but its
estimated rotation penalty did not increase consistently with
conditioning (Table [1](#tab:stage-a-rotation){reference-type="ref"
reference="tab:stage-a-rotation"}). The clearest penalty occurred on the
curved valley at $\kappa=100$, where the geometric mean gap after
rotation was $10^{0.586}\approx 3.9$ times the aligned gap. The paired
bend contrast indicates that bending increased the rotation penalty in
that setting. At $\kappa=10$, the estimated bend contrast went in the
opposite direction, but its interval spanned zero, leaving it unclear
whether bending reduced the penalty. The sensitivity to alignment is
consistent with the ablation study of @visser2025labcat, who found that
aligning the BO trust region with local principal components improved
performance particularly on low-dimensional unimodal benchmarks.

#### Budget sensitivity.

On eight additional seed blocks, in a straight valley setting,
checkpoints at 24, 48, 96 and 192 calls were examined. The positive
rotation effect for refined ARD at $\kappa=10$ in the main sample did
not recur at 96 calls in this smaller sample. At both 96 and 192 calls,
the intervals for the rotation effect of refined ARD relative to CMA-ES
included zero at every condition number. The check supports using 96
calls as a common short-budget comparison, while leaving open how the
curved-valley penalty changes with more evaluations.

#### Interpretation.

The rotation controls behaved as expected from their geometry. CMA-ES
can rotate its sampling covariance, so consistently rotating the search
state and draws preserves its behaviour [@hansen2016cma]. An isotropic
kernel likewise preserves similarity under rotation because it uses
Euclidean distances. An ARD kernel instead weights coordinate
differences separately and cannot represent a rotated anisotropic
distance using cross-coordinate terms [@rasmussen2006gpml]. This gives a
reason for ARD's sensitivity to orientation. The near-identical terminal
gaps of the pool-only BO variants, compared with the differences after
local refinement, show that the observed rotation effect also depended
on how the acquisition was optimised. This is consistent with benchmark
evidence that acquisition optimisation can materially affect BO
performance [@leriche2021coco].

The curved-valley comparison also shows why the condition number at the
optimum is an incomplete description of search difficulty. It is
unchanged by bending, yet the direction of the valley can vary over the
region sampled during search. Methods such as LABCAT explicitly adapt to
local directions on curved objectives, illustrating why this distinction
can matter for BO [@visser2025labcat]. The results identify a
geometry-dependent effect for the tested implementation and bend,
without establishing a general trend across condition numbers. These
noiseless results provide a reference for the next section, which
distinguishes whether noise prevents an optimiser from visiting good
points or from recognising and recommending them.

# Observation noise and evaluation budget {#sec:stage-b}

#### Experimental design.

Independent Gaussian noise was added to each objective evaluation, so
the optimiser observed $y(x)=f(x)+\sigma z$, with
$z\sim\mathcal{N}(0,1)$. The standard deviation $\sigma$ is expressed in
the normalised objective units defined earlier. The experiment crossed
$\sigma\in\{0,10^{-4},10^{-3},10^{-2}\}$ with all 12 earlier geometries.
Refined ARD BO, refined isotropic BO, full CMA-ES and uniform random
search were compared against one another over 24 seed blocks and 96
calls per run. BO and CMA-ES shared their first six locations and
observations. Corresponding calls in the two orientations used the same
noise draws within each block.

Both BO variants used the known noise variance in the GP and retained
the lower confidence bound based on latent uncertainty. At each
checkpoint, BO recommended the visited point with the lowest posterior
mean [@frazier2018bo]. CMA-ES ranked candidates by their observed
scores, and both CMA-ES and random search returned the point with the
lowest observation. The lowest-observation rule was also applied to the
same BO histories to assess how the recommendation rule affected the
result.

The budget extension retained both valley shapes and orientations at
$\kappa=100$ and used $\sigma\in\{0,10^{-3},10^{-2}\}$. The same search
strategies was continued to 192 calls and compared checkpoints at 96,
144 and 192 calls within each run. All additional calls followed the
original search policies.

#### Performance measures.

The clean objective was used to score the recommendation and the best
point visited. Their difference is the *selection loss*, the quality
lost by returning a different point from the best one found. It was
calculated within each run before aggregating across seeds. For the
budget comparison, the mean within-seed change was calculated from 96 to
192 calls and pointwise 95% intervals were obtained from 10,000
bootstrap resamples of seed blocks. Additionally, the fraction of
CMA-ES's three selected parents that would be replaced under the clean
ranking was measured, adapting a rank-based diagnostic
[@hansen2009uncertainty].

#### Effect of noise at 96 calls.

Both refined BO variants returned lower median gaps than CMA-ES across
the tested geometries. At high noise, however, their advantage in the
returned solution could be much larger than their advantage in search.
On the rotated curved valley at $\kappa=100$ and $\sigma=10^{-2}$,
isotropic BO and CMA-ES both visited points with median gaps close to
$10^{-4}$, yet isotropic BO returned a solution with a median gap about
3.6 times smaller. BO's search also deteriorated with noise. In the same
valley, isotropic BO's median best visited gap rose from
$1.92\times10^{-6}$ without noise to $1.02\times10^{-4}$ at the highest
noise level. Its posterior-mean recommendation improved on the
lowest-observation rule in the high-noise example, although the gains
were mixed at lower noise.

The evidence that noise changed the effect of rotation was weaker. For
refined ARD BO on the curved $\kappa=100$ valley, the estimated rotation
penalty for the recommendation gap decreased at all three nonzero noise
levels relative to the clean control. At $\sigma=10^{-2}$, the estimated
change was $-0.445$ on the $\log_{10}$ scale, with a 95% paired interval
of $[-0.940,\,0.104]$. The intervals at the other two noise levels also
included zero.

#### Effect of a larger budget.

At $\sigma=10^{-2}$, the mean paired recommendation gap improved from 96
to 192 calls for CMA-ES on both valley shapes and for both BO variants
on the curved valley. The intervals for BO on the straight valley
included zero. Both BO variants still returned lower median gaps than
CMA-ES in all four settings at 192 calls
(Figure [2](#fig:stage-b-budget){reference-type="ref"
reference="fig:stage-b-budget"}).

<figure id="fig:stage-b-budget" data-latex-placement="!htb">
<span class="image placeholder"
data-original-image-src="recommendation-budget-0.01.png"
data-original-image-title="" width="90%" height="36%"></span>
<figcaption>Median true recommendation gap at 96, 144 and 192 calls
across 24 seed blocks, for <span
class="math inline"><em>κ</em> = 100</span> and observation noise
standard deviation <span
class="math inline">10<sup>−2</sup></span>.</figcaption>
</figure>

::: {#tab:stage-b-budget}
+:--------------+:--------:+:---------:+:---------:+:---------:+
|               | Best visited gap     | Recommendation gap    |
+---------------+----------+-----------+-----------+-----------+
| Method        | 96 calls | 192 calls | 96 calls  | 192 calls |
+---------------+----------+-----------+-----------+-----------+
| Random search | 14.99    | 6.59      | 43.55     | 37.20     |
+---------------+----------+-----------+-----------+-----------+
| CMA-ES        | 1.01     | 0.45      | 31.08     | 22.69     |
+---------------+----------+-----------+-----------+-----------+
| ARD BO        | 2.20     | 1.32      | 7.55      | 7.35      |
+---------------+----------+-----------+-----------+-----------+
| Isotropic BO  | 1.02     | 0.53      | 8.64      | 4.48      |
+---------------+----------+-----------+-----------+-----------+

: Median true gaps on the rotated curved valley at $\kappa=100$ and
$\sigma=10^{-2}$. All entries are in units of $10^{-4}$ and lower is
better.
:::

Between 144 and 192 calls on the noisy rotated curved valley, the median
best visited gap fell for both CMA-ES and isotropic BO, while their
median recommendation gaps and selection losses increased. Within each
run, the best visited gap can only stay the same or decrease as more
points are evaluated. By 192 calls CMA-ES had a slightly lower median
best visited gap than isotropic BO, yet its median recommendation gap
was about five times larger
(Table [2](#tab:stage-b-budget){reference-type="ref"
reference="tab:stage-b-budget"}).

#### Interpretation.

In the later noisy generations, roughly 40--50% of CMA-ES's selected
parents differed from the clean top three, approaching the 50% reference
for independent rankings. The search distribution also remained wider
than in the clean runs. This is consistent with noise interfering with
late refinement [@hansen2009uncertainty]. CMA-ES nevertheless continued
to find good candidates, and its selection loss shows how much quality
was lost when it returned the point with the most favourable
observation.

Under high noise, BO queried farther from previous locations than in the
clean runs, even as its queries became more local with budget. In the
rotated cases, the exact best visited point was rarely among the five
lowest posterior means by call 192. Near-equivalent points can still
rank highly, so the recommendation gap measures the practical cost of
these ranking errors. Together, these results motivate testing whether
spending part of a fixed budget on reevaluation can improve final
selection enough to offset the reduction in search calls.

# Allocating evaluations under noise {#sec:stage-c}

An optimiser can visit better points as its budget increases while
returning a worse recommendation. Therefore, two ways to improve the
returned solution are tested: spending part of the budget on repeated
measurements and using the existing observations to make a better final
choice. Each policy had 96 objective calls, including initialisation.
Each round used 24 fresh seed blocks, paired across policies within the
round. The primary setting was the rotated curved valley at $\kappa=100$
and $\sigma=10^{-2}$. In the earlie section, isotropic BO and CMA-ES
visited points of comparable quality but returned recommendations with
substantially different gaps, making it a useful case for testing
changes to final selection. The geometric mean of within-seed
recommendation gap ratios, with pointwise 95% bootstrap intervals are
reported A ratio below one favours the first policy in a comparison.

#### Fixed allocation of repeats.

96 calls of continued search were compared, with two policies that spent
72 calls on search proposals and 24 on repeats. Finalist verification
used the first 72 calls of the continued-search run and allocated six
fresh measurements to each of four finalists. CMA-ES shortlisted points
by their observed scores, while BO used its GP posterior means. Both
returned the finalist with the lowest mean of the fresh measurements
[@boesel2003cleanup]. The late-repeat policy evaluated each of the last
24 proposals twice after 48 ordinary calls. CMA-ES ranked candidates by
the mean of each pair, and BO incorporated both observations into its
GP.

Neither repeat policy clearly improved on continued search in the
primary setting. Late repeats reduced some CMA-ES parent ranking errors,
but the policy completed 12 generations instead of 16 and had a worse
median best visited gap. Finalist verification increased ARD BO's paired
recommendation gap by a factor of $1.44\,[1.13,1.89]$. About 89% of its
mean selection loss arose before verification, from the gap between the
best shortlisted point and the best point visited.

#### Adaptive allocation and shortlist coverage.

For CMA-ES, repeats were capped at 12 before updates and the observed
third and fourth candidates were reevaluated when their scores differed
by at most $\sqrt{2}\sigma$ after generation eight. This targeted
uncertainty near the parent selection boundary [@hansen2009uncertainty].
After the last complete generation, any remaining calls alternated
between the two visited points with the lowest current average scores.
The policy completed a median of 14 generations, but showed no sustained
reduction in parent ranking errors at matched generations and no clear
recommendation gain.

For BO, the 72-call search history was retained. We compared two ways of
spending the remaining calls. First, we used draws from the joint GP
posterior to choose a shortlist covering more plausible near-best
points, then applied the original verification rule. This reduced mean
shortlist loss by $9.16\times10^{-4}$ for ARD BO and $3.28\times10^{-4}$
for isotropic BO, but did not reliably improve the returned solutions.
Better finalists still had to be distinguished using noisy measurements.
Second, we switched from the lower confidence bound (LCB) to a
finite-candidate knowledge-gradient (KG) approximation that could choose
either a repeat or a new location based on its expected effect on the
final recommendation [@scott2011kg]. This gave the lowest median
recommendation gap for isotropic BO in that round, although its paired
interval against continued LCB included one.

#### Final recommendation and switch timing.

On fresh seeds, an isotropic GP was fitted to each completed CMA-ES
history using the known noise variance. The visited point with the
lowest posterior mean was compared with the point returned by the
original recommendation rule [@frazier2018bo]. For isotropic BO, we
compared switching from LCB to KG after 60 or 72 calls with using LCB
throughout.

::: {#tab:stage-c-final}
  Comparison                                         Gap ratio \[95% interval\]   Seeds improved
  ------------------------------------------------- ---------------------------- ----------------
  Standard CMA-ES, GP vs. original recommendation       $0.44\,[0.27,0.72]$          $18/24$
  Cutoff CMA-ES, GP vs. original recommendation         $0.46\,[0.29,0.73]$          $15/24$
  Cutoff vs. standard CMA-ES, both using GP             $0.86\,[0.64,1.15]$           $9/24$
  Isotropic BO, KG after 72 calls vs. LCB               $0.74\,[0.40,1.27]$          $14/24$
  Isotropic BO, KG after 60 calls vs. LCB               $1.02\,[0.56,2.00]$          $12/24$

  : Paired true recommendation gap ratios at 96 calls on the rotated
  curved valley at $\kappa=100$ and $\sigma=10^{-2}$. Intervals are
  pointwise 95% bootstrap intervals across 24 fresh seed blocks.
:::

GP reranking gave the clearest improvement
(Table [3](#tab:stage-c-final){reference-type="ref"
reference="tab:stage-c-final"}). For standard CMA-ES, median selection
loss fell from $2.38\times10^{-3}$ to $1.70\times10^{-4}$. Both
recommendations came from the same visited set, so the gain was
attributable to selection and required no additional objective calls.
Under the same GP recommendation rule, cutoff repeats showed no clear
advantage over standard search. For BO, the earlier KG switch used a
median of 21 repeat calls, compared with 14.5 for the later switch,
while both returned almost identical median gaps. The BO results leave
the benefit of changing the acquisition uncertain.

#### Shared-history ablation.

We ran 24 fresh seed blocks (10100--10123) on the four valleys with
$\kappa=100$, at $\sigma\in\{0,10^{-3},10^{-2}\}$. Each isotropic BO
branch shared 72 LCB evaluations and had 24 calls left. The five
branches used continued LCB, alternation between new evaluations and
repeats (12 each), or a finite-candidate, one-step approximation to KG
[@frazier2018bo; @scott2011kg]. The three KG variants allowed both types
of action, new locations only, or repeats only. The criterion minimised
the expected smallest updated posterior mean over the visited locations
and the proposed point, using 24 quadrature nodes. It considered all
distinct visited locations and up to 24 new locations from the LCB
candidate procedure. All branches used the same GP specification,
refitting schedule and posterior-mean recommendation rule. The primary
contrast in this round was unrestricted versus new-only KG on the
rotated curved valley at $\sigma=10^{-2}$.

For visited set $V_B$ and recommendation $\hat{x}_B\in V_B$, we used
$$R_B=f(\hat{x}_B)-f^\star=S_B+L_B,\qquad
S_B=\min_{x\in V_B} f(x)-f^\star,\qquad
L_B=f(\hat{x}_B)-\min_{x\in V_B} f(x).$$ mean paired changes were
decomposed on the raw gap scale, so $\Delta R=\Delta S+\Delta L$ in
every block. The geometric means of within-block gap ratios were
reported, flooring gaps at $10^{-14}$ before taking logs. The new
intervals use 10,000 bootstrap resamples of complete seed blocks.
Intervals are pointwise 95% unless stated otherwise. The clean objective
scores outcomes and diagnostics only.

::: {#tab:mechanisms}
  Comparison                         Ratio \[95% interval\]   $\Delta S$   $\Delta L$   $\Delta R$
  -------------------------------- ------------------------ ------------ ------------ ------------
  KG unrestricted vs. new only          $0.76\,[0.46,1.29]$      $+0.16$      $-1.07$      $-0.91$
  KG unrestricted vs. LCB               $0.80\,[0.39,1.73]$      $-1.57$      $-0.34$      $-1.91$
  KG new only vs. LCB                   $1.06\,[0.51,2.13]$      $-1.73$      $+0.73$      $-1.01$
  KG repeats only vs. LCB               $0.73\,[0.37,1.26]$      $+0.64$      $-0.49$      $+0.15$
  Alternation vs. LCB                   $1.47\,[1.11,2.05]$      $+0.15$      $+2.87$      $+3.02$
  Standard CMA-ES, GP vs. native       $0.54\,[0.31,1.001]$          $0$     $-11.70$     $-11.70$
  Cutoff CMA-ES, GP vs. native          $0.43\,[0.29,0.64]$          $0$     $-15.57$     $-15.57$

  : Fresh seed blocks, rotated curved valley at $\sigma=10^{-2}$. Ratios
  below one favour the first policy. Raw changes are means in units of
  $10^{-4}$. GP reranking comparisons hold the history fixed.
:::

#### Candidate availability and selection.

At $\sigma=10^{-2}$, perfect selection from the first 72 BO points could
beat the call-96 LCB recommendation in 22, 22, 23 and 23 of 24 blocks on
the straight aligned, straight rotated, curved aligned and curved
rotated valleys, respectively. Restricting the oracle to the four lowest
posterior means reduced these counts to 10, 12, 12 and 11. Shortlist
omission accounted for 82--86% of the mean selection loss at call 72
across the four valleys. On the primary valley, for the two clean best
points in that shortlist, six independent measurements per point gave a
mean misordering probability of 49%. For a pair with true gap difference
$\Delta$, that probability is $\Phi(-\Delta\sqrt{6}/(\sqrt{2}\sigma))$.
Weighting by $\Delta$, the expected loss in this two-point comparison
fell from $5.74\times10^{-5}$ with one fresh measurement per point to
$5.52\times10^{-5}$ with six. The small expected loss alongside frequent
misorderings reflects near ties, while the shortlist omission cost was
substantially larger.

#### Fixed-budget shortlist size.

We held the 72-call histories fixed and allocated the remaining 24 calls
as 4 finalists with 6 repeats each, 8 with 3, or 12 with 2. The
shortlists were nested by the call-72 posterior mean. We simulated fresh
Gaussian measurement errors under the original noise model. Both
readouts used the same fresh observations and selected within the same
shortlist. One ranked fresh sample means. The other updated the call-72
GP with those means, using noise variance $\sigma^2/r$ for $r$ repeats
per finalist and holding the fitted hyperparameters fixed. We averaged
over 4,096 simulations within each seed block, then resampled the 24
blocks for pointwise intervals. These simulations estimate verification
performance conditional on the saved search histories.

::: {#tab:shortlist-size}
+:------------+----------:+------------:+-----------:+------------:+-----------:+
| Finalists   | Shortlist | Ranking loss             | Change vs. LCB           |
| $\times$    | loss      |                          |                          |
| repeats     |           |                          |                          |
+-------------+-----------+-------------+------------+-------------+------------+
|             |           | Fresh means | Updated GP | Fresh means | Updated GP |
+-------------+-----------+-------------+------------+-------------+------------+
| $4\times6$  | 7.47      | 1.81        | 1.57       | $+2.60$     | $+2.36$    |
+-------------+-----------+-------------+------------+-------------+------------+
| $8\times3$  | 6.33      | 3.01        | 2.58       | $+2.66$     | $+2.23$    |
+-------------+-----------+-------------+------------+-------------+------------+
| $12\times2$ | 5.66      | 3.56        | 3.22       | $+2.54$     | $+2.21$    |
+-------------+-----------+-------------+------------+-------------+------------+

: Verification on the rotated curved valley at $\sigma=10^{-2}$. Mean
gap units are $10^{-4}$. Ranking loss is the gap above the best
shortlisted point. Positive changes relative to continued LCB are worse.
:::

With four finalists, even the shortlist oracle had a worse mean gap than
continued LCB on the primary valley. Widening to 12 lowered shortlist
loss by $1.80\times10^{-4}$, while fresh-mean ranking loss rose by
$1.75\times10^{-4}$ (Table [5](#tab:shortlist-size){reference-type="ref"
reference="tab:shortlist-size"}). The net paired change for fresh means
was $-0.05\,[-0.80,0.72]$ in units of $10^{-4}$. All six verification
policies had higher mean gaps than continued LCB on each valley at high
noise, with some pointwise intervals spanning zero. On the primary
valley, the 12-finalist GP policy's mean excess gap was
$2.21\,[0.32,4.86]\times10^{-4}$.

#### Search and selection effects.

Unrestricted KG's point estimate against new-only KG combined lower
selection loss and a slightly higher best visited gap
(Table [4](#tab:mechanisms){reference-type="ref"
reference="tab:mechanisms"}). The mean change in selection loss was
$-1.07\,[-4.33,1.98]\times10^{-4}$. The primary ratio interval admits
gains and losses. Its 90% interval $[0.50,1.18]$ also extends outside an
illustrative $[0.80,1.25]$ similarity range. At this noise level, the
same contrast on each of the other three valleys also had a pointwise
interval spanning one. Relative to LCB, both KG policies that could
visit new points reduced the mean best visited gap, with pointwise
intervals for $\Delta S$ below zero. The repeats-only policy gave up
$6.38\times10^{-5}$ of mean search progress and recovered
$4.88\times10^{-5}$ through selection, so the mean recommendation gap
rose slightly even though the geometric mean ratio favoured repeats.
These summaries weight runs differently. The earlier transfer sample
gave a KG-to-LCB ratio of $0.56\,[0.37,0.83]$ in this setting. The fresh
estimate was closer to one and more uncertain.

At $\sigma=0$, unrestricted KG used a median of 19 repeats in the
rotated curved valley and had a ratio of $1.87\,[1.14,3.46]$ against
new-only KG. The implementation uses numerical jitter in the GP and does
not enforce exact interpolation or suppress repeats in the clean
control. A repeat supplies no new information about a deterministic
objective, so clean repeat choices reveal a property of the implemented
policy. The clean control does not isolate which numerical or modelling
choice drove them.

In the fresh CMA-ES runs, GP reranking reduced selection loss in both
search histories on average. For standard CMA-ES at high noise, the gap
ratio was 0.27 on the straight valleys and 0.54 on the curved valleys.
The paired rotations gave the same ratio within each shape. In the
primary valley, the fresh interval for standard search narrowly crossed
one, while the interval for cutoff search stayed below one. With the GP
readout held fixed, the cutoff-to-standard ratio was
$1.15\,[0.72,1.78]$. At matched generations 9--14, cutoff repeats
reduced the mean parent replacement fraction by 2.1 percentage points
($[-4.9,0.5]$), while using 84.5 new evaluations on average versus 96
for standard CMA-ES. The cutoff rule is a local heuristic. Full CMA-ES
uncertainty handling also uses rank changes to adapt evaluation effort
and search spread [@hansen2009uncertainty]. The budget trade-off is
related to the depth-versus-fidelity question studied by
@wanglu2026depth.

#### Supporting benchmark checks.

The comparison was repeated on two-dimensional BBOB f14 and f15 on
$[-5,5]^2$ (instances 201--224). Function f14 has increasing anisotropy
near the optimum [@cocof14], and f15 is highly multimodal [@cocof15].
Each instance's objective gap was scaled by its mean over 4,096 Sobol
reference points. At $\sigma=10^{-2}$, the median ratio of the clean
shortlist gap between its best two points to $\sigma$ was 0.003 on the
primary valley, 0.022 on f14 and 0.326 on f15. Continued LCB reduced the
mean best visited gap over calls 73--96 by $6.38\times10^{-5}$,
$2.86\times10^{-4}$ and $1.70\times10^{-3}$, respectively. The call-72
oracle could beat the continued recommendation in 23, 20 and 14 of 24
blocks, respectively. Thus the same global noise scale yielded different
local ranking margins and remaining opportunities for search.

The high-noise KG-to-LCB gap ratios were 1.03 on f14 and 1.04 on f15. To
compare each with the primary valley, we subtracted its within-seed log
policy effect from the benchmark's, then exponentiated the mean
difference. The resulting ratios of ratios were $1.29\,[0.60,2.61]$ and
$1.30\,[0.55,2.70]$, using Bonferroni-adjusted 97.5% intervals for the
two exploratory comparisons. Both intervals permit substantial
differences in either direction.

# Conclusion

The experiments distinguish four limitations. First, geometry can slow
discovery, as for refined ARD BO on the strongly conditioned rotated
curved valley. Second, a competitive point can be present but excluded
by the shortlist. Third, a wider shortlist can improve candidate
coverage while reducing the number of measurements available to
distinguish each finalist. Fourth, allocating calls to verification can
forgo useful search progress. The 192-call experiment also showed that
further search can improve the visited set even when the final
recommendation worsens.

For isotropic BO at 96 calls and high noise, continued LCB with a
posterior-mean recommendation is the better default among the fixed
allocations tested. The benefit of adaptive repeats under the KG
approximation remains uncertain. For CMA-ES, standard search with GP
reranking is a promising strategy when the noise variance is known. It
changes the recommendation using existing evaluations, and cutoff
repeats showed no clear added gain under the same readout. These
recommendations concern evaluation efficiency on the tested
two-dimensional objectives and do not compare wall-clock costs or
establish the best policy for other budgets or noise models.

#### Reproducibility.

The updated master notebook's `report` mode runs the main stages. In the
`full` profile, `new` regenerates the fresh seed blocks and the
shortlist follow-up. The `shortlist` mode reads the saved archive
without rerunning search. `transfer_original` reproduces the earlier
transfer sample. The report source includes the analysis script for the
complete-block checks, decompositions and cross-function intervals.

# AI declaration

Generative AI tools (ChatGPT and codex) were used for literature search,
coding and enhancing expressions.
