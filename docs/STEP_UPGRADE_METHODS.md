# STEP_UPGRADE: algorithm description with equations

Status: draft, 2026-10-18. Written as a methods supplement (in the style of
Chang et al. 2016, Supplemental Material S1) for later papers and posters.
It describes the code exactly. Where the code makes a choice the literature
does not specify, the choice is stated. Section 9 lists open questions found
while writing this document.

Code: `step/identification.py`, `step/identification_fast.py`,
`step/tracking.py` (v0.4.0b3), and the opt-in modules
`step/identification_regions.py`, `step/tracking_coverage.py`,
`step/lineage_events.py`, `step/chang2016.py`, `step/chang2016_hourly.py`.
Experimental numbers are in `docs/STEP_UPGRADE_SUPPLEMENT.md`.

## 1. Notation

| Symbol | Meaning |
|---|---|
| $\Omega$ | grid cells $s = (y, x)$ of the domain, spacing $\Delta x$ (km) |
| $t = 0, 1, \dots$ | frame index, spacing $\Delta t$ (h) |
| $r_t(s)$ | rain rate (mm/h), mean over the interval ending at frame $t$ |
| $V_t \subseteq \Omega$ | valid cells (finite data, optional mask, optional edge band) |
| $\theta$ | wet threshold (mm/h) |
| $\lVert s-u \rVert$ | Euclidean distance in cells |
| $\mathcal{D}_\rho(X)$ | disk dilation, $\{ s : \min_{u \in X} \lVert s-u\rVert \le \rho \}$ |
| $\mathcal{E}_\rho(X)$ | disk erosion, $\{ s : \{u : \lVert s-u\rVert \le \rho\} \subseteq X \}$ |
| $\mathrm{cc}(X)$ | 8-connected components of $X$ |
| $\lvert X \rvert$ | number of cells in $X$ |

The wet set at frame $t$ is

$$W_t = \{ s \in V_t : r_t(s) \ge \theta,\ r_t(s) > 0 \}.$$

The Chang et al. reference uses the strict form $r_t(s) > \theta$
(`strict=True`).

## 2. Objects and their attributes

An identification step returns disjoint objects $M_{t,1}, M_{t,2}, \dots
\subseteq W_t$. For an object $M$ at frame $t$:

$$A = \lvert M \rvert, \qquad P = \sum_{s \in M} r_t(s), \qquad \bar r = P / A,$$

$$c = \frac{\sum_{s \in M} r_t(s)\, s}{P} \quad (\text{if } P > 0; \text{ otherwise the unweighted mean of } M).$$

Here $A$ is in cells (area $A\,\Delta x^2$ km$^2$), $P$ is in mm/h summed over
cells, and $c$ is the intensity-weighted centroid. The bounding box
$\mathrm{bb}(M)$ is the half-open row/column range of $M$. Labels are
renumbered by the first cell of each object in row-major order, so they do
not depend on processing order.

## 3. Identification

### 3.1 Preset: almost-connected labelling (`identify`, `identify_disk`)

With bridge radius $\rho$ (cells):

$$G_t = \mathcal{D}_\rho(W_t) \cap V_t, \qquad \{M_{t,k}\} = \{ C \cap W_t : C \in \mathrm{cc}(G_t),\ C \cap W_t \ne \emptyset \}.$$

Objects with fewer than $n_{\min}$ cells are removed. Two wet cells belong
to the same object whenever a chain of wet cells with gaps of at most about
$2\rho$ connects them; no other condition applies. The morphology acts on
the binary mask, so labels cannot depend on label values (the original
STEP applied grey morphology to the label image).
`identify_disk` computes $\mathcal{D}_\rho$ with an exact Euclidean
distance transform, $\mathcal{D}_\rho(X) = \{ s : \mathrm{edt}(s; X)^2 \le
\rho^2 \}$, and gives identical labels.

Preset: $\theta = 1$ mm/h, $\rho = 9$ cells (36 km at 4 km), $n_{\min} = 1$.

### 3.2 Large/small regions (`identify_regions`, opt-in)

This follows S1, with the core radius $\rho_c$ and the bridge radius
$\rho_b$ allowed to differ ($\rho_b = \rho_c$ reproduces S1 and
`chang2016.identify_frame`).

1. Contiguous regions: $\{B_i\} = \mathrm{cc}(W_t)$.
2. Cores: $E = \mathcal{E}_{\rho_c}(W_t)$, with cells outside the domain
   counted as wet (`erosion_border='ignore'`, the behaviour of the original
   STEP). Region $i$ is *large* if $B_i \cap E \ne \emptyset$; the set of
   large regions is $\mathcal{L}$, and the rest are *small* ($\mathcal{S}$).
3. Large grouping (`large_join='almost'`): opening
   $O = \mathcal{D}_{\rho_c}(E) \cap W_t$, clusters
   $\{K_m\} = \mathrm{cc}(\mathcal{D}_{\rho_b}(O) \cap V_t)$. Build the
   bipartite graph with an edge $(i, m)$ whenever $B_i \cap O \cap K_m \ne
   \emptyset$; each connected component of this graph is one large segment
   $L_g$. With `large_join='region'`, each large region is its own segment.
4. Small attachment: each small region joins the large segment that its own
   dilation covers most,

   $$g(j) = \arg\max_g \lvert \mathcal{D}_{\rho_b}(B_j) \cap L_g \rvert, \qquad j \in \mathcal{S},$$

   with ties going to the lowest $g$. A small region is attached only if the
   maximum is positive. Because each small region is dilated on its own, it
   can never connect two large segments.
5. Unattached small regions $U$ (`small_join='almost'`): segments are
   $\mathrm{cc}(\mathcal{D}_{\rho_b}(\cup_{j \in U} B_j) \cap V_t)$
   restricted to $W_t$. With `small_join='region'`, each $B_j$ stays
   separate.

Edge band (`edge_cells` $=e$): $V_t \leftarrow V_t \cap \{ s : e \le y < H-e,\ e \le x < W-e \}$.
No rain is labelled in the band, and no dilation crosses it.

Recommended: $\theta = 1$ mm/h, $\rho_c = 3$, $\rho_b = 9$, $e = 15$ cells.

## 4. Lineage tracking (`Tracker`, preset policies)

Tracking is sequential. At frame $t$ the tracker keeps the *active*
terminals (the objects of frame $t-1$) and the *dormant* terminals (objects
last seen at $t-k$, $2 \le k \le G+1$, that were not continued). $G$ is
`max_gap`. Each terminal $p$ stores its object $M_p$, its time $t_p$ and,
when it continued a parent $p^-$, that parent's centroid and time.

### 4.1 Motion

$$v_p = \frac{c_p - c_{p^-}}{t_p - t_{p^-}} \ \text{ if } p \text{ continued a parent, else } v_p = 0,$$

$$k = t - t_p, \qquad \hat c_p = c_p + k\, v_p, \qquad \sigma_p = \mathrm{round}(k\, v_p),$$

where $\hat c_p$ is the predicted centroid and $\sigma_p$ the integer
advection shift. Children created by a split or merge carry no velocity.
With `velocity_reset='morphology'`, the velocity is also dropped when the
area ratio of the continuation lies outside $[1/2, 2]$ or its displacement
exceeds $\kappa$. The preset uses `'off'`.

### 4.2 Candidates

Let $R = \kappa k$, where $\kappa$ is `max_displacement` (cells per
frame), and $\tilde R = \max(R, 1)$. Child $j$ is a candidate of parent
$p$ if

$$\lVert c_j - \hat c_p \rVert \le R \ \lor\ \lVert c_j - c_p \rVert \le R \quad (\text{proximity}),$$

or if the bounding box of $j$ intersects that of $M_p$ or of $M_p + \sigma_p$
and $0 < \max(\pi_p, \pi_j, \pi^{a}_p, \pi^{a}_j)$ with that maximum
$\ge \varepsilon_{ev}$ (`event_overlap`), where the $\pi$ are defined next.

### 4.3 Overlap measures

For the raw masks and for the advected parent $M_p + \sigma_p$:

$$I = \lvert M_p \cap M_j \rvert, \quad \mathrm{IoU} = \frac{I}{\lvert M_p \cup M_j\rvert}, \quad \pi_p = \frac{I}{A_p}, \quad \pi_j = \frac{I}{A_j},$$

and $I^a, \mathrm{IoU}^a, \pi^a_p, \pi^a_j$ with $M_p$ replaced by
$M_p + \sigma_p$. The intensity similarity is

$$\varrho = \frac{\min(\bar r_p, \bar r_j)}{\max(\bar r_p, \bar r_j, 10^{-12})} \in [0, 1].$$

### 4.4 Score

Legacy score (used for gaps, $k \ge 2$, and with `score_policy='legacy'`):

$$S = 0.60 \max(\mathrm{IoU}, \mathrm{IoU}^a) + 0.25\, e^{-\lVert c_j - \hat c_p\rVert / \tilde R} + 0.15\, \varrho.$$

Coherent score (preset, `score_policy='coherent_adjacent'`, $k = 1$). It
compares a stationary hypothesis with a constant-velocity hypothesis:

$$S = \max\Big( 0.60\, \mathrm{IoU} + 0.25\, e^{-\lVert c_j - c_p\rVert/\tilde R},\ \ 0.60\, \mathrm{IoU}^a + 0.25\, e^{-\lVert c_j - \hat c_p\rVert/\tilde R} \Big) + 0.15\, \varrho.$$

$S \in [0, 1]$ is a matching score, not a probability. The weights are
fixed in the code.

### 4.5 Topology: splits, merges and strong overlaps

Under `adjacent_policy='overlap_ranked'` (preset) and `'overlap_balanced'`,
only raw (unadvected) overlap is used. A pair is an *event edge* if

$$I \ge \max(16, n_{ev}), \qquad \min(\pi_p, \pi_j) \ge 0.15, \qquad \max(\pi_p, \pi_j) \ge 0.5,$$

and a *strong edge* if $I \ge 16$ and $\min(\pi_p, \pi_j) \ge 0.30$.
Under `adjacent_policy='score'`, an event edge instead needs
$\max(\pi_p, \pi_j, \pi^a_p, \pi^a_j) \ge \varepsilon_{ev}$, an
intersection of at least $n_{ev}$ cells and (with
`event_policy='score_and_overlap'`) $S \ge \tau$; there are no strong edges.

Take the connected components of the graph formed by event and strong
edges. A component with one parent and several children is a **split**,
several parents and one child a **merge**, and several of each is
**complex**. One-to-one components are not events. In every event, each
child starts a new branch, and the families of all parents and children are
united (section 4.9). All parents and children of events are then removed
from further matching.

### 4.6 Continuations

The remaining pairs are linked in three passes, each using only parents and
children not linked yet.

1. **Strong overlap**: every remaining strong edge is a continuation.
2. **Assignment**: among the remaining candidates with $S \ge \tau$, choose
   the one-to-one matching that maximizes $\sum S$. Each parent may stay
   unmatched (an extra zero-score column per parent). This is solved with
   `linear_sum_assignment` separately on each connected component of the
   candidate graph.
3. **Mutual dominant overlap** (`overlap_ranked`): pair $(p, j)$ is linked if
   $I_{pj}$ is the largest intersection of $p$ over all children and of $j$
   over all parents, is larger than and at least twice the runner-up on both
   sides, and satisfies $I_{pj} \ge 16$ and $\min(\pi_p, \pi_j) \ge 0.15$.
   This lets weak-scoring but geometrically unambiguous pairs link.

A continuation gives the child the parent's branch and a velocity
(section 4.1).

### 4.7 Gaps

Dormant parents with $2 \le k \le G+1$, together with this frame's
unmatched active parents (which become dormant), are candidates for the
children left unlinked. Candidates are found as in section 4.2 with the
legacy score, and must satisfy $\min(A_p, A_j) \ge n_{gap}$
(`gap_min_pixels`).

* Conflict filter (`gap_conflict_policy='endpoint_overlap'`, preset): drop
  $(p, j)$ if some object $m$ at frame $t-1$ covers both endpoints,
  $\lvert M_m \cap M_p\rvert \ge \max(16, A_p/2)$ and
  $\lvert M_m \cap M_j\rvert \ge \max(16, A_j/2)$. In that case the frame in
  between was not empty, and the adjacent-frame logic already decided.
* Ambiguity: drop child $j$ if its two best gap scores differ by less than
  $\delta$ (`gap_ambiguity`).
* Assignment as in section 4.6, item 2, with threshold $\tau_{gap}$.

A gap link keeps the branch, and the edge records $k-1$ missing frames.

### 4.8 Identifiers

* *node*: one object at one frame.
* *branch*: a chain of continuation and gap edges. It ends at every split or
  merge, where each child starts a new branch.
* *family*: the union-find closure of branches joined by any split, merge or
  complex event. Families can only grow; the root is the smallest branch id.

### 4.9 Preset parameters (4 km, 1 h)

| Symbol | Name | Value | Physical meaning at 4 km / 1 h |
|---|---|---|---|
| $\theta$ | `threshold_mm_h` | 1 | 1 mm/h |
| $\rho$ | `bridge_radius_cells` | 9 | 36 km |
| $\tau$ | `tau` | 0.35 | |
| $\kappa$ | `max_displacement` | 20 cells/frame | 80 km/h |
| $G$ | `max_gap` | 1 | one missing hour |
| $\tau_{gap}$ | `gap_tau` | 0.45 | |
| $\delta$ | `gap_ambiguity` | 0.05 | |
| $\varepsilon_{ev}$ | `event_overlap` | 0.10 | |
| $n_{gap}$ | `gap_min_pixels` | 16 | 256 km$^2$ |
| | 16-cell minima (sections 4.5-4.7) | 16, hard-coded | 256 km$^2$ |

## 5. Opt-in tracking options (`CoverageTracker`)

All options default to the behaviour of section 4.

**Child coverage** (`overlap_metric='child_coverage'`): in the score
(section 4.4), $\mathrm{IoU} \to \pi_j$ and $\mathrm{IoU}^a \to \pi^a_j$.
A small child lying inside a large parent then scores
$0.60 + 0.25\,e^{-d/\tilde R} + 0.15\,\varrho$ instead of about
$0.25\,e^{-d/\tilde R} + 0.15\,\varrho$.
The overlap-coefficient variant (`'coverage'`) uses
$\max(\pi_p, \pi_j)$; it lets a large child inherit a small parent's branch
and is not recommended.

**Shared parents** (`share_parent=True`): after section 4, a child $j$ with
no incoming edge joins the branch of

$$p^\star = \arg\max_{p \in \text{frame } t-1} I_{pj} \quad \text{if} \quad \frac{I_{p^\star j}}{A_j} \ge \gamma \ \text{ and } \ I_{p^\star j} \ge n_{sh},$$

with $\gamma = 0.5$ and $n_{sh} = 1$ by default. Ties go to the earliest
parent. The edge type is `share`, and the child gets no velocity. This is
the S1 rule that several segments may have the same predecessor.

**Overlap rescue** (`rescue_overlap=True`, sensitivity option): before
sharing, pair each unlinked child $j$ with a parent $p$ that has no
successor, if $I_{pj}$ is the largest intersection of $j$ among such
parents and of $p$ among such children, and $I_{pj} \ge n_{rs}$. This is
done first with frame $t-1$ (edge `overlap_continue`) and then with
dormant objects from $t-2$ (edge `gap_overlap`). It remains one-to-one.

## 6. Events and the four-factor decomposition

### 6.1 Paper-style events (`paper_event_ids`)

Process nodes in time order. Let $\mathrm{pa}(n)$ be the parents of node $n$
over all edge types. Then

$$e(n) = \begin{cases} \text{new id} & \mathrm{pa}(n) = \emptyset \\ e(p^\star),\ p^\star = \arg\max_{p \in \mathrm{pa}(n)} (A_p, -\mathrm{id}_p) & \text{otherwise.} \end{cases}$$

So split children share their parent's event, and a merge continues the
largest parent, as in S1. In the data tested so far, the number of such
events equals the number of families (section 9, Q6).

### 6.2 Factors (Chang et al. 2016, Section 4b)

For any event labelling $e$, let $\mathcal{N}$ be the set of events and
$l(i)$ the number of frames in which event $i$ has at least one node. With
cell area $a = \Delta x^2$:

$$\text{Intensity} = \frac{\sum_n P_n}{\sum_n A_n}\ \text{(mm/h)}, \qquad \text{Size} = a\,\frac{\sum_n A_n}{\sum_i l(i)}\ (\text{km}^2),$$

$$\text{Duration} = \Delta t\,\frac{\sum_i l(i)}{\lvert \mathcal N\rvert}\ \text{(h)}, \qquad \text{Number} = \lvert \mathcal N\rvert,$$

$$\text{Intensity} \times \text{Size} \times \text{Duration} \times \text{Number} = a\, \Delta t \sum_n P_n \quad (\text{mm km}^2).$$

Intensity is independent of tracking. Size, Duration and Number depend on
tracking only through $\lvert \mathcal N\rvert$ and $\sum_i l(i)$. The units
are mm/h and km$^2$, never mm/h/km$^2$; 1 mm km$^2$ = 1000 m$^3$.

## 7. Reference and hourly paper trackers

### 7.1 Chang et al. (2016) S1 (`chang2016.track`)

Weights $w_n(s) = r(s)/P_n$. The similarity between segments $a$ and $b$ is

$$\mathrm{sim}(a, b) = \sum_{s \in M_a} \sum_{u \in M_b} w_a(s)\, w_b(u)\, e^{-\phi\, \Delta x \lVert s - u\rVert},$$

with $\phi = 1/120$ km$^{-1}$. A link from segment $j$ at $t$ to $p$ at
$t-1$ requires $\mathrm{sim} > \tau$ (0.05 in JJA, 0.01 in DJF), and either
$\Delta x\lVert c_j - c_p\rVert < 120$ km or the angle between
$c_j - c_p$ and $c_p - c_{p^-}$ is $< 120^\circ$. Among the qualifying
parents, the one with the most cells is chosen. The implementation decides
$\mathrm{sim} > \tau$ exactly, without summing the kernel for most pairs,
using certified bounds:

$$e^{-\phi \Delta x\, g_{ab}} \ \ge\ \mathrm{sim}(a, b)\ \ge\ e^{-\phi \Delta x\,(\lVert c_a - c_b\rVert + \bar d_a + \bar d_b)},$$

where $g_{ab}$ is the bounding-box gap and
$\bar d_n = \sum_s w_n(s) \lVert s - c_n \rVert$. The lower bound follows from
Jensen's inequality (the kernel is convex) and the triangle inequality.
Blocked versions of both bounds (32 and 8 cells) are tried before the exact
sum, which uses an untruncated FFT convolution for large pairs.

### 7.2 Hourly variant (`chang2016_hourly.track_hourly`)

Same similarity and event semantics. The motion rule is replaced by the
mask distance

$$d_M(j, p) = \Delta x \min_{s \in M_j,\, u \in M_p} \lVert s - u\rVert \le D \qquad (D = 40 \text{ km}),$$

and the parent with the largest $\mathrm{sim}$ is chosen (ties: smaller
$d_M$, then larger area). A child with no parent at $t-1$ may link to a
segment whose event ended at $t-1$, within $2D$.

## 8. Recommended configuration

Identification as in section 3.2 with $(\theta, \rho_c, \rho_b, e) = (1\
\text{mm/h}, 3, 9, 15)$. Tracking as in section 4 with the preset
parameters, plus `overlap_metric='child_coverage'` and
`share_parent=True`. Report factors with branch ids, family ids or
paper-style events (section 6.1).

## 9. Open questions for discussion

These came up while writing the equations. None is a coding error, but
each affects results or their interpretation.

* **Q1. Grid dependence of the fixed minima.** The 16-cell minima in
  sections 4.5-4.7 are hard-coded, while $\kappa$ and $\rho$ are
  configurable cells. The minima are 256 km$^2$ at 4 km but 2304 km$^2$ at
  12 km. Should they be given in km$^2$?
* **Q2. Uncalibrated score.** The weights (0.60, 0.25, 0.15) and $\tau =
  0.35$ are fixed by design, not fitted. In the distance term the length
  scale is the search radius $\tilde R = \kappa k$, so a gap ($k = 2$)
  halves the distance penalty, while $\tau_{gap} = 0.45$ is stricter.
  Is that intended?
* **Q3. Topology ignores advection.** Split/merge and strong-overlap edges
  use raw overlap only (section 4.5). Advection enters only the score. A
  fast-moving system that splits can therefore fail the topology tests even
  though its advected overlap is large.
* **Q4. Split criterion at 1 h.** Requiring $\min(\pi_p, \pi_j) \ge 0.15$
  and $\max \ge 0.5$ means that small fragments never form splits. In a
  simple symmetric test case (parent of radius 12 cells, two children of
  radius 7), no edge was formed at all. `share_parent` covers the small
  fragments, but the symmetric case may need a looser rule.
* **Q5. Erosion border.** With `erosion_border='ignore'`, cells outside the
  domain count as wet, so regions touching the border pass the large-region
  test more easily. This matches the original STEP. The edge band
  (section 3.2) makes the choice mostly irrelevant.
* **Q6. Families vs events.** On the tested data every merge re-joined
  objects from a common ancestor, so families equal paper-style events.
  This may not hold for long runs or other thresholds. Families never split
  and can only grow with record length.
* **Q7. Threshold.** $\theta = 1$ mm/h removes 16-25% of the rain, and the
  removed share changes with climate. The decomposition in section 6.2 is
  then of rain above 1 mm/h, not of total rain.
