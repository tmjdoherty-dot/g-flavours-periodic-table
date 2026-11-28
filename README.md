# The periodic table with g electron spin flavours

Author: Thomas Doherty, with computational and theory assistance from Claude (Anthropic)

## Summary

Electrons have two spin states, and every shell capacity in chemistry, 2(2l+1), quietly depends on that. This repository asks what the periodic table of free atoms looks like when electrons have **g** internal states (g = 1, 2, 3, …, ∞). It is a thought experiment. The aim is to separate what in the ground configurations comes from the flavour-blind direct Coulomb term and what comes from same-flavour terms: exchange, which in mean-field theory also cancels each electron's interaction with itself.

Rescaling turns the g-flavour atom into a one-parameter family. In per-flavour units the problem is the g = 1 atom with nuclear charge z = Z/g, and exchange enters multiplied by **λ = 1/g**:

- g → ∞ (λ = 0) is direct Coulomb only: the Hartree limit, which becomes exact as g → ∞ (self-interaction and correlation both vanish there);
- g = 2 (λ = ½) is real electrons;
- g = 1 (λ = 1) is spinless electrons.

Main results, from exchange-only LDA scans and then Hartree–Fock (HF) with exact exchange:

1. **Each block starts at the same z = Z/g for every g** (within ±0.85 per-flavour electrons), as Thomas–Fermi scaling predicts. Stronger exchange brings d and f shells in slightly earlier.
2. **The Cr/Cu-type and 4d-row exceptions of free atoms come from exchange between different orbitals, not from self-interaction removal.**
   - At g = 2, HF with maximal flavour polarisation (Hund) gets 11 of 12 real ground-state s occupancies in each d row, including Cr, Cu, Nb, Mo, Ru, Rh, Ag and Pd 4d¹⁰5s⁰.
   - With direct Coulomb only (self-interaction removed), the same atoms have no exceptions in the 3d row and only Pd in the 4d row (Result 6).
   - Configuration-average HF matches the lowest configuration averages of Wang et al. 2006 (Table 1) for 8/10 (3d; Fe and Co come out dⁿs² rather than dⁿ⁺¹s¹) and 10/10 (4d). As a configuration average Cr is d⁴s², as in that table; d⁵s¹ appears only for the high-spin (Hund) part, i.e. through the exchange stabilisation of the half-filled shell.
   - Across g, the 3d exception share in Hund mode is 0.80 / 0.30 / 0 / 0 / 0 / 0 / 0 for g = 1 / 2 / 3 / 4 / 5 / 6 / 8.
3. **The large-g limit keeps some 4d exceptions; real atoms without exchange do not.**
   - The 4d exception share falls 0.80 → 0.70 → 0.60 → 0.45 → 0.44 → 0.37 → 0.38 as g grows. It levels off near the direct-Coulomb (Hartree) drain at g → ∞ (0.35). The self-interaction-free direct-Coulomb model rises towards the same region from below (0.10 → 0.25).
   - The dial 1/g scales exchange and the self-interaction correction together. So the large-g limit is not the same as "the real atom without exchange", and the large-g 4d exceptions say little about real 4d atoms.
4. **With direct Coulomb only, 3d is not yet bound when 4s fills.** At λ = 0 the 3d level is an unbound continuum state until z ≈ 10.5, so 4p starts filling first. The head start of p before d grows smoothly with g in both LDA and HF. This is the d-orbital collapse after the ns shell fills (Connerade 1991; Chen, Li & Schwarz 2021), happening later when the same-flavour terms are switched off.
5. **The La/Ac step (one d electron before the f shell) depends on g.**
   - HF: absent at g = 1 (4f lower by 3.7 eV), present at g = 2 by only 0.08 eV, and strong for g ≥ 3 (3–9 eV).
   - At g = 2, HF matches 13 of 14 lanthanide ground configurations.
   - The real free La atom sits just on the d-first side of a boundary.
6. **Splitting the Coulomb terms at g = 2.** Self-interaction-free direct Coulomb → + exchange between orbitals → + flavour polarisation takes the d-block exception share from 0 → 0.20 → 0.30 (3d) and 0.10 → 0.60 → 0.70 (4d). Direct Coulomb *with* self-interaction (plain Hartree) does not even bind the outermost electron of a neutral atom at g = 2, so that split is ill-defined for real atoms.
7. **Cations follow the regular dⁿ s⁰ pattern (all valence electrons in d) only for small g.**
   - With one electron removed per flavour (M²⁺ at g = 2), HF gives d^n s⁰ for every 3d and 4d dication (Hund and configuration average), as observed and as already stated for free cations by Wang et al. 2006. The new part is what follows.
   - Direct Coulomb alone would leave s electrons in 40% (3d) and 30% (4d) of the dications, so the regular cation pattern also needs exchange.
   - From g = 4 on, cations increasingly keep s electrons (3d row, Hund: 0.30 / 0.57 / 0.65 at g = 4 / 6 / 8).
8. **Which exchange produces which exception (g = 2, configuration average, self-interaction removed).**
   - With direct Coulomb plus core–valence exchange alone, the 4d row already has Tc, Ru, Rh, Pd and Ag as d^(n+1)s¹ (5 of the 6 full-HF exceptions), while the 3d row has none.
   - The 3d row needs core–valence and within-shell (d–d) exchange together. That pair reproduces Table 1 of Wang et al. 2006 exactly (Fe–Cu as d^(n+1)s¹), and d–s exchange then pushes Fe and Co back to s².

9. **Core–valence exchange acts mainly through orbital size.** At g = 2 it shrinks ⟨r⟩ of (n−1)d by 23–24% on average and leaves ns almost unchanged (≤ 0.8%). Without it the 3d of Sc is as diffuse as 4s (3.88 vs 3.94 bohr; 1.63 bohr with it). About 70–75% of its effect on the d-vs-s preference comes from this relaxation. Size is necessary but not sufficient: at equal d/s size ratios, configurations with valence exchange move electrons from s to d far more often than without (Result 9).
10. **Whole-flavour locking does not grow with the number of flavours.** The extra stabilisation of whole-flavour d counts by flavour polarisation, relative to other counts, falls from 1.74 at g = 2 to 1.12 at g = 8 (3d row; 1.71 → 1.10 in the 4d row). Half- and whole-filled shells are most special at g = 2.
11. **Adding LDA correlation leaves the split intact.** In the 3d row, no combination without d–d exchange gives an exception, with or without correlation. Full HF plus LDA correlation matches the lowest configuration averages of Wang et al. 2006 for all of Ca–Cu.
12. **At g = 2 levels decide; with more flavours the averages do (3d row fully, 4d row partly).** The high-spin gain differs between competing configurations by 2.4× (3d) and 2.6× (4d) the gap between their configuration averages, so the lowest level, not the lowest average, picks the ground configuration (as Wang et al. 2006 state). The gain differences fall roughly as 1/g. In the 3d row the Hund and average ground configurations agree for every g ≥ 3. In the 4d row the averages stay nearly degenerate (6–8 mHa per flavour) and about a quarter of the row still differs at g = 8.
13. **A confining wall makes every g regular, at a radius fixed in per-flavour units.**
    - A soft spherical wall around the neutral atom drives the ground configuration of every d-block atom to the regular one, with all valence electrons in d (Cu-type d¹⁰s¹ and Zn-type d¹⁰s² aside). This holds for every g from 1 to 6, in Hund, average and self-interaction-free direct-Coulomb modes.
    - The radius where it happens is the same for every g in per-flavour units: about 6 for the 3d row and 8 for the 4d row, about the size of the free atom's outer s orbital (⟨r⟩ ≈ 7–8 per flavour for every g ≥ 2). At g = 2 that is 3 and 4 bohr (1.6 and 2.1 Å).
    - Exchange shifts the onset only by about one radius step. At g = 2, core–valence plus d–d exchange makes the 3d row regular at a larger radius, and d–s exchange delays it. Inside the wall the half-filled Hund configurations (Cr d⁵s¹, Mn d⁵s²) keep their s electrons longest.
    - Between free and confined there is a mixed stage. At g = 2 it is d^(n+1)s¹ for most of the row, like the metal-like pattern. At g ≥ 3 the s shell drains one electron at a time with no preferred count.

No earlier treatment of a g-flavour periodic table was found (searches on spinless, SU(N) and spin-degeneracy variants of the periodic table and on large-N atoms; see *Related work*).

## Scope and terminology
- **Free atoms only.** Every result concerns the ground configuration of an isolated atom in vacuum. Atoms bound in compounds and metals behave far more regularly. Cations in complexes are close to (n−1)d^(G−q) ns⁰, and metals are close to (n−1)d^(G−1) ns¹ throughout the row, whatever the free atom does (Chen, Li & Schwarz 2021; G is the group number and q the charge, not the flavour count g). Result 13 imitates the bound environment crudely with a confining wall. The free-atom irregularities studied here therefore have little direct relevance to chemistry. The question is where they come from.
- **Flavours, not spin.** The g states are an internal degeneracy with SU(g) symmetry, like the nuclear-spin flavours of fermionic Sr or Yb in cold-atom experiments, while the particles stay fermions. They are not electrons of spin (g−1)/2: real particles of integer spin would be bosons and obey different wave equations (Proca for spin 1, Rarita–Schwinger for spin 3/2).
- **Configurations, not states.** "Ground configuration" means the leading configuration of the lowest state. HF energies here are averages over a configuration or over its maximally flavour-polarised (high-spin) part, not individual angular-momentum-coupled states.
- **Madelung as bookkeeping.** "Exception" means a ground configuration that differs from "ns full before (n−1)d". This is used as a counting reference only. The (n+l, n) rule is not claimed to describe orbital energies; for the d block it does not (Wang, Qiu, Fang & Schwarz 2006; Chen, Li & Schwarz 2021).
- **How the Coulomb energy is split.** The Coulomb interaction has direct, exchange and correlation parts. The split between direct and exchange is not unique: the spurious self-interaction can be added to the direct term and subtracted in the exchange term. Here, "same-flavour terms" means everything multiplied by 1/g in the rescaled mean-field energy, including the self-interaction cancellation. Result 6 separates the two by removing self-interaction from the direct term first (self-interaction is not part of the original Hartree or Hartree–Fock equations; it appears only when the self-exchange term is added and subtracted). Result 8 splits the exchange further. Correlation is not included; at g → ∞ it vanishes, at g = 2 it does not.

## The one-dial rescaling

For g flavours, each spatial orbital holds g electrons. With s = g·r, energies divided by g², and z = Z/g, the Hamiltonian becomes

$$
\frac{H}{g^2} = \sum_i \left[ -\tfrac12 \nabla_{s_i}^2 - \frac{z}{s_i} \right] + \frac{1}{g} \sum_{i<j} \frac{1}{|s_i - s_j|}.
$$

In a mean-field treatment with per-flavour density ρ (integrating to z), the direct (Hartree) term is independent of g. Exchange acts only within a flavour and carries the factor 1/g. With Slater (LDA) exchange, the Kohn–Sham potential is

$$
v(s) = -\frac{z}{s} + v_H[\rho](s) - \lambda \left(\frac{6}{\pi}\right)^{1/3} \rho(s)^{1/3}, \qquad \lambda = \frac{1}{g}.
$$

So the whole family is one problem with one dial. Thomas–Fermi theory has no shells: it puts the first appearance of angular momentum l at Z_l ∝ g(2l+1)³, i.e. at the same z for every g. Any g-dependence of the filling order is therefore a shell or exchange effect.

## Methods

**LDA scans (`atom.py`, `g1_*.py`)**
- Spherical, non-relativistic atom with exchange-only LDA, on a logarithmic radial grid.
- Occupations found by direct minimisation over fractional shell occupations. Janak's theorem gives ∂E/∂f = ε, so projected-gradient steps are used. Partly filled shells share the Fermi level.
- The plain Fermi-filling version did not converge for transition metals (charge sloshing between 3d and 4s); the minimisation does.
- Scan: z = 0.25–60 for g = ∞, 4, 3, 2, 1, plus g = 2 with PZ81 correlation as a reference row.
- Two independent runs with different solver settings agree on the drain table exactly and on the onsets to within 0.16 in z (the onset table uses run 1).

**Hartree–Fock (`hf.py`)**
- Slater's average-of-configuration HF, generalised to g flavours. For shell a with q_a electrons:

$$
E = \sum_a q_a I_a + \frac1g \Big[ \sum_a \tfrac{q_a(q_a-1)}{2} E_{aa} + \sum_{a<b} q_a q_b E_{ab} \Big],
$$

$$
E_{aa} = \frac{2l+1}{g(2l+1)-1} \Big[ g F^0_{aa} - \sum_k \begin{pmatrix} l & k & l \\ 0&0&0 \end{pmatrix}^{2} F^k_{aa} \Big], \qquad
E_{ab} = F^0_{ab} - \frac1g \sum_k \begin{pmatrix} l_a & k & l_b \\ 0&0&0 \end{pmatrix}^{2} G^k_{ab}.
$$

  For g = 2 this is Slater's average-of-configuration energy.
- **Hund mode:** a single flavour-resolved configuration with maximal flavour polarisation; each shell fills flavour 1 to 2l+1 before flavour 2. At g = 2 this is the high-spin term, e.g. Cr ⁷S and Fe ⁵D.
- **Solver:** direct minimisation on the radial grid. Occupied–virtual steps use an Olsen (Jacobi–Davidson-type) correction with a tridiagonal LDA reference Hamiltonian, falling back to plain gradient steps when that direction is not downhill. Rotations between occupied shells of equal l and different occupation are handled explicitly. Steps use Löwdin re-orthonormalisation and backtracking on E.
- **Grid:** r_min = 10⁻⁹. Total energies on N = 1600 and 2400 points, extrapolated with the h² law.
- **Starting guess:** an LDA calculation with the configuration's own occupations held fixed. This matters for double-well f orbitals: with an empty-4f guess, Ce converged to the extended 4f solution, 0.14 Ha too high.
- **Modes:** configuration average, Hund, self-interaction-free direct Coulomb (`mode='sicfree'`: every exchange term dropped, pair counts q(q−1)/2, so no electron interacts with itself), and `mode='split'` (configuration average with core–core, core–valence, d–s and within-shell valence exchange switched separately; all on = configuration average, all off = 'sicfree', both checked to machine precision).
- **Correlation (optional):** PZ81 LDA correlation of the total density added to the HF energy, at a chosen strength (g = 2 only).
- **Confining wall (optional):** V(s) = V₀ / (1 + exp(−(s − R_w)/w)) added to the one-electron operator, in per-flavour units with V₀ = 0.5 and w = 0.5 and R_w varied (free, 12, 10, 8, 7, 6, 5). At g = 2 this is a wall of height 2 Ha and width 0.25 bohr at radius R_w/2 bohr. The starting orbitals include the wall. The wall runs use a single grid (N = 2000) without extrapolation, since only configuration differences are compared.
- **Step cap:** each preconditioned orbital step is capped at norm 0.5 (orbitals have norm 1), with convergence measured before the cap. Inside a wall, orbital energies can be positive and the preconditioner nearly singular. In free atoms the cap never acts.
- **Runs:** ~76,000 HF calculations on Modal, including a first wall scan that was repeated after the starting-orbital fix (`modal_jobs/hf_modal.py`).

## Validation

| check | result |
|---|---|
| LDA vs NIST LDA reference, Ne and Ar orbital energies | agree to ≤ 3.6 mHa |
| LDA Fe, Cu (fractional aufbau) | 3d⁶·⁹4s¹·¹, 3d¹⁰4s¹ (standard spherical-LDA result) |
| HF gradients vs finite differences | ratio 1 ± 10⁻⁶ |
| HF total energies (extrapolated) vs numerical-HF reference: Ne, Ar, Zn | +1·10⁻⁵, +3·10⁻⁵, +1·10⁻⁴ Ha |
| HF Cr 3d⁵4s¹ ⁷S, Fe 3d⁶4s² ⁵D (Hund mode) | +3.4·10⁻⁵, +4.4·10⁻⁵ Ha |
| Part A with two different starting guesses | identical ground states |
| Step cap on vs off (Cu, Fe ⁵D, Mo d⁵s¹; free atoms) | bit-identical energies and iteration counts |
| Wall runs: first attempt (wall missing from the starting orbitals) vs final | the first attempt stalled at the start in 40% of runs (now 1.8%, nearly all free g = 6 atoms). Where both converged, the g = 2 ground configurations agree |

## Results

### 1. Block onsets (LDA)
Persistent onset of the first d and f shells, in per-flavour z (the real element is Z = g·z). Madelung values are 10, 19, 28 and 44.

| shell | g = ∞ | 4 | 3 | 2 | 2 + corr | 1 |
|---|---|---|---|---|---|---|
| 3d | 10.58 | 10.33 | 10.25 | 10.17 | 10.5 | 9.67 |
| 4d | 19.67 | 19.42 | 19.33 | 19.17 | 19.5 | 18.92 |
| 4f | 29.0 | 28.75 | 28.58 | 28.33 | 28.5 | 27.58 |
| 5f | 45.2 | 44.83 | 44.67 | 44.4 | 44.5 | 43.67 |

### 2. Outer-s drain (LDA, fractional occupations)
Largest drop of the outer s shell inside each block, per flavour (1 = the s shell empties):

| block | g = ∞ | 4 | 3 | 2 | 1 |
|---|---|---|---|---|---|
| 3d (4s) | 0.00 | 0.17 | 0.31 | 0.58 | 1.00 |
| 4d (5s) | 0.35 | 0.77 | 0.92 | 1.00 | 1.00 |
| 4f + 5d (6s) | 0.43 | 0.83–0.87 | 1.00 | 1.00 | 1.00 |

Fractional HF scans (configuration average) give the same picture in the 4d row: 1.00 / 0.89 / 0.68 / 0.58 / 0.49 / 0.45 / 0.38 for g = 1 / 2 / 3 / 4 / 6 / 8 / 16. These scans are qualitative: about 45% of points end with residual gradients above 10⁻⁴, mostly in empty orbitals.

### 3. Ground configurations from HF (integer configurations)
For each element, every configuration [core] (n−1)dˣ nsʸ with y = 0…g was computed, and the lowest was taken. An exception is a ground state that is not "s full first".

Share of the d block with an exception:

| row | mode | g = 1 | 2 | 3 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|---|---|
| 3d | Hund | 0.80 | 0.30 | 0 | 0 | 0 | 0 | 0 |
| 3d | average | 0.80 | 0.20 | 0 | 0 | 0 | 0 | 0 |
| 4d | Hund | 0.80 | 0.70 | 0.60 | 0.45 | 0.44 | 0.37 | 0.38 |
| 4d | average | 0.80 | 0.60 | 0.53 | 0.45 | 0.40 | 0.40 | 0.38 |

Agreement with experiment at g = 2 (ground-state s occupancy):

| row | Hund | configuration average |
|---|---|---|
| 3d (K–Zn) | 11/12 (V wrong) | 10/12 (Cr, Ni wrong) |
| 4d (Rb–Cd) | 11/12 (Zr wrong; Pd 4d¹⁰5s⁰ right) | 9/12 (Nb, Tc, Rh wrong) |

Whether ground states lock at whole flavours (d counts that are multiples of 5) is tested in Result 10. The locking fades with g.

### 4. p before d at large g
- **LDA, λ = 0:** 4p fills to 0.49 per flavour (converged, residual 5·10⁻⁸) before 3d is bound. At λ = ¼ the 4p peak is 0.22; at λ = ½ there is none.
- **HF:** the p head start in the 4d row is 0.04 / 0.21 / 0.29 / 0.38 / 0.39 / 0.42 for g = 2 / 3 / 4 / 6 / 8 / 16.
- **Hysteresis:** up and down sweeps in z through the 3d onset agree exactly at λ = ¼ and ½. At λ = 0 they differ only at points that never converged, so the collapse cannot be called first-order.

### 5. The La/Ac step (HF, Hund mode)
Configurations [Xe] 6s^g 4fˣ 5dʸ with y ≤ 2:

| g | first element takes 5d? | E(4f¹) − E(5d¹) | elements with 5d occupied |
|---|---|---|---|
| 1 | no | −3.65 eV | 0/7 |
| 2 | yes | +0.08 eV | 2/14 (La, Gd; 13/14 match experiment, Ce missed) |
| 3 | yes | +8.62 eV | 15/21 |
| 4 | yes | +5.52 eV | 28/28 |
| 6 | yes | +3.26 eV | 40/42 |

Whole-flavour f counts (multiples of 7, a generalised Gd) are favoured above chance at g = 2–4 and less so at g = 6 (Result 10). At g ≥ 4 the ground state often sits at the y = 2 cap, so states with more 5d electrons may be lower.

### 6. Splitting the Coulomb terms (HF, neutral atoms)
Share of the d block with an exception, by model:

| row | model | g = 1 | 2 | 3 | 4 | 6 | 8 |
|---|---|---|---|---|---|---|---|
| 3d | direct Coulomb, self-interaction removed | 0.00 | 0.00 | 0.07 | 0.05 | 0.10 | 0.10 |
| 3d | + exchange (configuration average) | 0.80 | 0.20 | 0 | 0 | 0 | 0 |
| 3d | + flavour polarisation (Hund) | 0.80 | 0.30 | 0 | 0 | 0 | 0 |
| 4d | direct Coulomb, self-interaction removed | 0.00 | 0.10 | 0.13 | 0.25 | 0.20 | 0.25 |
| 4d | + exchange (configuration average) | 0.80 | 0.60 | 0.53 | 0.45 | 0.40 | 0.38 |
| 4d | + flavour polarisation (Hund) | 0.80 | 0.70 | 0.60 | 0.45 | 0.37 | 0.38 |

- At g = 2 the self-interaction-free direct-Coulomb model gets 10/12 (3d) and 6/12 (4d) real s occupancies. It puts every 3d atom at 4s², and every 4d atom at 5s² except Pd (5s¹).
- At small g exchange creates the exceptions. At g ≥ 3 in the 3d row it removes the few that the direct-Coulomb model has.
- 474 of 3176 direct-Coulomb runs did not converge and are excluded.

### 7. Cations (HF, one electron removed per flavour)
Share of the d block whose ground configuration keeps s electrons:

| row | model | g = 1 | 2 | 3 | 4 | 6 | 8 |
|---|---|---|---|---|---|---|---|
| 3d | direct Coulomb, self-interaction removed | 0.00 | 0.40 | 0.60 | 0.75 | 0.80 | 0.82 |
| 3d | configuration average | 0.00 | 0.00 | 0.00 | 0.20 | 0.40 | 0.53 |
| 3d | Hund | 0.00 | 0.00 | 0.00 | 0.30 | 0.57 | 0.65 |
| 4d | direct Coulomb, self-interaction removed | 0.00 | 0.30 | 0.47 | 0.45 | 0.47 | 0.50 |
| 4d | configuration average | 0.00 | 0.00 | 0.13 | 0.20 | 0.33 | 0.38 |
| 4d | Hund | 0.00 | 0.00 | 0.07 | 0.30 | 0.37 | 0.42 |

- At g = 2, every dication Sc²⁺–Zn²⁺ and Y²⁺–Cd²⁺ comes out d^n s⁰, as observed.
- All 8604 cation runs converged.

### 8. Exchange split (HF, g = 2, configuration average)
Direct Coulomb with self-interaction removed, core–core exchange on; then core–valence (CV), d–s, and within-shell valence (d–d; s–s exchange is zero for g = 2) exchange switched on or off. Exceptions are counted against "ns full first" (energy gain in eV); agreement is with Table 1 of Wang et al. 2006 (lowest configuration averages, groups 2–11).

| exchange on (besides core–core) | 3d exceptions | 3d vs Table 1 | 4d exceptions | 4d vs Table 1 |
|---|---|---|---|---|
| none | – | 6/10 | Rh, Pd | 4/10 |
| CV | – | 6/10 | Tc, Ru, Rh, Pd, Ag | 7/10 |
| d–s | – | 6/10 | Pd | 4/10 |
| within-shell | – | 6/10 | Ru, Rh, Pd | 6/10 |
| CV + d–s | – | 6/10 | Ru, Rh, Pd, Ag | 6/10 |
| CV + within-shell | Fe, Co, Ni, Cu (0.10–1.13) | **10/10** | Nb, Mo, Tc, Ru, Rh, Pd, Ag | 8/10 |
| d–s + within-shell | – | 6/10 | Ru, Rh, Pd | 5/10 |
| all (full HF) | Ni, Cu (0.13, 0.37) | 8/10 | Mo, Tc, Ru, Rh, Pd, Ag | **10/10** |

- In the 4d row, core–valence exchange alone does most of the work. It acts on the valence orbitals like part of the core potential, which fits the picture in which orbital size and core screening decide (compact (n−1)d vs diffuse ns, Wang et al. 2006).
- In the 3d row it is not enough, and the exceptions of the configuration averages also need the d–d exchange.
- 6 of 1024 runs did not converge and are excluded.

### 9. Orbital size (HF, g = 2 and across g)
Core–valence (CV) exchange switched on, starting from direct Coulomb (self-interaction removed) with core–core exchange, valence–valence exchange off. dE = E(d^(n+1)s¹) − E(dⁿs²). The CV effect on dE is split into the exchange energy at frozen orbitals (X) and orbital relaxation (R):

| row | mean shrink of ⟨r⟩ (n−1)d | mean shrink of ⟨r⟩ ns | CV effect on dE | X | R |
|---|---|---|---|---|---|
| 3d (Sc–Cu) | 24% | 0.1% | −2.03 eV | −0.51 eV | −1.53 eV |
| 4d (Y–Pd) | 23% | 0.7% | −1.36 eV | −0.40 eV | −0.97 eV |

Ag is excluded: its no-CV d¹⁰s¹ run converged to a different state. Pd is the one element where X dominates.

Across g (per-flavour units), the d/s size ratio ⟨r⟩_d/⟨r⟩_s of the s-full configuration and the share of the d block where moving one s electron to d lowers the energy:

| row | model | g = 1 | 2 | 3 | 4 | 6 | 8 |
|---|---|---|---|---|---|---|---|
| 3d | direct Coulomb only: ratio / share | 0.36 / 0 | 0.40 / 0 | 0.40 / 0 | 0.44 / 0 | 0.49 / 0 | 0.88 / 0 |
| 3d | configuration average: ratio / share | 0.33 / 1.00 | 0.35 / 0.22 | 0.37 / 0 | 0.40 / 0 | 0.67 / 0 | 0.78 / 0 |
| 4d | direct Coulomb only: ratio / share | 0.48 / 0 | 0.52 / 0.11 | 0.60 / 0.14 | 0.55 / 0.22 | 0.83 / 0.21 | 0.80 / 0.24 |
| 4d | configuration average: ratio / share | 0.43 / 1.00 | 0.46 / 0.67 | 0.48 / 0.57 | 0.51 / 0.47 | 0.73 / 0.41 | 0.82 / 0.38 |

The d orbital grows relative to s as g rises (exchange carries 1/g), and the share falls with it. But at equal ratio (~0.5 in the 4d row) the share is 0.11 without valence exchange and 0.47 with it.

### 10. Whole-flavour locking vs g (HF, Hund vs configuration average)
Flavour-polarisation gain (average minus Hund energy of the same configuration) at whole-flavour d counts divided by the gain at other counts:

| row | g = 2 | 3 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|
| 3d | 1.74 | 1.48 | 1.23 | 1.18 | 1.22 | 1.12 |
| 4d | 1.71 | 1.37 | 1.23 | 1.19 | 1.15 | 1.10 |

- The share of d-block ground states at whole counts is the same in Hund and average modes for the 3d row at g ≥ 3 (0.14–0.18), and ~0.1 higher in Hund mode for the 4d row at every g.
- In pair counting, every whole-flavour boundary costs the same 2l same-flavour partners at any g. With more flavours, though, every count is already mostly polarised, so whole counts become less special.

### 11. Exchange split with LDA correlation (HF + PZ81, g = 2)
Same split as Result 8, with PZ81 correlation of the total density at half and full strength (LDA overestimates atomic correlation roughly twofold: Ne −0.74 vs −0.39 Ha).

| exchange on (besides core–core) | 3d exceptions, corr ½ | 3d exceptions, corr 1 | 3d vs Table 1, corr 1 | 4d exceptions, corr 1 | 4d vs Table 1, corr 1 |
|---|---|---|---|---|---|
| none | – | – | 6/10 | Rh, Pd | 4/10 |
| CV | V | – | 6/10 | Mo, Tc, Ru, Rh, Pd, Ag | 8/10 |
| CV + d–s | – | – | 6/10 | Tc, Ru, Rh, Pd, Ag | 7/10 |
| CV + within-shell | Fe, Co, Ni, Cu | Mn, Fe, Co, Ni, Cu | 9/10 | Nb–Ag | 7/10 |
| all (full HF) | Co, Ni, Cu | Fe, Co, Ni, Cu | **10/10** | Nb–Ag | 8/10 |

11 of 3072 runs did not converge.

### 12. Levels vs configuration averages across g (HF, Hund vs average)
The J-level spreads of competing configurations overlap (Wang et al. 2006, Fig. 8), so the ground state can be set by the lowest level rather than by the lowest configuration average. In HF terms, the high-spin gain G = E_avg − E_Hund of a configuration plays the role of the spread. For each d-block atom, the gap between the two lowest configuration averages is compared with the difference in G between them (|ΔG|). A "flip" is an atom whose Hund ground configuration differs from its average ground configuration.

High-spin gains at g = 2, compared with the level positions in Wang et al. 2006 (Cr ⁷S more than 4 eV below the d⁵s¹ average; lowest d⁴s² level about 2.5 eV below its average):

| atom | dⁿ⁺²s⁰ | dⁿ⁺¹s¹ | dⁿs² |
|---|---|---|---|
| Cr | 2.67 eV | 5.84 eV | 3.64 eV |
| Mo | 2.27 eV | 4.92 eV | 2.75 eV |
| Mn | 1.48 eV | 3.86 eV | 6.54 eV |
| Ni | 0.00 eV | 0.09 eV | 0.77 eV |

HF overshoots the Cr gains by about 40%, as expected without correlation, but the ordering is right: the half-filled d⁵s¹ gains most.

Across g (per-flavour mHa, medians over the d block):

| row | g | flip share | gap between averages | \|ΔG\| | \|ΔG\| / gap |
|---|---|---|---|---|---|
| 3d | 2 | 0.33 | 5.9 | 12.4 | 2.35 |
| 3d | 3 | 0 | 13.7 | 7.6 | 0.49 |
| 3d | 4 | 0 | 16.1 | 5.1 | 0.29 |
| 3d | 5 | 0 | 17.0 | 3.8 | 0.20 |
| 3d | 6 | 0 | 17.4 | 3.0 | 0.16 |
| 3d | 8 | 0 | 17.7 | 2.1 | 0.11 |
| 4d | 2 | 0.44 | 6.6 | 10.4 | 2.55 |
| 4d | 3 | 0.14 | 6.2 | 5.9 | 0.81 |
| 4d | 4 | 0.26 | 6.1 | 3.3 | 0.39 |
| 4d | 5 | 0.21 | 7.0 | 2.7 | 0.43 |
| 4d | 6 | 0.21 | 7.5 | 2.5 | 0.42 |
| 4d | 8 | 0.26 | 8.0 | 1.6 | 0.22 |

- At g = 2 the high-spin gain differences are more than twice the gaps between averages, in both rows. The lowest level, not the lowest average, decides, as stated by Wang et al. 2006.
- The gain differences fall roughly as 1/g (g·|ΔG| goes 25 → 17 mHa from g = 2 to 8), since the flavour-polarisation energy is an exchange term.
- 3d row: the gaps between averages grow from 6 to 18 mHa and the ground configuration equals the average one from g = 3 on.
- 4d row: the gaps stay at 6–8 mHa, so a quarter of the row still flips at g = 8. The ratio is falling but convergence to the average configuration is incomplete in the range computed. The 4d configurations stay nearly degenerate at every g, which is also why the 4d exception share levels off rather than vanishing (Result 3).
- g = 1: one flavour, so Hund and average coincide.

### 13. Neutral atoms inside a confining wall (HF, all g)
A crude stand-in for a bound atom: the soft wall of *Methods* around the neutral atom, in the spirit of the He-cage calculation for Mn⁺ in Wang et al. 2006 (Fig. 10). Every configuration [core] (n−1)dˣ nsʸ was computed at each wall radius. "Regular" means the ground configuration has the fewest s electrons allowed: s⁰, or d¹⁰s¹ / d¹⁰s² for the Cu- and Zn-type atoms (d^(5g) plus 1 to g s electrons). An element counts only if every configuration converged.

Share of the d block with a regular ground configuration, by wall radius R_w (per flavour; real radius at g = 2 is R_w/2 bohr):

| row | mode | g | free | 12 | 10 | 8 | 7 | 6 | 5 |
|---|---|---|---|---|---|---|---|---|---|
| 3d | Hund | 2 | 0.20 | 0.20 | 0.20 | 0.50 | 0.50 | 1.00 | 1.00 |
| 3d | Hund | 3 | 0.07 | 0.07 | 0.07 | 0.13 | 0.27 | 0.87 | 1.00 |
| 3d | Hund | 4 | 0.06 | 0.06 | 0.05 | 0.10 | 0.20 | 0.80 | 1.00 |
| 3d | Hund | 6 | 0.12 | 0.05 | 0.04 | 0.08 | 0.27 | 0.92 | 1.00 |
| 3d | direct Coulomb only | 2 | 0.10 | 0.10 | 0.10 | 0.10 | 0.20 | 0.80 | 1.00 |
| 3d | direct Coulomb only | 4 | 0.08 | 0.06 | 0.05 | 0.05 | 0.20 | 0.70 | 1.00 |
| 3d | direct Coulomb only | 6 | 0.20 | 0.07 | 0.06 | 0.04 | 0.20 | 0.76 | 1.00 |
| 4d | Hund | 2 | 0.30 | 0.40 | 0.50 | 0.90 | 1.00 | 1.00 | 1.00 |
| 4d | Hund | 3 | 0.20 | 0.20 | 0.27 | 0.80 | 1.00 | 1.00 | 1.00 |
| 4d | Hund | 4 | 0.15 | 0.20 | 0.20 | 0.80 | 1.00 | 1.00 | 1.00 |
| 4d | Hund | 6 | 0.13 | 0.17 | 0.20 | 0.73 | 1.00 | 1.00 | 1.00 |
| 4d | direct Coulomb only | 2 | 0.10 | 0.10 | 0.20 | 0.80 | 0.90 | 1.00 | 1.00 |
| 4d | direct Coulomb only | 4 | 0.05 | 0.05 | 0.10 | 0.80 | 0.80 | 1.00 | 1.00 |
| 4d | direct Coulomb only | 6 | 0.04 | 0.03 | 0.03 | 0.73 | 0.80 | 1.00 | 1.00 |

"Direct Coulomb only" is the self-interaction-free model of Result 6. The configuration-average rows are similar (see `results/g2/W_summary.txt`). g = 1 is regular at every radius in Hund and average modes, since its free atoms are already dⁿs⁰.

Median radius from which an element stays regular (per flavour, d block without the Cu/Zn types):

| row | mode | g = 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|
| 3d | Hund | 6 | 6 | 6 | 6 |
| 3d | average | 7 | 6 | 6 | 6.5 |
| 3d | direct Coulomb only | 6 | 6 | 5 | – |
| 4d | Hund | 8 | 8 | 8 | 8 |
| 4d | average | 10 | 8 | 8 | 8 |
| 4d | direct Coulomb only | 8 | 8 | 8 | 8 |

For comparison, the mean ⟨r⟩ of the outer s orbital in the free atoms (s-full configuration, per flavour, Hund) is 6.9 / 7.3 / 7.5 / 7.7 (3d) and 7.6 / 8.0 / 8.2 / 8.4 (4d) for g = 2 / 3 / 4 / 6. The wall takes over at about the size of the free ns orbital. That size, like the block onsets of Result 1, hardly depends on g in per-flavour units, and so neither does the radius.

Ground s counts at g = 2, Sc–Zn and Y–Cd:

| R_w (real radius) | 3d Hund | 3d average | 3d direct Coulomb only | 4d Hund | 4d average |
|---|---|---|---|---|---|
| free | 2211222212 | 2222222112 | 2222222222 | 2111211012 | 2221110012 |
| 10 (5 bohr) | 1111222112 | 2221111112 | 2222222222 | 1101110012 | 1100000012 |
| 8 (4 bohr) | 0001211112 | 0111111112 | 2222222222 | 0001000012 | 0000000012 |
| 7 (3.5 bohr) | 0001111112 | 0000000012 | 0112222222 | 0000000012 | 0000000012 |
| 6 (3 bohr) | 0000000012 | 0000000012 | 0000001112 | 0000000012 | 0000000012 |

- In Hund mode the half-filled configurations resist longest: Cr keeps d⁵s¹ and Mn d⁵s² (then d⁶s¹) after their neighbours have gone to s⁰. At R_w = 8 (4 bohr), E(Mn d⁷) − E(Mn d⁵s²) has fallen from +9.6 eV (free) to +3.8 eV, while V is already at −4.5 eV. The wall shifts each s electron by roughly 3 eV here, the same order as the ~2 eV upshift of Mn⁺ 4s in a He cage at R ≈ 2 Å (Wang et al. 2006).
- In the average configurations the 3d row passes through d^(n+1)s¹ for Ti–Cu (R_w = 8), which is the metal-like pattern, before going to s⁰. At g ≥ 3 the mixed stage is not a single count: the s shell empties one electron at a time (y = g−1, g−2, …) as the wall closes in.

Exchange split inside the wall (g = 2, configuration average, as in Result 8). Number of the 10 d-block atoms with a regular ground configuration:

| exchange on (besides core–core) | 3d: free / 8 / 7 / 6 | 4d: free / 8 / 7 / 6 |
|---|---|---|
| none | 1 / 1 / 1 / 6 | 1 / 5 / 8 / 8 |
| CV | 1 / 2 / 7 / 10 | 2 / 10 / 10 / 10 |
| d–s | 1 / 1 / 1 / 3 | 1 / 1 / 8 / 8 |
| within-shell | 1 / 1 / 2 / 10 | 2 / 8 / 8 / 8 |
| CV + d–s | 1 / 1 / 4 / 9 | 2 / 10 / 10 / 10 |
| CV + within-shell | 2 / 10 / 10 / 10 | 5 / 10 / 10 / 10 |
| d–s + within-shell | 1 / 1 / 1 / 6 | 1 / 4 / 8 / 8 |
| all (full HF) | 2 / 3 / 10 / 10 | 4 / 10 / 10 / 10 |

- The same pair that produces the free-atom 3d exceptions (core–valence + d–d) makes the 3d row regular earliest. d–s exchange works against it, in free atoms and in the wall.
- In the 4d row, every combination with core–valence exchange is fully regular by R_w = 8. Without it, Pd stays d⁹s¹ and Ag d⁹s² down to R_w = 6, the smallest radius in this split (8/10). In the full direct-Coulomb model, without core–core exchange either, both become regular by R_w = 6.
- 314 of 17,892 wall runs (mostly free g = 6 atoms) and 1 of 2,048 split runs did not converge and are excluded.

## Limitations
- Non-relativistic throughout. The real 5d row has few exceptions (Pt, Au) because relativity contracts 6s; the model drains 6s fully. That gap is a measure of the relativistic effect, which is known (Pyykkö) but not quantified here.
- Correlation only as an LDA add-on (Result 11), not a correlated wave function. HF misses V, Zr and Ce at g = 2. Correlation is part of the Coulomb interaction; its g-dependence is not studied.
- Free atoms only (see *Scope*), apart from the model wall of Result 13. No statement is made about atoms in molecules or solids.
- Spherical averaging within shells (configuration average or flavour-resolved average); no individual multiplets beyond the Hund state.
- At λ = 0 the f-filling windows (z ≈ 30–34 and 46–50) never converge. Two fillings alternate, most likely the known dual collapsed/extended f solutions.
- The fractional HF scans are qualitative (see Result 2).
- The wall of Result 13 is a crude cage: spherical, much stiffer than the Pauli repulsion of real neighbours, and without bonding or charge transfer. It tests only whether confinement alone produces the regular configurations, not what a given compound does. Wall runs use one grid without extrapolation; the configuration energy differences discussed are 0.05–40 eV.

## Related work
- Thomas–Fermi onsets of l shells: Fermi (1928). Madelung failure at very large Z: Bjerg, Fournais, Hearnshaw, Solovej, [arXiv:2506.14227](https://arxiv.org/abs/2506.14227).
- Transition-metal configurations and the Madelung rule: W. H. E. Schwarz, *J. Chem. Educ.* **87**, 444 (2010).
- Free atoms versus bound atoms; d and s populations in metals and complexes: K. Chen, W.-L. Li, W. H. E. Schwarz, *Z. Naturforsch. B* (2021), [doi:10.1515/znb-2021-0139](https://doi.org/10.1515/znb-2021-0139).
- The so-called electron configurations of the transition metals, including the lowest configuration averages (Table 1) and configurational energies (Fig. 8) used here: S.-G. Wang, Y.-X. Qiu, H. Fang, W. H. E. Schwarz, *Chem. Eur. J.* **12**, 4101 (2006).
- d-orbital collapse after the ns shell fills (the mechanism behind Result 4): J.-P. Connerade, *J. Phys. B* **24**, L109 (1991).
- Fractional-occupation spherical HF/HFS for all atoms H–Og: S. Lehtola, *Phys. Rev. A* **101**, 012516 (2020), [arXiv:1908.02528](https://arxiv.org/abs/1908.02528).
- Exchange shaping orbital collapse in specific cases: term-dependent HF for the 4f collapse in Xe-like ions ([OSTI 5267753](https://www.osti.gov/biblio/5267753)); 5d collapse in the Yb⁺ isoelectronic sequence attributed mainly to exchange ([link](https://www.bohrium.com/paper-details/collapse-of-d-and-f-orbitals-in-the-isoelectronic-sequence-of-singly-ionized-ytterbium/811035570686918658-532)).
- Orbital collapse and dual solutions: Goeppert Mayer (1941); Tupitsyn et al., *Phys. Rev. A* **109**, 042807 (2024), [arXiv:2402.02609](https://arxiv.org/abs/2402.02609).
- Periodic table in other dimensions: [atom-ph/9512001](https://arxiv.org/abs/atom-ph/9512001).
- SU(N) fermions (N ≤ 10) exist in cold-atom experiments with alkaline-earth atoms, but not as electrons bound to a nucleus.

## Reproducing
Requirements: Python 3.11, NumPy, SciPy; Modal for the HF scans. Run everything from `code/`.

| step | command | time |
|---|---|---|
| LDA check | `python check_nist.py` | seconds |
| LDA scan | `python g1_scan.py` | ~13 min on 12 cores → `results/g1_parts/`, `results/g1_scan.npz` |
| LDA analysis | `python g1_analyse.py results/g1_scan_run2.npz` and `python g1_analyse2.py results/g1_scan_run2.npz` | seconds |
| LDA hysteresis | `python g1_hyst.py` | ~5 min |
| HF checks | `python hf_check.py` (reference energies), `python check_ce.py` (Ce with the fixed-occupation start) | ~1–2 min each |
| HF scans | `python -m modal run modal_jobs/hf_modal.py --parts A` (also `--parts A --guess fixed`, `--parts B`, `--parts C`, `--parts D`, `--parts E`, `--parts F`, `--parts G`, `--parts H`, `--parts WX`) | 3–55 min wall each on Modal (WX: ~24 min) |
| HF analysis | `python g2_analyse_A.py results/g2/A.json results/g2/A_fixed.json`, `python g2_analyse_B.py`, `python g2_analyse_C.py`, `python g2_analyse_DE.py`, `python g2_analyse_F.py`, `python g2_analyse_GH.py`, `python g2_lock.py`, `python g2_levels_vs_avg.py`, `python g2_analyse_W.py` | seconds |
| CV decomposition | `python g2_cv_decomp.py` | ~30 s on 12 cores |

The run data used for every table are included in `code/results/`. `notes/working_notes.md` is the running write-up made during the work.
