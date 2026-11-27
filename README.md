# g-flavours-periodic-table
The periodic table with g electron spin flavours

Author: Thomas Doherty, with computational and proof assistance from Claude (Anthropic)

## Summary

Electrons have two spin states, and every shell capacity in chemistry, 2(2l+1), quietly depends on that. This repository asks what the periodic table looks like when electrons have **g** internal states (g = 1, 2, 3, …, ∞). The aim is to separate what in the table comes from plain electrostatics from what comes from exchange.

Rescaling turns the g-flavour atom into a one-parameter family. In per-flavour units the problem is the g = 1 atom with nuclear charge z = Z/g, and exchange enters multiplied by **λ = 1/g**:

- g → ∞ (λ = 0) is pure electrostatics, the exact mean-field (Hartree) limit;
- g = 2 (λ = ½) is real chemistry;
- g = 1 (λ = 1) is spinless electrons.

Main results, from exchange-only LDA scans and then Hartree–Fock (HF) with exact exchange:

1. **Each block starts at the same z = Z/g for every g** (within ±0.85 per-flavour electrons), as Thomas–Fermi scaling predicts. Stronger exchange brings d and f shells in slightly earlier.
2. **The Cr/Cu-type exceptions in the 3d row are exchange effects.**
   - In HF with maximal flavour polarisation (Hund), the share of the 3d block with a non-Madelung ground configuration is 0.80 / 0.30 / 0 / 0 / 0 / 0 / 0 for g = 1 / 2 / 3 / 4 / 5 / 6 / 8.
   - At g = 2 the same HF gets 11 of 12 real ground-state s occupancies in each d row, including Cr, Cu, Nb, Mo, Ru, Rh, Ag and Pd 4d¹⁰5s⁰.
3. **About half of the 4d-row anomaly survives without exchange.**
   - The 4d exception share falls 0.80 → 0.70 → 0.60 → 0.45 → 0.44 → 0.37 → 0.38 as g grows. It levels off where the pure-electrostatics limit puts the 5s drain (0.35).
   - In the λ = 0 limit the 3d row is exactly Madelung, while the 4d row is not.
4. **Without exchange, 3d is not bound when 4s fills.** At λ = 0 the 3d level is an unbound continuum state until z ≈ 10.5, so 4p starts filling first. The head start of p before d grows smoothly with g in both LDA and HF.
5. **The La/Ac step (one d electron before the f shell) depends on g.**
   - HF: absent at g = 1 (4f lower by 3.7 eV), present at g = 2 by only 0.08 eV, and strong for g ≥ 3 (3–9 eV).
   - At g = 2, HF matches 13 of 14 lanthanide ground configurations.
   - Real chemistry sits just on the d-first side of a boundary.

No earlier treatment of a g-flavour periodic table was found (searches on spinless, SU(N) and spin-degeneracy variants of the periodic table and on large-N atoms; see *Related work*).

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
- **Runs:** ~18,000 HF calculations on Modal (`modal_jobs/hf_modal.py`).

## Validation

| check | result |
|---|---|
| LDA vs NIST LDA reference, Ne and Ar orbital energies | agree to ≤ 3.6 mHa |
| LDA Fe, Cu (fractional aufbau) | 3d⁶·⁹4s¹·¹, 3d¹⁰4s¹ (standard spherical-LDA result) |
| HF gradients vs finite differences | ratio 1 ± 10⁻⁶ |
| HF total energies (extrapolated) vs numerical-HF reference: Ne, Ar, Zn | +1·10⁻⁵, +3·10⁻⁵, +1·10⁻⁴ Ha |
| HF Cr 3d⁵4s¹ ⁷S, Fe 3d⁶4s² ⁵D (Hund mode) | +3.4·10⁻⁵, +4.4·10⁻⁵ Ha |
| Part A with two different starting guesses | identical ground states |

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

At large g, Hund-mode exceptions in the 4d row tend to sit at d counts that are multiples of 5, i.e. whole flavours filled (a generalised Cr/Pd stability): 56% vs 37% in the configuration average and 20% by chance. This is suggestive, not clean.

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

f counts lock at multiples of 7 (whole flavours, a generalised Gd). At g ≥ 4 the ground state often sits at the y = 2 cap, so states with more 5d electrons may be lower.

## Limitations
- Non-relativistic throughout. The real 5d row has few exceptions (Pt, Au) because relativity contracts 6s; the model drains 6s fully. That gap is a measure of the relativistic effect, which is known (Pyykkö) but not quantified here.
- No correlation, apart from the LDA reference row. HF misses V, Zr and Ce at g = 2.
- Spherical averaging within shells (configuration average or flavour-resolved average); no individual multiplets beyond the Hund state.
- At λ = 0 the f-filling windows (z ≈ 30–34 and 46–50) never converge. Two fillings alternate, most likely the known dual collapsed/extended f solutions.
- The fractional HF scans are qualitative (see Result 2).

## Related work
- Thomas–Fermi onsets of l shells: Fermi (1928). Madelung failure at very large Z: Bjerg, Fournais, Hearnshaw, Solovej, [arXiv:2506.14227](https://arxiv.org/abs/2506.14227).
- Transition-metal configurations and the Madelung rule: W. H. E. Schwarz, *J. Chem. Educ.* **87**, 444 (2010).
- Fractional-occupation spherical HF/HFS for all atoms H–Og: S. Lehtola, *Phys. Rev. A* **101**, 012516 (2020), [arXiv:1908.02528](https://arxiv.org/abs/1908.02528).
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
| HF scans | `python -m modal run modal_jobs/hf_modal.py --parts A` (also `--parts A --guess fixed`, `--parts B`, `--parts C`) | 35–55 min wall each on Modal |
| HF analysis | `python g2_analyse_A.py results/g2/A.json results/g2/A_fixed.json`, `python g2_analyse_B.py`, `python g2_analyse_C.py` | seconds |

The run data used for every table are included in `code/results/`. `notes/working_notes.md` is the running write-up made during the work.
