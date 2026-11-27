# The periodic table with g electron spin flavours

**Question.** Chemistry quietly assumes electrons have g = 2 spin states. What does the periodic table look like for g = 1, 3, 4, …, ∞? Which parts of the table are electrostatics and which are exchange?

**Novelty.** Nothing found (searches on spinless / SU(N) / spin-degeneracy periodic tables, large-N atoms). Neighbours: the periodic table in other dimensions (atom-ph/9512001), Madelung from Thomas–Fermi (Fermi 1928, Z_l ≈ 0.17(2l+1)³), and Madelung failing at very large Z (arXiv:2506.14227). SU(N) fermions exist in cold atoms (Sr, Yb, N ≤ 10), but not as electrons around a nucleus.

## The one-dial trick
Use per-flavour units: z = Z/g, s = g·r, energies / g². The g-flavour atom with Slater (LDA) exchange becomes the g = 1 problem with nuclear charge z, plus exchange multiplied by **λ = 1/g**:

- λ = 0 (g → ∞): pure Hartree electrostatics, the exact mean-field limit;
- λ = ½: real chemistry;
- λ = 1: spinless electrons.

Thomas–Fermi theory has no shells. In TF the onset of each new l falls at Z_l ∝ g(2l+1)³, i.e. at the same z for every g. So any g-dependence in the filling order comes from shell structure and exchange.

## Method
- Spherical, non-relativistic atom with exchange-only LDA (Slater, α = ⅔); PZ81 correlation added only for the real g = 2 reference row.
- Occupations minimised directly (Janak: ∂E/∂f = ε), so partly filled shells share the Fermi level ("fractional aufbau").
- Log grid with 1200 points.
- Check against the NIST LDA reference: Ne and Ar eigenvalues agree to ≤ 3.6 mHa; Fe comes out 3d⁶·⁹4s¹·¹ and Cu 3d¹⁰4s¹ (standard spherical LDA).
- Scan: z = 0.25–60 for g = ∞, 4, 3, 2, 1, plus g = 2 with correlation. Two independent runs with different solver settings agree (except the g = ∞ points noted below).
- Code: `atom.py`, `g1_scan.py`, `g1_analyse.py`, `g1_analyse2.py`, `g1_hyst.py`. Outputs: `results/`.

## Results

**1. Each new block starts at nearly the same z for every g.** Persistent onset of the first d and f shells (per-flavour z; Madelung values 10, 19, 28, 44):

| shell | g = ∞ | 4 | 3 | 2 | 2 + corr | 1 |
|---|---|---|---|---|---|---|
| 3d | 10.58 | 10.33 | 10.25 | 10.17 | 10.5 | 9.67 |
| 4d | 19.67 | 19.42 | 19.33 | 19.17 | 19.5 | 18.92 |
| 4f | 29.0 | 28.75 | 28.58 | 28.33 | 28.5 | 27.58 |
| 5f | 45.2 | 44.83 | 44.67 | 44.4 | 44.5 | 43.67 |

- The spread is ≤ 1.7 per-flavour electrons, as Thomas–Fermi scaling predicts.
- The shift is monotone: stronger exchange brings compact d/f orbitals in earlier.

**2. Draining of the outer s electrons into the d/f shell is set by exchange.** Largest drop of the outer s shell inside each block, per flavour (1 = the s shell empties):

| block | g = ∞ | 4 | 3 | 2 | 1 |
|---|---|---|---|---|---|
| 3d (4s) | 0.00 | 0.17 | 0.31 | 0.58 | 1.00 |
| 4d (5s) | 0.35 | 0.77 | 0.92 | 1.00 | 1.00 |
| 4f+5d (6s) | 0.43 | 0.83–0.87 | 1.00 | 1.00 | 1.00 |

- In the pure-electrostatics limit the 3d row is plain Madelung: 4s stays full and 3d fills on top.
- The Cr/Cu-type s→d transfer grows steadily with λ.
- At fixed λ it is stronger in deeper rows. That matches the real table, where the 4d row (Nb, Mo, Ru, Rh, Pd, Ag) has more exceptions than the 3d row (Cr, Cu). The real 5d row has fewer exceptions because of relativity, which this model leaves out.

**3. The La/Ac step (one d electron before the f shell) is the rule for g ≥ 3.**
- 5d starts at z = 28.17 for g = 3 and 4 (4f only at 28.6–28.8); 6d starts at 44.17, before 5f (44.7–44.8).
- At g = 2 (exchange-only LDA) the step is missing in both rows, although real La, Ac and Th show it.
- At g = 1 it appears in the actinides (6d 43.33 < 5f 43.67) but not in the lanthanides.

**4. At g = ∞ the next p shell starts before the d shell, because the d orbital isn't bound yet.**
- At λ = 0 the 3d level sits at +0.001 (an unbound box state; the box size was checked, rmax 120 vs 400) until z ≈ 10.5. So 4p fills to 0.49 per flavour first (converged, err 5·10⁻⁸). Then 3d drops in.
- At λ = ¼, 4p peaks at 0.22 before 3d takes over smoothly. At λ = ½ there is no 4p step.
- Exchange is what binds 3d in time for Madelung order.
- The same happens with 5p before 4d.
- Up/down sweeps at λ = ¼ and ½ agree exactly: no hysteresis. At λ = 0 the sweeps differ only at points that never converged (z = 10.67–11.17), so the collapse can't be called first-order.

## Hartree–Fock check (G2, run on Modal)
Average-of-configuration HF generalised to g flavours (`hf.py`): exact exchange, no self-interaction error. Two modes:
- **configuration average**;
- **Hund**: maximal flavour polarisation; at g = 2 this is the high-spin term, e.g. Cr ⁷S, Fe ⁵D.

Checks: Be, Ne, Ar, Zn, Cr ⁷S and Fe ⁵D reproduce the published HF energies to 4·10⁻⁵ to 1·10⁻⁴ Ha (two grids, h² extrapolation). Integer configurations (n−1)dˣ nsʸ were run for every allowed y, for g = 1–8, with two different starting guesses. The ground states are identical between the two guesses.

**Exceptions (ground state not "s full first") as a share of the d block, Hund mode:**

| row | g = 1 | 2 | 3 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|---|
| 3d | 0.80 | 0.30 | 0 | 0 | 0 | 0 | 0 |
| 4d | 0.80 | 0.70 | 0.60 | 0.45 | 0.44 | 0.37 | 0.38 |

- **At g = 2, HF (Hund) gets 11/12 real s occupancies in each row:**
  - 3d row: Cr and Cu right; V wrong.
  - 4d row: Nb, Mo, Ru, Rh, Ag right, and Pd 4d¹⁰5s⁰; Zr wrong.
  - The configuration average gets only 10/12 and 9/12, so high-spin exchange is what produces the real exceptions.
- **The 3d exceptions are pure exchange:** they vanish for every g ≥ 3.
- **The 4d exceptions only partly are:** they fall to ~0.38 of the block and stay there. LDA gave a drain of 0.35 at g = ∞ (no exchange at all), so about half of the 4d anomaly survives in pure electrostatics.
- **g = 1 (spinless):** HF also drains the outer s completely across the d block (0.80 = every element but the last). LDA's 100% drain was **not** a self-interaction artefact.
- **Locking at whole flavours:** in Hund mode, 56% of large-g 4d exceptions sit at d counts that are multiples of 5 (whole flavours filled), vs 37% for the configuration average and 20% by chance. That is a generalised Cr/Pd stability recurring g times per row. Suggestive, not clean.

**Fractional HF scans** (configuration average, qualitative: ~45% of points flagged, mostly residuals of empty orbitals):
- d onset moves later with g, as in LDA.
- The p shell's head start before d grows with g (4d row: 0.04 → 0.42 for g = 2 → 16; the pure-electrostatics limit is 0.49), so result 4 holds with exact exchange.
- 4d s drain 1.00 / 0.89 / 0.68 / 0.58 / 0.49 / 0.45 / 0.38 for g = 1–16.

**La/Ac step (lanthanide row, Hund):** the first element takes 5d before 4f:
- g = 1: no (4f wins by 3.7 eV);
- g = 2: yes, but only by 0.08 eV. Gd 4f⁷5d¹ is right, and 13/14 lanthanide ground configurations match experiment (Ce missed);
- g ≥ 3: yes, by 3–9 eV, with 5d occupied across most of the f block. f counts lock at multiples of 7.

So result 3 holds in HF, and the g = 2 case that exchange-only LDA missed comes out right. Caveat: only up to two 5d electrons were allowed.

Code: `hf.py`, `modal_jobs/hf_modal.py`, `g2_analyse_A.py`, `g2_analyse_B.py`, `g2_analyse_C.py`. Data: `results/g2/`.

## Status
- **New (nothing found):** the g-flavour table and the λ = 1/g dial; results 2–4, now confirmed with Hartree–Fock (3d exceptions vanish for g ≥ 3, 4d exceptions level off at ~0.38, La step for g ≥ 2).
- **Expected:** result 1 (Thomas–Fermi scaling).
- **Model limits:** exchange-only LDA, spherical, non-relativistic, fractional occupations. The real-chemistry row reproduces Madelung onsets to ±0.5 in z but not the individual integer exceptions; LDA's fractional 4s→3d transfer is a known feature.

## Open
- λ = 0 never converges where the f shells fill (z ≈ 30–34, 46–50): two fillings alternate. Most likely the known dual collapsed/extended solutions of double-well f orbitals (arXiv:2402.02609); HF in this work hit the same branch problem for Ce until the start had 4f occupied.
- Relativity: the non-relativistic model drains the 6s fully in the 5d row, while real 5d has few exceptions (Pt, Au). The size of that gap is a measure of the relativistic 6s contraction. That is known physics (Pyykkö), but a scalar-relativistic rerun would put a number on it inside this framework.
- HF puts real chemistry (g = 2) just 0.08 eV on the d-first side of the La step, so g = 2 sits right next to the boundary.
