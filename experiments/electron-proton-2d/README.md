# electron-proton-2d

Electron wavepacket near a fixed proton, 2D time-dependent Schrödinger equation. Question: can the electron get captured into a bound state?

## Run

Open `index.html` in a browser (no build, no dependencies). Drag on the field to place the electron and give it a momentum. Presets cover both stages.

`node check.mjs` runs the physics headless and checks the results below (~1.5 min).

## Model

- Atomic units (ħ = mₑ = e = 1). Proton fixed at the origin, soft-core potential V(r) = −1/√(r² + a²), a = 0.5 a₀.
- 256×256 grid, 96×96 a₀ box (Δx = 0.375 a₀), Δt = 0.05 a.u. Split-operator FFT: e^{−iVΔt/2} e^{−iTΔt} e^{−iVΔt/2}.
- Complex absorbing potential ∝ s⁴ in a 22 a₀ edge band, so outgoing waves leave instead of wrapping around.
- Bound states from the radial equation (one symmetric tridiagonal eigenproblem per |m|, Sturm bisection + inverse iteration), interpolated onto the grid. Only states with ≥ 99 % of their probability inside r < 26 a₀ are tracked: n ≤ 3 (6 levels). P(bound) is the sum of |⟨φₙₘ|ψ⟩|² over those.

### Absorber

The first version (64 a₀ box, 6 a₀ quadratic absorber) reflected 2.5 % at k = 1 and 15–30 % for slow waves, which showed up as a fake "remnant" near the proton. Reflection is now measured against the exact free Gaussian (excess probability left inside):

| k (a.u.) | 0.15 | 0.3 | 0.6 | 1 | 2 | 3 |
|---|---|---|---|---|---|---|
| old | 0.34 | 0.30 | 0.15 | 2.5e-2 | 5e-6 | 4e-6 |
| new | 0.09 | 0.035 | 1.0e-3 | 1.3e-6 | 3.5e-6 | 1.5e-4 |

Very slow waves (wavelength longer than the absorber) are still partly reflected; fixing that needs a much bigger box or a different boundary method.

### Energy sink: radiation reaction (stage 2)

An accelerating charge radiates (Larmor, P = 2a²/3c³) and feels the recoil force F = τ da/dt with τ = 2/3c³ (Abraham–Lorentz in Landau–Lifshitz form). Applied semiclassically: a = ⟨−∇V⟩ (mean Coulomb acceleration, Ehrenfest), F acts as a uniform force (extra potential −F·r) on the whole wavefunction.

The real τ = 2.6·10⁻⁷ a.u. is far too small to watch (hydrogen 2p lifetime ≈ 1.6 ns ≈ 7·10⁷ a.u.), so it is multiplied by a boost:

| boost | τ (a.u.) | regime |
|---|---|---|
| 10⁶ | 0.26 | τ ≪ orbital times; work done by F ≈ −(Larmor energy) |
| 10⁸ | 26 | LL assumption violated; F still removes energy but only qualitatively |

### Classical ensemble

A second view runs N classical point electrons in the same potential (velocity Verlet, 4 substeps per quantum step), sampled from the current ψ on demand:

- position ~ |ψ|²
- momentum ~ Gaussian with the conditional mean and covariance of the Wigner function at that position: mean ∇S (phase gradient), covariance −¼∇∇ ln|ψ|² (clamped ≥ 0, since near nodes it can go negative / diverge).

Exact for Gaussian states; in general it reproduces ⟨p⟩, ⟨p²⟩ and therefore ⟨H⟩, but not interference. Fractions (remaining, bound = E < 0) are scaled by the quantum norm at sampling time so they plot on the same axes as the quantum curves.

With radiation on, each classical electron loses energy at its own Larmor rate τa², taken from its kinetic energy. That is not the Landau–Lifshitz force used for the quantum mean field (LL is unstable for point particles at large τ: near the core τ·|∇∇V| ≫ 1/Δt), only the same loss rate. A classical bound electron keeps radiating until it sits at the bottom of the well (E = −1/a = −2 Ha) — the classical collapse of the atom; the quantum ground state does not radiate.

Limits of the mean-field picture: an excited eigenstate alone has ⟨a⟩ = 0 and never decays (real atoms do, via vacuum fluctuations); only superpositions radiate. And the wavefunction stays one coherent object instead of branching into "captured + photon" and "free electron".

## Checks (`node check.mjs`)

- Radial eigenvalues vs. energy of the same state on the 2D grid: max |ΔE| ≈ 1·10⁻⁴ Ha.
- Absorber reflection at k = 0.6 / 1 / 2 below 3·10⁻³ / 10⁻⁴ / 10⁻⁴.
- No sink, fly-by (E ≈ +0.49 Ha): P(bound) ≈ 2·10⁻⁶ throughout.
- No sink, packet released next to the proton: P(bound) = 0.8695 ± 0.0001 over t = 0…80 a.u.
- Ehrenfest on the grid: Δ⟨p⟩ = ∫⟨−∇V⟩dt to 4 digits (checked once, not in `check.mjs`).
- Sink, boost 10⁶, ground + 2p superposition: relaxes to P(1s) = 0.992; work −0.2216 Ha vs radiated 0.2221 Ha.
- Sink, boost 10⁷, ground state alone: radiates nothing (⟨a⟩ = 0 by symmetry).
- Sink, boost 10⁸, slow fly-by: P(bound) 0.03 → 0.59. Same result at Δt = 0.025 and 0.0125 (converged in time step).
- Ensemble sampled from the Gaussian packet: ⟨p⟩ = (1.002, 0.000), var p = 0.040 (exact 1/4σ² = 0.04).
- Ensemble sampled from the ground state: classical ⟨E⟩ = −0.6545 vs quantum −0.6555 Ha.
- Classical energy conservation without sink: max |ΔE| = 4·10⁻⁴ Ha over t = 40 (20 000 electrons).
- Classical ground-state ensemble with radiation (boost 10⁶): ⟨E⟩ drops to −1.6 Ha within t = 40, far below the quantum ground state.

## Findings

1. **No sink → no capture.** With a time-independent Hamiltonian each eigenstate only acquires a phase, so bound-state populations are constants of motion. What is bound at t = 0 stays bound; the continuum part scatters and leaves.
2. **The remnant after a fly-by is free, slow electron.** Positive energy, overlap with bound states ~10⁻⁶. It is the slow tail of the packet's momentum distribution and drains away roughly as a power law in t, not exponentially (×~10 per doubling of t).
3. **Radiation reaction gives capture — but only when it is unphysically strong.** At boost 10⁶ (where the model is valid) a single pass captures nothing measurable; at 10⁸ a slow electron ends ~60 % bound, a head-on 13 eV electron ~40 %. Matches the real picture: radiative recombination is a rare process per encounter. Slow electrons are captured much more easily.
4. **Display note:** the brightness used to auto-scale to the current peak, which made a 10⁻⁴ remnant look as bright as the original packet. Default is now a fixed scale (initial peak); auto is a toggle.
5. **Captured population lands in excited levels first and cascades down**; a ground + 2p superposition relaxes to the ground state, with the emitted energy matching the 2p–1s gap.

## Next

1. ~~Stage 1: plain TDSE, show that P(bound) is conserved.~~
2. ~~Stage 2: energy sink (semiclassical radiation reaction).~~
3. Stage 3: two electrons in 1D + proton (2D configuration space), three-body recombination — energy carried away by the second electron, no fudge factor.
