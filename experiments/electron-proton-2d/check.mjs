// Headless sanity checks for the physics in index.html.
// Run: node check.mjs   (takes ~1–2 min)
import { readFileSync } from 'node:fs';

const html = readFileSync(new URL('./index.html', import.meta.url), 'utf8');
const code = html.match(/<script id="physics">([\s\S]*?)<\/script>/)[1];
const { Sim, Basis, Ensemble } = new Function(`${code}; return { Sim, Basis, Ensemble };`)();

let failures = 0;
function expect(label, ok, detail) {
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${label}  ${detail}`);
  if (!ok) failures++;
}
const fmt = (v) => v.toExponential(1);

const sim = new Sim();
const basis = new Basis(sim);
const lv = (n, m) => basis.levels.find((l) => l.n === n && l.m === m);
console.log(`grid ${sim.N}² over ${sim.L} a₀, ${basis.levels.length} bound levels tracked\n`);

// 1. Radial eigenvalues vs. energy of the same state on the 2D grid.
let maxDiff = 0;
for (const l of basis.levels) {
  sim.setSuperposition([{ level: l, m: l.m, c: [1, 0] }]);
  const m = sim.measure();
  maxDiff = Math.max(maxDiff, Math.abs(m.E / m.norm - l.E));
}
expect('bound-state energies, radial vs 2D grid', maxDiff < 3e-4, `max |ΔE| = ${fmt(maxDiff)} Ha`);

// 2. Absorber: free packet into the edge; excess interior probability vs exact free solution.
function reflection(k, sigma = 3) {
  const s = new Sim({ soft: 1e12 });
  const inner = s.L / 2 - s.absorb, x0 = inner - 10;
  s.setPacket({ x0, y0: 0, kx: k, ky: 0, sigma });
  s.step(Math.ceil(((10 + 2 * s.absorb + 4 * sigma) / k + 10) / s.dt));
  const t = s.t, st = sigma * Math.sqrt(1 + (t / (2 * sigma * sigma)) ** 2), dA = s.dx ** 2;
  let num = 0, ex = 0;
  for (let j = 0; j < s.N; j++) for (let i = 0; i < s.N; i++) {
    const x = s.xs[i], y = s.xs[j];
    if (Math.abs(x) > inner || Math.abs(y) > inner) continue;
    const idx = j * s.N + i;
    num += (s.re[idx] ** 2 + s.im[idx] ** 2) * dA;
    ex += Math.exp(-((x - x0 - k * t) ** 2 + y * y) / (2 * st * st)) / (2 * Math.PI * st * st) * dA;
  }
  return num - ex;
}
for (const [k, tol] of [[0.6, 3e-3], [1, 1e-4], [2, 1e-4]]) {
  const r = reflection(k);
  expect(`absorber reflection k=${k}`, r < tol, `${fmt(r)} (tolerance ${fmt(tol)})`);
}

// 3. No radiation: P(bound) is a constant of motion.
function pRange(init, steps) {
  init();
  let lo = Infinity, hi = -Infinity;
  for (let s = 0; s <= steps; s += 200) {
    if (s) sim.step(200);
    const p = basis.project(sim.re, sim.im).total;
    lo = Math.min(lo, p); hi = Math.max(hi, p);
  }
  return [lo, hi];
}
sim.boost = 0;
let [lo, hi] = pRange(() => sim.setPacket({ x0: -18, y0: 4, kx: 1, ky: 0, sigma: 2.5 }), 1600);
expect('fly-by, no sink: P(bound) stays 0', hi < 1e-4, `range ${fmt(lo)} .. ${fmt(hi)}`);
[lo, hi] = pRange(() => sim.setPacket({ x0: -3, y0: 0, kx: 0, ky: 0, sigma: 1.5 }), 1600);
expect('released near proton, no sink: P(bound) constant', hi - lo < 2e-3, `range ${lo.toFixed(5)} .. ${hi.toFixed(5)}`);

// 4. Radiation reaction.
sim.boost = 1e6;
sim.setSuperposition([{ level: lv(1, 0), m: 0, c: [1, 0] }, { level: lv(2, 1), m: 1, c: [1, 0] }]);
sim.step(8000);
let pops = basis.project(sim.re, sim.im).pops;
expect('ground + 2p, boost 1e6: relaxes to ground state', pops[0] > 0.98, `P(1s) = ${pops[0].toFixed(4)} at t = ${sim.t.toFixed(0)}`);
expect('  energy balance (work ≈ −radiated)', Math.abs(sim.rrWork + sim.Erad) < 0.05 * sim.Erad,
  `work ${sim.rrWork.toFixed(4)}, radiated ${sim.Erad.toFixed(4)} Ha`);

sim.boost = 1e7;
sim.setSuperposition([{ level: lv(1, 0), m: 0, c: [1, 0] }]);
sim.step(400);
expect('ground state, boost 1e7: does not radiate', sim.Erad < 1e-8, `radiated ${fmt(sim.Erad)} Ha`);

sim.boost = 1e8;
sim.setPacket({ x0: -16, y0: 3, kx: 0.4, ky: 0, sigma: 3 });
const p0 = basis.project(sim.re, sim.im).total;
sim.step(2400);
const p1 = basis.project(sim.re, sim.im).total;
expect('slow fly-by, boost 1e8: capture', p1 > p0 + 0.3, `P(bound) ${p0.toFixed(3)} → ${p1.toFixed(3)}`);
expect('  energy removed', sim.rrWork < -0.1, `work ${sim.rrWork.toFixed(4)} Ha (radiated ${sim.Erad.toFixed(4)}; τ = 26 a.u., LL balance not expected)`);

// 5. Classical ensemble sampled from ψ.
let seed = 12345;
const rng = () => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
sim.boost = 0;
const ens = new Ensemble(sim);
sim.setPacket({ x0: -18, y0: 4, kx: 1, ky: 0, sigma: 2.5 });
ens.sample(20000, rng);
const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const vr = (a) => { const m = mean(a); return mean(a.map((v) => (v - m) ** 2)); };
const [mvx, mvy, vvx, vvy, vx_] = [mean(ens.vx), mean(ens.vy), vr(ens.vx), vr(ens.vy), vr(ens.x)];
expect('sampling a Gaussian packet: momentum mean and variance',
  Math.abs(mvx - 1) < 0.01 && Math.abs(mvy) < 0.01 && Math.abs(vvx - 0.04) < 0.003 && Math.abs(vvy - 0.04) < 0.003,
  `⟨p⟩ = (${mvx.toFixed(3)}, ${mvy.toFixed(3)}), var = (${vvx.toFixed(4)}, ${vvy.toFixed(4)}) vs 0.04; var x = ${vx_.toFixed(2)} vs 6.25`);

sim.setSuperposition([{ level: lv(1, 0), m: 0, c: [1, 0] }]);
ens.sample(20000, rng);
const Ecl = mean(Array.from({ length: ens.n }, (_, k) => ens.energy(k)));
expect('sampling the ground state: classical ⟨E⟩ ≈ quantum ⟨H⟩', Math.abs(Ecl - lv(1, 0).E) < 0.02,
  `${Ecl.toFixed(4)} vs ${lv(1, 0).E.toFixed(4)} Ha`);

const E0s = Array.from({ length: ens.n }, (_, k) => ens.energy(k));
ens.step(800);
let maxDE = 0;
for (let k = 0; k < ens.n; k++) if (ens.alive[k]) maxDE = Math.max(maxDE, Math.abs(ens.energy(k) - E0s[k]));
expect('classical energy conservation (no sink, t = 40)', maxDE < 1e-3, `max |ΔE| = ${fmt(maxDE)} Ha`);

ens.step(800, 1e6 * 2 / (3 * 137.035999 ** 3));
const st = ens.stats();
expect('classical ground-state ensemble with radiation (boost 1e6) collapses below the quantum ground state',
  st.Emean < lv(1, 0).E - 0.2, `⟨E⟩ = ${st.Emean.toFixed(3)} Ha (quantum ground ${lv(1, 0).E.toFixed(3)}, well bottom −2)`);

console.log(failures ? `\n${failures} check(s) failed` : '\nall checks passed');
process.exit(failures ? 1 : 0);
