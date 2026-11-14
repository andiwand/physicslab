// Headless sanity checks for the physics in index.html.
// Run: node check.mjs   (takes ~1–2 min)
import { readFileSync } from 'node:fs';

const html = readFileSync(new URL('./index.html', import.meta.url), 'utf8');
const code = html.match(/<script id="physics">([\s\S]*?)<\/script>/)[1];
const { Sim, Basis } = new Function(`${code}; return { Sim, Basis };`)();

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

console.log(failures ? `\n${failures} check(s) failed` : '\nall checks passed');
process.exit(failures ? 1 : 0);
