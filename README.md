# morphogenesis

> **Form from formlessness. Turing patterns for agent development.**

[![crates.io](https://img.shields.io/crates/v/morphogenesis.svg)](https://crates.io/crates/morphogenesis)
[![docs.rs](https://docs.rs/morphogenesis/badge.svg)](https://docs.rs/morphogenesis)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Rust library implementing reaction-diffusion systems, Turing instability analysis, symmetry breaking detection, and morphogen signaling for self-organizing systems. Models how homogeneous populations of agents spontaneously develop structure — roles, specializations, and boundaries — through local interactions alone.

---

## Table of Contents

- [What is Morphogenesis?](#what-is-morphogenesis)
- [Why Does This Matter?](#why-does-this-matter)
- [Architecture](#architecture)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
- [Mathematical Background](#mathematical-background)
- [Installation](#installation)
- [Related Crates](#related-crates)
- [License](#license)

---

## What is Morphogenesis?

In 1952, Alan Turing published "The Chemical Basis of Morphogenesis" — showing how two chemicals (morphogens) diffusing at different rates and reacting with each other can spontaneously break symmetry and create patterns. Spots, stripes, spirals — all from initially uniform conditions.

This library applies the same principle to **agent development**:

- How do agents in a uniform pool spontaneously specialize into roles?
- How do behavioral patterns emerge from initially identical agents?
- How do boundaries form between agent teams without central control?

The answer: **reaction-diffusion dynamics**. Agents exchange signals locally, react to their neighbors' states, and patterns emerge from the interplay of activation and inhibition spreading at different rates.

```
Time 0:  ○ ○ ○ ○ ○ ○ ○ ○ ○ ○    (uniform — all agents identical)
Time 10: ○ ● ○ ○ ● ○ ○ ● ○ ○    (perturbation begins)
Time 50: ● ● ○ ○ ○ ● ● ○ ○ ○    (pattern forming)
Time 200:● ● ● ○ ○ ● ● ● ○ ●    (stable Turing pattern — roles fixed)
```

## Why Does This Matter?

**For swarm robotics**: Turing patterns explain how identical robots can spontaneously divide into specialized teams without any central coordinator or pre-assigned roles.

**For multi-agent systems**: Morphogen signaling provides a biologically-inspired alternative to explicit role assignment. Agents self-organize based on local interactions alone.

**For development AI**: Understanding how structure emerges from homogeneity is fundamental to building agents that develop, adapt, and grow rather than being hard-coded.

**For systems biology**: The Gray-Scott model implemented here directly models real chemical morphogenesis — the same math that creates patterns on animal coats and fish scales.

## Architecture

```
morphogenesis
│
├── reaction_diffusion          ← Gray-Scott model
│   ├── ReactionDiffusion1D         Two-morphogen system (activator u, inhibitor v)
│   ├── step()                      Single time step (Euler method)
│   ├── run(n)                      Run n steps
│   ├── laplacian_1d()              Spatial second derivative
│   ├── amplitude()                 Pattern strength measure
│   └── is_stable()                 Convergence check
│
├── turing_analysis             ← Linear stability analysis
│   ├── is_turing_unstable()        Can patterns form?
│   ├── critical_wavenumber()       Most unstable spatial mode
│   ├── dispersion_relation()       Growth rate vs wavenumber
│   ├── fastest_growing_mode()      Dominant pattern wavelength
│   ├── in_turing_space()           Parameter space check
│   └── stability_margin()          Distance from instability boundary
│
├── symmetry                    ← Pattern formation detection
│   ├── is_symmetry_broken()        Has the uniform state broken?
│   ├── symmetry_breaking_measure() How far from uniform
│   ├── order_parameter()           Long-range order metric
│   ├── is_lr_symmetric()          Left-right symmetry check
│   └── breaking_speed()            Rate of pattern formation
│
└── morphogen                   ← Individual morphogen fields
    ├── Morphogen                   Named chemical field
    ├── diffuse()                   Spatial diffusion step
    ├── decay()                     Natural decay
    ├── produce()                   Source production
    ├── gradient()                  Concentration gradient at position
    ├── threshold()                 Binary thresholding
    └── interpolate()               Continuous field lookup
```

## Quick Start

```rust
use morphogenesis::{
    ReactionDiffusion1D,
    turing_analysis::is_turing_unstable,
    symmetry::is_symmetry_broken,
};

// Create a Gray-Scott reaction-diffusion system
// Parameters: du (activator diffusion), dv (inhibitor diffusion),
//             dt (time step), dx (spatial step), f (feed), k (kill)
let mut rd = ReactionDiffusion1D::new(100, 0.16, 0.08, 1.0, 1.0, 0.035, 0.065);

// Perturb the uniform state (necessary to seed pattern formation)
rd.init_perturbation();
// Or use random initialization: rd.init_random(42);

// Run the simulation for 1000 steps
rd.run(1000);

// Check the results
let pattern = rd.pattern();        // The activator field u
let amp = rd.amplitude();          // Pattern strength
let stable = rd.is_stable();       // Has it converged?

println!("Pattern amplitude: {:.4}", amp);
println!("Stable: {}", stable);
println!("Symmetry broken: {}", is_symmetry_broken(pattern, 0.01));

// Analyze: can these parameters produce Turing patterns?
let turing_ok = is_turing_unstable(0.16, 0.08, -1.0, 1.0, 2.0, -1.5);
println!("Turing unstable: {}", turing_ok);
```

## API Reference

### ReactionDiffusion1D

| Method | Returns | Description |
|--------|---------|-------------|
| `new(size, du, dv, dt, dx, f, k)` | `Self` | Create Gray-Scott system |
| `init_perturbation()` | `()` | Seed center with activator |
| `init_random(seed)` | `()` | Random seed placement |
| `step()` | `()` | Single Euler time step |
| `run(n)` | `()` | Run n steps |
| `pattern()` | `&[f64]` | Activator field u |
| `amplitude()` | `f64` | Pattern strength (std dev of u) |
| `mean_u()` / `mean_v()` | `f64` | Mean concentration |
| `is_stable()` | `bool` | Has system converged? |

### Turing Analysis

| Function | Returns | Description |
|----------|---------|-------------|
| `is_turing_unstable(fu, fv, gu, gv, du, dv)` | `bool` | Can patterns form? |
| `critical_wavenumber(fu, fv, gu, gv, du, dv)` | `f64` | Most unstable mode |
| `dispersion_relation(k, fu, fv, gu, gv, du, dv)` | `f64` | Growth rate at wavenumber k |
| `fastest_growing_mode(du, dv, fu, fv, gu, gv)` | `f64` | Wavenumber with max growth |
| `in_turing_space(fu, fv, gu, gv, du, dv)` | `bool` | Parameters in Turing space? |
| `stability_margin(fu, fv, gu, gv, du, dv)` | `f64` | Distance from instability |
| `homogeneous_stable(fu, fv, gu, gv)` | `bool` | Is uniform state stable? |

### Symmetry & Pattern Detection

| Function | Returns | Description |
|----------|---------|-------------|
| `is_symmetry_broken(pattern, threshold)` | `bool` | Has pattern emerged? |
| `symmetry_breaking_measure(pattern)` | `f64` | Distance from uniform |
| `order_parameter(pattern)` | `f64` | Long-range order |
| `is_lr_symmetric(pattern, tol)` | `bool` | Left-right symmetry |
| `breaking_speed(patterns)` | `Vec<f64>` | Rate of pattern change |
| `dominant_wavelength(pattern, dx)` | `f64` | Primary spatial frequency |

### Morphogen

| Method | Returns | Description |
|--------|---------|-------------|
| `new(name, size, diffusion, decay, production)` | `Self` | Create named morphogen field |
| `set_uniform(value)` | `()` | Set all cells to value |
| `add_source(position, amount)` | `()` | Add point source |
| `diffuse(dx, dt)` | `()` | Diffusion step |
| `decay(dt)` | `()` | Natural decay |
| `produce(dt)` | `()` | Auto-production |
| `gradient(pos, dx)` | `f64` | Spatial gradient |
| `threshold(t)` | `Vec<bool>` | Binary field |

## Mathematical Background

### Gray-Scott Model

The Gray-Scott reaction-diffusion equations:

```
∂u/∂t = Du ∇²u − uv² + f(1 − u)
∂v/∂t = Dv ∇²v + uv² − (f + k)v
```

Where:
- **u** = activator concentration (short-range activation)
- **v** = inhibitor concentration (long-range inhibition)
- **Du** < **Dv** = activator diffuses slower than inhibitor (critical for pattern formation)
- **f** = feed rate (replenishment of u)
- **k** = kill rate (removal of v)

The nonlinear term **uv²** models autocatalysis: the inhibitor v is produced when both u and v are present, creating a positive feedback loop.

### Turing Instability

A parameter set (Du, Dv, f, k) is **Turing unstable** when:
1. The homogeneous (uniform) state is stable to uniform perturbations
2. The homogeneous state is unstable to spatially varying perturbations

This requires **differential diffusion**: Du < Dv. The inhibitor spreads faster than the activator, creating local activation surrounded by lateral inhibition — the recipe for patterns.

### Dispersion Relation

The growth rate σ(k) of a perturbation with wavenumber k:

```
σ(k) = −[J11 + J22 + (Du + Dv)k²]/2 ± √{[J11 − J22 + (Du − Dv)k²]²/4 − det(J − k²D)}
```

Where J is the Jacobian of the reaction kinetics and D = diag(Du, Dv). Turing instability requires σ(k) > 0 for some k > 0.

### Laplacian (1D, periodic)

```
∇²u_i = (u_{i-1} − 2u_i + u_{i+1}) / Δx²
```

With periodic boundary conditions: u_{-1} = u_{N-1} and u_N = u_0.

## Installation

```bash
cargo add morphogenesis
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
morphogenesis = "0.1"
```

## Related Crates

Part of the **SuperInstance Exocortex** ecosystem:

- **[markov-blanket](https://github.com/SuperInstance/markov-blanket)** — Statistical boundary between agent and world
- **[free-energy](https://github.com/SuperInstance/free-energy)** — Variational free energy computation
- **[active-inference](https://github.com/SuperInstance/active-inference)** — Action as surprise minimization
- **[signal-transduction](https://github.com/SuperInstance/signal-transduction)** — Biological signal cascading
- **[graph-homology](https://github.com/SuperInstance/graph-homology)** — Topological structure of graphs

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project.
