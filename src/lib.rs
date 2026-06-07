//! Turing pattern formation for self-organizing systems.
//!
//! Implements reaction-diffusion dynamics, Turing instability analysis,
//! wavelength computation, symmetry breaking, and morphogen signaling.

// ── Module: reaction_diffusion ───────────────────────────────────────────────

pub mod reaction_diffusion {
    /// A 1D reaction-diffusion system with two morphogens (activator u, inhibitor v).
    #[derive(Clone, Debug)]
    pub struct ReactionDiffusion1D {
        pub u: Vec<f64>,
        pub v: Vec<f64>,
        pub du: f64,  // diffusion rate of activator
        pub dv: f64,  // diffusion rate of inhibitor
        pub dt: f64,  // time step
        pub dx: f64,  // spatial step
        pub f: f64,   // feed rate
        pub k: f64,   // kill rate
    }

    impl ReactionDiffusion1D {
        pub fn new(size: usize, du: f64, dv: f64, dt: f64, dx: f64, f: f64, k: f64) -> Self {
            ReactionDiffusion1D {
                u: vec![1.0; size],
                v: vec![0.0; size],
                du, dv, dt, dx, f, k,
            }
        }

        /// Initialize with a perturbation at the center.
        pub fn init_perturbation(&mut self) {
            let mid = self.u.len() / 2;
            let radius = self.u.len() / 10;
            for i in mid.saturating_sub(radius)..=(mid + radius).min(self.u.len() - 1) {
                self.u[i] = 0.5;
                self.v[i] = 0.25;
            }
        }

        /// Initialize with random perturbation.
        pub fn init_random(&mut self, seed: u64) {
            let mut state = seed;
            for i in 0..self.u.len() {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let rand_val = ((state >> 33) as f64) / (1u64 << 31) as f64;
                if rand_val < 0.1 {
                    self.u[i] = 0.5 + rand_val;
                    self.v[i] = 0.25 + rand_val * 0.5;
                }
            }
        }

        /// Compute Laplacian of a 1D field with periodic boundary.
        pub fn laplacian_1d(field: &[f64], dx: f64) -> Vec<f64> {
            let n = field.len();
            let mut lap = vec![0.0; n];
            for i in 0..n {
                let left = if i == 0 { field[n - 1] } else { field[i - 1] };
                let right = if i == n - 1 { field[0] } else { field[i + 1] };
                lap[i] = (left - 2.0 * field[i] + right) / (dx * dx);
            }
            lap
        }

        /// Step forward in time (Gray-Scott model).
        pub fn step(&mut self) {
            let lap_u = Self::laplacian_1d(&self.u, self.dx);
            let lap_v = Self::laplacian_1d(&self.v, self.dx);
            let n = self.u.len();
            let mut new_u = vec![0.0; n];
            let mut new_v = vec![0.0; n];
            for i in 0..n {
                let uvv = self.u[i] * self.v[i] * self.v[i];
                new_u[i] = self.u[i] + self.dt * (self.du * lap_u[i] - uvv + self.f * (1.0 - self.u[i]));
                new_v[i] = self.v[i] + self.dt * (self.dv * lap_v[i] + uvv - (self.f + self.k) * self.v[i]);
                new_u[i] = new_u[i].max(0.0);
                new_v[i] = new_v[i].max(0.0);
            }
            self.u = new_u;
            self.v = new_v;
        }

        /// Run multiple steps.
        pub fn run(&mut self, steps: usize) {
            for _ in 0..steps {
                self.step();
            }
        }

        /// Get the current pattern (u field).
        pub fn pattern(&self) -> &[f64] {
            &self.u
        }

        /// Get pattern amplitude (max - min of u).
        pub fn amplitude(&self) -> f64 {
            let max = self.u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let min = self.u.iter().cloned().fold(f64::INFINITY, f64::min);
            max - min
        }

        /// Get pattern mean.
        pub fn mean_u(&self) -> f64 {
            self.u.iter().sum::<f64>() / self.u.len() as f64
        }

        /// Get pattern mean for v.
        pub fn mean_v(&self) -> f64 {
            self.v.iter().sum::<f64>() / self.v.len() as f64
        }

        /// Check if NaN has appeared (numerical instability).
        pub fn is_stable(&self) -> bool {
            self.u.iter().all(|x| x.is_finite()) && self.v.iter().all(|x| x.is_finite())
        }

        /// Size of the domain.
        pub fn size(&self) -> usize {
            self.u.len()
        }
    }
}

// ── Module: turing_instability ───────────────────────────────────────────────

pub mod turing_instability {
    /// Check if a reaction-diffusion system exhibits Turing instability.
    /// Conditions: (1) stable without diffusion, (2) unstable with diffusion.
    pub fn is_turing_unstable(
        fu: f64, fv: f64, gu: f64, gv: f64,
        du: f64, dv: f64,
    ) -> bool {
        // Condition 1: stable without diffusion
        let trace = fu + gv;
        let det = fu * gv - fv * gu;
        let stable_no_diff = trace < 0.0 && det > 0.0;

        // Condition 2: exists k such that det(J - k²D) < 0
        // where D = diag(du, dv)
        // This happens when du * gv + dv * fu > 0 (with du != dv typically)
        let can_destabilize = if du != dv {
            du * gv + dv * fu > 0.0
        } else {
            false
        };

        stable_no_diff && can_destabilize
    }

    /// Compute the critical wavenumber for Turing instability.
    pub fn critical_wavenumber(fu: f64, _fv: f64, _gu: f64, gv: f64, du: f64, dv: f64) -> f64 {
        let numerator = du * gv + dv * fu;
        let denominator = 2.0 * du * dv;
        if denominator > 0.0 && numerator > 0.0 {
            (numerator / denominator).sqrt()
        } else {
            0.0
        }
    }

    /// Compute the Turing space (range of parameters leading to patterns).
    pub fn in_turing_space(
        fu: f64, fv: f64, gu: f64, gv: f64,
        du: f64, dv: f64,
    ) -> bool {
        is_turing_unstable(fu, fv, gu, gv, du, dv)
    }

    /// Compute dispersion relation: σ(k²) for a given wavenumber.
    pub fn dispersion_relation(
        k_sq: f64, fu: f64, fv: f64, gu: f64, gv: f64,
        du: f64, dv: f64,
    ) -> f64 {
        let trace = (fu + gv) - k_sq * (du + dv);
        let det = (fu - k_sq * du) * (gv - k_sq * dv) - fv * gu;
        // Growth rate is the larger eigenvalue
        let discriminant = trace * trace - 4.0 * det;
        if discriminant >= 0.0 {
            (trace + discriminant.sqrt()) / 2.0
        } else {
            trace / 2.0
        }
    }

    /// Find the wavenumber with maximum growth rate.
    pub fn fastest_growing_mode(
        fu: f64, fv: f64, gu: f64, gv: f64,
        du: f64, dv: f64,
    ) -> f64 {
        let mut best_k: f64 = 0.0;
        let mut best_sigma: f64 = f64::NEG_INFINITY;
        for i in 0..1000 {
            let k = i as f64 * 0.1;
            let sigma = dispersion_relation(k * k, fu, fv, gu, gv, du, dv);
            if sigma > best_sigma {
                best_sigma = sigma;
                best_k = k;
            }
        }
        best_k
    }

    /// Compute stability margin (how far from instability boundary).
    pub fn stability_margin(
        fu: f64, fv: f64, gu: f64, gv: f64,
        du: f64, dv: f64,
    ) -> f64 {
        let kc = critical_wavenumber(fu, fv, gu, gv, du, dv);
        if kc > 0.0 {
            dispersion_relation(kc * kc, fu, fv, gu, gv, du, dv)
        } else {
            0.0
        }
    }

    /// Check if homogeneous steady state is stable (without diffusion).
    pub fn homogeneous_stable(fu: f64, fv: f64, gu: f64, gv: f64) -> bool {
        let trace = fu + gv;
        let det = fu * gv - fv * gu;
        trace < 0.0 && det > 0.0
    }
}

// ── Module: wavelength ───────────────────────────────────────────────────────

pub mod wavelength {
    /// Compute wavelength from wavenumber: λ = 2π/k.
    pub fn from_wavenumber(k: f64) -> f64 {
        if k > 0.0 { 2.0 * std::f64::consts::PI / k } else { f64::MAX }
    }

    /// Compute wavenumber from wavelength: k = 2π/λ.
    pub fn from_wavelength(wavelength: f64) -> f64 {
        if wavelength > 0.0 { 2.0 * std::f64::consts::PI / wavelength } else { 0.0 }
    }

    /// Estimate dominant wavelength from a 1D pattern using autocorrelation.
    pub fn dominant_wavelength(pattern: &[f64], dx: f64) -> f64 {
        let n = pattern.len();
        if n < 4 { return 0.0; }
        let mean = pattern.iter().sum::<f64>() / n as f64;
        let centered: Vec<f64> = pattern.iter().map(|&x| x - mean).collect();

        // Find first peak in autocorrelation
        let variance: f64 = centered.iter().map(|x| x * x).sum();
        if variance == 0.0 { return 0.0; }

        let mut best_lag = 0;
        let mut best_corr = f64::NEG_INFINITY;
        for lag in 1..n / 2 {
            let mut corr = 0.0;
            for i in 0..(n - lag) {
                corr += centered[i] * centered[i + lag];
            }
            corr /= variance;
            if corr > best_corr {
                best_corr = corr;
                best_lag = lag;
            }
        }
        best_lag as f64 * dx
    }

    /// Compute number of wavelengths that fit in a domain.
    pub fn wavelengths_in_domain(domain_length: f64, wavelength: f64) -> f64 {
        if wavelength > 0.0 { domain_length / wavelength } else { 0.0 }
    }

    /// Compute wavelength from Turing parameters.
    pub fn turing_wavelength(du: f64, dv: f64, fu: f64, fv: f64, gu: f64, gv: f64) -> f64 {
        let kc = super::turing_instability::critical_wavenumber(fu, fv, gu, gv, du, dv);
        from_wavenumber(kc)
    }

    /// Compute pattern wavelength range (band of unstable modes).
    pub fn unstable_band(
        fu: f64, fv: f64, gu: f64, gv: f64,
        du: f64, dv: f64,
    ) -> (f64, f64) {
        let mut k_min: f64 = f64::MAX;
        let mut k_max: f64 = 0.0;
        for i in 1..10000 {
            let k = i as f64 * 0.01;
            let sigma = super::turing_instability::dispersion_relation(
                k * k, fu, fv, gu, gv, du, dv
            );
            if sigma > 0.0 {
                k_min = k_min.min(k);
                k_max = k_max.max(k);
            }
        }
        (from_wavenumber(k_max), from_wavenumber(k_min))
    }
}

// ── Module: symmetry_breaking ────────────────────────────────────────────────

pub mod symmetry_breaking {
    /// Detect if a pattern has broken symmetry (non-uniform).
    pub fn is_symmetry_broken(pattern: &[f64], threshold: f64) -> bool {
        let mean = pattern.iter().sum::<f64>() / pattern.len() as f64;
        let variance: f64 = pattern.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>() / pattern.len() as f64;
        variance > threshold
    }

    /// Compute symmetry breaking measure (normalized variance).
    pub fn symmetry_breaking_measure(pattern: &[f64]) -> f64 {
        let mean = pattern.iter().sum::<f64>() / pattern.len() as f64;
        if mean.abs() < 1e-10 { return 0.0; }
        let variance: f64 = pattern.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>() / pattern.len() as f64;
        variance.sqrt() / mean.abs()
    }

    /// Find the time step at which symmetry breaks.
    pub fn symmetry_breaking_time(
        patterns: &[Vec<f64>],
        threshold: f64,
    ) -> Option<usize> {
        patterns.iter().position(|p| is_symmetry_broken(p, threshold))
    }

    /// Compute symmetry order parameter (0 = symmetric, 1 = fully broken).
    pub fn order_parameter(pattern: &[f64]) -> f64 {
        let n = pattern.len();
        if n < 2 { return 0.0; }
        let mean = pattern.iter().sum::<f64>() / n as f64;
        let sum_sq: f64 = pattern.iter().map(|&x| (x - mean).powi(2)).sum();
        let sum_abs: f64 = pattern.iter().map(|&x| (x - mean).abs()).sum();
        if sum_abs == 0.0 { return 0.0; }
        (sum_sq / n as f64).sqrt() / (sum_abs / n as f64)
    }

    /// Check for left-right symmetry.
    pub fn is_lr_symmetric(pattern: &[f64], tolerance: f64) -> bool {
        let n = pattern.len();
        for i in 0..n / 2 {
            if (pattern[i] - pattern[n - 1 - i]).abs() > tolerance {
                return false;
            }
        }
        true
    }

    /// Compute the symmetry breaking speed (rate of change of variance).
    pub fn breaking_speed(patterns: &[Vec<f64>]) -> Vec<f64> {
        patterns.windows(2).map(|w| {
            let m0 = symmetry_breaking_measure(&w[0]);
            let m1 = symmetry_breaking_measure(&w[1]);
            m1 - m0
        }).collect()
    }
}

// ── Module: morphogen ────────────────────────────────────────────────────────

pub mod morphogen {
    /// A morphogen with concentration and properties.
    #[derive(Clone, Debug)]
    pub struct Morphogen {
        pub name: String,
        pub concentration: Vec<f64>,
        pub diffusion_rate: f64,
        pub decay_rate: f64,
        pub production_rate: f64,
    }

    impl Morphogen {
        pub fn new(name: &str, size: usize, diffusion: f64, decay: f64, production: f64) -> Self {
            Morphogen {
                name: name.to_string(),
                concentration: vec![0.0; size],
                diffusion_rate: diffusion,
                decay_rate: decay,
                production_rate: production,
            }
        }

        /// Set uniform concentration.
        pub fn set_uniform(&mut self, value: f64) {
            for c in &mut self.concentration {
                *c = value;
            }
        }

        /// Add a source at a position.
        pub fn add_source(&mut self, position: usize, amount: f64) {
            if position < self.concentration.len() {
                self.concentration[position] += amount;
            }
        }

        /// Diffuse one step (1D).
        pub fn diffuse(&mut self, dx: f64, dt: f64) {
            let n = self.concentration.len();
            let lap = super::reaction_diffusion::ReactionDiffusion1D::laplacian_1d(
                &self.concentration, dx
            );
            for i in 0..n {
                self.concentration[i] += dt * self.diffusion_rate * lap[i];
                self.concentration[i] = self.concentration[i].max(0.0);
            }
        }

        /// Apply decay.
        pub fn decay(&mut self, dt: f64) {
            for c in &mut self.concentration {
                *c *= (1.0 - self.decay_rate * dt).max(0.0);
            }
        }

        /// Apply production.
        pub fn produce(&mut self, dt: f64) {
            for c in &mut self.concentration {
                *c += self.production_rate * dt;
            }
        }

        /// Get total amount.
        pub fn total(&self) -> f64 {
            self.concentration.iter().sum()
        }

        /// Get mean concentration.
        pub fn mean(&self) -> f64 {
            self.total() / self.concentration.len() as f64
        }

        /// Get gradient at a position.
        pub fn gradient(&self, pos: usize, dx: f64) -> f64 {
            if pos == 0 || pos >= self.concentration.len() - 1 {
                return 0.0;
            }
            (self.concentration[pos + 1] - self.concentration[pos - 1]) / (2.0 * dx)
        }

        /// Interpolate concentration at a continuous position.
        pub fn interpolate(&self, pos: f64) -> f64 {
            let idx = pos.floor() as usize;
            let frac = pos - pos.floor();
            if idx >= self.concentration.len() - 1 {
                return self.concentration.last().copied().unwrap_or(0.0);
            }
            self.concentration[idx] * (1.0 - frac) + self.concentration[idx + 1] * frac
        }

        /// Threshold the concentration.
        pub fn threshold(&self, threshold: f64) -> Vec<bool> {
            self.concentration.iter().map(|&c| c > threshold).collect()
        }

        /// Size of the concentration field.
        pub fn size(&self) -> usize {
            self.concentration.len()
        }
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Reaction-diffusion tests ──

    #[test]
    fn test_rd_creation() {
        let rd = reaction_diffusion::ReactionDiffusion1D::new(100, 0.2, 0.1, 1.0, 1.0, 0.04, 0.06);
        assert_eq!(rd.size(), 100);
    }

    #[test]
    fn test_rd_init_perturbation() {
        let mut rd = reaction_diffusion::ReactionDiffusion1D::new(100, 0.2, 0.1, 1.0, 1.0, 0.04, 0.06);
        rd.init_perturbation();
        assert!(rd.u[50] < 1.0);
        assert!(rd.v[50] > 0.0);
    }

    #[test]
    fn test_laplacian_uniform() {
        let field = vec![1.0; 10];
        let lap = reaction_diffusion::ReactionDiffusion1D::laplacian_1d(&field, 1.0);
        for l in &lap {
            assert!(l.abs() < 1e-10);
        }
    }

    #[test]
    fn test_laplacian_peak() {
        let mut field = vec![0.0; 11];
        field[5] = 1.0;
        let lap = reaction_diffusion::ReactionDiffusion1D::laplacian_1d(&field, 1.0);
        assert!(lap[5] < 0.0); // Peak should have negative Laplacian
    }

    #[test]
    fn test_rd_step() {
        let mut rd = reaction_diffusion::ReactionDiffusion1D::new(50, 0.2, 0.1, 0.1, 0.5, 0.04, 0.06);
        rd.init_perturbation();
        let u_before = rd.u.clone();
        rd.step();
        // Something should have changed
        assert_ne!(rd.u, u_before);
    }

    #[test]
    fn test_rd_stability() {
        let mut rd = reaction_diffusion::ReactionDiffusion1D::new(50, 0.2, 0.1, 0.5, 0.5, 0.04, 0.06);
        rd.init_perturbation();
        rd.run(100);
        assert!(rd.is_stable());
    }

    #[test]
    fn test_rd_amplitude() {
        let mut rd = reaction_diffusion::ReactionDiffusion1D::new(50, 0.2, 0.1, 1.0, 1.0, 0.04, 0.06);
        rd.init_perturbation();
        assert!(rd.amplitude() >= 0.0);
    }

    #[test]
    fn test_rd_mean() {
        let rd = reaction_diffusion::ReactionDiffusion1D::new(100, 0.2, 0.1, 1.0, 1.0, 0.04, 0.06);
        assert!((rd.mean_u() - 1.0).abs() < 1e-10);
        assert!((rd.mean_v()).abs() < 1e-10);
    }

    #[test]
    fn test_rd_init_random() {
        let mut rd = reaction_diffusion::ReactionDiffusion1D::new(100, 0.2, 0.1, 1.0, 1.0, 0.04, 0.06);
        rd.init_random(42);
        // Some values should differ from initial
        assert!(rd.v.iter().any(|&x| x > 0.0));
    }

    #[test]
    fn test_rd_run_multiple() {
        let mut rd = reaction_diffusion::ReactionDiffusion1D::new(50, 0.2, 0.1, 0.5, 0.5, 0.04, 0.06);
        rd.init_perturbation();
        rd.run(50);
        assert!(rd.is_stable());
    }

    // ── Turing instability tests ──

    #[test]
    fn test_turing_unstable_case() {
        // Classic activator-inhibitor parameters
        let unstable = turing_instability::is_turing_unstable(
            1.0, -1.0, 2.0, -1.5, 1.0, 10.0
        );
        assert!(unstable);
    }

    #[test]
    fn test_not_turing_unstable_same_diffusion() {
        let unstable = turing_instability::is_turing_unstable(
            -1.0, 1.0, -1.0, -1.0, 1.0, 1.0
        );
        assert!(!unstable);
    }

    #[test]
    fn test_critical_wavenumber() {
        let kc = turing_instability::critical_wavenumber(1.0, -1.0, 2.0, -1.5, 1.0, 10.0);
        assert!(kc > 0.0);
    }

    #[test]
    fn test_dispersion_relation() {
        let sigma = turing_instability::dispersion_relation(
            1.0, 1.0, -1.0, 2.0, -1.5, 1.0, 10.0
        );
        // Should be finite
        assert!(sigma.is_finite());
    }

    #[test]
    fn test_fastest_growing_mode() {
        let k = turing_instability::fastest_growing_mode(1.0, -1.0, 2.0, -1.5, 1.0, 10.0);
        assert!(k > 0.0);
    }

    #[test]
    fn test_stability_margin() {
        let margin = turing_instability::stability_margin(1.0, -1.0, 2.0, -1.5, 1.0, 10.0);
        // Should be positive for unstable system
        assert!(margin > 0.0);
    }

    #[test]
    fn test_homogeneous_stable() {
        let stable = turing_instability::homogeneous_stable(-2.0, 1.0, -1.0, -1.5);
        let trace = -2.0 + (-1.5);
        let det = (-2.0) * (-1.5) - 1.0 * (-1.0);
        assert_eq!(stable, trace < 0.0 && det > 0.0);
    }

    #[test]
    fn test_in_turing_space() {
        let in_space = turing_instability::in_turing_space(1.0, -1.0, 2.0, -1.5, 1.0, 10.0);
        assert!(in_space);
    }

    // ── Wavelength tests ──

    #[test]
    fn test_wavelength_from_wavenumber() {
        let lambda = wavelength::from_wavenumber(1.0);
        assert!((lambda - 2.0 * std::f64::consts::PI).abs() < 1e-10);
    }

    #[test]
    fn test_wavenumber_from_wavelength() {
        let k = wavelength::from_wavelength(2.0 * std::f64::consts::PI);
        assert!((k - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_roundtrip_wavelength() {
        let k = 2.5;
        let lambda = wavelength::from_wavenumber(k);
        let k2 = wavelength::from_wavelength(lambda);
        assert!((k - k2).abs() < 1e-10);
    }

    #[test]
    fn test_dominant_wavelength_uniform() {
        let pattern = vec![1.0; 50];
        let wl = wavelength::dominant_wavelength(&pattern, 1.0);
        assert_eq!(wl, 0.0);
    }

    #[test]
    fn test_wavelengths_in_domain() {
        let n = wavelength::wavelengths_in_domain(10.0, 2.0);
        assert!((n - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_turing_wavelength() {
        let wl = wavelength::turing_wavelength(1.0, 10.0, 1.0, -1.0, 2.0, -1.5);
        assert!(wl > 0.0);
    }

    #[test]
    fn test_unstable_band() {
        let (wl_max, wl_min) = wavelength::unstable_band(1.0, -1.0, 2.0, -1.5, 1.0, 10.0);
        // Should return valid range
        assert!(wl_min >= wl_max); // Note: higher k = lower wavelength
    }

    // ── Symmetry breaking tests ──

    #[test]
    fn test_symmetry_not_broken() {
        let pattern = vec![1.0; 50];
        assert!(!symmetry_breaking::is_symmetry_broken(&pattern, 0.01));
    }

    #[test]
    fn test_symmetry_broken() {
        let mut pattern = vec![1.0; 50];
        pattern[25] = 5.0;
        assert!(symmetry_breaking::is_symmetry_broken(&pattern, 0.01));
    }

    #[test]
    fn test_symmetry_breaking_measure_uniform() {
        let pattern = vec![1.0; 50];
        let m = symmetry_breaking::symmetry_breaking_measure(&pattern);
        assert!(m.abs() < 1e-10);
    }

    #[test]
    fn test_symmetry_breaking_measure_patterned() {
        let pattern: Vec<f64> = (0..50).map(|i| if i < 25 { 0.5 } else { 1.5 }).collect();
        let m = symmetry_breaking::symmetry_breaking_measure(&pattern);
        assert!(m > 0.0);
    }

    #[test]
    fn test_symmetry_breaking_time() {
        let uniform = vec![1.0; 50];
        let mut patterned = vec![1.0; 50];
        patterned[25] = 5.0;
        let patterns = vec![uniform.clone(), uniform.clone(), patterned];
        let t = symmetry_breaking::symmetry_breaking_time(&patterns, 0.01);
        assert_eq!(t, Some(2));
    }

    #[test]
    fn test_order_parameter_uniform() {
        let pattern = vec![1.0; 50];
        let op = symmetry_breaking::order_parameter(&pattern);
        assert!(op.abs() < 1e-10);
    }

    #[test]
    fn test_lr_symmetric() {
        let pattern = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        assert!(symmetry_breaking::is_lr_symmetric(&pattern, 0.01));
    }

    #[test]
    fn test_not_lr_symmetric() {
        let pattern = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(!symmetry_breaking::is_lr_symmetric(&pattern, 0.01));
    }

    #[test]
    fn test_breaking_speed() {
        let uniform = vec![1.0; 50];
        let mut p1 = vec![1.0; 50];
        p1[25] = 1.5;
        let mut p2 = vec![1.0; 50];
        p2[25] = 2.0;
        let speeds = symmetry_breaking::breaking_speed(&[uniform, p1, p2]);
        assert_eq!(speeds.len(), 2);
        assert!(speeds[1] > speeds[0] || speeds[0] > 0.0);
    }

    // ── Morphogen tests ──

    #[test]
    fn test_morphogen_creation() {
        let m = morphogen::Morphogen::new("activator", 50, 0.1, 0.01, 0.05);
        assert_eq!(m.size(), 50);
        assert_eq!(m.name, "activator");
    }

    #[test]
    fn test_set_uniform() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 0.05);
        m.set_uniform(2.0);
        assert!(m.concentration.iter().all(|&c| (c - 2.0).abs() < 1e-10));
    }

    #[test]
    fn test_add_source() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 0.05);
        m.add_source(5, 3.0);
        assert!((m.concentration[5] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_diffuse() {
        let mut m = morphogen::Morphogen::new("test", 20, 0.1, 0.01, 0.05);
        m.add_source(10, 10.0);
        m.diffuse(1.0, 0.1);
        // Should spread
        assert!(m.concentration[9] > 0.0);
        assert!(m.concentration[11] > 0.0);
    }

    #[test]
    fn test_decay() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.5, 0.05);
        m.set_uniform(2.0);
        m.decay(1.0);
        assert!(m.concentration[0] < 2.0);
    }

    #[test]
    fn test_produce() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 1.0);
        m.produce(1.0);
        assert!((m.concentration[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_total() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 0.05);
        m.set_uniform(3.0);
        assert!((m.total() - 30.0).abs() < 1e-10);
    }

    #[test]
    fn test_mean() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 0.05);
        m.set_uniform(5.0);
        assert!((m.mean() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_gradient() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 0.05);
        for i in 0..10 {
            m.concentration[i] = i as f64;
        }
        let g = m.gradient(5, 1.0);
        assert!((g - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate() {
        let mut m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 0.05);
        m.concentration[3] = 3.0;
        m.concentration[4] = 5.0;
        let val = m.interpolate(3.5);
        assert!((val - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_threshold() {
        let mut m = morphogen::Morphogen::new("test", 5, 0.1, 0.01, 0.05);
        m.concentration = vec![0.1, 0.3, 0.5, 0.7, 0.9];
        let t = m.threshold(0.5);
        assert_eq!(t, vec![false, false, false, true, true]);
    }

    #[test]
    fn test_gradient_boundary() {
        let m = morphogen::Morphogen::new("test", 10, 0.1, 0.01, 0.05);
        assert_eq!(m.gradient(0, 1.0), 0.0);
        assert_eq!(m.gradient(9, 1.0), 0.0);
    }

    #[test]
    fn test_interpolate_boundary() {
        let mut m = morphogen::Morphogen::new("test", 5, 0.1, 0.01, 0.05);
        m.concentration[4] = 7.0;
        let val = m.interpolate(10.0);
        assert!((val - 7.0).abs() < 1e-10);
    }
}
