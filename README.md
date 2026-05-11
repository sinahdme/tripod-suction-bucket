# Tripod Suction Bucket Foundation Simulator

A MATLAB framework for analyzing **tripod (3-bucket) suction bucket foundations** for offshore wind turbines under multidirectional lateral-moment loading.

The model uses a Winkler spring approach with coupled **p-y**, **t-z**, and **q-z** soil-structure interaction curves to compute foundation capacity, validated against centrifuge test data from Kim et al. (2014).

> **Note:** The full validated framework and parametric study results are currently under journal review. This repository demonstrates the modeling approach and code structure.

---

## Quick Start

```matlab
% Open MATLAB, navigate to this folder, then run:
THREE_POD_MAINCODE
```

This generates **H-theta** (lateral load vs. rotation) and **M-theta** (moment vs. rotation) curves for a tripod suction bucket foundation with representative soil parameters.

---

## File Reference

### Main Script

| File | Description |
|------|-------------|
| `THREE_POD_MAINCODE.m` | Main simulation entry point. Sets up tripod geometry (buckets at 0, 120, 240 deg), soil parameters, bucket mesh, and runs the iterative solver to compute the lateral-moment capacity envelope. |

### Validation & Polar Sweep

| File | Description |
|------|-------------|
| `THREE_POD_KIM2014_VALIDATION_LD123.m` | Validates the model against Kim et al. (2014) centrifuge tests (D=6.5 m, L=8.0 m, L/D=1.23, Dr=82.7%, 70g). Uses `pycurve_LD123` with extrapolated beta coefficients and group interaction p-multiplier. Target: M_y ~ 93 MN-m at 0.6 deg rotation. |
| `THREE_POD_POLAR.m` | Sweeps loading direction 0-360 deg at fixed eccentricity to generate the multidirectional moment capacity envelope. Higher resolution mesh (n=32, n_l=32). |
| `THREE_POD_POLAR_ver02.m` | Polar sweep version 2 with refined discretization. |
| `THREE_POD_POLAR_ver03.m` | Polar sweep version 3, latest iteration. |

### p-y Curve Functions (Lateral Springs)

| File | Signature | Description |
|------|-----------|-------------|
| `pycurve.m` | `[P1] = pycurve(fi, gamma, d, z, do_b, l_b, n)` | Standard p-y spring curve. Computes passive/active Rankine coefficients (Kp, Ka, K0), uses shape coefficients beta1-beta4 from `BETA13` and `BETA24`, and returns lateral resistance via hyperbolic tangent formulation: `P = P_R * (beta1*tanh(beta2*d) + beta3*tanh(beta4*d))`. |
| `pycurve2.m` | `[P1] = pycurve2(fi, gamma, d, z, do_b, l_b, n, n_l, A_c, m_c)` | Alternative p-y curve using oedometric modulus (`E_oed = m_c * 100`) for soil stiffness variation with depth. |
| `pycurve_LD123.m` | `[P1] = pycurve_LD123(fi, gamma, d, z, do_b, l_b, n)` | Extrapolated p-y curve for L/D=1.23 geometry (Kim 2014). Hardcoded beta coefficients: beta1=1.5228, beta2=101.96, beta3=0.6396, beta4=25.71. |

### t-z Curve Functions (Vertical Skin Friction Springs)

| File | Signature | Description |
|------|-----------|-------------|
| `tzcurve.m` | `[tz_in, tz_out] = tzcurve(P_i_out, fi, gamma, do_b, di_b, z, disp_mm, l_b, n_l, n, fi_crit, delta, D_r, e_max, e_min)` | Vertical skin friction spring model. Implements critical-state soil mechanics with Bolton dilatancy. Handles inside and outside bucket skirt friction, with mudline partial-embedding logic for nodes crossing z=0. |
| `tzcurve_deadlaod.m` | *(similar to tzcurve)* | t-z curve variant for computing skin friction under static dead-load settlement (before rotation begins). |

### q-z Curve Functions (Tip & Lid Bearing Springs)

| File | Signature | Description |
|------|-----------|-------------|
| `qzcurve2.m` | `[qb1] = qzcurve2(fi_end, tip_end, tip_end_fixed, gamma, n, thickness, Di, Do)` | Skirt-tip bearing resistance. Linear Nq interpolation as f(phi), normalized displacement ratio W/D, iterative solver for coupled spring stiffness. |
| `qzcurve3.m` | *(variant)* | Alternative q-z formulation (intermediate version). |
| `qzcurve4.m` | `[qb2] = qzcurve4(fi, Z_lid, r_tip_end, r_tip_end_fixed, gamma, n, n_q, R_split, Di, D_r, S_gamma)` | Lid bearing resistance with radial discretization (`R_split` array). Elastic-plastic coupling (K_e_qz, K_p_qz). More detailed than `qzcurve2`. |

### Force Coefficient Functions

| File | Signature | Description |
|------|-----------|-------------|
| `coeffinder.m` | `[KKx, KKy, P_i_out] = coeffinder(endd, endd_fixed, z, l_b, n_l, fi, n, do_b, gamma, Q, A_c, m_c, translation_vector, th1, f_g)` | Evaluates p-y spring forces at each node. Calls `pycurve`, extracts displacements/angles from the deformation gradient, and applies the group interaction p-multiplier (f_g). Returns horizontal force components KKx, KKy and lateral pressure P_i_out for t-z coupling. |
| `coeffinder_LD123.m` | *(same as coeffinder)* | Force coefficient finder using `pycurve_LD123` for Kim 2014 calibrated beta coefficients. |
| `coeffinder_deadload.m` | *(variant)* | Computes p-y forces affecting t-z skin friction during dead-load settlement (before rotation). |

### Shape Coefficient Functions

| File | Signature | Description |
|------|-----------|-------------|
| `BETA13.m` | `[beta_1, beta_3] = BETA13(a1, a2, b1, b2, fi, L)` | Computes shape coefficients beta1 and beta3 for the p-y curve by solving a quadratic: `p = [1, -(a1*fi/L+a2), (b1*fi/L+b2)]`. Calibrated from Knudsen et al. (2013). |
| `BETA24.m` | `[beta_2, beta_4] = BETA24(c1, c2, c3, d1, d2, d3, fi, L)` | Computes shape coefficients beta2 and beta4 with higher-order `(fi/L)^2` terms. |

### Geometry & Transformation Functions

| File | Signature | Description |
|------|-----------|-------------|
| `transformation.m` | `transformedPoint = transformation(P, Q, rotInput)` | 3D point rotation about center Q. Supports legacy Euler angles `[phi, theta, psi]` or direct 3x3 rotation matrix input. Core geometric transformation for computing bucket positions after tilt. |
| `rodrigues_rotation.m` | `R = rodrigues_rotation(dd, alpha_deg)` | Builds a rotation matrix for tilt magnitude `dd` (rad) in direction `alpha` (deg) using the Rodrigues formula. Modern replacement for Euler-angle approach. |
| `points.m` | *(complex signature)* | Generates bucket node positions: circumferential (`n`) x axial (`n_l`) mesh of skirt nodes, pile tip nodes, and lid nodes for each bucket. Core geometry builder. |
| `points_deadload.m` | *(variant)* | Generates node positions for dead-load (vertical-only) analysis. |
| `relative_points.m` | *(complex signature)* | Computes relative displacement between body-frame and fixed-frame node positions for each bucket. Differentiates fixed-frame and body-frame node sets. |
| `LocalFrameFunction.m` | `[new_position] = LocalFrameFunction(translation_vector, Theta)` | Defines the local coordinate frame for each bucket via translation and z-axis rotation. Used by `coeffinder` and `points` to position buckets. |

### Soil Mechanics Functions

| File | Signature | Description |
|------|-----------|-------------|
| `Rankine.m` | `[PR] = Rankine(fi, gamma, d, z, do_b, l_b, n, n_l)` | Computes Rankine earth pressure (passive/active) as reference for p-y curve. |
| `Zvalues.m` | -- | Calculates the embedded depth of each spring node in the soil. |
| `fi_finder.m` | -- | Determines peak friction angle as a function of relative density (Dr) using Bolton-type formula. |
| `apdefiner.m` | -- | Calculates coefficient alpha_p to account for the effect of lateral soil resistance on interface normal stress. Called within `tzcurve`. |
| `apdefiner_general.m` | `ap_out = apdefiner_general(ap, KK, n, n_l)` | Generalized alpha_p mirror: transfers ap from passive to active (tension) side using actual p-y force distribution rather than hardcoded axes. Direction-independent. |

### Group Interaction Functions

| File | Signature | Description |
|------|-----------|-------------|
| `detect_fg_per_bucket.m` | `fg_vec = detect_fg_per_bucket(s_end_cells, s_fixed_cells, adj_matrix, SD_ratio, delta_ref, LD_ratio, config)` | Computes displacement-dependent group p-multiplier (f_g) per bucket. Ramps from 1.0 (no reduction) to f_g_min based on differential compression. Uses enriched Barari et al. (2023) database with 2D interpolation for S/D and L/D. |
| `apply_group_pmult.m` | `[KKx_out, KKy_out, fg_applied] = apply_group_pmult(...)` | Applies the p-multiplier to p-y forces. Detects compression vs. tension buckets based on vertical displacement, then applies f_g only to compression buckets that have adjacent compression neighbors (passive wedge overlap). |
| `barari_database.m` | -- | Enriched Barari et al. (2023) lookup table for group efficiency (E_cM) as a function of L/D and S/D ratios. |

### Visualization

| File | Description |
|------|-------------|
| `tripod_visualization.m` | 3D visualization of tripod bucket configuration showing untilted vs. tilted states, compression (red/orange) vs. tension (blue) buckets. |
| `animate_3d_heatmaps.m` | Generates animated GIFs (`py_3d_animation.gif`, `tz_3d_animation.gif`) of p-y and t-z force distributions across bucket surfaces during tilt. Auto-detects tripod (3 buckets) vs. tetrapod (4 buckets). |

---

## Model Architecture

```
THREE_POD_MAINCODE.m
├── points.m / points_deadload.m        (generate bucket node mesh)
│   └── LocalFrameFunction.m            (local coordinate frame per bucket)
├── transformation.m / rodrigues_rotation.m  (rotate nodes for tilt)
├── relative_points.m                   (compute body-frame vs fixed-frame displacements)
│
├── [Iterative solver loop]
│   ├── coeffinder.m                    (evaluate p-y spring forces)
│   │   ├── pycurve.m                   (p-y curve model)
│   │   │   ├── BETA13.m / BETA24.m     (shape coefficients)
│   │   │   └── Rankine.m              (earth pressure reference)
│   │   └── fi_finder.m                (peak friction angle)
│   │
│   ├── tzcurve.m                       (t-z vertical skin friction)
│   │   └── apdefiner.m / apdefiner_general.m  (lateral-to-normal stress coupling)
│   │
│   ├── qzcurve2.m                      (skirt-tip bearing)
│   ├── qzcurve4.m                      (lid bearing with radial discretization)
│   │
│   └── detect_fg_per_bucket.m          (group interaction p-multiplier)
│       ├── apply_group_pmult.m         (apply f_g to compression buckets)
│       └── barari_database.m           (efficiency lookup table)
│
└── [Output: H-theta, M-theta curves]
```

---

## Tripod Configuration

```
         Bucket 1 (0 deg)
            *
           / \
          /   \
         /     \
        /   Q   \         Q = rotation center (centroid)
       /         \
      *-----------*
  Bucket 2       Bucket 3
  (120 deg)      (240 deg)
```

- **3 buckets** at 120 deg spacing on a circle of radius `R_tripod = S_cc / sqrt(3)`
- **Symmetry:** 120-deg periodicity in moment capacity envelope

---

## Default Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Bucket diameter | D | 10.0 | m |
| Skirt length | L | 14.0 | m |
| Wall thickness | t | 0.030 | m |
| Center-to-center spacing | S_cc | 25.5 | m |
| Tower height (eccentricity) | e | 100.0 | m |
| Submerged unit weight | gamma' | 9.5 | kN/m3 |
| Relative density | Dr | 0.70 | - |
| Critical friction angle | phi_cv | 33.0 | deg |

---

## Validation

The model has been validated against:

- **Kim et al. (2014)** -- Centrifuge tests on tripod suction bucket foundations at 70g (D=6.5 m, L=8.0 m, L/D=1.23, Dr=82.7%). See `THREE_POD_KIM2014_VALIDATION_LD123.m`.

---

## Requirements

- MATLAB R2018b or later (uses basic matrix operations, no special toolboxes required)

---

## Author

**Sina Hadadi**
Ph.D. Candidate -- Wind Energy (Offshore Foundations)
Kunsan National University, South Korea

Email: sinahdme@gmail.com
GitHub: [github.com/sinahdme](https://github.com/sinahdme)

---

## License

This project is shared for academic and research purposes. Please cite appropriately if used in publications.
"# tripod-suction-bucket" 
"# tripod-suction-bucket" 
"# tripod-suction-bucket" 
