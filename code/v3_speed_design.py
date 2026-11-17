"""Recalculate E2 sweep identifiability and propagate uncertainty through v3."""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


DESIGNS = {
    'six_speeds': (5, 10, 20, 50, 100, 200),
    'four_speeds': (5, 10, 20, 50),
    'min_20': (20, 50, 100, 200),
}
VS_SCENARIOS = (5, 10, 20)  # mm/s
REPEATS_PER_SPEED = 3
RUN_MU_SD = 0.01
SIMULATED_DATASETS = 1000
SPAN_SPEEDS = (10, 20)


def steady_mu(speed, mu_s, mu_k, log_vs):
    fraction = np.exp(-(speed / np.exp(log_vs)) ** 2)
    return mu_k + (mu_s - mu_k) * fraction


def steady_jacobian(speed, mu_s, mu_k, log_vs):
    a = (speed / np.exp(log_vs)) ** 2
    fraction = np.exp(-a)
    return np.column_stack((-np.expm1(-a), 2 * (mu_s - mu_k) * a * fraction))


def fit(speed, measured_mu, mu_s, nominal_mu_k):
    result = least_squares(
        lambda p: steady_mu(speed, mu_s, *p) - measured_mu,
        [nominal_mu_k, np.log(10.0)],
        jac=lambda p: steady_jacobian(speed, mu_s, *p),
        bounds=([0.0, np.log(1.0)], [mu_s, np.log(100.0)]),
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.x


def assess(output, estimates, mu_k, span, save_rows):
    mu_s = estimates['center']
    mean_sd = RUN_MU_SD / np.sqrt(REPEATS_PER_SPEED)
    rng = np.random.default_rng(20261112)
    designs, recovery, ratios, predictions = [], [], [], []
    for vs in VS_SCENARIOS:
        params = np.array([mu_k, np.log(vs)])
        all_speeds = np.array(DESIGNS['six_speeds'])
        steady = steady_mu(all_speeds, mu_s, *params)
        reference = steady[all_speeds == 50].item()
        for speed, coefficient in zip(all_speeds, steady):
            ratios.append(dict(true_vs_mm_s=vs, speed_mm_s=speed, mu_steady=coefficient,
                               m1_ratio_to_50=1.0, m2_ratio_to_50=coefficient / reference))

        for name, speed_values in DESIGNS.items():
            speeds = np.array(speed_values)
            expected = steady_mu(speeds, mu_s, *params)
            jac = steady_jacobian(speeds, mu_s, *params)
            step = 1e-4 * np.eye(2)
            finite_diff = np.column_stack([
                (steady_mu(speeds, mu_s, *(params + h)) - steady_mu(speeds, mu_s, *(params - h)))
                / (2 * h.max()) for h in step
            ])
            np.testing.assert_allclose(jac, finite_diff, atol=2e-8, rtol=2e-5)
            inverse = np.linalg.pinv(jac)
            covariance = mean_sd**2 * (inverse @ inverse.T)
            mu_k_se, log_vs_se = np.sqrt(covariance.diagonal())
            designs.append(dict(design=name, speeds_mm_s=' '.join(map(str, speed_values)),
                                true_vs_mm_s=vs, repeats_per_speed=REPEATS_PER_SPEED,
                                assumed_run_mu_sd=RUN_MU_SD, jac_condition=np.linalg.cond(jac),
                                mu_k_se=mu_k_se, vs_local_se_mm_s=vs * log_vs_se,
                                ln_vs_se=log_vs_se))
            if name != 'four_speeds':
                continue

            np.testing.assert_allclose(fit(speeds, expected, mu_s, mu_k), params, atol=1e-6, rtol=0)
            synthetic = expected + rng.normal(0, mean_sd, (SIMULATED_DATASETS, len(speeds)))
            fitted = np.array([fit(speeds, sample, mu_s, mu_k) for sample in synthetic])
            lower_e1 = fit(speeds, expected, estimates['lower'], mu_k)
            upper_e1 = fit(speeds, expected, estimates['upper'], mu_k)
            for index, label, transform in ((0, 'mu_k', float), (1, 'vs_mm_s', np.exp)):
                low, median, high = np.quantile(fitted[:, index], [0.025, 0.5, 0.975])
                recovery.append(dict(true_vs_mm_s=vs, parameter=label, datasets=SIMULATED_DATASETS,
                                     truth=transform(params[index]), q025=transform(low),
                                     median=transform(median), q975=transform(high),
                                     fit_at_e1_lower=transform(lower_e1[index]),
                                     fit_at_e1_upper=transform(upper_e1[index])))

            for speed in SPAN_SPEEDS:
                nominal = span(speed, mu_s, params)
                gradient = np.array([(span(speed, mu_s, params + h) - span(speed, mu_s, params - h))
                                     / (2 * h.max()) for h in step])
                # Refit E2's exact mean curve for each E1 endpoint before rerunning v3.
                span_lower_e1 = span(speed, estimates['lower'], lower_e1)
                span_upper_e1 = span(speed, estimates['upper'], upper_e1)
                predictions.append(dict(true_vs_mm_s=vs, board_speed_mm_s=speed, span_um=nominal,
                                        span_se_um=np.sqrt(gradient @ covariance @ gradient),
                                        span_at_e1_lower_um=span_lower_e1,
                                        span_at_e1_upper_um=span_upper_e1))

    output.mkdir(parents=True, exist_ok=True)
    for filename, rows in [('e2_designs.csv', designs), ('e2_synthetic_recovery.csv', recovery),
                           ('e2_steady_ratios.csv', ratios), ('e2_span_uncertainty.csv', predictions)]:
        save_rows(output / filename, rows)
