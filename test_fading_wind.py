#!/usr/bin/env python3
"""Synthetic recovery test for the strict full-correlation estimator."""

from __future__ import annotations

import unittest

import numpy as np

import fading_wind


class FadingWindTest(unittest.TestCase):
    def test_recovers_advected_gaussian_pattern(self) -> None:
        generator = np.random.default_rng(4)
        sample_interval_s = 0.25
        time_s = np.arange(0.0, 15.0 * 60.0, sample_interval_s)
        mode_count = 1000
        wavevector = generator.normal(size=(mode_count, 2)) / 85.0
        phase = generator.uniform(0.0, 2.0 * np.pi, mode_count)
        amplitude = generator.normal(size=mode_count) / np.sqrt(mode_count)
        pattern_velocity_ms = np.asarray((60.0, 40.0))
        channels = []
        for antenna_position in fading_wind.ANTENNA_EN_M:
            argument = (
                (
                    antenna_position[None, :]
                    - time_s[:, None] * pattern_velocity_ms[None, :]
                )
                @ wavevector.T
                + phase
            )
            gaussian_field = np.sqrt(2.0) * np.cos(argument) @ amplitude
            channels.append(np.exp(0.25 * gaussian_field))
        result = fading_wind.fit_robust_fading_window(
            np.stack(channels, axis=1),
            sample_interval_s,
        )
        self.assertTrue(result["valid"])
        np.testing.assert_allclose(
            (
                result["zonal_wind_ms"],
                result["meridional_wind_ms"],
            ),
            0.5 * pattern_velocity_ms,
            atol=8.0,
        )


if __name__ == "__main__":
    unittest.main()
