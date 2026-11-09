"""Run the PHM Studio 19-feature example with `python -m examples.studio_features`."""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.fft import irfft, rfft

from .vibration_features import CATALOG_FEATURE_LABELS, catalogue_features
from .vibration_spectrum import (
    DEFAULT_HIGH_FREQUENCY_START_HZ,
    spectral_features,
    validated_signal,
)


# Keep catalogue order, interface labels and CSV headers aligned with Studio.
FEATURES = (
    ("acc_raw_rms", "加速度 RMS（原始）", "加速度 RMS（原始）"),
    ("acc_rms", "加速度 RMS（去均值）", "加速度有效值"),
    ("vel_rms", "速度 RMS", "振动烈度"),
    ("acc_abs_peak", "绝对峰值", "绝对峰值"),
    ("acc_peak_to_peak", "峰峰值", "峰峰值"),
    ("acc_dominant_hz", "主峰频率", "主峰频率 (Hz)"),
    ("acc_spectral_centroid_hz", "谱质心", "谱质心 (Hz)"),
    ("acc_high_band_rms", "高频段 RMS", "高频段 RMS"),
    ("acc_high_band_energy_pct", "高频能量占比", "高频能量占比 (%)"),
    ("acc_kurt_FF", "全频段峭度", "全频段峭度"),
    ("acc_crest_factor_raw", "峰值因子（原始 RMS）", "峰值因子（原始 RMS）"),
    ("acc_peak_factor", "峰值因子（去均值 RMS）", "峰值因子"),
    ("acc_impulse_factor", "脉冲因子", "冲击因子"),
    ("acc_envelope_cv_pct", "包络变异系数", "包络变异系数 (%)"),
    ("acc_modulation_peak_hz", "包络调制主频", "包络调制主频 (Hz)"),
    ("acc_instantaneous_freq_std_hz", "瞬时频率标准差", "瞬时频率标准差 (Hz)"),
    ("acc_phase_residual_rms_rad", "相位残差 RMS", "相位残差 RMS (rad)"),
    ("acc_harmonic_energy_pct", "谐波能量占比", "谐波能量占比 (%)"),
    ("acc_thd_pct", "总谐波失真 THD", "总谐波失真 THD (%)"),
)


def _studio_fft_filter(samples, fs, low, high):
    """Match Studio's legacy FFT mask, truncated inverse and length scaling."""
    spectrum = rfft(samples)
    n = samples.size
    first_bin = max(1, int(np.round(low / fs * n)))
    end_bin = int(np.round(high / fs * n)) + 1
    spectrum[:first_bin] = 0
    filtered = irfft(spectrum[:end_bin])
    return filtered * filtered.size / n


def _legacy_acceleration(samples, fs):
    signal = samples
    if fs > 25_600:
        try:
            signal = _studio_fft_filter(samples - np.mean(samples), fs, 0, 12_800)
        except ValueError:
            # A very short record can leave too few bins for the inverse FFT.
            pass

    rms = np.std(signal)
    peak = np.max(np.abs(signal))
    abs_mean = np.mean(np.abs(signal))
    kurtosis = np.sum((signal - np.mean(signal)) ** 4) / rms**4 / signal.size if rms else None
    return {
        "acc_rms": rms,
        "acc_kurt_FF": kurtosis,
        "acc_peak_factor": peak / rms if rms else 0,
        "acc_impulse_factor": peak / abs_mean if abs_mean else 0,
    }


def _velocity_rms(velocity, fs):
    low = max(10, fs / velocity.size)
    high = min(1000, fs / 2 - 1)
    filtered = _studio_fft_filter(velocity, fs, low, high) if high > low else velocity
    rms = np.std(filtered)
    if np.isfinite(rms) and rms > 0:
        return rms

    rms = np.std(velocity)
    return rms if np.isfinite(rms) and rms > 0 else 0.001


def calculate_features(
    acceleration,
    fs,
    *,
    velocity=None,
    velocity_fs=None,
    high_frequency_start_hz=DEFAULT_HIGH_FREQUENCY_START_HZ,
    modulation_carrier_hz=None,
    modulation_half_band_hz=None,
    harmonic_fundamental_hz=None,
):
    """Return a feature-ID mapping for finite, one-dimensional input channels.

    Leave a parameter unset to skip its parameterized group. Supply frequency
    settings in Hz, following the ranges and signal conditions in the guide.
    """
    samples = validated_signal(acceleration, fs)
    params = {
        "modulation_carrier_hz": modulation_carrier_hz,
        "modulation_half_band_hz": modulation_half_band_hz,
        "harmonic_fundamental_hz": harmonic_fundamental_hz,
    }
    values = catalogue_features(samples, fs, set(CATALOG_FEATURE_LABELS), params)
    values.update(spectral_features(samples, fs, high_frequency_start_hz))
    values.update(_legacy_acceleration(samples, fs))
    if velocity is not None:
        speed_fs = fs if velocity_fs is None else velocity_fs
        values["vel_rms"] = _velocity_rms(validated_signal(velocity, speed_fs), speed_fs)

    ordered = {feature_id: values.get(feature_id) for feature_id, _, _ in FEATURES}
    return {key: float(value) if value is not None and np.isfinite(value) else None
            for key, value in ordered.items()}


def demo_features():
    """Synthetic two-channel measurement with known spectral and modulation components."""
    fs = 10_240
    t = np.arange(2 * fs) / fs
    acceleration = (
        0.4
        + (1 + 0.3 * np.cos(2 * np.pi * 80 * t)) * np.sin(2 * np.pi * 800 * t)
        + 0.1 * np.sin(2 * np.pi * 1_600 * t)
        + 0.05 * np.sin(2 * np.pi * 2_700 * t)
    )
    velocity = 0.8 * np.sin(2 * np.pi * 250 * t)
    return calculate_features(
        acceleration,
        fs,
        velocity=velocity,
        high_frequency_start_hz=2_000,
        modulation_carrier_hz=800,
        modulation_half_band_hz=200,
        harmonic_fundamental_hz=800,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="also write one row with the Studio feature headers")
    args = parser.parse_args()
    values = demo_features()
    for feature_id, label, _ in FEATURES:
        value = values[feature_id]
        display = f"{value:.8g}" if value is not None else "—"
        print(f"{feature_id:34} {label}  {display}")

    if args.csv:
        with args.csv.open("w", encoding="utf-8-sig", newline="") as output:
            writer = csv.writer(output)
            writer.writerow([header for _, _, header in FEATURES])
            writer.writerow([values[feature_id] for feature_id, _, _ in FEATURES])


if __name__ == "__main__":
    main()
