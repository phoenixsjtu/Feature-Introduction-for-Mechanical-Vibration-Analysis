"""Spectral measurements for uniformly sampled vibration signals."""

import numpy as np
from scipy.signal import welch


DEFAULT_HIGH_FREQUENCY_START_HZ = 2000
DEFAULT_WELCH_WINDOW = 2048
SPECTRAL_FEATURE_LABELS = {
    "acc_dominant_hz": "主峰频率 (Hz)",
    "acc_spectral_centroid_hz": "谱质心 (Hz)",
    "acc_high_band_rms": "高频段 RMS",
    "acc_high_band_energy_pct": "高频能量占比 (%)",
}
HIGH_BAND_FEATURE_IDS = {"acc_high_band_rms", "acc_high_band_energy_pct"}


def validated_signal(samples, fs):
    signal = np.asarray(samples, dtype=float)
    if signal.ndim != 1 or signal.size < 2 or not np.isfinite(signal).all():
        raise ValueError("请使用至少含两个有限数值的单通道信号。")
    if not np.isfinite(fs) or fs <= 0:
        raise ValueError("采样率必须是有限的正数。")
    return signal


def welch_density(samples, fs, window_size=DEFAULT_WELCH_WINDOW):
    signal = validated_signal(samples, fs)
    segment_size = min(signal.size, max(2, int(window_size)))
    frequencies, density = welch(
        signal,
        fs=fs,
        window="hann",
        nperseg=segment_size,
        noverlap=segment_size // 2,
        detrend="constant",
        scaling="density",
    )
    if not np.isfinite(density).all():
        raise ValueError("信号幅值超出功率谱的计算范围。")
    return frequencies, density, segment_size


def spectral_features(samples, fs, high_frequency_start):
    values = dict.fromkeys(SPECTRAL_FEATURE_LABELS)
    try:
        frequencies, density, segment_size = welch_density(samples, fs)
    except ValueError:
        return values

    # The one-sided PSD bin sum times its width estimates mean square.
    total_density = float(density.sum())
    if total_density > 0:
        if np.any(density[1:] > 0):
            values["acc_dominant_hz"] = float(frequencies[1 + np.argmax(density[1:])])
        values["acc_spectral_centroid_hz"] = float(frequencies @ density / total_density)

    high_bins = frequencies >= high_frequency_start
    if high_frequency_start < fs / 2 and high_bins.any():
        high_density = float(density[high_bins].sum())
        values["acc_high_band_rms"] = float(np.sqrt(high_density * (fs / segment_size)))
        if total_density > 0:
            values["acc_high_band_energy_pct"] = 100 * high_density / total_density
    return values
