"""Scalar features for the vibration feature catalogue."""

import numpy as np
from scipy.signal import butter, hilbert, periodogram, sosfiltfilt

from .vibration_spectrum import validated_signal, welch_density


TIME_FEATURE_LABELS = {
    "acc_raw_rms": "加速度 RMS（原始）",
    "acc_abs_peak": "绝对峰值",
    "acc_peak_to_peak": "峰峰值",
    "acc_crest_factor_raw": "峰值因子（原始 RMS）",
}
MODULATION_FEATURE_LABELS = {
    "acc_envelope_cv_pct": "包络变异系数 (%)",
    "acc_modulation_peak_hz": "包络调制主频 (Hz)",
    "acc_instantaneous_freq_std_hz": "瞬时频率标准差 (Hz)",
    "acc_phase_residual_rms_rad": "相位残差 RMS (rad)",
}
HARMONIC_FEATURE_LABELS = {
    "acc_harmonic_energy_pct": "谐波能量占比 (%)",
    "acc_thd_pct": "总谐波失真 THD (%)",
}
CATALOG_FEATURE_LABELS = TIME_FEATURE_LABELS | MODULATION_FEATURE_LABELS | HARMONIC_FEATURE_LABELS


def _time_features(signal):
    peak = float(np.max(np.abs(signal)))
    rms = peak * float(np.sqrt(np.mean((signal / peak) ** 2))) if peak else 0.0
    return {
        "acc_raw_rms": rms,
        "acc_abs_peak": peak,
        "acc_peak_to_peak": float(signal.max()) - float(signal.min()),
        "acc_crest_factor_raw": peak / rms if rms else None,
    }


def _modulation_features(signal, fs, carrier, half_bandwidth):
    if not 0 < carrier - half_bandwidth < carrier + half_bandwidth < fs / 2:
        return {}

    # Leave three inverse-bandwidth intervals at each end for filter transients.
    trim = int(np.ceil(3 * fs / half_bandwidth))
    if signal.size - 2 * trim < 32:
        return {}

    sos = butter(4, [carrier - half_bandwidth, carrier + half_bandwidth], "bandpass", fs=fs, output="sos")
    analytic = hilbert(sosfiltfilt(sos, signal))[trim:-trim]
    envelope = np.abs(analytic)
    mean_envelope = float(np.mean(envelope))
    if mean_envelope <= np.max(np.abs(signal)) * 1e-6 or mean_envelope == 0:
        return {}

    freq, density, _ = welch_density(envelope, fs)
    envelope_cv = float(np.std(envelope)) / mean_envelope
    # Modulation is resolved within the retained carrier sideband bandwidth.
    sideband = (freq > 0) & (freq <= half_bandwidth)
    peak_frequency = None
    if envelope_cv > 1e-6 and sideband.any() and np.any(density[sideband] > 0):
        peak_frequency = float(freq[sideband][np.argmax(density[sideband])])

    # Phase is undefined when the analytic amplitude approaches zero.
    freq_std = phase_rms = None
    if float(envelope.min()) > mean_envelope * 1e-3:
        phase = np.unwrap(np.angle(analytic))
        x = np.arange(phase.size) - (phase.size - 1) / 2
        centered_phase = phase - phase.mean()
        residual = centered_phase - x * ((x @ centered_phase) / (x @ x))
        freq_std = float(np.std(np.diff(phase) * fs / (2 * np.pi)))
        phase_rms = float(np.sqrt(np.mean(residual**2)))

    return {
        "acc_envelope_cv_pct": 100 * envelope_cv,
        "acc_modulation_peak_hz": peak_frequency,
        "acc_instantaneous_freq_std_hz": freq_std,
        "acc_phase_residual_rms_rad": phase_rms,
    }


def _harmonic_features(signal, fs, fundamental):
    resolution = fs / signal.size
    half_width = 2 * resolution
    harmonics = [order * fundamental for order in range(2, 6) if order * fundamental + half_width <= fs / 2]
    if fundamental < 5 * resolution or not harmonics:
        return {}

    freq, density = periodogram(signal, fs=fs, window="hann", detrend="constant", scaling="density")
    total = float(density.sum())
    base = float(density[np.abs(freq - fundamental) <= half_width].sum())
    if total <= 0 or base <= total * 1e-12:
        return {}

    harmonic = sum(float(density[np.abs(freq - center) <= half_width].sum()) for center in harmonics)
    return {
        "acc_harmonic_energy_pct": 100 * harmonic / total,
        "acc_thd_pct": 100 * np.sqrt(harmonic / base),
    }


def catalogue_features(samples, fs, selected, params):
    values = dict.fromkeys(selected.intersection(CATALOG_FEATURE_LABELS))
    try:
        signal = validated_signal(np.atleast_1d(np.squeeze(samples)), fs)
    except ValueError:
        return values

    measured = {}
    if values.keys() & TIME_FEATURE_LABELS.keys():
        measured.update(_time_features(signal))
    if values.keys() & MODULATION_FEATURE_LABELS.keys():
        carrier = params.get("modulation_carrier_hz")
        half_bandwidth = params.get("modulation_half_band_hz")
        if carrier is not None and half_bandwidth is not None:
            measured.update(_modulation_features(signal, fs, carrier, half_bandwidth))
    if values.keys() & HARMONIC_FEATURE_LABELS.keys():
        fundamental = params.get("harmonic_fundamental_hz")
        if fundamental is not None:
            measured.update(_harmonic_features(signal, fs, fundamental))

    for key in values:
        value = measured.get(key)
        if value is not None and np.isfinite(value):
            values[key] = float(value)
    return values
