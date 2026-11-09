# PHM Studio 振动特征

本文与 [PHM Studio](https://phm-studio.onrender.com/card2/) 特征工程的 19 项计算对齐。表中区分界面名称、特征 ID 和 CSV 列名。用于数据处理或跨版本比对时，优先使用特征 ID。

## 输入与结果

- 加速度和速度均使用单通道、等间隔、有限的数值序列。频率相关计算使用各通道的实际采样率。设原始加速度为 $x$，原始速度为 $v$，加速度采样率为 $f_s$，记录长度为 $N$。
- 记 $\overline{q}$ 为序列 $q$ 的算术平均值，$\sigma(q)=\sqrt{\overline{(q-\overline{q})^2}}$ 为总体标准差，$M(q)=\max_i\lvert q_i\rvert$ 为绝对峰值。
- 传统加速度统计量使用序列 $u$。当 $f_s\leq25\,600$ Hz 时，$u=x$。当 $f_s>25\,600$ Hz 时，先去均值，再做截至 12,800 Hz 的 FFT 截断重建。截断后频点不足以做逆变换时，回退到原始样本。原始 RMS、绝对峰值、峰峰值、Welch、调制和谐波特征直接使用原始数据。
- 条件不足的结果在界面显示为空、在 CSV 中留空。表中列出的零值和速度 RMS 的兼容值是现有计算约定。RMS 和峰值的单位继承输入信号，频率为 Hz，比例型特征按表中约定输出无量纲值或百分数。

## 能量特征（9 项）

设 $(f_k,P_k)$ 是原始加速度的单边 Welch 功率谱密度，$T=\sum_k P_k$，$H=\sum_{k:f_k\ge f_h} P_k$。每段长度 $L=\min(N,2048)$，采用周期 Hann 窗、每段去均值、$\lfloor L/2\rfloor$ 样本重叠，按密度标度取段均值。频率间隔为 $\Delta f=f_s/L$。高频起点 $f_h$ 默认 2,000 Hz，界面允许 1 至 100,000 Hz，并且要严格低于 Nyquist 频率 $f_s/2$ 才能得到高频结果。

| 界面名称 | 特征 ID | CSV 列名 | 计算公式 | 适用条件 |
| --- | --- | --- | --- | --- |
| 加速度 RMS（原始） | `acc_raw_rms` | 加速度 RMS（原始） | $\sqrt{\overline{x^2}}$ | 用于衡量包含直流分量的整体幅值。 |
| 加速度 RMS（去均值） | `acc_rms` | 加速度有效值 | $\sigma(u)$ | 用于衡量去均值后的波动幅值。高采样率时遵循上文的带限规则。 |
| 速度 RMS | `vel_rms` | 振动烈度 | $\sigma(v_{\mathrm{band}})$ | 需要显式的速度通道。滤波和回退规则见下文。 |
| 绝对峰值 | `acc_abs_peak` | 绝对峰值 | $M(x)$ | 用于查看原始记录中绝对幅度最大的样本，对离群点和直流偏置敏感。 |
| 峰峰值 | `acc_peak_to_peak` | 峰峰值 | $\max_i x_i-\min_i x_i$ | 用于衡量记录的幅度跨度，对离群点敏感。 |
| 主峰频率 | `acc_dominant_hz` | 主峰频率 (Hz) | $f_{\arg\max_{k>0} P_k}$ | 需要有正功率的非直流频点。输出为 Welch 频点。 |
| 谱质心 | `acc_spectral_centroid_hz` | 谱质心 (Hz) | $\sum_k f_kP_k/T$ | $T>0$。计算包含直流频点。 |
| 高频段 RMS | `acc_high_band_rms` | 高频段 RMS | $\sqrt{H\,\Delta f}$ | $f_h<f_s/2$ 且至少有一个频点落入高频段。带内功率为零时结果为零。 |
| 高频能量占比 | `acc_high_band_energy_pct` | 高频能量占比 (%) | $100H/T$ | 满足高频段 RMS 的频段条件，且 $T>0$。 |

**速度 RMS 的频带。** 设速度采样率为 $f_v$，速度记录长度为 $N_v$。下限为 $\max(10,f_v/N_v)$ Hz，上限为 $\min(1000,f_v/2-1)$ Hz。上限大于下限时，对速度执行 FFT 截断重建后取总体标准差。频带不可用、滤波失败或结果非正时，回退到原始速度的总体标准差。如果该值也非正或非有限，Studio 返回兼容值 `0.001`。这是占位值，不代表有效测量。

**传统 FFT 截断重建。** 对长度为 $n$、采样率为 $f$ 的序列做实数 FFT。保留从 $\max(1,\operatorname{round}(f_{\mathrm{low}}n/f))$ 到 $\operatorname{round}(f_{\mathrm{high}}n/f)$ 的频率下标（两端含，超出原频谱范围的下标自然忽略）。将低于起点的频点置零，截去终点之后的频点。对截断频谱使用默认长度做实数逆 FFT，并将结果乘以输出长度与原长度之比。

## 冲击特征（4 项）

| 界面名称 | 特征 ID | CSV 列名 | 计算公式 | 适用条件 |
| --- | --- | --- | --- | --- |
| 全频段峭度 | `acc_kurt_FF` | 全频段峭度 | $\overline{(u-\overline{u})^4}/\sigma(u)^4$ | $\sigma(u)>0$。采用 Pearson 峭度，正态分布的理论值为 3。 |
| 峰值因子（原始 RMS） | `acc_crest_factor_raw` | 峰值因子（原始 RMS） | $M(x)/\sqrt{\overline{x^2}}$ | 原始 RMS 大于零。峰值和 RMS 都包含直流偏置。 |
| 峰值因子（去均值 RMS） | `acc_peak_factor` | 峰值因子 | $M(u)/\sigma(u)$ | 传统口径只在分母内去均值。分母为零时返回 0。 |
| 脉冲因子 | `acc_impulse_factor` | 冲击因子 | $M(u)/\overline{\lvert u\rvert}$ | 使用平均绝对幅值作为分母。分母为零时返回 0。 |

## 调制特征（4 项）

给定载频 $f_c$、半带宽 $B$，要求 $0<f_c-B<f_c+B<f_s/2$。使用 SciPy 的四阶 Butterworth 带通设计，采用 SOS 形式做前后向零相位滤波，再用 Hilbert 变换构造解析信号 $z$。两端各舍弃 $\lceil3f_s/B\rceil$ 个样本后，至少保留 32 个样本。定义包络 $A_i=\lvert z_i\rvert$、展开相位 $\phi_i=\operatorname{unwrap}(\arg z_i)$。包络平均值须大于原始绝对峰值的 $10^{-6}$ 倍且为正。选择主要由单个载波及其边带占据的频段，便于解释瞬时频率和相位。

| 界面名称 | 特征 ID | CSV 列名 | 计算公式 | 适用条件 |
| --- | --- | --- | --- | --- |
| 包络变异系数 | `acc_envelope_cv_pct` | 包络变异系数 (%) | $100\,\sigma(A)/\overline{A}$ | 满足共同的载频、长度和包络条件。 |
| 包络调制主频 | `acc_modulation_peak_hz` | 包络调制主频 (Hz) | $f_{\arg\max_{k:0<f_k\le B}P_k^{(A)}}$，$P_k^{(A)}$ 采用上述 Welch 口径 | 包络变异系数（乘 100 前）大于 $10^{-6}$，且 $(0,B]$ 内至少有一个正功率频点。 |
| 瞬时频率标准差 | `acc_instantaneous_freq_std_hz` | 瞬时频率标准差 (Hz) | $\sigma\left(f_s(\phi_{i+1}-\phi_i)/(2\pi)\right)$ | 包络最小值须大于平均包络的 $10^{-3}$ 倍，以避开近零幅值处的相位不定。 |
| 相位残差 RMS | `acc_phase_residual_rms_rad` | 相位残差 RMS (rad) | $\sqrt{\overline{(\phi_i-\widehat\phi_i)^2}}$，$\widehat\phi_i=a+bi$ 为最小二乘直线 | 适用条件同瞬时频率标准差。 |

## 谐波特征（2 项）

给定已知基频 $f_0$。对整段原始信号去均值并加周期 Hann 窗，求单边周期图功率谱密度 $Q_k$。频率间隔为 $\delta f=f_s/N$。设 $E_\mathrm{tot}=\delta f\sum_k Q_k$，$E_j=\delta f\sum_{k:\lvert f_k-j f_0\rvert\le 2\delta f}Q_k$。要求 $f_0\ge 5\delta f$，$E_\mathrm{tot}>0$，$E_1>10^{-12}E_\mathrm{tot}$，且至少有一个 2 至 5 阶谐波的积分窗口完整处于 Nyquist 频率以内。定义 $E_h=\sum_j E_j$，只加总满足 $j f_0+2\delta f\le f_s/2$ 的 2 至 5 阶谐波。

| 界面名称 | 特征 ID | CSV 列名 | 计算公式 | 适用条件 |
| --- | --- | --- | --- | --- |
| 谐波能量占比 | `acc_harmonic_energy_pct` | 谐波能量占比 (%) | $100E_h/E_\mathrm{tot}$ | 满足上文基频分辨率和带宽条件。适合基频稳定且已知的记录。 |
| 总谐波失真 THD | `acc_thd_pct` | 总谐波失真 THD (%) | $100\sqrt{E_h/E_1}$ | 取可解析的 2 至 5 阶谐波与基频的 RMS 比。 |

载频、半带宽和基频在 Studio 界面的输入范围为 0.01 至 100,000 Hz。与采样率、记录长度有关的条件仍须满足。跨文件比较时，保持采样率、窗长和频带设置一致，并记录设备工况。

## 运行示例

在仓库根目录执行。示例生成两秒加速度和速度数据，无需下载数据集。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r examples/requirements.txt
python -m examples.studio_features --csv studio-features.csv
```

终端打印 19 个特征的 ID、界面名称和数值。CSV 使用 Studio 的特征列名。这个合成信号的主峰约为 800 Hz，包络调制主频为 80 Hz，THD 约为 10%。[示例源码](../examples/studio_features.py)

处理自己的数组时，可以调用 `calculate_features(acceleration, fs, velocity=velocity, velocity_fs=velocity_fs, ...)`。省略速度通道时，`vel_rms` 返回 `None`。缺少载频、半带宽或基频时，相应的参数化特征返回 `None`。调用示例中同一测量记录只传入一条加速度和一条可选速度信号。
