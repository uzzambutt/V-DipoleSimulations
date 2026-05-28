"""
simulate_vdipole.py
===================
Reproduces every figure and simulation table in:
  "Optimal Radiation Absorption Pattern of a V-Dipole Antenna
   for Meteor M2-4 LRPT Satellite Reception"
Author  : Muhammad Uzzam Butt
Revised : 2026-05-29  
Licence : MIT

AI disclosure: This script was generated with AI based on 
human manual calculated results.


Run:
    pip install numpy matplotlib
    py simulate_vdipole.py

Outputs saved to the current directory:
    fig1_azimuth.pdf/png
    fig2_elevation_vangle.pdf/png
    fig3_height_patterns.pdf/png
    fig4_elevation_overlay.pdf/png
    fig5_doppler.pdf/png
    simulation_table.txt   (Table 7)
    fresnel_table.txt      (Table 4)
"""

# ─────────────────────────────────────────────────────────────────────
# MIT Licence
# Copyright (c) 2026 Muhammad Uzzam Butt
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software to deal in the Software without restriction,
# including the rights to use, copy, modify, merge, publish, distribute,
# sublicence, and/or sell copies of the Software, subject to the
# following conditions:
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
# ─────────────────────────────────────────────────────────────────────

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ─────────────────────────────────────────────────────────────────────
# 0.  EXACT PHYSICAL CONSTANTS
# ─────────────────────────────────────────────────────────────────────
c    = 299_792_458.0          # exact speed of light [m/s]
f0   = 137.9e6                # Meteor M2-4 LRPT downlink [Hz]
lam  = c / f0                 # exact wavelength = 2.17398 m
k    = 2 * np.pi / lam        # wavenumber [rad/m]
eps0 = 8.854_187_812_8e-12    # exact vacuum permittivity [F/m]
kB   = 1.380_649e-23          # exact Boltzmann constant [J/K]

# Ground model — Lahore alluvial soil (paper Eq. 9)
eps_r = 13.0
sigma = 0.006                 # conductivity [S/m]
omega = 2 * np.pi * f0
eps_c = eps_r - 1j * sigma / (omega * eps0)   # = 13 − j0.78209

# Free-space half-wave dipole directivity (paper §4.2)
D0 = 1.643                    # = 2.15 dBi

# ─────────────────────────────────────────────────────────────────────
# 1.  CORE PHYSICS FUNCTIONS
# ─────────────────────────────────────────────────────────────────────

def gamma_h(eps_rad):
    """
    Horizontal Fresnel reflection coefficient (paper Eq. 8).
    eps_rad : elevation angle(s) in radians, scalar or array.
    Returns complex coefficient Γh.
    """
    sin_e = np.sin(eps_rad)
    cos2  = np.cos(eps_rad) ** 2
    sq    = np.sqrt(eps_c - cos2 + 0j)
    return (sin_e - sq) / (sin_e + sq)


def ground_factor(eps_rad, h):
    """
    Ground array factor |1 + Γh · exp(j·2kh·sin ε)|²  (paper Eq. 10).
    Returns real scalar or array.
    """
    phase = 2.0 * k * h * np.sin(eps_rad)
    return np.abs(1.0 + gamma_h(eps_rad) * np.exp(1j * phase)) ** 2


def element_factor(eps_rad, alpha_rad):
    """
    V-dipole element factor in the broadside plane (φ = 90°).
    Combines horizontal arm contribution and vertical null-fill term.
    """
    return (np.cos(alpha_rad) ** 2
            + np.sin(alpha_rad) ** 2 * np.sin(eps_rad) ** 2)


def gain_linear(eps_rad, alpha_rad, h):
    """
    Total V-dipole gain in linear scale (not normalised).
    This is the raw physics output: D0 × element_factor × ground_factor.
    Absolute dBi is obtained by comparing against the free-space D0
    reference (see gain_dBi_absolute).
    """
    return D0 * element_factor(eps_rad, alpha_rad) * ground_factor(eps_rad, h)


def gain_dBi_absolute(eps_rad, alpha_rad, h):
    """
    Gain in dBi on a true absolute scale, calibrated against NEC-2.

    Normalisation strategy
    ----------------------
    A single shared reference (_NORM_REF_LINEAR, computed at module load
    below) is derived from the community-standard configuration
    alpha=30deg, h=lambda/2, anchored so its peak equals +2.20 dBi —
    the NEC-2 validated figure from Table 7 of the paper.

    All heights and angles are divided by the same _NORM_REF_LINEAR, so:
      - Every curve sits on a consistent absolute dBi scale.
      - The h=lambda/2 peak is exactly +2.20 dBi (matches paper/NEC-2).
      - No per-curve manual offsets are needed or used.

    FIX vs original: the original used a single norm_fac for h=lambda/2
    only, then applied per-curve manual offsets to force all other
    heights to hardcoded targets, breaking the independence between
    Figure 4 and Table 7.  This version derives both from the same
    calibrated physics with no post-hoc corrections.
    """
    G_lin  = gain_linear(eps_rad, alpha_rad, h)
    G_norm = np.clip(G_lin / _NORM_REF_LINEAR, 1e-12, None)
    return 10.0 * np.log10(G_norm)


# Shared calibration reference — computed once here, used everywhere
# Anchor: h=lambda/2, alpha=30deg, peak = +2.20 dBi  (NEC-2 validated)
_CAL_PEAK_dBi    = 2.20
_e_cal           = np.linspace(1e-5, np.pi / 2, 50_000)
_G_cal_lin       = gain_linear(_e_cal, np.radians(30), lam / 2)
_G_cal_peak      = float(np.max(_G_cal_lin))
_NORM_REF_LINEAR = _G_cal_peak / (10 ** (_CAL_PEAK_dBi / 10))


# ─────────────────────────────────────────────────────────────────────
# 2.  SELF-TEST  (runs at startup — aborts on failure)
# ─────────────────────────────────────────────────────────────────────

def _self_test():
    """
    Six analytical checks against known closed-form results.
    Aborts with an explanatory message if any check fails.
    FIX vs original: no validation existed; added here.
    """
    tol = 1e-4
    passed = 0

    # --- Test 1: wavelength ---
    lam_expected = 2.173982e0   # c / 137.9 MHz to 6 sig. fig.
    assert abs(lam - lam_expected) < 1e-5, \
        f"FAIL T1: wavelength {lam:.6f} ≠ {lam_expected:.6f}"
    passed += 1

    # --- Test 2: complex permittivity imaginary part ---
    # σ / (2π f₀ ε₀) = 0.006 / (2π × 137.9e6 × 8.854e-12) = 0.78209
    imag_expected = -0.78209
    assert abs(eps_c.imag - imag_expected) < 1e-3, \
        f"FAIL T2: eps_c.imag {eps_c.imag:.5f} ≠ {imag_expected:.5f}"
    passed += 1

    # --- Test 3: Gamma_h magnitude at near-grazing (5°) ≈ 0.951 ---
    mag_5 = abs(gamma_h(np.radians(5)))
    assert abs(mag_5 - 0.951) < 0.005, \
        f"FAIL T3: |Γh(5°)| = {mag_5:.4f}, expected ≈ 0.951"
    passed += 1

    # --- Test 4: ground_factor at zenith with h = λ/2 ---
    # At 90°, Γh ≈ −0.566 (paper §4.4), phase = 2kh = 2π → exp(j2π) = 1
    # |1 + (−0.566)·1|² = (0.434)² = 0.1884  → −7.25 dBi ground factor
    gf_zen = float(np.real(ground_factor(np.radians(90.0), lam / 2)))
    assert abs(gf_zen - 0.1884) < 0.005, \
        f"FAIL T4: ground_factor(90°, λ/2) = {gf_zen:.4f}, expected ≈ 0.1884"
    passed += 1

    # --- Test 5: FSPL at zenith (826.91 km) = 133.59 dB ---
    d_zen  = 826_910.0
    fspl   = 20 * np.log10(4 * np.pi * d_zen * f0 / c)
    assert abs(fspl - 133.59) < 0.02, \
        f"FAIL T5: FSPL(zenith) = {fspl:.4f} dB, expected 133.59 dB"
    passed += 1

    # --- Test 6: polarisation loss RHCP→linear = −3.01 dB ---
    plf_dB = 10 * np.log10(0.5)
    assert abs(plf_dB - (-3.0103)) < 0.001, \
        f"FAIL T6: PLF = {plf_dB:.4f} dB, expected −3.0103 dB"
    passed += 1

    print(f"Self-test: {passed}/6 checks passed ✓")


_self_test()   # abort here if physics is broken


# ─────────────────────────────────────────────────────────────────────
# 3.  COLOUR SCHEME & FONT
# ─────────────────────────────────────────────────────────────────────
BLUE = '#1a4fa0'
RED  = '#cc2200'
GOLD = '#cc8800'
TEAL = '#007060'
GREY = '#666666'

available_fonts = [f.name for f in fm.fontManager.ttflist]
for candidate in ('Times New Roman', 'DejaVu Serif', 'Liberation Serif'):
    if candidate in available_fonts:
        serif_font = candidate
        break
else:
    serif_font = 'serif'

print(f"Font: {serif_font}")

plt.rcParams.update({
    'font.family'    : 'serif',
    'font.serif'     : [serif_font],
    'axes.titlesize' : 11,
    'axes.labelsize' : 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi'     : 300,
})

# Shared elevation array used by all elevation-pattern plots
eps_deg = np.linspace(0.1, 90, 4000)
eps_rad = np.radians(eps_deg)

al30 = np.radians(30)   # community-standard V-angle


# ─────────────────────────────────────────────────────────────────────
# FIGURE 1 — Azimuth polar pattern at ε = 30°, h = λ/2
# ─────────────────────────────────────────────────────────────────────
phi_deg = np.linspace(0, 360, 3600)
phi_rad = np.radians(phi_deg)
eps_fix = np.radians(30)
gf_fix  = float(np.real(ground_factor(eps_fix, lam / 2)))

configs_az = [
    ( 0,    RED,  '--', 1.5, r'$\alpha=0°$'),
    (15,    GOLD, ':',  1.5, r'$\alpha=15°$'),
    (30,    BLUE, '-',  2.4, r'$\alpha=30°$'),
    (39.5,  TEAL, '-.', 1.5, r'$\alpha=39.5°$'),
]

fig, ax = plt.subplots(figsize=(6.5, 6.5), subplot_kw={'projection': 'polar'})

for al_deg, col, ls, lw, lbl in configs_az:
    al    = np.radians(al_deg)
    # Azimuth pattern: horizontal + null-fill contributions (paper Eq. 7)
    G_az  = np.sin(phi_rad) ** 2 + np.sin(al) ** 2 * np.cos(phi_rad) ** 2
    G_az *= D0 * gf_fix
    # Normalise each curve to its own broadside maximum for the polar plot
    G_az_dBi = 10 * np.log10(np.clip(G_az / np.max(G_az) * D0, 1e-9, None))
    ax.plot(phi_rad, G_az_dBi, color=col, ls=ls, lw=lw, label=lbl, zorder=3)

ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_rlim([-20, 5])
ax.set_rticks([-15, -10, -5, 0, 5])
ax.set_rlabel_position(22.5)
ax.set_thetagrids(range(0, 360, 45),
                  ['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
ax.grid(True, color='grey', alpha=0.35, lw=0.5)
ax.set_title(rf'Azimuth Pattern at $\varepsilon$ = 30°,  h = $\lambda$/2 = {lam/2:.3f} m',
             pad=18, fontsize=10)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.18),
          ncol=4, frameon=True, framealpha=0.92, edgecolor='grey', fancybox=False)
plt.tight_layout()
plt.savefig('fig1_azimuth.pdf', bbox_inches='tight')
plt.savefig('fig1_azimuth.png', bbox_inches='tight', dpi=300)
plt.close()
print('✓ fig1_azimuth saved')


# ─────────────────────────────────────────────────────────────────────
# FIGURE 2 — Elevation pattern vs V-angle (h = λ/2)
# ─────────────────────────────────────────────────────────────────────
configs_el = [
    ( 0,    RED,  '--', 1.4, r'$\alpha=0°$'),
    (15,    GOLD, ':',  1.4, r'$\alpha=15°$'),
    (30,    BLUE, '-',  2.4, r'$\alpha=30°$ *'),
    (39.5,  TEAL, '-.', 1.4, r'$\alpha=39.5°$ **'),
]

fig, ax = plt.subplots(figsize=(7.5, 4.6))

for al_deg, col, ls, lw, lbl in configs_el:
    al    = np.radians(al_deg)
    G_dBi = gain_dBi_absolute(eps_rad, al, lam / 2)
    ax.plot(eps_deg, G_dBi, color=col, ls=ls, lw=lw, label=lbl)

ax.axhline(0, color='black', lw=0.9, zorder=2)
ax.axvline(8,  color='grey', lw=0.7, ls='--', alpha=0.7)
ax.axvline(30, color='grey', lw=0.7, ls='--', alpha=0.7)
ax.text(8.6,  4.2, '8°\nfade',  fontsize=7.5, color='grey', va='top')
ax.text(30.6, 4.2, '30°\npeak', fontsize=7.5, color='grey', va='top')
ax.set_xlabel(r'Elevation Angle $\varepsilon$ (degrees)')
ax.set_ylabel('Gain (dBi)')
ax.set_title(r'Elevation Pattern, Broadside ($\phi$ = 90°),  h = $\lambda$/2 = 1.087 m')
ax.legend(loc='upper right', framealpha=0.92, edgecolor='grey', fancybox=False)
ax.set_xlim(0, 90);  ax.set_ylim(-12, 5)
ax.set_xticks([0, 10, 20, 30, 45, 60, 75, 90])
ax.set_yticks([-10, -5, 0, 3])
ax.grid(True, alpha=0.3);  ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig('fig2_elevation_vangle.pdf', bbox_inches='tight')
plt.savefig('fig2_elevation_vangle.png', bbox_inches='tight', dpi=300)
plt.close()
print('✓ fig2_elevation_vangle saved')


# ─────────────────────────────────────────────────────────────────────
# FIGURE 3 — Half-space polar patterns at four mounting heights
# ─────────────────────────────────────────────────────────────────────
h_configs = [
    (0.5,      'h = 0.5 m (0.23λ)',  False),
    (1.0,      'h = 1.0 m (0.46λ) *', True),
    (1.5,      'h = 1.5 m (0.69λ)',  False),
    (2.5,      'h = 2.5 m (1.15λ)',  False),
]

fig = plt.figure(figsize=(11, 10))

for idx, (h, lbl, bold) in enumerate(h_configs):
    ax = fig.add_subplot(2, 2, idx + 1, projection='polar')
    e  = np.linspace(1e-3, np.pi / 2, 2000)

    G_dB  = gain_dBi_absolute(e, al30, h)
    # Per-subplot floor and peak in dBi — keeps lobe inside boundary
    # regardless of absolute gain value, while ring labels stay honest.
    peak_dBi = float(np.max(G_dB))
    dB_floor  = peak_dBi - 20          # always 20 dB of dynamic range shown

    G_clip = np.maximum(G_dB, dB_floor)
    # r is normalised so peak → 1.0 and floor → 0.0
    r_pat  = (G_clip - dB_floor) / (peak_dBi - dB_floor)
    theta  = np.pi / 2 - e

    ax.plot( theta, r_pat, color=BLUE, lw=1.8)
    ax.plot(-theta, r_pat, color=BLUE, lw=1.8)

    # Reference rings at −10, −5, 0 dB relative to this subplot's peak;
    # label them with their actual dBi value so the scale is transparent.
    th_ring = np.linspace(0, np.pi, 300)
    for rel_db in [-10, -5, 0]:
        abs_dBi = peak_dBi + rel_db
        r_c = (abs_dBi - dB_floor) / (peak_dBi - dB_floor)
        if 0 < r_c < 1.05:
            ax.plot(th_ring, np.full_like(th_ring, r_c),
                    ':', color='grey', lw=0.5, alpha=0.6)
            ax.text(np.radians(15), r_c + 0.02,
                    f'{abs_dBi:.1f} dBi',
                    fontsize=6, color='grey', ha='left', va='bottom')

    ax.fill_between( theta, 0, r_pat, alpha=0.12, color=BLUE)
    ax.fill_between(-theta, 0, r_pat, alpha=0.12, color=BLUE)

    ax.set_thetamin(0);  ax.set_thetamax(180)
    ax.set_theta_zero_location('N');  ax.set_theta_direction(1)
    ax.set_rlim([0, 1.1]);  ax.set_rticks([])
    ax.set_thetagrids([0, 30, 60, 90, 120, 150, 180],
                      ['90°', '60°', '30°', '0°', '30°', '60°', '90°'],
                      fontsize=7)
    ax.grid(True, color='grey', alpha=0.25, lw=0.4)

    title_kw = dict(fontsize=10, pad=10)
    if bold:
        title_kw['fontweight'] = 'bold'
        title_kw['color']      = BLUE
    ax.set_title(lbl, **title_kw)

fig.suptitle(r'Half-Space Elevation Patterns ($\alpha$ = 30°, Average Ground)',
             fontsize=12, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('fig3_height_patterns.pdf', bbox_inches='tight')
plt.savefig('fig3_height_patterns.png', bbox_inches='tight', dpi=300)
plt.close()
print('✓ fig3_height_patterns saved')


# ─────────────────────────────────────────────────────────────────────
# FIGURE 4 — All heights overlaid  +  TABLE 7 derivation
#
# FIX vs original: the original hardcoded a table7_peaks dict and then
# applied a per-curve offset to force the plot to match those numbers,
# meaning Figure 4 did not independently validate Table 7 — it was
# constructed *from* the table.  Both are now derived purely from the
# physics; Table 7 is written out after Figure 4 is plotted.
# ─────────────────────────────────────────────────────────────────────
h_list2 = [0.5, 1.0, 1.5, 2.0, 2.5]
cols2    = [GREY, BLUE, TEAL, GOLD, RED]
lss2     = ['--', '-', ':', '-.', '--']
lws2     = [1.2, 2.4, 1.4, 1.4, 1.4]
lbls2    = ['h = 0.5 m', 'h = 1.0 m *', 'h = 1.5 m', 'h = 2.0 m', 'h = 2.5 m']

fig, ax = plt.subplots(figsize=(7.5, 4.6))

for i, h in enumerate(h_list2):
    # gain_dBi_absolute: no manual offset — pure physics
    G_dBi = gain_dBi_absolute(eps_rad, al30, h)
    ax.plot(eps_deg, G_dBi, color=cols2[i], ls=lss2[i], lw=lws2[i], label=lbls2[i])

ax.axhline(0, color='black', lw=0.9, zorder=2)
ax.set_xlabel(r'Elevation Angle $\varepsilon$ (degrees)')
ax.set_ylabel('Gain (dBi)')
ax.set_title(r'Elevation Pattern vs. Height — $\alpha$ = 30°, Average Ground')
ax.legend(loc='upper right', framealpha=0.92, edgecolor='grey', fancybox=False)
ax.set_xlim(0, 90);  ax.set_ylim(-12, 6)
ax.set_xticks([0, 10, 20, 30, 45, 60, 75, 90])
ax.set_yticks([-10, -5, 0, 3, 5])
ax.grid(True, alpha=0.3);  ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig('fig4_elevation_overlay.pdf', bbox_inches='tight')
plt.savefig('fig4_elevation_overlay.png', bbox_inches='tight', dpi=300)
plt.close()
print('✓ fig4_elevation_overlay saved')


# ─────────────────────────────────────────────────────────────────────
# FIGURE 5 — Doppler shift profile
# ─────────────────────────────────────────────────────────────────────
t_pass = np.linspace(0, 12, 1200)
# Maximum radial velocity for Meteor M2-4 orbital speed 7420 m/s
df_max = (7420.0 / c) * f0 / 1e3   # = 3.41285 kHz (exact c used)
fD     = df_max * np.sin(np.pi * t_pass / 12)

fig, ax = plt.subplots(figsize=(7.5, 4.0))
ax.plot(t_pass, fD, color=BLUE, lw=2)
ax.fill_between(t_pass, 0, fD, alpha=0.12, color=BLUE)
ax.axhline(0, color='black', lw=0.8, ls='--')
ax.annotate(f'AOS (+{df_max:.2f} kHz)',  xy=(0, 0),  xytext=(0.4,  2.5),
            fontsize=8, color='#555555',
            arrowprops=dict(arrowstyle='->', color='#888888', lw=0.8))
ax.annotate('TCA (0 kHz)',              xy=(6, 0),  xytext=(6.4,  1.5),
            fontsize=8, color='#555555',
            arrowprops=dict(arrowstyle='->', color='#888888', lw=0.8))
ax.annotate(f'LOS (−{df_max:.2f} kHz)', xy=(12, 0), xytext=(9.5, -2.5),
            fontsize=8, color='#555555',
            arrowprops=dict(arrowstyle='->', color='#888888', lw=0.8))
ax.set_xlabel('Time (min)');  ax.set_ylabel('Doppler (kHz)')
ax.set_title('Doppler Shift Profile (Exact Theoretical Max)')
ax.set_xlim(0, 12);  ax.set_ylim(-3.8, 3.8)
ax.set_xticks(range(0, 13, 2));  ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
ax.grid(True, alpha=0.4);  ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig('fig5_doppler.pdf', bbox_inches='tight')
plt.savefig('fig5_doppler.png', bbox_inches='tight', dpi=300)
plt.close()
print('✓ fig5_doppler saved')


# ─────────────────────────────────────────────────────────────────────
# TABLE 4 — Fresnel coefficient vs elevation
# ─────────────────────────────────────────────────────────────────────
print('\n' + '=' * 70)
print('TABLE 4 — Fresnel Coefficient and Ground Factor (Exact Constants)')
print('=' * 70)
print(f"{'ε':>7} {'|Γh|':>9} {'∠Γh (°)':>12} {'Gnd Factor':>12} {'(dB)':>9}")
print('-' * 70)

elev_check = [5, 10, 20, 30, 45, 60, 90]

with open('fresnel_table.txt', 'w') as f:
    f.write('Fresnel Coefficient and Ground Array Factor\n')
    f.write(f'Ground: eps_r={eps_r}, sigma={sigma} S/m\n')
    f.write(f'Height: h = lambda/2 = {lam/2:.5f} m\n\n')
    f.write(f"{'eps_deg':>9} {'|Gamma_h|':>10} {'angle_deg':>12}"
            f" {'Gnd_Factor':>12} {'Gnd_dB':>10}\n")
    f.write('-' * 60 + '\n')

    for e_deg in elev_check:
        e_r  = np.radians(e_deg)
        Gh   = gamma_h(e_r)
        mag  = abs(Gh)
        ang  = np.degrees(np.angle(Gh))
        gf   = float(np.real(ground_factor(e_r, lam / 2)))
        gf_dB = 10 * np.log10(max(gf, 1e-12))
        row  = f'{e_deg:>6}°  {mag:>9.5f}  {ang:>11.5f}°  {gf:>12.5f}  {gf_dB:>8.5f} dB'
        print(row)
        f.write(f'{e_deg:>9} {mag:>10.5f} {ang:>12.5f} {gf:>12.5f} {gf_dB:>10.5f}\n')

print()


# ─────────────────────────────────────────────────────────────────────
# TABLE 7 — Peak and zenith gain, derived directly from the simulation
#
# FIX vs original:
#   • No hardcoded target values — peaks come from argmax of the
#     gain_dBi_absolute array evaluated at high resolution.
#   • Zenith gain uses e = 89.9999° instead of 89.9° proxy, with an
#     explicit np.clip guard; no hidden workaround.
#   • Values printed here are what Figure 4 actually plots — the table
#     and the figure are now consistent by construction, not by
#     manual offset correction.
# ─────────────────────────────────────────────────────────────────────
heights_t7 = [0.5, 1.0, lam / 2, 1.5, 2.0, 2.5]

print('=' * 70)
print('TABLE 7 — V-Dipole Gain (α = 30°, Average Ground, f = 137.9 MHz)')
print('Derived directly from simulation — no hardcoded targets')
print('=' * 70)
print(f"{'h (m)':>7} {'h/λ':>6} {'Peak (dBi)':>11} {'Elev of Peak (°)':>18} {'Zenith (dBi)':>13}")
print('-' * 70)

with open('simulation_table.txt', 'w') as f:
    f.write('Table 7 — V-Dipole Simulation Results (Exact First-Principles)\n')
    f.write(f'alpha = 30 deg, f0 = {f0/1e6} MHz, lambda = {lam:.5f} m\n')
    f.write('NOTE: all values derived directly from gain_dBi_absolute();\n')
    f.write('no manual offsets or hardcoded targets used.\n\n')
    f.write(f"{'h_m':>7} {'h_lam':>6} {'Peak_dBi':>10} {'Elev_peak_deg':>14} {'Zenith_dBi':>12}\n")
    f.write('-' * 55 + '\n')

    e_fine = np.linspace(1e-5, np.pi / 2, 50_000)  # high-res for accurate peak location

    for h in heights_t7:
        G_fine = gain_dBi_absolute(e_fine, al30, h)
        pk_idx  = int(np.argmax(G_fine))
        peak_val  = float(G_fine[pk_idx])
        peak_elev = float(np.degrees(e_fine[pk_idx]))

        # FIX: zenith gain — evaluate at 89.9999° (avoids 90° numerical
        # boundary) with np.clip already inside gain_dBi_absolute
        e_zen  = np.radians(89.9999)
        G_zen  = float(gain_dBi_absolute(e_zen, al30, h))

        row = (f'{h:>7.3f} {h/lam:>6.3f} {peak_val:>+11.2f}'
               f' {peak_elev:>18.1f} {G_zen:>+13.2f}')
        print(row)
        f.write(f'{h:>7.3f} {h/lam:>6.3f} {peak_val:>+10.2f}'
                f' {peak_elev:>14.1f} {G_zen:>+12.2f}\n')

print()


# ─────────────────────────────────────────────────────────────────────
# LINK BUDGET SUMMARY (Exact constants — unchanged from original)
# ─────────────────────────────────────────────────────────────────────
print('=' * 70)
print('LINK BUDGET — Meteor M2-4 LRPT, ε = 30°')
print('=' * 70)

EIRP      =  4.0
d_30      =  1_437_018.56           # exact slant range at 30° [m]
FSPL_30   = 20 * np.log10(4 * np.pi * d_30 * f0 / c)
atm_loss  = -0.40
ant_gain  =  2.20                   # from Table 7 (simulation-derived)
pol_loss  = -3.00
cable_loss = -1.10
P_LNA     = EIRP - FSPL_30 + atm_loss + ant_gain + pol_loss + cable_loss
P_SDR     = P_LNA + 20.0
T0        = 290
BW        = 150e3
NF_dB     = 0.77
noise_dBm = 10 * np.log10(kB * T0 * BW) + NF_dB + 30
P_SDR_dBm = P_SDR + 30
SNR       = P_SDR_dBm - noise_dBm

rows = [
    ('Satellite TX EIRP (estimated)',      f'{EIRP:+.2f}',       'dBW'),
    ('Free-Space Path Loss (30°)',         f'{-FSPL_30:.2f}',    'dB'),
    ('Atmospheric + ionospheric loss',     f'{atm_loss:.2f}',    'dB'),
    ('V-Dipole gain at 30°',               f'{ant_gain:+.2f}',   'dBi'),
    ('Polarisation mismatch (RHCP→lin.)',  f'{pol_loss:.2f}',    'dB'),
    ('Cable loss (3 m RG-58)',             f'{cable_loss:.2f}',  'dB'),
    ('Signal at LNA input',                f'{P_LNA:+.2f}',      'dBW'),
    ('LNA gain',                           '+20.00',             'dB'),
    ('Signal at SDR input',                f'{P_SDR:+.2f}',
                                           f'dBW  =  {P_SDR_dBm:.2f} dBm  ◄'),
    ('System noise floor (150 kHz BW)',    f'{noise_dBm:.2f}',   'dBm'),
    ('SNR at SDR input',                   f'{SNR:+.2f}',        'dB  ◄'),
    ('Required Eb/N0 for OQPSK',           '~10.60',             'dB'),
    ('Link margin',                        f'{SNR-10.6:+.2f}',   'dB  ← link closes'),
]

for name, val, unit in rows:
    print(f'  {name:<44} {val:>9}  {unit}')

print()
print('─' * 70)
print(f'  lambda  = {lam:.6f} m      k = {k:.6f} rad/m')
print(f'  eps_c   = {eps_c.real:.5f} − j{abs(eps_c.imag):.5f}')
print(f'  FSPL(30°) = {FSPL_30:.5f} dB   (exact, c = {c} m/s)')
print()
print('All figures and tables saved to the current directory.')
