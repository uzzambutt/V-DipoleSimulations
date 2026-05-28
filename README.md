# V-DipoleSimulation Output

**Research Paper**: https://zenodo.org/records/20424227 

```
E:\ResearchPaperSimulations>py sum1.py
Using font: Times New Roman
✓  fig1_azimuth  saved
✓  fig2_elevation_vangle  saved
✓  fig3_height_patterns  saved
✓  fig4_elevation_overlay  saved
✓  fig5_doppler  saved

======================================================================
TABLE 4 — Fresnel Coefficient and Ground Factor (Exact Constants)
======================================================================
 ε (deg)      |Γh|     ∠Γh (deg)    Ground Factor   Factor (dB)
----------------------------------------------------------------------
       5°   0.95101     179.90633°        0.27892    -5.54523 dB
      10°   0.90479     179.81354°        0.97823    -0.09561 dB
      20°   0.82133     179.63406°        2.56348     4.08830 dB
      30°   0.75034     179.46793°        3.06362     4.86235 dB
      45°   0.66708     179.25508°        1.81692     2.59336 dB
      60°   0.61007     179.09662°        0.57386    -2.41193 dB
      90°   0.56623     178.96692°        0.18834    -7.25051 dB

======================================================================
TABLE 7 — Simulated Peak Gain (Exact Constants)
======================================================================
  h (m)     h/λ   Peak (dBi)    Elev of Peak (°)   Zenith (dBi)
----------------------------------------------------------------------
  0.500   0.230        +2.10                90.0          +2.10
  1.000   0.460        +2.20                32.4          -6.44
  1.087   0.500        +2.20                29.5          -9.01
  1.500   0.690        +2.30                20.8          +1.48
  2.000   0.920        +2.50                15.5          -3.08
  2.500   1.150        +2.60                12.3          +0.46

======================================================================
LINK BUDGET — Meteor M2-4 LRPT, ε = 30° (Exact Constants)
======================================================================
  Satellite TX EIRP (estimated)                    +4.00  dBW
  Free-Space Path Loss (30°)                     -138.39  dB
  Atmospheric + ionospheric loss                   -0.40  dB
  V-Dipole gain at 30°                             +2.20  dBi
  Polarisation mismatch (RHCP→lin.)                -3.00  dB
  Cable loss (3 m RG-58)                           -1.10  dB
  Signal at LNA input                            -136.69  dBW
  LNA gain                                        +20.00  dB
  Signal at SDR input                            -116.69  dBW  = -86.69 dBm ◄
  System noise floor (150 kHz BW)                -121.44  dBm
  SNR at SDR input                                +34.76  dB ◄
  Required Eb/N0 for OQPSK                        ~10.60  dB
  Link margin                                     +24.16  dB  ← link closes ◄

All figures and tables saved to current directory.
  lambda = 2.17398 m    k = 2.89017 rad/m
  eps_c  = 13.00000 - j0.78209
  Exact FSPL at 30° = 138.38832 dB
```
