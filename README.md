# V-DipoleSimulation Output

**Research Paper**: https://zenodo.org/records/20424227 

```
E:\ResearchPaperSimulations>py simulation.py
Self-test: 6/6 checks passed ✓
Font: Times New Roman
✓ fig1_azimuth saved
✓ fig2_elevation_vangle saved
✓ fig3_height_patterns saved
✓ fig4_elevation_overlay saved
✓ fig5_doppler saved

======================================================================
TABLE 4 — Fresnel Coefficient and Ground Factor (Exact Constants)
======================================================================
      ε      |Γh|      ∠Γh (°)   Gnd Factor      (dB)
----------------------------------------------------------------------
     5°    0.95101    179.90633°       0.27892  -5.54523 dB
    10°    0.90479    179.81354°       0.97823  -0.09561 dB
    20°    0.82133    179.63406°       2.56348   4.08830 dB
    30°    0.75034    179.46793°       3.06362   4.86235 dB
    45°    0.66708    179.25508°       1.81692   2.59336 dB
    60°    0.61007    179.09662°       0.57386  -2.41193 dB
    90°    0.56623    178.96692°       0.18834  -7.25051 dB

======================================================================
TABLE 7 — V-Dipole Gain (α = 30°, Average Ground, f = 137.9 MHz)
Derived directly from simulation — no hardcoded targets
======================================================================
  h (m)    h/λ  Peak (dBi)   Elev of Peak (°)  Zenith (dBi)
----------------------------------------------------------------------
  0.500  0.230       +2.06               90.0         +2.06
  1.000  0.460       +2.17               32.4         -6.47
  1.087  0.500       +2.20               29.5         -9.01
  1.500  0.690       +2.34               20.8         +1.52
  2.000  0.920       +2.46               15.5         -3.12
  2.500  1.150       +2.55               12.3         +0.41

======================================================================
LINK BUDGET — Meteor M2-4 LRPT, ε = 30°
======================================================================
  Satellite TX EIRP (estimated)                    +4.00  dBW
  Free-Space Path Loss (30°)                     -138.39  dB
  Atmospheric + ionospheric loss                   -0.40  dB
  V-Dipole gain at 30°                             +2.20  dBi
  Polarisation mismatch (RHCP→lin.)                -3.00  dB
  Cable loss (3 m RG-58)                           -1.10  dB
  Signal at LNA input                            -136.69  dBW
  LNA gain                                        +20.00  dB
  Signal at SDR input                            -116.69  dBW  =  -86.69 dBm  ◄
  System noise floor (150 kHz BW)                -121.44  dBm
  SNR at SDR input                                +34.76  dB  ◄
  Required Eb/N0 for OQPSK                        ~10.60  dB
  Link margin                                     +24.16  dB  ← link closes

──────────────────────────────────────────────────────────────────────
  lambda  = 2.173984 m      k = 2.890170 rad/m
  eps_c   = 13.00000 − j0.78209
  FSPL(30°) = 138.38832 dB   (exact, c = 299792458.0 m/s)

All figures and tables saved to the current directory.
```
