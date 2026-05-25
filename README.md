# V-DipoleSimulation Output
```
======================================================================
TABLE 4 — Fresnel Coefficient and Ground Factor  (εr=13, σ=0.006 S/m)
======================================================================
 ε (deg)      |Γh|     ∠Γh (deg)    Ground Factor   Factor (dB)
----------------------------------------------------------------------
       5°     0.951         179.9°          0.279        -5.5 dB
      10°     0.905         179.8°          0.978        -0.1 dB
      20°     0.821         179.6°          2.563         4.1 dB
      30°     0.750         179.5°          3.064         4.9 dB
      45°     0.667         179.3°          1.817         2.6 dB
      60°     0.610         179.1°          0.574        -2.4 dB
      90°     0.566         179.0°          0.188        -7.3 dB

======================================================================
TABLE 7 — Simulated Peak Gain  (α=30°, Average Ground, f=137.9 MHz)
======================================================================
  h (m)     h/λ   Peak (dBi)    Elev of Peak (°)   Zenith (dBi)
----------------------------------------------------------------------
   0.50   0.230         +2.1                  90           +2.1
   1.00   0.460         +2.2                  32           -6.5
   1.09   0.501         +2.2                  29           -9.0
   1.50   0.690         +2.3                  21           +1.5
   2.00   0.920         +2.5                  15           -3.1
   2.50   1.150         +2.6                  12           +0.4

======================================================================
LINK BUDGET — Meteor M2-4 LRPT,  ε = 30°,  With LNA at feed-point
======================================================================
  Satellite TX EIRP (estimated)                     +4.0  dBW
  Free-Space Path Loss (30°)                      -139.6  dB
  Atmospheric + ionospheric loss                    -0.4  dB
  V-Dipole gain at 30°                              +2.2  dBi
  Polarisation mismatch (RHCP→lin.)                 -3.0  dB
  Cable loss (3 m RG-58)                            -1.1  dB
  Signal at LNA input                             -137.9  dBW
  LNA gain                                         +20.0  dB
  Signal at SDR input                             -117.9  dBW  = -87.9 dBm ◄
  System noise floor (150 kHz BW)                 -121.4  dBm
  SNR at SDR input                                 +33.6  dB ◄
  Required Eb/N0 for OQPSK                         ~10.6  dB
  Link margin                                      +23.0  dB  ← link closes ◄

All figures and tables saved to current directory.
  lambda = 2.1740 m    k = 2.8901 rad/m
  eps_c  = 13.00 - j0.7821
  FSPL at 30° = 139.6 dB
  ```
