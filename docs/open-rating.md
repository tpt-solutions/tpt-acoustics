# TPT Isolation Index (TII) v1 and v2

Status: **project definition, not a standard.** TII is not ASTM E413, not STC,
and not any other published rating. Do not report a TII value as an STC,
OITC or any other standard's result, or as compliance with one.

## Purpose

A single-number summary of a partition's transmission loss, defined openly so
that anyone can reproduce it, check it, and argue with it. It is for
comparing designs and for teaching. It is not for compliance or contracts.

## Inputs

- A transmission-loss spectrum `(frequency_hz, tl_db)` at the 16 bands
  125 Hz to 4 kHz in 1/3-octave steps (the band centres listed in
  `STC_BAND_HZ`). Each band must be present, matched to a 1 % tolerance.

## Reference contour

A straight line on a log-frequency axis:

    contour(f) = R + 5 · log2(f / 500 Hz)     [dB]

- `R` is the rating, the contour value at 500 Hz.
- The slope, 5 dB per octave, is a **design choice** of this project, close
  to the diffuse-field mass-law slope above coincidence.

## Procedure

For each integer `R` from 150 down to −60, compute at each band the deficit
`d_i = contour_i(R) − TL_i`, keeping only positive parts. The rating is the
first `R` (highest) where

- the total deficit `Σ d_i` is at most **32 dB**, and
- the largest single deficit `max d_i` is at most **8 dB**.

Deficits grow with `R`, so the first feasible value from the top is the rating.

## Validation (what the tests check)

These use no standard data.

1. **Line spectrum.** If TL is exactly on the TII line at level `a` (at 500 Hz),
   then `R = a + 2`. Each of the 16 bands can sit 2 dB above the contour before
   the total deficit passes 32 dB. Tested for `a` = 20, 30, 41.
2. **Mass law.** The diffuse mass law (`mass_law_diffuse_tl`) gives a rating
   step of 5 to 7 dB when the surface mass doubles from 5 to 10 kg/m². The
   mass law adds about 6 dB per doubling, against a 5 dB/octave reference, and
   the integer rating rounds the difference. The step observed is 5 dB, at the
   lower edge of the tested range, so the check is deliberately loose.
3. **Shift.** Raising every TL by 3 dB raises the rating by exactly 3.
4. **Missing band.** Any spectrum without all 16 bands is rejected.

These checks test the procedure and its physical behaviour. They do **not**
show that TII agrees with any standard or with measured performance.

## Version 2 (comparison variant)

v2 uses the same line contour and the same procedure, on 20 bands: the 16
STC bands plus 63, 80, 100 Hz and 5 kHz. The total-deficit limit is 2 dB per
band (40 dB), so a line spectrum still gives `level + 2`. The single-band limit
stays at 8 dB.

v2 exists to show low-frequency sensitivity. Below 125 Hz, the diffuse mass
law and a 5 dB/octave line are weak models of real partitions, because
stiffness and panel resonances dominate there. Treat v2 values as comparison
figures, and keep v1 as the headline value.

## Known limits

- It uses the 16 STC bands only, so it says nothing about low frequencies
  below 125 Hz, where mass-law and coincidence effects dominate many real
  partitions.
- It is a single number, so two partitions with the same rating can differ
  widely in performance at particular frequencies.
- The slope and limits are choices. Other choices give other ratings.
  Changing them changes the definition, so it needs a new version number.

## Provenance and licensing

- The method is written here from first principles and from the general
  description of single-number ratings that is widely taught. No text, tables
  or limits were copied from any published standard.
- The deficit limits, 32 dB and 8 dB, are the same numbers that appear in
  secondary descriptions of the STC procedure. They are numerical parameters,
  not expression copied from a standard. If you want the definition to be
  fully independent of that rule, change them and bump the version.
- The built-in E413-shaped contour in `tpt_acou_vibro::transmission` is a
  separate thing. Its values have **not** been checked against the published
  standard, and it is labelled as such in the code and reports.

## Reproducing a result

    cargo test -p tpt-acou-vibro tii_

The partition study (`tpt-acou study templates/partition-study.toml`) prints
the TII value next to the STC-style rating, with its label.
