# Validation record — 2026-11-13

## Default Euler orientation and crystal-system display

- Updated the default Bunge ZXZ angles to (269.178546842, 27.3683882123, 93.4990822005) degrees, including the GUI hint, help text, README and example orientation CSV. Existing explicitly loaded orientations remain unchanged.
- CIF metadata now displays the crystal system and short Hermann–Mauguin symbol beneath the lattice parameters, including Unicode screw-axis subscripts and overbars. Unresolved groups remain explicitly unresolved; invalid input clears the previous description.
- All 43 tests and source GUI smoke checks passed. Added checks covering all seven crystal systems, screw-axis versus separate-axis notation, both rhombohedral settings and unresolved groups. Independently checked the formatted symbols against Gemmi short symbols for all 230 reference groups (normalizing H to conventional R).
- Source and rebuilt standalone smoke checks passed: 206 Laue contributions and 165 harmonic groups with the new orientation. Visually checked the metadata and Euler controls. Earlier validation entries below describe previous defaults.

## Nickel startup crystal

- Default example is now examples/Ni.cif, block Ni, space group 225. Updated the README and portable-settings expectation.
- All 42 tests passed. Source and rebuilt standalone smoke tests confirmed the default selection and Ni calculations: 18 powder lines, 62 monochromatic grain events, 212 Laue contributions and 179 harmonic groups. No commit or push performed.

## Tungsten defaults and additional detector presets

- Default Laue source is W at 120 kV, with the existing 5-120 keV energy window. Explicitly edited or loaded source settings remain preserved across mode switches.
- Added user-specified CCD_4x_3.362um (2032 columns x 2032 rows, 3.362 um), CCD_quadro_24um (2040 x 2040, 24 um), and Varian_127um (1456 columns x 1840 rows, 127 um). Pixel pitches apply to both axes; camera magnification is already included in the supplied effective pitch.
- All 42 tests passed, including extended detector dimension/pitch checks. Source and rebuilt standalone smoke tests passed with all 12 presets, including preset application/beam centring and the new Laue source defaults. Updated README. No commit or push performed.

## Bunge Euler inputs/exports and 5-120 keV defaults

- Replaced the single nine-value U field with three labelled degree entries under Single EulerZXZ (Bunge), and changed orientation file import/export to phi1_deg,Phi_deg,phi2_deg. Internal diffraction matrices remain crystal-to-sample U=Rz(phi1) Rx(Phi) Rz(phi2), the transpose of the passive Bunge matrix. The bundled examples/orientations.csv now uses Euler triples.
- Default angles are (301.891410, 5.389322, 57.994617) degrees, matching the previously supplied near-[001] cubic IPF-Z orientation within 2e-8 in matrix elements. Default Laue energy range is 5-120 keV; explicitly edited/loaded values remain preserved across mode switches.
- Legacy single_u JSON settings convert to Euler angles when single_euler is absent. Legacy nine-column orientation files are rejected with guidance; no silent reinterpretation. Laue provenance exports now include Euler triples and an explicit orientation convention alongside internal U matrices.
- All 42 tests passed in ImageD11_py3p12, including independent elementary-rotation composition, 1000 random Euler/matrix round trips, exact/near singular Phi=0/180 cases, CSV round trips/errors, default IPF-Z direction and legacy single-U settings migration.
- Source and rebuilt executable smoke tests passed, including three-entry/settings synchronization, mode-switch energy defaults, mono/Laue calculations and energy hover. Updated default case: 229 Laue contributions, 195 harmonic groups. Monochromatic reference remains 16 powder lines/61 events.
- Visually inspected refreshed monochromatic and Laue GUI screenshots in docs/images. Updated README and application help. Executable size is 24,235,814 bytes. No commit or push performed.

## Fixed-omega Laue mode and source spectra

- Added energy-independent signed-HKL enumeration, fixed-omega analytical wavelengths and projection through the existing tilted detector. Cached Laue enumeration keys include CIF path, size, nanosecond modification time, block and spacing bound.
- Added single-U input, random/imported orientations, separate mono/Laue controls, energy band, source preview, energy/weight colouring, harmonic-aware hover and double-click details, and per-HKL exports with embedded source provenance and actual orientation matrices.
- All 37 tests passed in ImageD11_py3p12. New checks cover analytical wavelength/position/harmonic energy ratios, same-position weight sums, nonzero-omega tilted-detector agreement with the monochromatic solver, band/sign/angular rejection, cancellation, invalid inputs, CSV validation/interpolation, reference-source support, cache reuse and full export content.
- Generated seven reference tube spectra with isolated build-only SpekPy 2.5.4: W 40/80/120 kV, Cr and Cu 40/50 kV. Fixed conditions and per-keV units are recorded in every CSV and spectra/README.md. SpekPy and its numerical dependencies are not application imports or executable dependencies.
- Source and standalone GUI smoke tests passed: original 16 powder lines/61 mono grain events, Laue 62 contributions/46 harmonic groups, seven bundled spectra, nine detectors, 177 CIFs. Checks include switching tables, disabled scan controls, energy hover, empty Laue results, selection highlighting and PNG rendering.
- Visually inspected the real Laue GUI screenshot docs/images/laue_gui.png, captured with the user-approved Pillow method. It shows a single Fe bcc grain, omega 13 degrees, 5-80 keV, W 80 kV reference spectrum, and the energy table/legend.
- Rebuilt dist/XDiffCalculator.exe, 24,232,783 bytes, and passed its expanded smoke test. No commit or push performed.
- Scientific boundary: spectral weights are S_E*F_squared*P, not quantitative Laue intensity. No Laue integration Jacobian, anomalous scattering, absorption, finite line/rocking profile or detector response. Uniform markers are geometry-only. Exact harmonics are grouped per grain; different-grain overlaps are not merged. Powder white-beam integration and continuous-rotation white-beam exposures remain outside this implementation.

## Logo and README GUI preview

- Created assets/logo.png using the built-in image-generation tool; the prompt is recorded in assets/logo_prompt.txt. Converted it to a 48-pixel Tk icon and a multi-resolution Windows ICO, preserving transparency.
- Added the logo to the application header/window and PyInstaller build. Rebuilt dist/XDiffCalculator.exe (24,132,337 bytes).
- All 30 existing tests passed in ImageD11_py3p12. The rebuilt executable smoke test passed with 16 powder lines, 61 grain events, 9 detector presets and 175 example CIFs.
- Native computer-use screenshot capture timed out. With the user's approval, Pillow captured the actual Tk demo window instead. Visually reviewed docs/images/gui.png (1400 x 845 pixels): Fe bcc, 40 keV, EIGER2 X CdTe 4M, 150 mm distance and zero tilts. The screenshot includes the new logo, input controls, intensity scale and reflection table.
- README references both images using repository-relative paths. No commit or push was performed.

## Detector tilts and GitHub documentation

- Added Tilt X/Y/Z inputs in degrees, default zero. Right-handed fixed lab rotations are applied X then Y then Z, R=Rz Ry Rx, about the direct-beam point (D,0,0); D remains the beam-axis distance.
- Added tilted ray-plane projection/inverse pixel geometry, corner-based angular coverage and clipped cone projections for powder rings. Rejects invalid, edge-on/back-facing and partly behind-sample detector planes.
- Grain polarization/Lorentz weights use outgoing laboratory azimuth rather than tilted pixel offsets. Powder hover uses inverse lab geometry. Old JSON settings without tilt fields reset all tilts to zero.
- All 30 tests passed in ImageD11_py3p12: new tests cover analytical rotation signs/order, zero-tilt compatibility, fixed beam point, 100 sampled pixel/ray round trips, angular coverage, clipping, ring-cone consistency and intensity invariance for common events.
- Tilted GUI/hover/PNG and old-settings checks passed. A PNG with X=12°, Y=-18°, Z=23° was visually inspected. Source and standalone candidate smoke checks passed including tilted calculation and PNG encoding.
- README.md was rewritten for GitHub with setup, usage, conventions, limits and tests. GITHUB_SETUP.md provides manual repository creation and release steps. No git init, commit, remote change or push was performed.
- The new executable was initially tested at build/release_candidate/XDiffCalculator.exe while the existing dist executable was running. After the user closed the application, dist/XDiffCalculator.exe was replaced, its SHA256 verified against the tested candidate, and the final executable passed the tilted-geometry/PNG smoke check. The temporary candidate executable was removed.

## PNG snapshot and taller initial plot

- Replaced SVG export with Export snapshot as png. The new snapshot.py renderer uses the displayed canvas items in draw order, including masks, text and legend, to render at 300 DPI. It preserves zoom, pan, frame, thresholds, selected-ring thickness and preview decimation. It excludes temporary selection rectangles and other UI controls.
- Added Pillow 12.0.0 to the isolated application environment and pinned runtime requirements; no change to the conda environment. Uses Windows system fonts and does not need a screen grab, Ghostscript or Matplotlib.
- Increased target window height (screen-aware, up to 1000 pixels), plot pane requested height to 700 and set initial bottom table height to 190. On the validation screen the default canvas increased from 853x329 to 853x406 pixels, giving a 2663x1268 PNG at 300 DPI.
- All 25 calculation tests passed in ImageD11_py3p12. Export checks verified the saved PNG dimensions and 300-DPI metadata, and that zoom changes the snapshot without resetting the view. Full-view PNG previews were visually inspected for labels, borders, ring colours and intensity legend.
- Source and rebuilt standalone smoke tests rendered and encoded a PNG using bundled Pillow and system fonts, then reopened it to verify dimensions/DPI. Temporary QA preview files were removed; only the current executable remains under dist.

## Live cursor readout

- Added a fixed-height readout below the plot for zero-based fractional detector column/row coordinates, nearest rendered signal intensity in arbitrary units/percent and ring HKL. Nearby grain spots also report grain ID and HKL. The six-screen-pixel hover tolerance stays constant when zooming.
- Coordinates use the inverse of the current drawing transform. Nearby ring segments and spot items are queried through Tk's spatial display list; no diffraction calculation or full redraw is performed on ordinary mouse motion.
- Ring intensity uses local polarization at the cursor azimuth. Grey grain-mode ring guides contribute no intensity. Away from features the schematic signal is zero; outside-detector readout is unavailable. This is not a rasterized detector-count image and does not sum overlapping features.
- Hidden-GUI checks passed for background coordinates, known ring index/HKL/intensity, zoom alignment, outside/leave behavior, grain spot readout and frame filtering. The fixed-height readout avoids hover-induced resize loops.
- All 25 existing tests passed in ImageD11_py3p12; source and rebuilt standalone smoke checks passed. Visual layout was not independently inspected.

## Plot border fix

- Clipping masks previously covered the bottom/right border strokes. The visible plot frame is now drawn after all four masks, using the detector/viewport intersection so zoomed borders cannot spill into the legend.
- Source and rebuilt standalone smoke checks passed. Hidden-GUI checks verified frame layering after all four masks and bounds within the viewport at full fit and zoomed/panned views. No independent visual screenshot inspection was performed.

## Adjustable panes and reciprocal-spacing update

- Added a horizontal workspace split (left controls versus plot/table) and a vertical plot/table split, both with draggable dividers and minimum pane sizes. The left form follows its available width; the bottom table expands in height and has horizontal/vertical scrolling.
- Replaced the metadata Treeview with aligned ttk labels using the same background as the surrounding controls, without borders or table lines. Headings are a (Å), b (Å), c (Å), α (°), β (°), γ (°), sgno.
- Added 1/d [Å⁻¹] immediately after d [Å] in the HKL table. Both reflection/ring CSV exports include inv_d_A_inverse computed from the unrounded d spacing.
- All 25 calculation tests passed in ImageD11_py3p12. Source and standalone smoke checks passed with 175 bundled examples.
- An invisibly mapped GUI test moved both actual sashes and verified canvas width/height changes, retained zoom/intensity values, label headings, reciprocal column order/values and d*(1/d)=1 in both CSV exports.
- Layout behavior was checked programmatically; no independent screenshot inspection was performed.

## Metadata table and rectangle navigation update

- SAMPLE now uses a two-row table headed a, b, c, α, β, γ, sgno; units are indicated beneath. Invalid metadata clears all values.
- Default end angle is 2 degrees and frame width is 0.25 degrees (eight frames). The broad-scan workflow test now explicitly requests 180 degrees because a narrow scan may legitimately have no events for a given small set of orientations.
- Replaced zoom-in/out buttons with one-shot rectangle selection. After a valid rectangle, left-drag returns to panning; a simple click still inspects spots. Mouse-wheel cursor zoom and Restore view remain. No right-button navigation bindings exist.
- All 25 tests passed in ImageD11_py3p12. Controlled hidden-GUI checks verified table headings/values, rectangle zoom in both directions, automatic return to pan mode, left-drag versus click behavior, wheel zoom, restore, unchanged intensity normalization, and absent right-button bindings.
- Source and rebuilt standalone smoke checks passed with 175 bundled CIFs, 16 powder lines and 61 events for the 20-grain default narrow scan. Frame 4 displays [1,1.25) degrees. Only dist/XDiffCalculator.exe is distributed.
- Navigation was verified programmatically; visual layout was not independently inspected in this update.

## Crystal metadata and plot navigation

- SAMPLE displays a, b, c in Å; alpha, beta, gamma in degrees; and sgno. CIF/block changes, direct path edits and loaded settings refresh metadata without computing reflections. Metadata reading does not require scattering sites; unresolved space groups remain explicitly unresolved.
- Added zoom-in/out and restore buttons, cursor-anchored mouse-wheel zoom and right-drag panning. Detector drawing is clipped to a fixed viewport; legend and intensity normalization stay fixed. New calculations restore the full detector; frame changes preserve navigation.
- SVG export intentionally remains full-detector, independent of screen zoom, and its button/README explicitly say so.
- All 25 tests passed in ImageD11_py3p12, including known Fe/Al space-group numbers and multi-block non-cubic metadata without atom sites.
- Source and rebuilt standalone hidden-Tk checks passed for selection metadata, stale-value clearing on invalid block, anchored zoom, zoom-out, panning, restore transformation and unchanged intensity normalization.
- Bundled all 175 current example CIFs unchanged. Discovery/bundling does not certify that every user-supplied CIF is valid. Only dist/XDiffCalculator.exe remains as the distribution artifact (15,312,396 bytes).
- GUI behavior was verified programmatically; visual layout was not independently inspected during this update.

## Relative-intensity update

- Random seed defaults to 200; powder-line indices start at zero. FCC Al starts with index 0 for (1,1,1). Directions renamed Multiplicities.
- Added optional signed linear polarization degree (default 0.95; blank disables P), azimuth-dependent P, powder-profile Lorentz and omega-rotation Lorentz factors. Raw F² is preserved and corrected intensities are exported separately.
- Added a linear 0-100% colour scale, using one normalization over the full grain scan. Powder rings show azimuth-dependent colours; the table uses full-ring mean P. Both GUI and SVG use the same colour function. Values are relative, not photon counts.
- All 24 tests passed in ImageD11_py3p12: includes Thomson polarization limits/average, analytical powder and rotation Lorentz cases, tangent rejection, weight products, preserved event geometry, FCC first index and colour endpoints.
- Source and standalone hidden-Tk smoke checks passed: intensity-scale items, zero-based table index, frame normalization and scan controls. Updated 20-grain seed-200 smoke calculation produced 16 powder lines and 4818 events.
- Temporary-file export checks passed for powder and grain modes: SVG parses with intensity scale and matching spot count; CSV indices, intensity products and full event counts match the calculation.
- No new runtime dependencies. Visual layout was not independently inspected; no experimental calibration/validation was performed. Near-tangent rotation intensities remain limited by the ideal zero-width reflection model, as documented in README.md.

## Scan clarity and single executable update

- Defaults: EIGER2 X CdTe 4M, distance 150 mm, maximum 2theta 40 degrees, matching dimensions/pitches and centred beam.
- Scan controls explain half-open angular exposures and the rotation convention, and are disabled for powder mode. All-frames view disables the frame selector; individual-frame view displays the angular interval and clamps frame indices to the calculated range.
- All 20 existing tests passed in ImageD11_py3p12 with Gemmi from the isolated application environment.
- Source and rebuilt single-file hidden-Tk smoke checks passed, including scan/control-state and frame-12 interval checks. The new geometry produced 16 lines and 4838 detected events for 20 grains; all 9 presets and 91 example filenames were available.
- The build script now produces only dist/XDiffCalculator.exe. Removed the obsolete dist/portable and dist/XDiffCalculator folder builds after the new executable passed validation.
- GUI controls were exercised programmatically; visual layout and experimental accuracy were not independently validated in this update. Earlier build paths and counts below are historical.

## Detector/example update

- Added EIGER2 X CdTe 4M (2068 x 2162) and 16M (4148 x 4362), both 75 um pitch, using manufacturer specifications linked in `detectors.json`.
- Added five 2048 x 2048 imaging-camera configurations with effective pitches exactly as supplied by the user: 0.55, 1.1, 1.47, 2.2 and 3.25 um.
- Added automatic example selection and kept all user CIF files unchanged.
- 20 tests passed in `ImageD11_py3p12`; updated source GUI smoke test passed and exercised all 9 presets, including beam-centre reset and both pixel pitches. It discovered 91 example CIFs.
- Read-only CIF loading audit: 85 / 91 files loaded. This checks parsing/structure loading, not experimental correctness or reference diffraction intensities.
- Six files need attention before they can be simulated:
  - `Al2O3.cif`: expected `data_` block header at the start.
  - `BaTiO3_cubic.cif`: parse error, line 3.
  - `H2O-Ice-VI.cif`: unsupported atomic scattering factor for site `Wat1`.
  - `TiCu-B2-3.14-221_adjust.cif`: parse error, line 14.
  - `Tin_tetragonal.cif`: parse error, line 10.
  - `UO2.cif`: unterminated quoted string, line 471.
- The original validation record below describes the initial build; updated build checks are recorded after it.

## Initial delivery

- Current single-file application: `dist/XDiffCalculator.exe` (15,122,981 bytes).
- Current portable application: `dist/portable/XDiffCalculator/XDiffCalculator.exe`; distribute its entire containing folder.
- Earlier `dist/XDiffCalculator/` is the initial QA build and is not the final delivery.
- Runtime: Python 3.12.13, Gemmi 0.7.5, Tcl/Tk 8.6.15. Packaging: PyInstaller 6.22.3.
- No NumPy/SciPy/Matplotlib/Qt/ImageD11 runtime imports; Gemmi is the sole third-party application dependency.

## Completed checks

- 17 deterministic unit tests passed in `ImageD11_py3p12`, loading Gemmi from the isolated project environment. The ImageD11 environment was not modified.
- Analytical BCC/FCC structure factors confirm correct handling of special-position occupancies; powder multiplicity is not applied to individual grain spots.
- Randomized diffraction roots satisfy elastic scattering and analytical detector ring radii.
- Orientation matrices, convention signs, detector row direction, half-open omega boundaries and reproducible isotropic sampling are covered.
- Bundled example paths are serialized independently of one-file temporary extraction directories.
- Source GUI smoke test passed.
- Final single-file executable smoke test passed (`singlefile_smoke.json`).
- Final portable executable smoke test passed (`portable_smoke.json`).
- A copy of the single executable, alone in `build/standalone_check`, passed with PATH restricted to Windows directories (`isolated_smoke.json`).
- Each executable check created a hidden Tk window, loaded bundled inputs, simulated 20 random BCC Fe grains, produced 29 powder lines and 13,444 detector events, and exercised integrated/frame rendering and the reflection table.

## Validation boundaries

- Native desktop screenshot capture repeatedly timed out; no successful visual screenshot review is claimed. Programmatic Tk rendering checks did pass.
- Tests were run on this Windows machine, not a separate clean Windows installation or another architecture.
- Experimental intensity accuracy is not established. Intensity corrections, detector gaps/response, finite rocking profiles and rasterized photon-count images are outside this version.
- Python 2.7 runtime is not supported by the selected modern standalone dependencies.
- No git repository initialization, commits or pushes were performed.

## Updated executable checks

- Rebuilt both `dist/XDiffCalculator.exe` and `dist/portable/XDiffCalculator/XDiffCalculator.exe` with the new presets and all 91 example CIFs.
- Single executable size: 15,197,044 bytes.
- Both packaged smoke tests passed: 9 presets, 91 examples, 29 baseline Fe powder lines and 13,444 baseline random-grain events. See `singlefile_smoke.json` and `portable_smoke.json`.
- Existing initial-QA application windows were not closed or modified. Restart using the current executable path to see the new controls and presets.
