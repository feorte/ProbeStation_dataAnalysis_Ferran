# ProbeStation data analysis

C++/ROOT analysis pipeline for the electrical characterisation of silicon pixel
sensors measured with a probe station. The code reads IV and CV ASCII files,
stores every measurement in ROOT `TTree`s, extracts the main sensor parameters,
and produces per-channel plots and 2D wafer maps.

The workflow is:

```
raw_data/*.txt
      |
      v
read_store.cpp  ->  stored_data/stored_data_<sensor>_<type>.root
      |
      v
analyse_data.cpp  ->  results/<type>_<sensor>.root
      |
      v
plot_means_groups*.cpp  ->  all_sensors_results*.root
```

## Repository contents

| Path | Description |
| --- | --- |
| `read_store.cpp` | Reads a raw probe-station file and stores it in two ROOT trees. |
| `analyse_data.cpp` | Analyses one stored sensor file (CV or IV) and writes plots and maps. |
| `plot_means_groups.cpp` | Summary plots of the means for sensors 09, 10, 13, 14, 15, 16, 17, 18 and 20. |
| `plot_means_groups_CSIS.cpp` | Summary plots for the CSIS2025 sensors 001-008. |
| `raw_data/` | Input measurement files. |
| `stored_data/` | ROOT files produced by `read_store.cpp`. |
| `results/` | Per-sensor ROOT output produced by `analyse_data.cpp`. |
| `mappingCSIS.txt` | LabVIEW-to-sensor channel mapping used for the CSIS sensors. |
| `script_sensors` | Batch script that runs the analysis over a list of sensors. |
| `alternative_code_versions/` | Earlier versions of the code, kept for reference only. |

## Requirements

- ROOT 6 with C++17 support (`std::filesystem` is used).
- A C++ compiler supported by ROOT.
- macOS or Linux. The examples below assume that `root` is in `PATH`.

## Input data format

Each raw file starts with a 30-line header followed by tab-separated columns.
Two header fields are used:

- Line 1: measurement type, `IV` or `CV`.
- Line 16: `Custom channels`, either `None` (all 264 channels) or a
  comma-separated list of channels.

The data columns are:

- **IV**: voltage, channel, current, current error, total current, actual
  voltage, time, temperature, humidity, followed by eight unused columns.
- **CV**: voltage, channel, `Cs`, `Cs` error, total current, actual voltage,
  time, temperature, humidity, `Cp`, `Cp` error, impedance, impedance error,
  phase, phase error, uncorrected `Cs`, uncorrected `Cp`.

Units in the files are V, nA and pF.

## Usage

### 1. Store raw measurements

`read_store.cpp` creates `stored_data/stored_data_<sensor>_<type>.root`. The
default arguments read the CSIS2025_008 IV file. To process another file, pass
the path and the sensor identifier:

```bash
root -l -q 'read_store("raw_data/CALICE_6in_256ch_CSIS2025_008_IV_CV_20250523_0_IV.txt", "CSIS_008")'
root -l -q 'read_store("raw_data/CALICE_6in_256ch_CSIS2025_008_IV_CV_20250523_0_CV.txt", "CSIS_008")'
```

The output file contains two trees:

- `raw_measurements`: one entry per channel and voltage, with scalar branches.
- `analysis`: one entry per voltage, with one value per channel stored in
  `std::vector` branches.

The `analysis` tree is the one used by `analyse_data.cpp`.

### 2. Analyse one sensor

```bash
root -l -q 'analyse_data.cpp("CSIS_008", "CV")'
root -l -q 'analyse_data.cpp("CSIS_008", "IV")'
```

The first argument is the sensor identifier used in the `stored_data` file
name, and the second is `IV` or `CV`. Results are written to:

- `results/CV_<sensor>.root` for CV.
- `results/IV_<sensor>.root` for IV.

CV analysis requires `stored_data/systemCapacitanceCV.root`. This file contains
the system capacitance measured without a sensor and is subtracted from every
CV curve.

### 3. Run the batch script

`script_sensors` calls `analyse_data.cpp` for a list of sensors and both
measurement types:

```bash
./script_sensors
```

Edit the `sensors` array in the script to change the list. Note that the script
assumes the corresponding `stored_data` files already exist.

### 4. Summary plots

```bash
root -l -q 'plot_means_groups.cpp'
root -l -q 'plot_means_groups_CSIS.cpp'
```

These macros use the mean values extracted from the individual analyses and
write `all_sensors_results.root` or `all_sensors_results_CSIS.root`.

## Analysis details

### Channel mapping and sensor geometry

The CSIS sensors use a non-sequential LabVIEW channel numbering. The mapping
is implemented in `kCsisChannelMap` in `analyse_data.cpp` and is the same as
`mappingCSIS.txt`.

For the maps, channel `ch` (1-256) is placed at:

```
x = (ch - 1) % 16 + 1
y = (ch - 1) / 16 + 1
```

The 16x16 grid represents the 256 pixels. Channels 257-264 are extra channels
used in the per-channel histograms and in the mean calculations; they are not
shown on the 16x16 maps.

Channels on the outer row/column are classified as **border** channels. All
other channels are **inner** channels. Several plots compare these two groups.

### CV: depletion voltage, capacitance and donor density

For every channel:

1. The system capacitance is subtracted and its error is propagated.
2. The curve is fitted in log-log scale:
   - a straight line (`pol1`) for the rising part, from the second point to
     `ln(35 V)`;
   - a constant (`pol0`) for the plateau, from `ln(60 V)` to the last point.
3. The depletion voltage is the intersection of both fits:
   `V_dep = exp((p0_plateau - p0_rise) / p1_rise)`.
4. The full-depletion capacitance is `exp(p0_plateau)`.
5. The donor density is obtained from a linear fit of `1/C^2` versus voltage:

```
N_don = 2 / (q * eps_Si * A^2 * slope)
```

where `q` is the electron charge, `eps_Si` is the silicon permittivity and `A`
is the pad area (`5.5 x 5.5 mm^2`).

### IV: leakage current

For every channel the leakage current at **120 V** is extracted. The value is
stored per channel and used to build the 2D current map and the current
distribution histograms. Gaussian fits are applied to the inner and border
populations.

### Invalid channels

The code marks a channel as invalid when its curve fails basic quality checks
(for example, a negative final capacitance in CV, or a current above a
threshold in IV). The neighbours of an invalid channel are also flagged in CV.
A few sensors have additional hard-coded bad channels.

The main validation thresholds are:

- CV: last capacitance below 0 pF or first capacitance below 4 pF.
- IV: current at 120 V above 0.15 nA at the border, above 0.06 nA in the
  inner region, or below 0 nA. The four corner channels use a 0.3 nA
  threshold.

These conditions and the sensor-specific corrections are intentionally kept in
the code because they depend on the probe-station setup.

## Output contents

`results/CV_<sensor>.root` contains:

- `CV_graphs/`: one `TGraphErrors` per channel.
- `Depletion_voltage/`: one canvas per channel with the log-log curve and fits.
- `Donnor_density/`: one canvas per channel with the `1/C^2` fit.
- Histograms and canvases with the per-channel values and the 2D maps of
  depletion voltage, donor density and full-depletion capacitance.

`results/IV_<sensor>.root` contains:

- `IV_graphs/`: one `TGraphErrors` per channel.
- Per-channel current histogram and the 2D current map at 120 V.

## Configuration notes

The code keeps several values hard-coded because they are part of the current
measurement campaign:

- IV map voltage: `120 V`.
- CV fit ranges: `35 V` and `60 V`.
- Pad area: `5.5 x 5.5 mm^2`.
- System-capacitance corrections for sensors 15, 16, 17 and 18.
- Sensor-specific invalid-channel lists.

The CSIS channel mapping is currently enabled for every sensor (`csis = true`
in both analysis functions). It should be set to `false` for sensors whose
stored channels already use the physical sensor numbering.

The summary input in `plot_means_groups_CSIS.cpp` contained the typo
`1.27645e+111`; it has been corrected to `1.27645e+11`.
