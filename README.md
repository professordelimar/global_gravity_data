# Global Gravity Anomaly Modeling

Python workflows for interpolating global geophysical grids, calculating approximate gravity-related quantities, and extracting geographic subsets.

**Maintainer:** [professordelimar](https://github.com/professordelimar)

> **Scientific caution:** This repository implements exploratory approximations, not a validated geodetic gravity-reduction pipeline. Verify the provenance, reference systems, units, and physical meaning of the supplied GDF columns before interpreting results. The scripts have not been independently validated against reference solutions.

## Repository contents

| File | Purpose | Input structure |
| --- | --- | --- |
| `grav.gdf` | Gravity input grid | longitude, latitude, auxiliary height, gravity |
| `etopo.gdf` | Elevation/bathymetry input grid | longitude, latitude, elevation |
| `geoide.gdf` | Geoid undulation input grid | longitude, latitude, undulation |
| `gravity.py` | Global interpolation and approximate gravity computations | reads all three GDF files |
| `masked.py` | Geographic subset extraction from the computed CSV | reads `gravity_results.txt` |

The provided GDF files are whitespace-delimited, with no header. Longitude and latitude are in degrees. The original files contain a regular grid sampled at approximately **2° spacing** in both directions. Their upstream model, epoch, datum, and licensing terms are **not identified within this archive**; do not attribute them to a particular data provider without external verification. `grav.gdf` contains two numerical values after the coordinates; the script treats the first as an auxiliary height and the second as gravity. Its precise physical meaning must be checked against the original data documentation.

## Requirements

- Python 3.9+
- NumPy
- SciPy
- tqdm

Install:

```bash
python -m pip install numpy scipy tqdm
```

## Running the global calculation

Place all five files in the same directory, then run:

```bash
python gravity.py
```

The script creates `gravity_results.txt` (comma-separated, with a header) in the same directory. It evaluates a global target grid from -180° to +180° longitude and -90° to +90° latitude at `STEP = 0.05` degrees, yielding **25,930,801 target points**. Interpolation to 0.05° does **not** add spatial detail to the original ~2° source grids. Output generation can take substantial time and disk space, particularly because terrain corrections sample surrounding points at every target location. Processing is organized in `CHUNK_ROWS = 8` latitude rows per batch; the complete result is written progressively.

### Configuration (`gravity.py`)

| Setting | Default | Meaning |
| --- | --- | --- |
| `STEP` | `0.05` | Target latitude/longitude increment (degrees) |
| `CHUNK_ROWS` | `8` | Number of target latitude rows per processing batch |
| `GRAV_KIND` | `absolute` | Required gravity input interpretation; other choices are rejected |
| `GRAV_UNIT` | `mGal` | Input gravity unit (`mGal` or `m/s2`) |
| `GEOID_UNIT` | `m` | Documented geoid unit; not programmatically used for conversion |
| `DENSITY` | `2670.0` | Assumed crustal density (kg/m³) |
| `MANTLE_DENSITY` | `3300.0` | Assumed mantle density (kg/m³) |
| `ROOT_DEPTH` | `30000.0` | Defined but unused in the current calculation (m) |
| `TERRAIN_RADIUS_KM` | `20.0` | Maximum terrain-correction sampling radius (km) |
| `TERRAIN_SAMPLING_KM` | `2.0` | Radial sampling interval (km) |
| `G` | `6.67430e-11` | Gravitational constant (m³ kg⁻¹ s⁻²) |
| `MPS2_TO_MGAL` | `1e5` | Conversion factor from m/s² to mGal |

## Mathematical formulation

The following equations document **what the code computes**, rather than claiming rigorous equivalence to all conventionally defined geodetic anomalies. Unless stated otherwise, gravity-related outputs are in **mGal** and heights are in **m**.

### 1. Interpolation and heights

Each input field is interpolated on the target latitude/longitude grid using SciPy `RegularGridInterpolator`, with linear interpolation by default. Out-of-bounds evaluations return `NaN`.

Let:

- `H` = interpolated ETOPO elevation/bathymetry (m), used as an approximate orthometric height;
- `N` = interpolated geoid undulation (m);
- `h = H + N` = approximate ellipsoidal height (m);
- `g` = interpolated gravity (mGal), assuming the configured unit interpretation;
- `H_land = max(H, 0)` = nonnegative elevation used in slab and root approximations.

No vertical-datum reconciliation or coordinate-reference transformation is performed. Negative elevations remain negative for the free-air and anomaly calculations, but are set to zero for the slab and root terms.

### 2. Normal gravity

The script implements a Somigliana-style GRS80 approximation at latitude `φ`:

```text
s² = sin²(φ)
γ₀(φ) = 9.7803267715 × (1 + 0.001931851353 × s²)
         / sqrt(1 - 0.00669438002290 × s²)                [m/s²]
γ(φ,z) = [γ₀(φ) - (3.0877×10⁻⁶ - 4.3×10⁻⁹ × s²) × z
          + 7.2×10⁻¹³ × z²] × 10⁵                         [mGal]
```

Here `z` is height in meters. The height correction is approximate and should not be treated as a fully rigorous normal-gravity continuation.

### 3. Gravity anomaly (script definition)

```text
Δg_script = g - γ(φ,N) + 0.3086 × H
```

The coefficient `0.3086` has units mGal/m. This is the implementation's approximate geoid-referenced reduction; it is not automatically equivalent to a rigorous geodetic gravity anomaly.

### 4. Gravity disturbance (script definition)

```text
δg_script = g - γ(φ,H+N)
```

Both gravity values are compared at an approximate common ellipsoidal height. Interpretation depends on the actual observation height and the reference frame of the gravity source.

### 5. Free-air anomaly (approximate)

```text
Δg_FA = g - γ(φ,0) + 0.3086 × H
```

This expression applies a constant vertical-gradient approximation and does not include all geodetic or geophysical refinements.

### 6. Simple Bouguer anomaly (planar slab)

```text
C_B = 2π G ρ H_land × 10⁵
Δg_B = Δg_FA - C_B
```

Here `ρ = 2670 kg/m³`; `G` is the gravitational constant. The slab correction is applied only to positive elevations; this implementation does not model water or bathymetric density contrasts explicitly.

### 7. Approximate local terrain correction

The script samples 10 radial distances, `r = 2, 4, ..., 20 km`, and 16 azimuths per radius. For each sample, it interpolates terrain height `H_sample`, then computes `ΔH = H_sample - H`. Its contribution is:

```text
A(r) = (2000 m)² × (2π/16) × (r_km / 2)
C_T ≈ Σ [G ρ A(r) × (1/R - 1/sqrt(R² + ΔH²)) × 10⁵]
R = 1000 × r_km                                              [m]
Δg_CB = Δg_B + C_T
```

The summation runs over all radial/azimuth samples. The result is nonnegative for finite heights, depends on this particular ring-area discretization, and is truncated at 20 km. This is a **heuristic terrain correction**, not a standard high-accuracy terrain integration. Sampling may extend beyond the GDF grid at boundaries, producing `NaN` values.

### 8. Approximate Airy isostatic term

The script estimates a local compensating-root thickness:

```text
T_root = ρ H_land / (ρ_m - ρ)
C_root = -2π G (ρ_m - ρ) T_root × 10⁵
Δg_iso_script = Δg_CB - C_root
```

where `ρ_m = 3300 kg/m³`. **Important:** This is a local infinite-slab approximation. `ROOT_DEPTH` is not used, and the gravitational effect of a root at finite depth is not integrated. Therefore `isostatica` should **not** be described as a rigorous isostatic anomaly. Under these formulas, the root term algebraically cancels the simple Bouguer slab correction for positive elevations, leaving `Δg_iso_script = Δg_FA + C_T` (up to floating-point rounding).

## Output schema

The generated `gravity_results.txt` has the following comma-separated columns:

| Column | Description | Units |
| --- | --- | --- |
| `longitude` | Longitude | degrees |
| `latitude` | Latitude | degrees |
| `etopo` | Interpolated elevation/bathymetry `H` | m (assumed) |
| `orto` | Auxiliary third column from `grav.gdf` | source-dependent; verify |
| `geoide` | Geoid undulation `N` | m (assumed) |
| `gravidade` | Input gravity `g`, after optional conversion | mGal |
| `disturbio` | Script gravity disturbance | mGal |
| `anomalia` | Script gravity anomaly | mGal |
| `freeair` | Approximate free-air anomaly | mGal |
| `bouguer` | Simple Bouguer anomaly | mGal |
| `bouguer_completa` | Approximate complete Bouguer anomaly | mGal |
| `isostatica` | Approximate local isostatic quantity | mGal |

The source code writes numerical results with rounded decimal formatting; precision in the file should not be confused with accuracy.

## Geographic subset extraction (`masked.py`)

The extraction script reads `gravity_results.txt`, snaps user-specified geographic limits to the closest existing longitude and latitude coordinates, and writes matching records to `output.txt`. The default bounding box is:

```python
LONG_MIN = -52.1
LONG_MAX = -51.9
LAT_MIN = -3.5
LAT_MAX = -3.2
```

Run:

```bash
python masked.py
```

**Known compatibility issue in the supplied code:** `masked.py` attempts to find and exclude a column named `orto_grav`, but `gravity.py` writes `orto`. Without modification, this produces a `ValueError`. To keep the `orto` column, remove the exclusion logic; to drop it, change:

```python
indice_excluir = cabecalho.index("orto_grav")
```

to:

```python
indice_excluir = cabecalho.index("orto")
```

The script streams rows from disk but scans the entire input four times to find nearest coordinate values, then once more to extract rows. It does not mask a raster or reproject coordinates. Its selected bounds are inclusive.

## Reproducibility and limitations

1. **Data provenance unknown:** the archive does not identify the original geophysical model(s), license(s), geodetic datum, epoch, uncertainty, or gravity measurement definition. Verify these before scientific publication or redistribution.
2. **Input grid vs. output grid:** ~2° source sampling is interpolated to 0.05°; no new high-frequency geophysical information is generated.
3. **Absolute gravity assumption:** `gravity.py` explicitly rejects `GRAV_KIND='anomaly'` and `'disturbance'`; incorrect interpretation of `grav.gdf` invalidates derived values.
4. **Height compatibility:** the code assumes `H + N` is meaningful, without checking vertical reference consistency or whether input gravity is observed at that height.
5. **Simplified corrections:** free-air, slab, terrain, and isostatic terms omit important physical refinements and should be validated before research use.
6. **Global edge effects:** terrain sampling near ±180° and ±90° clips coordinates rather than wrapping longitude or using a pole-aware integration, and may sample outside the input domain.
7. **Runtime and output size:** the global 0.05° run can be computationally expensive and produce a large CSV. Consider testing on a smaller domain first.
8. **No automated tests:** the archive contains no test suite, reference solutions, or documented validation results.

## Suggested repository workflow

```bash
git init
git branch -M main
git remote add origin https://github.com/professordelimar/global-gravity-anomaly-modeling.git
git add grav.gdf etopo.gdf geoide.gdf gravity.py masked.py README.md
git commit -m "Initial commit: global gravity datasets and modeling workflows"
git push -u origin main
```

Create the empty repository on GitHub first. If cloning an existing repository, do not run `git init` or re-add the remote. Do not commit generated `gravity_results.txt` or `output.txt` unless explicitly needed; they can be very large.

## Attribution and licensing

Project scripts are maintained under the [professordelimar](https://github.com/professordelimar) account. No license file or data-provider attribution was included in the supplied archive. **Do not assume the GDF datasets are public domain or redistributable**; establish their source and redistribution rights before making the repository public. A software license should be added only after the maintainer chooses its terms; it would not automatically license third-party input datasets.
