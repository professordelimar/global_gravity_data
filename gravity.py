from pathlib import Path
import time
from tqdm import tqdm
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

BASE = Path(__file__).resolve().parent
STEP = 0.05
CHUNK_ROWS = 8
GRAV_KIND = 'absolute'  # 'absolute' (mGal), 'anomaly' (mGal), 'disturbance' (mGal)
GRAV_UNIT = 'mGal'  # 'mGal' or 'm/s2'
GEOID_UNIT = 'm'
DENSITY = 2670.0  # kg/m3, crust
MANTLE_DENSITY = 3300.0
ROOT_DEPTH = 30000.0  # m
TERRAIN_RADIUS_KM = 20.0
TERRAIN_SAMPLING_KM = 2.0
G = 6.67430e-11
MPS2_TO_MGAL = 1e5


def read_gdf(path, n_values=1):
    """Read longitude, latitude and n_values columns without longitude conversion."""
    records = []
    with path.open(encoding='utf-8-sig', errors='replace') as f:
        for line in f:
            parts = line.strip().replace(',', ' ').split()
            if len(parts) < 2 + n_values:
                continue
            try:
                values = tuple(map(float, parts[:2 + n_values]))
            except ValueError:
                continue
            if np.isfinite(values).all():
                records.append(values)
    if not records:
        raise ValueError(f'No valid records found: {path.name}')
    a = np.asarray(records, dtype=float)
    lon, lat = a[:, 0], a[:, 1]
    if np.any((lon < -180) | (lon > 180)) or np.any((lat < -90) | (lat > 90)):
        raise ValueError(f'{path.name}: coordinates outside the expected bounds')
    lo, la = np.unique(lon), np.unique(lat)
    grids = np.full((n_values, len(la), len(lo)), np.nan)
    ii, jj = np.searchsorted(la, lat), np.searchsorted(lo, lon)
    for k in range(n_values):
        grids[k, ii, jj] = a[:, 2 + k]
    if np.isnan(grids).any():
        raise ValueError(f'{path.name}: incomplete or irregular grid')
    return la, lo, grids


def interpolators(name, n_values=1):
    la, lo, grids = read_gdf(BASE / name, n_values)
    return [RegularGridInterpolator((la, lo), grids[k], bounds_error=False,
                                    fill_value=np.nan) for k in range(n_values)]


def normal_gravity(lat_deg, h_m=0):
    """Somigliana/GRS80 with first-order height correction; mGal."""
    phi = np.deg2rad(lat_deg)
    s2 = np.sin(phi) ** 2
    gamma0 = 9.7803267715 * (1 + 0.001931851353 * s2) / np.sqrt(1 - 0.00669438002290 * s2)
    return (gamma0 - (3.0877e-6 - 4.3e-9 * s2) * h_m + 7.2e-13 * h_m**2) * MPS2_TO_MGAL


def terrain_correction(h, lat, lon, interp):
    """Positive terrain correction (mGal), approximate mass ring integration.
    Integration is truncated at TERRAIN_RADIUS_KM; this is not a global correction.
    """
    r = np.arange(TERRAIN_SAMPLING_KM, TERRAIN_RADIUS_KM + 0.01, TERRAIN_SAMPLING_KM)
    az = np.linspace(0, 2*np.pi, 16, endpoint=False)
    corr = np.zeros_like(h)
    for rr in r:
        for angle in az:
            dy = rr * np.cos(angle) / 111.195
            dx = rr * np.sin(angle) / (111.195 * np.maximum(np.cos(np.deg2rad(lat)), 0.02))
            pts = np.column_stack((np.clip(lat + dy, -90, 90), np.clip(lon + dx, -180, 180)))
            ht = interp(pts)
            dh = ht - h
            # Approximate attraction integral for relative relief: always positive
            R = rr * 1000
            area = (TERRAIN_SAMPLING_KM * 1000)**2 * (2*np.pi/16) * (rr / TERRAIN_SAMPLING_KM)
            corr += G * DENSITY * area * (1/R - 1/np.sqrt(R**2 + dh**2)) * MPS2_TO_MGAL
    return corr


def main():
    started = time.perf_counter()
    print('[1/3] Reading GDF grids...', flush=True)
    for filename in ('grav.gdf', 'etopo.gdf', 'geoide.gdf'):
        if not (BASE / filename).is_file():
            raise FileNotFoundError(f'Missing input file: {BASE / filename}')
    grav_orto, grav = interpolators('grav.gdf', 2)
    topo = interpolators('etopo.gdf')[0]
    geoid = interpolators('geoide.gdf')[0]
    lons = np.round(np.arange(-180, 180 + STEP/2, STEP), 8)
    lats = np.round(np.arange(-90, 90 + STEP/2, STEP), 8)
    output = BASE / 'gravity_results.txt'
    if GRAV_KIND not in ('absolute', 'anomaly', 'disturbance'):
        raise ValueError('Invalid GRAV_KIND')
    if GRAV_UNIT not in ('mGal', 'm/s2'):
        raise ValueError('Invalid GRAV_UNIT')
    if GRAV_KIND != 'absolute':
        raise ValueError('Absolute gravity is required to calculate gravity anomaly and disturbance separately (GRAV_KIND=absolute).')
    total_points = len(lons) * len(lats)
    print(f'[2/3] Interpolating and calculating {total_points:,} points in blocks...', flush=True)
    print('Calculation order: gravity anomaly, gravity disturbance, free-air, simple Bouguer, complete Bouguer, isostatic.', flush=True)
    with output.open('w', encoding='utf-8', buffering=1024*1024) as f:
        f.write('longitude,latitude,etopo,orto,geoide,gravidade,disturbio,anomalia,freeair,bouguer,bouguer_completa,isostatica\n')
        with tqdm(total=total_points, desc='Processing', unit='points', unit_scale=True, dynamic_ncols=True) as progress:
            for start in range(0, len(lats), CHUNK_ROWS):
                ys = lats[start:start+CHUNK_ROWS]
                xx, yy = np.meshgrid(lons, ys)
                lon, lat = xx.ravel(), yy.ravel()
                pts = np.column_stack((lat, lon))
                H = topo(pts)  # approximate orthometric elevation / bathymetry
                N = geoid(pts)
                altura_terreno = np.maximum(H, 0)
                g = grav(pts)
                H_grav = grav_orto(pts)  # third column of grav.gdf, used for output orto
                if GRAV_UNIT == 'm/s2': g = g * MPS2_TO_MGAL
                h_ell = H + N
                gamma_surface = normal_gravity(lat, h_ell)
                gamma_ellipsoid = normal_gravity(lat)
                # 1 Gravity anomaly: observed gravity minus normal gravity
                # reduced to the geoid using an approximate vertical gradient.
                gamma_geoid = normal_gravity(lat, N)
                anomaly = g - gamma_geoid + 0.3086 * H
                # 2 Gravity disturbance: difference at the same spatial point.
                disturbance = g - gamma_surface
                # 3 Free-air anomaly: approximate free-air reduction.
                freeair = g - gamma_ellipsoid + 0.3086 * H
                # 4 Simple Bouguer anomaly: planar slab for land elevations only.
                simple_bouguer = freeair - 2*np.pi*G*DENSITY*altura_terreno*MPS2_TO_MGAL
                # 5 Complete Bouguer anomaly: approximate local terrain correction.
                terrain = terrain_correction(H, lat, lon, topo)
                complete_bouguer = simple_bouguer + terrain
                # 6 Isostasy: local Airy model, crustal root below topography.
                # Infinite slab approximation; does NOT integrate the gravitational effect
                # of roots at depth and is not a rigorous isostatic anomaly.
                root_thickness = DENSITY * altura_terreno / (MANTLE_DENSITY - DENSITY)
                root_effect = -2*np.pi*G*(MANTLE_DENSITY-DENSITY)*root_thickness*MPS2_TO_MGAL
                isostatic = complete_bouguer - root_effect
                columns = np.column_stack((lon, lat, H, H_grav, N, g, disturbance,
                                           anomaly, freeair, simple_bouguer, complete_bouguer, isostatic))
                np.savetxt(f, columns, delimiter=',', fmt=['%.3f', '%.3f', '%.3f', '%.1f', '%.5f', '%.4f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f'])
                progress.update(len(lon))
    print('[3/3] Output saved:', output, flush=True)
    print(f'Total elapsed time: {(time.perf_counter() - started) / 60:.2f} minutes', flush=True)
    #print('WARNING: terrain and isostasy are exploratory approximations; verify GDF metadata before scientific use.')


if __name__ == '__main__':
    main()