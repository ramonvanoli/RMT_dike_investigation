import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.path import Path
from matplotlib.ticker import FuncFormatter

# This script reads a GeoTIFF DEM (.tif) and plots a 3D elevation surface.
# Edit the paths/settings below, then press Run. The chosen topography is exported to a .xyz file.
# There are some settings that can be changed where holes or spikes in the .tif file can be filled or removed.

plt.style.use("ggplot")

try:
    import rasterio
    from rasterio.transform import xy
    from pyproj import Transformer
except ImportError as exc:
    raise ImportError(
        "Install rasterio and pyproj in this environment: pip install rasterio pyproj"
    ) from exc

# ── user settings ─────────────────────────────────────────────────────────────
TIF_PATH = r"C:\Users\Baar\OneDrive - Vanoli AG\HS_25\Master Thesis\Inversion\3D\Dike3\AHN5_M_103000_430000.TIF"   # input GeoTIFF
DOWNSAMPLE = 2                          # use 4–10 for large rasters (keeps every Nth pixel)
OUTPUT = None                           # e.g. r"dem_3d.png" to save; None = only show plot
EXPORT_XYZ = r"C:/Users/Baar/OneDrive - Vanoli AG/HS_25/Master Thesis/Inversion/3D/Dike3/topography_dike1_FINAL.xyz" # lon lat elev; None = skip

# AOI for export/plot — choose one of the two options below:
SELECT_AOI_ON_MAP = True   # True: click 4 corners on a map; False: use AOI dict
OVERVIEW_DOWNSAMPLE = 8    # resolution of the click-map only (higher = faster)

# DEM repair (hole fill / spikes). Ignored when RAW_FROM_TIF is True.
# ── RAW vs repair ─────────────────────────────────────────────────────────────
# True  = use TIFF elevations as-is inside the AOI (recommended baseline).
#         No hole fill, no spike removal. Interior nodata stays as in the TIFF.
# False = allow optional repair below (FILL_HOLES / REMOVE_SPIKES).
RAW_FROM_TIF = True

# Only used when RAW_FROM_TIF is False:
FILL_HOLES = False          # fill small true-nodata gaps inside the AOI
MAX_HOLE_CELLS = 2         # skip connected voids larger than this
REMOVE_SPIKES = False       # remove isolated junk pixels
SPIKE_THRESHOLD_M = 1.5     # |elev − local median| above this → candidate
SPIKE_WINDOW = 3
MAX_SPIKE_CELLS = 5
SPIKE_ABS_MARGIN_M = 1.0

# M3D needs a full rectangular grid. Cells outside the clicked quad have no TIFF
# value; if True they are nearest-neighbour filled from the interior edge only.
# Interior TIFF nodata is NEVER filled by this step (unlike the old export behaviour).
FILL_EXPORT_CORNERS = True

# Clip in the TIF's native map coordinates (RD metres for AHN). None = full tile.
# Used only when SELECT_AOI_ON_MAP is False.
AOI = dict(min_x=102900, max_x=103000, min_y=418200, max_y=418400)
# AOI = None
# ───────────────────────────────────────────────────────────────────────────────


def read_dem(tif_path, aoi=None):
    from rasterio.windows import from_bounds

    with rasterio.open(tif_path) as src:
        crs = src.crs
        if aoi is None:
            window = None
            transform = src.transform
        else:
            window = from_bounds(
                aoi["min_x"],
                aoi["min_y"],
                aoi["max_x"],
                aoi["max_y"],
                transform=src.transform,
            )
            transform = src.window_transform(window)

        elev = src.read(1, window=window).astype(np.float64)
        nodata = src.nodata
        if nodata is not None:
            elev = np.where(elev == nodata, np.nan, elev)
        elev = np.where(elev <= -1e30, np.nan, elev)

        rows, cols = np.indices(elev.shape)
        xs, ys = xy(transform, rows.ravel(), cols.ravel(), offset="center")
        x = np.asarray(xs, dtype=np.float64).reshape(elev.shape)
        y = np.asarray(ys, dtype=np.float64).reshape(elev.shape)

    return x, y, elev, crs


def map_to_lonlat(x, y, crs):
    to_wgs84 = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    return to_wgs84.transform(x, y)


def lonlat_to_local_metres(lon, lat):
    """Azimuthal equidistant projection centred on the AOI (metres, true local scale)."""
    lon0 = float(np.nanmean(lon))
    lat0 = float(np.nanmean(lat))
    local_crs = f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +ellps=WGS84 +units=m +no_defs"
    to_local = Transformer.from_crs("EPSG:4326", local_crs, always_xy=True)
    east, north = to_local.transform(lon, lat)
    return east, north, lon0, lat0, local_crs


def print_raster_info(tif_path):
    with rasterio.open(tif_path) as src:
        b = src.bounds
        print(f"CRS: {src.crs}")
        print(
            f"Native bounds — min_x: {b.left}, min_y: {b.bottom}, "
            f"max_x: {b.right}, max_y: {b.top}"
        )
        print(f"Pixel size: {src.res[0]} m, grid size: {src.width} x {src.height}")

        corners_x = [b.left, b.right, b.left, b.right]
        corners_y = [b.bottom, b.bottom, b.top, b.top]
        lon, lat = map_to_lonlat(corners_x, corners_y, src.crs)
        print(
            f"Lon/lat bounds — min_lon: {min(lon):.6f}, max_lon: {max(lon):.6f}, "
            f"min_lat: {min(lat):.6f}, max_lat: {max(lat):.6f}"
        )


def downsample_grid(x, y, elev, factor):
    if factor <= 1:
        return x, y, elev
    return x[::factor, ::factor], y[::factor, ::factor], elev[::factor, ::factor]


def fill_holes(
    elev,
    inside_mask=None,
    remove_spikes=False,
    spike_threshold_m=1.5,
    window=3,
    max_hole_cells=50,
    max_spike_cells=5,
    spike_abs_margin_m=1.0,
):
    """Fill true nodata (and optionally clear junk spikes) inside the AOI only.

    Important design choices so the DEM stays essentially unchanged:
    - Valid elevations are never smoothed or overwritten en masse.
    - Spike removal (when enabled) is for AHN-scale relief (~ metres, not tens):
        * absolute outliers outside the AOI [p1, p99] ± margin
        * local |elev − median| > threshold, but ONLY in small isolated blobs
          (<= *max_spike_cells*). Continuous dike slopes are left alone.
    - Only small connected nodata voids (<= *max_hole_cells*) are filled with
      nearest-neighbour. Large voids stay NaN.

    Cells outside *inside_mask* (beyond the clicked quad) are never modified.
    """
    from scipy.interpolate import griddata
    from scipy.ndimage import label, median_filter

    elev = elev.astype(np.float64, copy=True)
    if inside_mask is None:
        inside_mask = np.ones(elev.shape, dtype=bool)
    else:
        inside_mask = inside_mask.astype(bool, copy=False)

    finite = np.isfinite(elev)
    if not finite.any():
        return elev

    # --- optional junk spikes: mark as NaN, then fill like other holes ----------
    n_spikes = 0
    if remove_spikes:
        interior = finite & inside_mask
        z_in = elev[interior]
        p1, p99 = np.percentile(z_in, [1.0, 99.0])
        abs_lo, abs_hi = p1 - spike_abs_margin_m, p99 + spike_abs_margin_m
        spikes_abs = interior & ((elev < abs_lo) | (elev > abs_hi))

        filler = float(np.median(z_in))
        tmp = np.where(finite, elev, filler)
        med = median_filter(tmp, size=window, mode="nearest")
        spikes_local = interior & (np.abs(elev - med) > spike_threshold_m)

        # Keep only small connected local-anomaly blobs (single bad pixels /
        # tiny clusters). Large connected "anomalies" are real landforms.
        labeled_loc, n_loc = label(spikes_local)
        spikes_iso = np.zeros_like(spikes_local, dtype=bool)
        for cid in range(1, n_loc + 1):
            component = labeled_loc == cid
            if int(component.sum()) <= max_spike_cells:
                spikes_iso |= component

        spikes = spikes_abs | spikes_iso
        n_spikes = int(spikes.sum())
        elev[spikes] = np.nan
        print(
            f"  Spike mask: {int(spikes_abs.sum())} absolute outliers "
            f"(outside [{abs_lo:.2f}, {abs_hi:.2f}] m), "
            f"{int(spikes_iso.sum())} isolated local outliers "
            f"(>{spike_threshold_m:g} m from {window}×{window} median, "
            f"≤{max_spike_cells} cells)."
        )

    # --- only fill small interior nodata components ----------------------------
    holes = ~np.isfinite(elev) & inside_mask
    n_holes = int(holes.sum())
    n_filled = 0
    n_skipped_large = 0
    if n_holes:
        labeled, n_comp = label(holes)
        fillable = np.zeros_like(holes, dtype=bool)
        for cid in range(1, n_comp + 1):
            component = labeled == cid
            size = int(component.sum())
            if size <= max_hole_cells:
                fillable |= component
            else:
                n_skipped_large += size

        n_to_fill = int(fillable.sum())
        if n_to_fill:
            valid = np.isfinite(elev)
            rows, cols = np.indices(elev.shape)
            filled = griddata(
                np.column_stack([rows[valid], cols[valid]]),
                elev[valid],
                np.column_stack([rows[fillable], cols[fillable]]),
                method="nearest",
            )
            elev[fillable] = filled
            n_filled = int(np.sum(np.isfinite(filled)))

    n_outside = int(np.sum(~inside_mask & ~np.isfinite(elev)))
    print(
        f"Hole repair: removed {n_spikes} spikes, filled {n_filled}/{n_holes} "
        f"interior nodata cells "
        f"(skipped {n_skipped_large} in voids larger than {max_hole_cells} cells; "
        f"{n_outside} exterior cells left masked)."
    )
    return elev


def lonlat_corners_to_native_bbox(corners_lonlat, crs):
    """Axis-aligned native bbox that contains all four lon/lat corners."""
    to_native = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    xs, ys = to_native.transform(corners_lonlat[:, 0], corners_lonlat[:, 1])
    return {
        "min_x": float(min(xs)),
        "max_x": float(max(xs)),
        "min_y": float(min(ys)),
        "max_y": float(max(ys)),
    }


def mask_elev_to_quad(lon, lat, elev, quad_lonlat):
    """Keep elevation only inside the lon/lat quadrilateral.

    Returns (elev_masked, inside_mask) where *inside_mask* is True for every
    cell geometrically inside the quad (even if the TIFF value is nodata).
    """
    path = Path(quad_lonlat)
    inside = path.contains_points(np.column_stack([lon.ravel(), lat.ravel()]))
    inside = inside.reshape(elev.shape)
    return np.where(inside, elev, np.nan), inside


def read_dem_quad(tif_path, quad_lonlat):
    """Read DEM bbox around a 4-corner selection and mask outside the quad."""
    with rasterio.open(tif_path) as src:
        crs = src.crs
    aoi = lonlat_corners_to_native_bbox(quad_lonlat, crs)
    x, y, elev, crs = read_dem(tif_path, aoi=aoi)
    lon, lat = map_to_lonlat(x, y, crs)
    elev, inside_mask = mask_elev_to_quad(lon, lat, elev, quad_lonlat)
    return x, y, elev, crs, inside_mask


def crop_to_valid_bbox(lon, lat, elev, *extra):
    """Crop to the smallest grid rectangle that still contains all valid cells."""
    valid = np.isfinite(elev)
    if not np.any(valid):
        raise ValueError("No valid elevation values inside the selected area.")

    row_ok = np.any(valid, axis=1)
    col_ok = np.any(valid, axis=0)
    r0, r1 = int(np.where(row_ok)[0][0]), int(np.where(row_ok)[0][-1] + 1)
    c0, c1 = int(np.where(col_ok)[0][0]), int(np.where(col_ok)[0][-1] + 1)
    cropped = (
        lon[r0:r1, c0:c1],
        lat[r0:r1, c0:c1],
        elev[r0:r1, c0:c1],
    )
    if extra:
        cropped += tuple(arr[r0:r1, c0:c1] for arr in extra)
    return cropped


def select_aoi_on_map(tif_path, overview_downsample=8):
    """Show a lon/lat map of the full tile; return 4 corners after clicks."""
    print(
        "Opening map — click 4 corners of the export area in order "
        "(clockwise or counter-clockwise around the rectangle)."
    )
    x, y, elev, crs = read_dem(tif_path, aoi=None)
    x, y, elev = downsample_grid(x, y, elev, overview_downsample)
    lon, lat = map_to_lonlat(x, y, crs)

    valid = np.isfinite(elev)
    if not np.any(valid):
        raise ValueError("No valid elevation values in the raster.")
    elev_plot = np.where(valid, elev, np.nanmin(elev[valid]))

    fig, ax = plt.subplots(figsize=(11, 8))
    cf = ax.contourf(lon, lat, elev_plot, cmap="terrain")
    fig.colorbar(cf, ax=ax, shrink=0.85, label="Elevation (m)")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Click 4 corners of the export area (in order around the rectangle)")
    plt.tight_layout()

    pts = plt.ginput(4, timeout=0)
    if len(pts) < 4:
        plt.close(fig)
        raise ValueError("Export cancelled — select four corners on the map.")

    quad_lonlat = np.asarray(pts, dtype=np.float64)
    ax.add_patch(
        Polygon(
            quad_lonlat,
            closed=True,
            fill=False,
            edgecolor="red",
            linewidth=2.5,
        )
    )
    for i, (px, py) in enumerate(quad_lonlat):
        ax.plot(px, py, "ro", ms=8)
        ax.annotate(
            str(i + 1),
            (px, py),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=10,
            color="red",
            fontweight="bold",
        )
    ax.set_title("Selected export area — close window to continue")
    plt.show()
    plt.close(fig)

    aoi = lonlat_corners_to_native_bbox(quad_lonlat, crs)
    print("Selected corners (lon, lat):")
    for i, (lo, la) in enumerate(quad_lonlat, start=1):
        print(f"  {i}: lon={lo:.6f}, lat={la:.6f}")
    print(
        f"Bounding native AOI — min_x: {aoi['min_x']:.1f}, max_x: {aoi['max_x']:.1f}, "
        f"min_y: {aoi['min_y']:.1f}, max_y: {aoi['max_y']:.1f}"
    )
    return quad_lonlat


def export_xyz(path, lon, lat, elev, interior_mask=None, fill_export_corners=True):
    """Write lon, lat, elevation (m) for M3D xyzt import (regular grid, sorted).

    M3D does: topo_lon = unique(A(:,1)); topo_lat = unique(A(:,2));
              topo_z = reshape(A(:,3), n_lon, n_lat)';
    So rows must be every (lon, lat) pair, lat slow / lon fast, sorted uniques.

    By default only *exterior* cells (outside the clicked quad) are filled for
    the rectangular M3D grid. Interior TIFF elevations and interior nodata are
    left unchanged.
    """
    n_lat, n_lon = elev.shape
    lon_vec = lon[0, :]
    lat_vec = lat[:, 0]

    if interior_mask is None:
        interior_mask = np.isfinite(elev)
    else:
        interior_mask = interior_mask.astype(bool, copy=False)

    z_grid = elev.astype(np.float64, copy=True)
    rows, cols = np.indices(elev.shape)
    valid_src = interior_mask & np.isfinite(elev)
    if not valid_src.any():
        raise ValueError("No valid interior elevation values to export.")

    exterior_nan = (~interior_mask) & ~np.isfinite(z_grid)
    interior_nan = interior_mask & ~np.isfinite(z_grid)
    n_interior_nan = int(interior_nan.sum())
    n_exterior_filled = 0

    if fill_export_corners and exterior_nan.any():
        from scipy.interpolate import griddata

        filled = griddata(
            np.column_stack([rows[valid_src], cols[valid_src]]),
            elev[valid_src],
            np.column_stack([rows[exterior_nan], cols[exterior_nan]]),
            method="nearest",
        )
        z_grid[exterior_nan] = filled
        n_exterior_filled = int(np.sum(np.isfinite(filled)))
        print(
            f"Export: filled {n_exterior_filled} exterior (outside-quad) cells "
            f"for rectangular M3D grid; interior TIFF values untouched."
        )
    elif exterior_nan.any():
        print(
            f"Export: leaving {int(exterior_nan.sum())} exterior cells as NaN "
            f"(FILL_EXPORT_CORNERS=False)."
        )

    if n_interior_nan:
        print(
            f"Export: {n_interior_nan} interior nodata cells left as in the TIFF "
            f"(not filled)."
        )

    z_flat = z_grid.ravel()
    lon_2d, lat_2d = np.meshgrid(lon_vec, lat_vec)
    data = np.column_stack([lon_2d.ravel(), lat_2d.ravel(), z_flat])
    np.savetxt(path, data, fmt="%.8f %.8f %.3f")
    print(
        f"Exported {len(data)} points ({n_lon} lon x {n_lat} lat) to {path}\n"
        f"  In M3D choose menu option 'xyzt file' (not Geotiff)."
    )


# ── load and plot ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print_raster_info(TIF_PATH)

    if SELECT_AOI_ON_MAP:
        quad_lonlat = select_aoi_on_map(TIF_PATH, overview_downsample=OVERVIEW_DOWNSAMPLE)
        x, y, elev, crs, interior_mask = read_dem_quad(TIF_PATH, quad_lonlat)
    else:
        x, y, elev, crs = read_dem(TIF_PATH, aoi=AOI)
        interior_mask = np.ones(elev.shape, dtype=bool)
    x, y, elev = downsample_grid(x, y, elev, DOWNSAMPLE)
    interior_mask = interior_mask[::DOWNSAMPLE, ::DOWNSAMPLE] if DOWNSAMPLE > 1 else interior_mask

    if not RAW_FROM_TIF and FILL_HOLES:
        elev = fill_holes(
            elev,
            inside_mask=interior_mask,
            remove_spikes=REMOVE_SPIKES,
            spike_threshold_m=SPIKE_THRESHOLD_M,
            window=SPIKE_WINDOW,
            max_hole_cells=MAX_HOLE_CELLS,
            max_spike_cells=MAX_SPIKE_CELLS,
            spike_abs_margin_m=SPIKE_ABS_MARGIN_M,
        )
    elif RAW_FROM_TIF:
        print(
            "RAW_FROM_TIF=True — no hole fill / spike removal; "
            "AOI elevations are as in the TIFF."
        )

    lon, lat = map_to_lonlat(x, y, crs)

    if EXPORT_XYZ:
        lon, lat, elev, interior_mask = crop_to_valid_bbox(
            lon, lat, elev, interior_mask
        )
        export_xyz(
            EXPORT_XYZ,
            lon,
            lat,
            elev,
            interior_mask=interior_mask,
            fill_export_corners=FILL_EXPORT_CORNERS,
        )

    nrows, ncols = elev.shape
    cell_x = float(np.nanmedian(np.abs(np.diff(lon[0, :])))) if ncols > 1 else 0.0
    cell_y = float(np.nanmedian(np.abs(np.diff(lat[:, 0])))) if nrows > 1 else 0.0
    print(f"Plot grid: {ncols} x {nrows} pixels  |  cell size: {cell_x:.6f}° x {cell_y:.6f}°")

    valid = np.isfinite(elev)
    if not np.any(valid):
        raise ValueError("No valid elevation values in the raster.")

    z_plot = np.where(valid, elev, np.nanmin(elev[valid]))

    east, north, _, _, local_crs = lonlat_to_local_metres(lon, lat)
    to_wgs84 = Transformer.from_crs(local_crs, "EPSG:4326", always_xy=True)

    # Local metres for true-scale 3D; axis ticks show WGS84 longitude / latitude.
    east0, north0 = np.nanmin(east), np.nanmin(north)
    z0 = np.nanmin(z_plot)
    x_plot = east - east0
    y_plot = north - north0
    z_plot = z_plot - z0
    north_mid = np.nanmean(north)

    x_range = np.nanmax(x_plot) - np.nanmin(x_plot)
    y_range = np.nanmax(y_plot) - np.nanmin(y_plot)
    z_range = np.nanmax(z_plot) - np.nanmin(z_plot)

    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")
    # Matplotlib 3.x defaults to rcount/ccount=50 and silently downsamples large grids.
    # Use every pixel so the 0.5 m AHN resolution is preserved (may be slow if > ~500²).
    surf = ax.plot_surface(
        x_plot,
        y_plot,
        z_plot,
        rcount=nrows,
        ccount=ncols,
        cmap="terrain",
        linewidth=0,
        antialiased=True,
        shade=True,
    )
    fig.colorbar(surf, ax=ax, shrink=0.55, pad=0.08, label="Elevation (m)")

    ax.set_box_aspect((max(x_range, 1e-9), max(y_range, 1e-9), max(z_range, 1e-9)))
    ax.set_xlim(0, x_range)
    ax.set_ylim(0, y_range)
    ax.set_zlim(0, z_range)

    def _fmt_lon(east_plot, _pos):
        lon, _lat = to_wgs84.transform(east_plot + east0, north_mid)
        return f"{lon:.5f}°"

    def _fmt_lat(north_plot, _pos):
        _lon, lat = to_wgs84.transform(east0 + x_range / 2, north_plot + north0)
        return f"{lat:.5f}°"

    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_lon))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_lat))

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_zlabel(f"Elevation (m, +{z0:.1f})")
    ax.set_title(f"3D elevation (WGS84, true scale) — {TIF_PATH}")
    ax.view_init(elev=35, azim=-60)

    plt.tight_layout()
    if OUTPUT:
        fig.savefig(OUTPUT, dpi=200, bbox_inches="tight")
        print(f"Saved: {OUTPUT}")
    plt.show()
