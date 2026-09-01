import os
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# This script plots topography (from TIF_to_DEM.py xyz export) with the MT survey layout.
# Each survey line has its own coordinates file in this format:
#   ID         Latitude              Longitude            Elevation
#   1201   +51.7504280933333    +4.63354051833333    1.225

plt.style.use("ggplot")

# ── user settings ─────────────────────────────────────────────────────────────
XYZ_PATH = r"C:/Users/Baar/OneDrive - Vanoli AG/HS_25/Master Thesis/Inversion/3D/Dike3/topography_dike1_forfigure.xyz"

# One coordinates file per survey line: (line name, path)
SURVEY_LINES = [
    ("Line 2", r"C:\Users\Baar\OneDrive - Vanoli AG\HS_25\Master Thesis\Data_Processed\EDI_standard\Dike1\Line2\coordinates_D1_L2.txt"),
    ("Line 3", r"C:\Users\Baar\OneDrive - Vanoli AG\HS_25\Master Thesis\Data_Processed\EDI_standard\Dike1\Line3\coordinates_D1_L3.txt"),
    ("Line 4", r"C:\Users\Baar\OneDrive - Vanoli AG\HS_25\Master Thesis\Data_Processed\EDI_standard\Dike1\Line4\coordinates_D1_L4.txt"),


]

PLOT_DOWNSAMPLE = 1  # use 2–4 if the plot feels slow
LAT_MAX = 51.75 + 0.00125  # northern map limit (°)
CONNECT_STATIONS = True  # draw a line through stations in file order
LABEL_STATIONS = True    # annotate station IDs

# Publication-style text sizes
FONT_TITLE = 16
FONT_LABEL = 14
FONT_TICK = 12
FONT_LEGEND = 12
FONT_STATION = 11
FONT_CBAR = 14
SAVE_DPI = 600  # high-quality PNG export
# ───────────────────────────────────────────────────────────────────────────────

LINE_COLORS = ["crimson", "royalblue", "darkorange", "seagreen", "purple", "saddlebrown"]


def _output_path_for(xyz_path: str | Path, ext: str = "svg") -> Path:
    xyz = Path(xyz_path)
    return xyz.with_name(f"{xyz.stem}_surveylayout.{ext}")


def read_xyz(path, downsample=1):
    """Read lon lat elev grid written by TIF_to_DEM.py."""
    size_mb = os.path.getsize(path) / 1e6
    print(f"Loading topography: {path}")
    print(f"  file size: {size_mb:.1f} MB")
    if size_mb > 500:
        raise ValueError(
            f"Topography file is {size_mb:.0f} MB — too large to plot. "
            "Re-export with TIF_to_DEM.py using an AOI and DOWNSAMPLE "
            "(expected size is well under 10 MB for your survey area)."
        )

    t0 = time.time()
    data = np.loadtxt(path)
    print(f"  loaded {len(data):,} points in {time.time() - t0:.1f} s")
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"Expected 3 columns (lon lat elev) in {path}")

    lon_col = data[:, 0]
    row_breaks = np.where(np.diff(lon_col) < 0)[0]
    if len(row_breaks) == 0:
        raise ValueError(
            f"XYZ file is not a regular grid. Re-export with TIF_to_DEM.py."
        )
    n_lon = int(row_breaks[0] + 1)
    if len(data) % n_lon != 0:
        raise ValueError(
            f"XYZ file is not a regular grid. Re-export with TIF_to_DEM.py."
        )
    n_lat = len(data) // n_lon
    lon_vec = data[:n_lon, 0]
    lat_vec = data[::n_lon, 1]
    elev_2d = data[:, 2].reshape(n_lat, n_lon)
    if downsample > 1:
        elev_2d = elev_2d[::downsample, ::downsample]
        lon_vec = lon_vec[::downsample]
        lat_vec = lat_vec[::downsample]
        print(f"  plot grid after downsample: {len(lon_vec)} x {len(lat_vec)}")

    lon_2d, lat_2d = np.meshgrid(lon_vec, lat_vec)
    return lon_2d, lat_2d, elev_2d


def read_coords(path):
    """Read station ID, latitude, longitude, elevation (whitespace-separated, header)."""
    stations = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.lower().startswith("id"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            stations.append(
                {
                    "id": parts[0],
                    "lat": float(parts[1]),
                    "lon": float(parts[2]),
                    "elev": float(parts[3]),
                }
            )
    if not stations:
        raise ValueError(f"No stations found in {path}")
    return stations


def plot_topography_survey(
    xyz_path,
    survey_lines,
    output=None,
    connect_stations=True,
    label_stations=True,
    downsample=1,
):
    lon, lat, elev = read_xyz(xyz_path, downsample=downsample)

    lat_max = LAT_MAX
    row_mask = lat[:, 0] <= lat_max
    if not np.any(row_mask):
        raise ValueError(f"No topography rows at or below latitude {lat_max:.5f}°")
    lon = lon[row_mask, :]
    lat = lat[row_mask, :]
    elev = elev[row_mask, :]

    if not np.any(np.isfinite(elev)):
        raise ValueError("No valid elevation values in the topography file.")

    print("Drawing map...")
    fig, ax = plt.subplots(figsize=(11, 8))
    levels = np.linspace(np.nanmin(elev), np.nanmax(elev), 25)
    cf = ax.contourf(lon, lat, elev, levels=levels, cmap="terrain", extend="both")
    ax.contour(lon, lat, elev, levels=levels[::2], colors="k", linewidths=0.3, alpha=0.35)

    for i, (line_name, coords_path) in enumerate(survey_lines):
        stations = read_coords(coords_path)
        color = LINE_COLORS[i % len(LINE_COLORS)]

        st_lon = np.array([s["lon"] for s in stations])
        st_lat = np.array([s["lat"] for s in stations])
        st_ids = [s["id"] for s in stations]

        if connect_stations:
            ax.plot(
                st_lon,
                st_lat,
                "-o",
                color="white",
                linewidth=2.5,
                markersize=0,
                zorder=3 + 2 * i,
            )
            ax.plot(
                st_lon,
                st_lat,
                "-o",
                color=color,
                linewidth=1.4,
                markersize=5,
                label=line_name,
                zorder=4 + 2 * i,
            )
        else:
            ax.scatter(
                st_lon,
                st_lat,
                s=55,
                c=color,
                edgecolors="black",
                linewidths=0.8,
                label=line_name,
                zorder=5 + i,
            )

        if label_stations:
            for sid, x, y in zip(st_ids, st_lon, st_lat):
                ax.annotate(
                    sid,
                    (x, y),
                    textcoords="offset points",
                    xytext=(4, 4),
                    fontsize=FONT_STATION,
                    color="black",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75),
                    zorder=10 + i,
                )

    ax.set_xlabel("Longitude (°)", fontsize=FONT_LABEL)
    ax.set_ylabel("Latitude (°)", fontsize=FONT_LABEL)
    ax.set_title("Topography and MT survey layout", fontsize=FONT_TITLE)
    ax.set_ylim(lat.min(), lat_max)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(axis="both", labelsize=FONT_TICK)

    lat_fmt = ScalarFormatter(useOffset=False)
    lat_fmt.set_scientific(False)
    ax.yaxis.set_major_formatter(lat_fmt)
    ax.legend(loc="upper right", fontsize=FONT_LEGEND)

    cbar = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label("Elevation (m)", fontsize=FONT_CBAR)
    cbar.ax.tick_params(labelsize=FONT_TICK)

    plt.tight_layout()

    out_svg = Path(output) if output else _output_path_for(xyz_path, "svg")
    out_png = out_svg.with_suffix(".png")

    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    fig.savefig(
        out_png,
        format="png",
        dpi=SAVE_DPI,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
        edgecolor="none",
    )
    print(f"Saved SVG: {out_svg.resolve()}")
    print(f"Saved PNG: {out_png.resolve()} ({SAVE_DPI} dpi)")
    print("Opening plot window...")
    plt.show()


if __name__ == "__main__":
    plot_topography_survey(
        XYZ_PATH,
        SURVEY_LINES,
        connect_stations=CONNECT_STATIONS,
        label_stations=LABEL_STATIONS,
        downsample=PLOT_DOWNSAMPLE,
    )
