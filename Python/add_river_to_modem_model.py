#!/usr/bin/env python3
"""
Add a river (or other body) to a ModEM 3-D resistivity model by drawing a
polygon on a map that shows the model resistivity with an AHN GeoTIFF overlay.

Edit the settings at the top, then run this file (Run / F5).

1. Left-click to place polygon vertices (3+ points)
2. Right-click or click \"Finish polygon\" to close the outline
3. Set layer range and resistivity, then \"Apply + Save\"

The TIF (EPSG:7415 / RD New) is registered to ModEM local coordinates using
a survey origin estimated from the data file (or set manually below).
"""

from __future__ import annotations

import sys
import traceback
import warnings
from pathlib import Path
from tkinter import (
    BOTH,
    BOTTOM,
    LEFT,
    TOP,
    X,
    BooleanVar,
    Button,
    Checkbutton,
    Entry,
    Frame,
    Label,
    Scale,
    StringVar,
    Tk,
    messagebox,
)

import matplotlib

matplotlib.use("TkAgg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Polygon

from modem_model import ModemModel, load_modem_model, save_modem_model

try:
    import rasterio
    from pyproj import Transformer
except ImportError as exc:
    raise ImportError(
        "Install rasterio and pyproj: pip install rasterio pyproj"
    ) from exc

# --- user settings: edit these, then click Run ---
MASTER_THESIS = Path(r"\\wsl.localhost\Ubuntu\home\ramonvanoli\MasterThesis")
MODEL_PATH = MASTER_THESIS / "Forward/Dike1/2mgrid/FINAL_WITH_RIVER/Modelfile_with_river.model"
OUTPUT_PATH = MASTER_THESIS / "Forward/Dike1/2mgrid/FINAL_WITH_RIVER/Modelfile_with_river_for_plotting.model"
DATA_PATH = MASTER_THESIS / "Forward/Dike1/2mgrid/FINAL_WITH_RIVER/Datafile_raw.data"

TIF_PATH = Path(r"C:\Users\Baar\OneDrive - Vanoli AG\HS_25\Master Thesis\Inversion\3D\Dike1\AHN5_M_102000_418000.TIF")

# Manual survey origin in RD metres (EPSG:7415). If None, estimated from DATA_PATH.
# Model: X = northing = RD_north - ORIGIN_RD_NORTH
#        Y = easting  = RD_east  - ORIGIN_RD_EAST
ORIGIN_RD_EAST: float | None = None
ORIGIN_RD_NORTH: float | None = None

START_LAYER = 1
DEFAULT_LAYER_FROM = 1
DEFAULT_LAYER_TO = 10
DEFAULT_RESISTIVITY = 28.0  # ohm-m for cells inside the river polygon
TIF_OVERLAY_ALPHA = 0.55
TIF_DOWNSAMPLE = 4  # higher = faster overlay (keeps every Nth pixel)
SHOW_TIF_BY_DEFAULT = True
# -------------------------------------------------


def _output_path_for(model_path: Path) -> Path:
    if OUTPUT_PATH is not None:
        return Path(OUTPUT_PATH).resolve()
    return model_path.resolve().with_name(
        f"{model_path.stem}_with_river{model_path.suffix}"
    )


def _jet_like_cmap():
    base = plt.cm.jet(np.linspace(0, 1, 24))
    keep = [0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    return mcolors.ListedColormap(base[keep][::-1])


def estimate_rd_origin_from_data(data_path: Path) -> tuple[float, float]:
    """
    Estimate RD easting/northing of ModEM (0, 0) from station Lat/Lon and X/Y.

    ModEM X = northing [m], Y = easting [m] in the local model frame.
    """
    stations: dict[str, tuple[float, float, float, float]] = {}
    for line in Path(data_path).read_text(errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith(">"):
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            lat = float(parts[2])
            lon = float(parts[3])
            x = float(parts[4])
            y = float(parts[5])
        except ValueError:
            continue
        key = parts[1].split("_")[0]
        stations[key] = (lat, lon, x, y)

    if len(stations) < 2:
        raise ValueError(f"Need station Lat/Lon/X/Y in {data_path}")

    lats = np.array([v[0] for v in stations.values()])
    lons = np.array([v[1] for v in stations.values()])
    xs = np.array([v[2] for v in stations.values()])
    ys = np.array([v[3] for v in stations.values()])

    to_rd = Transformer.from_crs("EPSG:4326", "EPSG:7415", always_xy=True)
    rdx, rdy = to_rd.transform(lons, lats)
    east0 = float(np.mean(np.asarray(rdx) - ys))
    north0 = float(np.mean(np.asarray(rdy) - xs))
    return east0, north0


def load_tif_in_model_km(
    tif_path: Path,
    origin_rd_east: float,
    origin_rd_north: float,
    *,
    downsample: int = 4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (easting_km, northing_km, elev) grids for pcolormesh overlay.

    easting_km / northing_km are 2-D arrays of cell corners/centres matching
    the ModEM plot axes (Y easting, X northing), in kilometres.
    """
    stride = max(1, int(downsample))
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Setting the shape on a NumPy array has been deprecated",
            category=DeprecationWarning,
        )
        with rasterio.open(tif_path) as src:
            elev = np.asarray(src.read(1), dtype=np.float64)
            nodata = src.nodata
            transform = src.transform
            crs = src.crs

    if nodata is not None:
        elev = np.where(elev == nodata, np.nan, elev)
    elev = np.where(elev <= -1e30, np.nan, elev)

    elev = elev[::stride, ::stride]
    # Affine for downsampled pixel centres
    a, b, c, d, e, f = transform.a, transform.b, transform.c, transform.d, transform.e, transform.f
    rows, cols = np.indices(elev.shape)
    rd_east = c + (cols + 0.5) * a * stride + (rows + 0.5) * b * stride
    rd_north = f + (cols + 0.5) * d * stride + (rows + 0.5) * e * stride

    # If TIF is not already RD, reproject (AHN5 is EPSG:7415)
    if crs is not None and crs.to_epsg() not in (7415, 28992):
        to_rd = Transformer.from_crs(crs, "EPSG:7415", always_xy=True)
        rd_east, rd_north = to_rd.transform(rd_east, rd_north)

    model_east_m = np.asarray(rd_east) - origin_rd_east
    model_north_m = np.asarray(rd_north) - origin_rd_north
    easting_km = model_east_m / 1000.0
    northing_km = model_north_m / 1000.0
    return easting_km, northing_km, elev


def _finite_resistivity_bounds(model: ModemModel) -> tuple[float, float]:
    finite = model.rho[np.isfinite(model.rho) & (model.rho > 0) & (model.rho < 1e15)]
    if finite.size == 0:
        return 1.0, 1000.0
    return float(np.min(finite)), float(np.max(finite))


def plot_layer(ax, model, layer, vmin, vmax, mesh, title):
    data = model.layer_slice(layer)
    log_data = np.log10(data)

    y_km = model.y_edges / 1000.0
    x_km = model.x_edges / 1000.0
    yy, xx = np.meshgrid(y_km, x_km)

    if mesh is None:
        mesh = ax.pcolormesh(
            yy,
            xx,
            log_data,
            shading="flat",
            cmap=_jet_like_cmap(),
            norm=plt.Normalize(vmin=np.log10(vmin), vmax=np.log10(vmax)),
            zorder=1,
        )
    else:
        mesh.set_array(log_data.ravel())
        mesh.set_clim(np.log10(vmin), np.log10(vmax))

    z0, z1 = model.layer_depth_range_km(layer)
    title.set_text(
        f"Layer {layer + 1}/{model.nz} | depth {z0:.3f} to {z1:.3f} km"
    )
    ax.set_xlabel("Easting Y (km)")
    ax.set_ylabel("Northing X (km)")
    ax.set_aspect("equal", adjustable="box")
    return mesh


def interactive_view(
    model: ModemModel,
    *,
    start_layer: int = 1,
    output_path: Path | None = None,
    tif_easting_km: np.ndarray | None = None,
    tif_northing_km: np.ndarray | None = None,
    tif_elev: np.ndarray | None = None,
) -> None:
    start_layer = int(np.clip(start_layer, 1, model.nz))
    layer_idx = start_layer - 1
    vmin, vmax = _finite_resistivity_bounds(model)

    root = Tk()
    root.title("Add river to ModEM model")
    root.minsize(950, 750)
    root.geometry("1100x900")

    panel = Frame(root, padx=10, pady=8, relief="ridge", bd=1)
    panel.pack(side=BOTTOM, fill=X)

    fig = Figure(figsize=(9, 6.5), dpi=100)
    ax = fig.add_subplot(111)
    title = ax.set_title("")
    mesh = plot_layer(ax, model, layer_idx, vmin, vmax, None, title)

    tif_mesh = None
    if tif_easting_km is not None and tif_elev is not None:
        elev_plot = np.ma.masked_invalid(tif_elev)
        tif_mesh = ax.pcolormesh(
            tif_easting_km,
            tif_northing_km,
            elev_plot,
            shading="auto",
            cmap="terrain",
            alpha=TIF_OVERLAY_ALPHA if SHOW_TIF_BY_DEFAULT else 0.0,
            zorder=2,
        )

    # Clip view to model extent (TIF tile is larger)
    ax.set_xlim(model.y_edges[0] / 1000.0, model.y_edges[-1] / 1000.0)
    ax.set_ylim(model.x_edges[0] / 1000.0, model.x_edges[-1] / 1000.0)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    tick_vals = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
    tick_vals = [v for v in tick_vals if vmin <= v <= vmax] or [vmin, vmax]
    cbar.set_ticks(np.log10(tick_vals))
    cbar.set_ticklabels([str(v) for v in tick_vals])
    cbar.set_label("Resistivity (ohm-m)")
    fig.suptitle(
        f"{model.nx}x{model.ny}x{model.nz}  |  polygon river edit  |  "
        f"saves to: {output_path.name}"
    )

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True, padx=8, pady=8)

    status_var = StringVar(
        value="Step 1: left-click polygon vertices (≥3). Right-click or Finish polygon to close."
    )
    Label(
        panel,
        textvariable=status_var,
        anchor="w",
        justify="left",
        bg="#eef3ff",
        fg="#111",
        font=("Segoe UI", 10),
        padx=8,
        pady=6,
        relief="groove",
    ).pack(fill=X, pady=(0, 6))

    slider_row = Frame(panel)
    slider_row.pack(fill=X, pady=4)

    layer_label_var = StringVar(value=f"View layer: {start_layer} / {model.nz}")
    Label(
        slider_row, textvariable=layer_label_var, font=("Segoe UI", 10), width=18, anchor="w"
    ).pack(side=LEFT, padx=(0, 8))
    layer_scale = Scale(
        slider_row,
        from_=1,
        to=model.nz,
        orient="horizontal",
        length=480,
        resolution=1,
        showvalue=True,
    )
    layer_scale.set(start_layer)
    layer_scale.pack(side=LEFT, fill=X, expand=True)

    show_tif_var = BooleanVar(value=SHOW_TIF_BY_DEFAULT and tif_mesh is not None)
    Checkbutton(
        slider_row,
        text="Show TIF overlay",
        variable=show_tif_var,
        font=("Segoe UI", 10),
    ).pack(side=LEFT, padx=(12, 0))

    edit_row = Frame(panel)
    edit_row.pack(fill=X, pady=4)

    Label(edit_row, text="Edit layers — from:", font=("Segoe UI", 10)).pack(side=LEFT)
    entry_from = Entry(edit_row, width=6, font=("Segoe UI", 11), relief="solid", bd=1)
    entry_from.insert(0, str(DEFAULT_LAYER_FROM))
    entry_from.pack(side=LEFT, padx=(6, 16))

    Label(edit_row, text="to:", font=("Segoe UI", 10)).pack(side=LEFT)
    entry_to = Entry(edit_row, width=6, font=("Segoe UI", 11), relief="solid", bd=1)
    entry_to.insert(0, str(DEFAULT_LAYER_TO))
    entry_to.pack(side=LEFT, padx=(6, 16))

    Label(edit_row, text="Rho (ohm-m):", font=("Segoe UI", 10)).pack(side=LEFT)
    entry_rho = Entry(edit_row, width=10, font=("Segoe UI", 11), relief="solid", bd=1)
    entry_rho.insert(0, str(DEFAULT_RESISTIVITY))
    entry_rho.pack(side=LEFT, padx=(6, 0))

    state = {
        "points": [],
        "markers": [],
        "polygon": None,
        "closed": False,
        "dirty": False,
    }

    def set_status(msg: str) -> None:
        status_var.set(msg)
        print(msg)
        root.update_idletasks()

    def redraw_current_layer() -> None:
        plot_layer(ax, model, layer_idx, vmin, vmax, mesh, title)
        canvas.draw_idle()

    def on_layer_scale(_value) -> None:
        nonlocal layer_idx
        layer_idx = int(float(layer_scale.get())) - 1
        layer_label_var.set(f"View layer: {layer_idx + 1} / {model.nz}")
        redraw_current_layer()

    layer_scale.configure(command=on_layer_scale)

    def on_tif_toggle() -> None:
        if tif_mesh is None:
            return
        tif_mesh.set_alpha(TIF_OVERLAY_ALPHA if show_tif_var.get() else 0.0)
        canvas.draw_idle()

    show_tif_var.trace_add("write", lambda *_: on_tif_toggle())

    def clear_picks() -> None:
        state["points"] = []
        state["closed"] = False
        for marker in state["markers"]:
            marker.remove()
        state["markers"] = []
        if state["polygon"] is not None:
            state["polygon"].remove()
            state["polygon"] = None
        canvas.draw_idle()

    def update_polygon_outline(*, closed: bool | None = None) -> None:
        if closed is not None:
            state["closed"] = closed
        if state["polygon"] is not None:
            state["polygon"].remove()
            state["polygon"] = None
        if len(state["points"]) < 2:
            return
        pts = state["points"]
        is_closed = state["closed"] and len(pts) >= 3
        outline = pts + [pts[0]] if is_closed else pts
        state["polygon"] = ax.add_patch(
            Polygon(
                outline,
                closed=False,
                linewidth=2.2,
                edgecolor="cyan",
                facecolor="cyan",
                alpha=0.25 if is_closed else 0.12,
                zorder=5,
            )
        )

    def finish_polygon() -> None:
        if len(state["points"]) < 3:
            messagebox.showwarning(
                "Need more points",
                f"A polygon needs at least 3 vertices.\nYou have {len(state['points'])}.",
            )
            return
        update_polygon_outline(closed=True)
        canvas.draw_idle()
        set_status(
            f"Polygon closed ({len(state['points'])} vertices). "
            "Edit layers/rho, then click Apply + Save."
        )

    def on_map_click(event) -> None:
        if event.inaxes is not ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        # Right-click finishes the polygon
        if event.button == 3:
            finish_polygon()
            return

        if event.button != 1:
            return

        # Starting a new polygon after a closed one
        if state["closed"]:
            clear_picks()

        y_km, x_km = float(event.xdata), float(event.ydata)
        state["points"].append((y_km, x_km))
        marker = ax.plot(
            y_km, x_km, "wo", markeredgecolor="black", markersize=7, zorder=6
        )[0]
        state["markers"].append(marker)
        update_polygon_outline(closed=False)
        canvas.draw_idle()

        n = len(state["points"])
        if n < 3:
            set_status(f"Point {n} placed. Need at least {3 - n} more, then Finish.")
        else:
            set_status(
                f"{n} vertices. Right-click or Finish polygon to close, "
                "or keep adding points."
            )

    canvas.mpl_connect("button_press_event", on_map_click)

    def parse_int(text: str, name: str) -> int:
        value = int(float(text.strip()))
        if value < 1 or value > model.nz:
            raise ValueError(f"{name} must be between 1 and {model.nz}")
        return value

    def on_apply() -> None:
        set_status("Working: apply + save ...")
        try:
            if not state["closed"] or len(state["points"]) < 3:
                messagebox.showwarning(
                    "Finish the polygon",
                    "Place ≥3 points and click Finish polygon (or right-click) first.",
                )
                set_status("Finish the polygon before Apply + Save.")
                return

            layer_from = parse_int(entry_from.get(), "Layer from")
            layer_to = parse_int(entry_to.get(), "Layer to")
            resistivity = float(entry_rho.get().strip())
            if resistivity <= 0:
                raise ValueError("Resistivity must be positive")

            vertices = np.array(state["points"], dtype=float)
            in_polygon = int(model.mask_cells_in_polygon_km(vertices).sum())

            changed = model.apply_resistivity_polygon(
                vertices, layer_from, layer_to, resistivity
            )
            if changed == 0:
                messagebox.showwarning(
                    "No cells updated",
                    f"{in_polygon} cells inside your polygon, but layers "
                    f"{layer_from}-{layer_to} are air there.\n\n"
                    "Try a different layer range.",
                )
                set_status("No earth cells in that layer range.")
                return

            state["dirty"] = True
            redraw_current_layer()

            set_status(f"Saving to {output_path} ...")
            root.update_idletasks()
            saved = save_modem_model(model, output_path)
            state["dirty"] = False

            msg = (
                f"Updated {changed} cells inside polygon ({in_polygon} columns).\n"
                f"Rho = {resistivity} ohm-m, layers {layer_from}-{layer_to}.\n\n"
                f"Saved to:\n{saved}"
            )
            set_status(f"Done. Saved {changed} cells to {saved.name}")
            messagebox.showinfo("Success", msg)
        except Exception as exc:
            set_status(f"Error: {exc}")
            messagebox.showerror("Failed", f"{exc}\n\n{traceback.format_exc()}")

    def on_clear() -> None:
        clear_picks()
        set_status("Points cleared. Left-click new polygon vertices.")

    btn_row = Frame(panel)
    btn_row.pack(fill=X, pady=(8, 0))

    Button(
        btn_row,
        text="Finish polygon",
        command=finish_polygon,
        font=("Segoe UI", 10, "bold"),
        bg="#1f6feb",
        fg="white",
        activebackground="#1558c0",
        activeforeground="white",
        padx=12,
        pady=8,
        cursor="hand2",
    ).pack(side=LEFT, padx=(0, 8))

    Button(
        btn_row,
        text="Apply + Save",
        command=on_apply,
        font=("Segoe UI", 11, "bold"),
        bg="#2e9b3c",
        fg="white",
        activebackground="#248a31",
        activeforeground="white",
        padx=16,
        pady=8,
        cursor="hand2",
    ).pack(side=LEFT, padx=(0, 8))

    Button(
        btn_row,
        text="Clear points",
        command=on_clear,
        font=("Segoe UI", 10),
        padx=12,
        pady=8,
    ).pack(side=LEFT)

    def on_close() -> None:
        if state["dirty"]:
            try:
                save_modem_model(model, output_path)
                print(f"Saved on close: {output_path}")
            except Exception as exc:
                print(f"Save on close failed: {exc}", file=sys.stderr)
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.bind("<Return>", lambda _e: on_apply())
    root.mainloop()


def main() -> int:
    model_path = Path(MODEL_PATH)
    if not model_path.is_file():
        print(f"Model file not found: {model_path}", file=sys.stderr)
        return 1
    if not TIF_PATH.is_file():
        print(f"TIF file not found: {TIF_PATH}", file=sys.stderr)
        return 1

    if ORIGIN_RD_EAST is not None and ORIGIN_RD_NORTH is not None:
        east0, north0 = float(ORIGIN_RD_EAST), float(ORIGIN_RD_NORTH)
        print(f"Using manual RD origin: east={east0:.3f}, north={north0:.3f}")
    elif DATA_PATH is not None and Path(DATA_PATH).is_file():
        east0, north0 = estimate_rd_origin_from_data(Path(DATA_PATH))
        print(
            f"Estimated RD origin from {Path(DATA_PATH).name}: "
            f"east={east0:.3f}, north={north0:.3f}"
        )
    else:
        print(
            "Set ORIGIN_RD_EAST / ORIGIN_RD_NORTH or provide DATA_PATH "
            "with station Lat/Lon.",
            file=sys.stderr,
        )
        return 1

    print(f"Loading TIF overlay: {TIF_PATH}")
    tif_e, tif_n, tif_z = load_tif_in_model_km(
        TIF_PATH, east0, north0, downsample=TIF_DOWNSAMPLE
    )
    print(
        f"  TIF in model km — easting [{np.nanmin(tif_e):.4f}, {np.nanmax(tif_e):.4f}], "
        f"northing [{np.nanmin(tif_n):.4f}, {np.nanmax(tif_n):.4f}]"
    )

    output_path = _output_path_for(model_path)
    model = load_modem_model(model_path)
    print(f"Loaded: {model_path.resolve()}")
    print(f"Output: {output_path}")
    interactive_view(
        model,
        start_layer=START_LAYER,
        output_path=output_path,
        tif_easting_km=tif_e,
        tif_northing_km=tif_n,
        tif_elev=tif_z,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
