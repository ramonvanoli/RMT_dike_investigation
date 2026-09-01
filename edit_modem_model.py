#!/usr/bin/env python3
"""
This script is used to remove topography from ModEM 3D resistivity models.
It allows the user to click 4 points and the topography outside the polygon is flattened to the target elevation.
Interactive ModEM model editor: clip topography to a clicked rectangle.

Edit the settings at the top, then run this file (Run / F5).

1. Click 4 corners on the map (same coordinates as view_modem_model.py)
2. Everything OUTSIDE the polygon is flattened to the target elevation
3. Click the green "Apply + Save" button
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path
from tkinter import (
    BOTH,
    BOTTOM,
    LEFT,
    TOP,
    X,
    Button,
    Entry,
    Frame,
    Label,
    OptionMenu,
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

# --- user settings: edit these, then click Run ---
MASTER_THESIS = Path(r"\\wsl.localhost\Ubuntu\home\ramonvanoli\MasterThesis")
MODEL_PATH = MASTER_THESIS / "Data/Dike1/FINAL2/Modelfile.model"
OUTPUT_PATH = MASTER_THESIS / "Data/Dike1/FINAL2/Modelfile_CUT.model"
START_LAYER = 1
DEFAULT_TARGET_ELEVATION_M = 0.0
NUM_PICK_POINTS = 4
# -------------------------------------------------


def _output_path_for(model_path: Path) -> Path:
    if OUTPUT_PATH is not None:
        return Path(OUTPUT_PATH).resolve()
    return model_path.resolve().with_name(f"{model_path.stem}_flat_topo{model_path.suffix}")


def _jet_like_cmap():
    base = plt.cm.jet(np.linspace(0, 1, 24))
    keep = [0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    return mcolors.ListedColormap(base[keep][::-1])


def _finite_resistivity_bounds(model: ModemModel) -> tuple[float, float]:
    finite = model.rho[np.isfinite(model.rho) & (model.rho > 0) & (model.rho < 1e15)]
    if finite.size == 0:
        return 1.0, 1000.0
    return float(np.min(finite)), float(np.max(finite))


def plot_resistivity_layer(ax, model, layer, vmin, vmax, mesh, title):
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
        )
    else:
        mesh.set_array(log_data.ravel())
        mesh.set_clim(np.log10(vmin), np.log10(vmax))

    z0, z1 = model.layer_depth_range_km(layer)
    title.set_text(
        f"Resistivity — layer {layer + 1}/{model.nz} | depth {z0:.3f} to {z1:.3f} km"
    )
    ax.set_xlabel("Easting Y (km)")
    ax.set_ylabel("Northing X (km)")
    ax.set_aspect("equal", adjustable="box")
    return mesh


def plot_topography(ax, model, mesh, title, keep_mask=None):
    elev = model.topography_surface_elevation_m()
    y_km = model.y_edges / 1000.0
    x_km = model.x_edges / 1000.0
    yy, xx = np.meshgrid(y_km, x_km)

    vmin = float(np.nanmin(elev))
    vmax = float(np.nanmax(elev))
    if vmin == vmax:
        vmax = vmin + 1.0

    if mesh is None:
        mesh = ax.pcolormesh(
            yy,
            xx,
            elev,
            shading="flat",
            cmap="terrain",
            vmin=vmin,
            vmax=vmax,
        )
    else:
        mesh.set_array(elev.ravel())
        mesh.set_clim(vmin, vmax)

    title.set_text(
        f"Topography surface | {vmin:.2f} to {vmax:.2f} m a.s.l."
    )
    ax.set_xlabel("Easting Y (km)")
    ax.set_ylabel("Northing X (km)")
    ax.set_aspect("equal", adjustable="box")

    if keep_mask is not None and keep_mask.any():
        outside = np.ma.masked_where(keep_mask, elev)
        ax.pcolormesh(
            yy,
            xx,
            outside,
            shading="flat",
            cmap=mcolors.ListedColormap([(1.0, 0.2, 0.2, 0.35)]),
            vmin=0,
            vmax=1,
        )

    return mesh


def interactive_edit(
    model: ModemModel,
    start_layer: int = 1,
    output_path: Path | None = None,
) -> None:
    start_layer = int(np.clip(start_layer, 1, model.nz))
    layer_idx = start_layer - 1
    vmin, vmax = _finite_resistivity_bounds(model)

    root = Tk()
    root.title("ModEM topography clip editor")
    root.minsize(900, 700)
    root.geometry("1050x820")

    panel = Frame(root, padx=10, pady=8, relief="ridge", bd=1)
    panel.pack(side=BOTTOM, fill=X)

    fig = Figure(figsize=(9, 6), dpi=100)
    ax = fig.add_subplot(111)
    title = ax.set_title("")
    view_mode = StringVar(value="Topography")
    plot_state = {"mesh": None, "cbar": None}
    plot_state["mesh"] = plot_topography(ax, model, None, title)
    plot_state["cbar"] = fig.colorbar(plot_state["mesh"], ax=ax, pad=0.02)
    plot_state["cbar"].set_label("Elevation (m a.s.l.)")
    fig.suptitle(
        f"{model.nx}x{model.ny}x{model.nz} cells  |  "
        f"outside polygon -> flat {DEFAULT_TARGET_ELEVATION_M:g} m  |  "
        f"saves to: {output_path.name}"
    )

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True, padx=8, pady=8)

    status_var = StringVar(
        value=f"Step 1: click {NUM_PICK_POINTS} corners around the area to KEEP."
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

    view_row = Frame(panel)
    view_row.pack(fill=X, pady=4)
    Label(view_row, text="Display:", font=("Segoe UI", 10)).pack(side=LEFT)
    OptionMenu(view_row, view_mode, "Topography", "Resistivity").pack(side=LEFT, padx=8)

    slider_row = Frame(panel)
    slider_row.pack(fill=X, pady=4)
    layer_label_var = StringVar(value=f"View layer: {start_layer} / {model.nz}")
    Label(
        slider_row,
        textvariable=layer_label_var,
        font=("Segoe UI", 10),
        width=18,
        anchor="w",
    ).pack(side=LEFT, padx=(0, 8))
    layer_scale = Scale(
        slider_row,
        from_=1,
        to=model.nz,
        orient="horizontal",
        length=500,
        resolution=1,
        showvalue=True,
    )
    layer_scale.set(start_layer)
    layer_scale.pack(side=LEFT, fill=X, expand=True)

    edit_row = Frame(panel)
    edit_row.pack(fill=X, pady=4)
    Label(
        edit_row,
        text="Flatten outside polygon to elevation (m a.s.l.):",
        font=("Segoe UI", 10),
    ).pack(side=LEFT)
    entry_elev = Entry(edit_row, width=10, font=("Segoe UI", 11), relief="solid", bd=1)
    entry_elev.insert(0, str(DEFAULT_TARGET_ELEVATION_M))
    entry_elev.pack(side=LEFT, padx=(6, 0))

    state = {
        "points": [],
        "polygon": None,
        "dirty": False,
        "keep_mask": None,
    }

    def set_status(msg: str) -> None:
        status_var.set(msg)
        print(msg)
        root.update_idletasks()

    def current_keep_mask():
        if len(state["points"]) < NUM_PICK_POINTS:
            return None
        return model.mask_cells_in_polygon_km(np.array(state["points"], dtype=float))

    def redraw() -> None:
        mode = view_mode.get()
        keep = current_keep_mask()
        state["keep_mask"] = keep

        if plot_state["cbar"] is not None:
            plot_state["cbar"].remove()
            plot_state["cbar"] = None

        ax.clear()
        title_obj = ax.set_title("")
        if mode == "Resistivity":
            plot_state["mesh"] = plot_resistivity_layer(
                ax, model, layer_idx, vmin, vmax, None, title_obj
            )
            plot_state["cbar"] = fig.colorbar(plot_state["mesh"], ax=ax, pad=0.02)
            tick_vals = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
            tick_vals = [v for v in tick_vals if vmin <= v <= vmax] or [vmin, vmax]
            plot_state["cbar"].set_ticks(np.log10(tick_vals))
            plot_state["cbar"].set_ticklabels([str(v) for v in tick_vals])
            plot_state["cbar"].set_label("Resistivity (ohm-m)")
        else:
            plot_state["mesh"] = plot_topography(
                ax, model, None, title_obj, keep_mask=keep
            )
            plot_state["cbar"] = fig.colorbar(plot_state["mesh"], ax=ax, pad=0.02)
            plot_state["cbar"].set_label("Elevation (m a.s.l.)")

        if keep is not None and len(state["points"]) >= NUM_PICK_POINTS:
            pts = state["points"] + [state["points"][0]]
            state["polygon"] = ax.add_patch(
                Polygon(
                    pts,
                    closed=True,
                    linewidth=2.5,
                    edgecolor="lime",
                    facecolor="lime",
                    alpha=0.18,
                )
            )

        for y_km, x_km in state["points"]:
            ax.plot(
                y_km, x_km, "wo", markeredgecolor="black", markersize=8, zorder=5
            )
        canvas.draw_idle()

    def on_view_change(*_args) -> None:
        redraw()

    view_mode.trace_add("write", on_view_change)

    def on_layer_scale(_value) -> None:
        nonlocal layer_idx
        layer_idx = int(float(layer_scale.get())) - 1
        layer_label_var.set(f"View layer: {layer_idx + 1} / {model.nz}")
        if view_mode.get() == "Resistivity":
            redraw()

    layer_scale.configure(command=on_layer_scale)

    def clear_picks() -> None:
        state["points"] = []
        state["keep_mask"] = None
        redraw()

    def on_map_click(event) -> None:
        if event.inaxes is not ax or event.button != 1:
            return
        if event.xdata is None or event.ydata is None:
            return

        if len(state["points"]) >= NUM_PICK_POINTS:
            clear_picks()

        y_km, x_km = float(event.xdata), float(event.ydata)
        state["points"].append((y_km, x_km))

        n = len(state["points"])
        if n < NUM_PICK_POINTS:
            set_status(f"Step 1: point {n}/{NUM_PICK_POINTS} placed. Click next corner.")
        else:
            keep = int(current_keep_mask().sum())
            outside = model.nx * model.ny - keep
            set_status(
                f"Step 2: keep {keep} columns inside, flatten {outside} outside. "
                "Click Apply + Save."
            )
        redraw()

    canvas.mpl_connect("button_press_event", on_map_click)

    def on_apply() -> None:
        set_status("Working: flatten topography outside polygon ...")
        try:
            if len(state["points"]) < NUM_PICK_POINTS:
                messagebox.showwarning(
                    "Need 4 points",
                    f"Click {NUM_PICK_POINTS} corners on the map first.\n"
                    f"You have placed {len(state['points'])}.",
                )
                set_status(
                    f"Need {NUM_PICK_POINTS} points on the map (have {len(state['points'])})."
                )
                return

            target_elev = float(entry_elev.get().strip())
            vertices = np.array(state["points"], dtype=float)
            keep = int(model.mask_cells_in_polygon_km(vertices).sum())
            outside = model.nx * model.ny - keep

            changed = model.flatten_topography_outside_polygon(
                vertices,
                target_elevation_m=target_elev,
            )
            if changed == 0:
                messagebox.showinfo(
                    "No changes",
                    "Topography outside the polygon already matches the target elevation.",
                )
                set_status("No cells needed updating.")
                return

            state["dirty"] = True
            redraw()

            set_status(f"Saving to {output_path} ...")
            root.update_idletasks()
            saved = save_modem_model(model, output_path)
            state["dirty"] = False

            elev = model.topography_surface_elevation_m()
            inside = model.mask_cells_in_polygon_km(vertices)
            msg = (
                f"Flattened topography outside polygon to {target_elev:g} m a.s.l.\n"
                f"Updated {changed} cells.\n"
                f"Kept {keep} columns inside, flattened {outside} outside.\n"
                f"Inside elevation range: "
                f"{elev[inside].min():.2f} to {elev[inside].max():.2f} m\n"
                f"Outside elevation range: "
                f"{elev[~inside].min():.2f} to {elev[~inside].max():.2f} m\n\n"
                f"Saved to:\n{saved}"
            )
            set_status(f"Done. Saved {changed} cells to {saved.name}")
            messagebox.showinfo("Success", msg)
        except Exception as exc:
            set_status(f"Error: {exc}")
            messagebox.showerror("Failed", f"{exc}\n\n{traceback.format_exc()}")

    def on_clear() -> None:
        clear_picks()
        set_status(f"Points cleared. Click {NUM_PICK_POINTS} corners on the map.")

    btn_row = Frame(panel)
    btn_row.pack(fill=X, pady=(8, 0))

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

    output_path = _output_path_for(model_path)
    model = load_modem_model(model_path)
    print(f"Loaded: {model_path.resolve()}")
    print(f"Output: {output_path}")
    interactive_edit(model, start_layer=START_LAYER, output_path=output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
