from pathlib import Path
from datetime import datetime
import csv
import json
import sys
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import rasterio
from rasterio.features import shapes
from rasterio.plot import plotting_extent
from rasterio.transform import Affine
from rasterio.warp import transform_geom
from rasterio.windows import from_bounds
from scipy import ndimage as ndi
from shapely.affinity import affine_transform
from shapely.geometry import box, mapping, shape
from shapely.ops import unary_union
from shapely.prepared import prep

WORKING_DIR = Path.cwd()
DATA_DIR = WORKING_DIR.parent / "Data"
OUTPUT_DIR = Path.home() / "YPL_outputs" / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
VERTICAL_SCALE = 1.0
GROUND_PERCENTILE = 10
MIN_HEIGHT = 2.0
GREEN_LIMIT = 0.03
MIN_FACE_AREA = 3.0
DIRECTION_BIN = 45.0
EDGE_DISTANCE = 0.6
MAIN_FACE_AREA = 10.0
MIN_MAIN_SLOPE = 25.0
MAX_MAIN_SLOPE = 55.0
HEIGHT_DEVIATION = 0.3
BLUE_CONTRAST = 0.26
PANEL_LENGTH = 1.762
PANEL_WIDTH = 1.134
PANEL_GAP = 0.02
EDGE_CLEARANCE = 0.35
OBSTRUCTION_CLEARANCE = 0.35
MIN_ARRAY_PANELS = 2
print("Python:", sys.executable)
print("Working folder:", WORKING_DIR)
print("Output folder:", OUTPUT_DIR)

with rasterio.open(DATA_DIR / "rgb.tif") as source:
    rgb = source.read([1, 2, 3]).transpose(1, 2, 0)
    rgb_crs, rgb_transform = (source.crs, source.transform)
with rasterio.open(DATA_DIR / "dsm_full.tif") as source:
    dsm = source.read(1, masked=True).filled(np.nan).astype(float) * VERTICAL_SCALE
    profile = source.profile.copy()
    extent = plotting_extent(source)
transform = profile["transform"]
pixel_size = transform.a
pixel_area = abs(transform.a * transform.e)
if not (
    rgb_crs == profile["crs"]
    and rgb.shape[:2] == dsm.shape
    and rgb_transform.almost_equals(transform)
    and (transform.a > 0)
    and (transform.e < 0)
    and (transform.b == 0)
    and (transform.d == 0)
    and np.isclose(transform.a, -transform.e)
):
    raise ValueError("Expected aligned, north-up rasters with square pixels.")
with rasterio.open(DATA_DIR / "dsm_clipped.tif") as source:
    clip = source.read(1, masked=True)
    if source.crs != profile["crs"] or not np.allclose(source.res, (pixel_size, pixel_size)):
        raise ValueError("The clipped DSM must share the full DSM grid.")
    window = from_bounds(*source.bounds, transform=transform)
offsets = [window.col_off, window.row_off]
if not np.allclose(offsets, np.round(offsets), atol=1e-06, rtol=0):
    raise ValueError("The property clip is shifted off the DSM grid.")
row_start, col_start = (round(window.row_off), round(window.col_off))
if (
    row_start < 0
    or col_start < 0
    or row_start + clip.shape[0] > dsm.shape[0]
    or (col_start + clip.shape[1] > dsm.shape[1])
):
    raise ValueError("The property clip extends outside the full DSM.")
property_mask = np.zeros(dsm.shape, dtype=bool)
property_mask[row_start : row_start + clip.shape[0], col_start : col_start + clip.shape[1]] = (
    ~np.ma.getmaskarray(clip) & np.isfinite(clip.data)
)
valid = np.isfinite(dsm)
property_mask = ndi.binary_fill_holes(property_mask)
if not (property_mask & valid).any():
    raise ValueError("The property contains no valid elevations.")
rows, cols = np.indices(dsm.shape)
x = (cols + 0.5) * pixel_size
y = -(rows + 0.5) * pixel_size
print("Inputs loaded:", dsm.shape, profile["crs"])

view_bounds = [extent[0], extent[1], extent[2], extent[3]]


def map_view(ax, values, title, cmap="viridis", limits=(None, None), alpha=1):
    if values.ndim == 3:
        image = ax.imshow(values, extent=extent)
    else:
        image = ax.imshow(
            values, extent=extent, cmap=cmap, vmin=limits[0], vmax=limits[1], alpha=alpha
        )
    ax.set(
        xlim=view_bounds[:2],
        ylim=view_bounds[2:],
        title=title,
        xlabel="Easting (m)",
        ylabel="Northing (m)",
    )
    ax.set_aspect("equal")
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(axis="x", labelrotation=45)
    return image


def save_figure(fig, filename, dpi=170):
    output_path = OUTPUT_DIR / filename
    with output_path.open("wb") as output:
        fig.savefig(output, format="png", dpi=dpi)
    plt.show()
    plt.close(fig)


def polygon_parts(geometry):
    if geometry.is_empty:
        return []
    return [geometry] if geometry.geom_type == "Polygon" else list(geometry.geoms)


def draw_polygon(ax, geometry, **style):
    for part in polygon_parts(geometry):
        coordinates = np.asarray(part.exterior.coords)
        ax.fill(transform.c + coordinates[:, 0], transform.f + coordinates[:, 1], **style)


fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")
map_view(axes[0], rgb, "RGB input")
image = map_view(axes[1], np.ma.masked_invalid(dsm), "DSM input", "terrain")
fig.colorbar(image, ax=axes[1], shrink=0.65, label="Elevation (metres assumed)")
save_figure(fig, "step_01_inputs.png")

red, green, blue = rgb.astype(float).transpose(2, 0, 1)
green_index = (2 * green - red - blue) / (red + green + blue + 1)
ground_estimate = float(np.percentile(dsm[property_mask & valid], GROUND_PERCENTILE))
roof_candidate = (
    property_mask & valid & (dsm > ground_estimate + MIN_HEIGHT) & (green_index < GREEN_LIMIT)
)
roof_candidate = ndi.binary_opening(roof_candidate, iterations=2)
roof_candidate = ndi.binary_closing(roof_candidate, iterations=2) & property_mask
objects, object_count = ndi.label(roof_candidate)
if object_count == 0:
    raise ValueError("No roof region found. Review the height and colour thresholds.")
object_sizes = np.bincount(objects.ravel())
object_sizes[0] = 0
initial_roof = ndi.binary_fill_holes(objects == object_sizes.argmax()) & property_mask
print(f"Estimated local ground elevation: {ground_estimate:.2f}")
print(f"Initial roof plan area: {initial_roof.sum() * pixel_area:.2f} m²")

plot_rows, plot_cols = np.where(initial_roof)
view_bounds = [
    transform.c + (plot_cols.min() - 15) * pixel_size,
    transform.c + (plot_cols.max() + 15) * pixel_size,
    transform.f - (plot_rows.max() + 15) * pixel_size,
    transform.f - (plot_rows.min() - 15) * pixel_size,
]
fig, axes = plt.subplots(1, 3, figsize=(15, 6), layout="constrained")
image = map_view(
    axes[0],
    np.where(property_mask, dsm - ground_estimate, np.nan),
    "Height above estimated ground",
    "terrain",
)
fig.colorbar(image, ax=axes[0], shrink=0.65, label="Metres assumed")
image = map_view(axes[1], green_index, "Green index", "RdYlGn", (-0.1, 0.15))
fig.colorbar(image, ax=axes[1], shrink=0.65)
map_view(axes[2], rgb, "Initial roof outline")
axes[2].contour(
    transform.c + x, transform.f + y, initial_roof, levels=[0.5], colors=["lime"], linewidths=1.3
)
save_figure(fig, "step_02_roof_detection.png")

elevations = dsm[property_mask & valid]
percentiles = np.arange(101)
levels = np.percentile(elevations, percentiles)
low_pixels = property_mask & valid & (dsm <= ground_estimate)
fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
axes[0].plot(percentiles, levels)
axes[0].scatter(GROUND_PERCENTILE, ground_estimate, color="red", zorder=3)
axes[0].axvline(GROUND_PERCENTILE, color="red", linestyle="--", alpha=0.5)
axes[0].axhline(ground_estimate, color="red", linestyle="--", alpha=0.5)
axes[0].set(
    xlabel="Position among readings (%)",
    ylabel="Elevation (metres assumed)",
    title=f"Chosen reference: {ground_estimate:.2f} m",
)
axes[0].grid(alpha=0.2)
axes[1].imshow(rgb, extent=extent)
axes[1].imshow(
    np.ma.masked_where(~low_pixels, np.ones(dsm.shape)),
    extent=extent,
    cmap="Blues",
    vmin=0,
    vmax=1,
    alpha=0.65,
)
axes[1].set_title(f"Blue: elevations at or below {ground_estimate:.2f} m")
axes[1].set_axis_off()
save_figure(fig, "step_02_ground_reference.png")

filled = dsm.copy()
if not valid.all():
    nearest = ndi.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    filled[~valid] = dsm[tuple(nearest[:, ~valid])]
filtered = ndi.gaussian_filter(filled, sigma=1)
dz_dy, dz_dx = np.gradient(filtered, -pixel_size, pixel_size)
slope = np.degrees(np.arctan(np.hypot(dz_dx, dz_dy)))
aspect = np.degrees(np.arctan2(-dz_dx, -dz_dy)) % 360
aspect[slope < 1] = np.nan
print("Slope and downslope aspect calculated.")

fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")
for ax, values, title, cmap, limits in [
    (axes[0], slope, "Pixel slope", "viridis", (0, 65)),
    (axes[1], aspect, "Downslope aspect from grid north", "hsv", (0, 360)),
]:
    image = map_view(ax, np.where(initial_roof, values, np.nan), title, cmap, limits)
    fig.colorbar(image, ax=ax, shrink=0.65, label="Degrees")
save_figure(fig, "step_03_slope_aspect.png")

directions = np.floor((np.nan_to_num(aspect) + DIRECTION_BIN / 2) % 360 / DIRECTION_BIN).astype(int)
face_seeds = np.zeros(dsm.shape, dtype=np.int16)
face_count = 0
for direction in range(int(360 / DIRECTION_BIN)):
    stable = initial_roof & (directions == direction) & (slope > 5) & (slope < 65) & valid
    stable = ndi.binary_opening(stable, iterations=2)
    regions, region_count = ndi.label(stable)
    region_sizes = np.bincount(regions.ravel())
    for region_id in range(1, region_count + 1):
        if region_sizes[region_id] * pixel_area >= MIN_FACE_AREA:
            face_count += 1
            face_seeds[regions == region_id] = face_count
if face_count == 0:
    raise ValueError("No stable roof faces found.")
distance, nearest = ndi.distance_transform_edt(face_seeds == 0, return_indices=True)
roof_mask = (
    ndi.binary_fill_holes(initial_roof & (distance * pixel_size <= EDGE_DISTANCE)) & property_mask
)
face_labels = np.where(roof_mask, face_seeds[tuple(nearest)], 0).astype(np.int16)
print("Initial roof faces:", face_count)
print(f"Refined roof plan area: {roof_mask.sum() * pixel_area:.2f} m²")

fig, axes = plt.subplots(1, 3, figsize=(16, 6), layout="constrained")
image = map_view(
    axes[0], np.ma.masked_where(~initial_roof, directions), "Direction groups", "tab10", (-0.5, 7.5)
)
bar = fig.colorbar(image, ax=axes[0], shrink=0.65, ticks=range(8))
bar.ax.set_yticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
for ax, values, title in [
    (axes[1], face_seeds, "Stable roof areas"),
    (axes[2], face_labels, f"Initial roof faces: {face_count}"),
]:
    map_view(ax, rgb, title)
    map_view(
        ax,
        np.ma.masked_where(values == 0, values),
        title,
        "tab20",
        (0.5, face_count + 0.5),
        alpha=0.65,
    )
for face_id in range(1, face_count + 1):
    centre_row, centre_col = ndi.center_of_mass(face_seeds == face_id)
    axes[2].text(
        transform.c + (centre_col + 0.5) * pixel_size,
        transform.f - (centre_row + 0.5) * pixel_size,
        str(face_id),
        ha="center",
        va="center",
        fontsize=7,
        bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 1},
    )
save_figure(fig, "step_04_connected_faces.png")

def fit_roof_plane(points):
    design = np.column_stack((points[:, 0], points[:, 1], np.ones(len(points))))
    keep = np.ones(len(points), dtype=bool)
    for _ in range(5):
        if np.linalg.matrix_rank(design[keep]) < 3:
            raise ValueError("A roof face has degenerate XY coordinates.")
        coefficients = np.linalg.lstsq(design[keep], points[keep, 2], rcond=None)[0]
        residuals = points[:, 2] - design @ coefficients
        centre = np.median(residuals[keep])
        spread = 1.4826 * np.median(np.abs(residuals[keep] - centre))
        keep = np.abs(residuals - centre) <= max(0.15, 3 * spread)
    coefficients = np.linalg.lstsq(design[keep], points[keep, 2], rcond=None)[0]
    residuals = points[:, 2] - design @ coefficients
    rmse = float(np.sqrt(np.mean(residuals[keep] ** 2)))
    return (coefficients, rmse, int(keep.sum()))


planes = []
initial_surface = np.full(dsm.shape, np.nan)
for face_id in range(1, face_count + 1):
    fit_pixels = face_seeds == face_id
    points = np.column_stack((x[fit_pixels], y[fit_pixels], dsm[fit_pixels]))
    (a, b, c), rmse, inlier_count = fit_roof_plane(points)
    pixels = face_labels == face_id
    initial_surface[pixels] = a * x[pixels] + b * y[pixels] + c
    face_slope = float(np.degrees(np.arctan(np.hypot(a, b))))
    azimuth = float(np.degrees(np.arctan2(-a, -b)) % 360) if face_slope >= 1 else None
    planes.append(
        {
            "face_id": face_id,
            "a": float(a),
            "b": float(b),
            "c": float(c),
            "plan_area_m2": float(pixels.sum() * pixel_area),
            "slope_deg": face_slope,
            "downslope_azimuth_deg": azimuth,
            "seed_fit_rmse_m": rmse,
            "seed_inliers": inlier_count,
            "status": "initial_plane",
        }
    )
    print(
        f"Face {face_id:2d} | {pixels.sum() * pixel_area:6.2f} m² | slope {face_slope:5.1f}° | fit RMSE {rmse:.3f}"
    )
residual = np.where(roof_mask & valid, dsm - initial_surface, np.nan)

height_limits = (float(np.nanmin(dsm[roof_mask])), float(np.nanmax(dsm[roof_mask])))
fig, axes = plt.subplots(1, 3, figsize=(15, 6), layout="constrained")
for ax, values, title in [
    (axes[0], dsm, "Measured roof elevations"),
    (axes[1], initial_surface, "Elevations predicted by fitted planes"),
]:
    image = map_view(ax, np.where(roof_mask, values, np.nan), title, "terrain", height_limits)
    fig.colorbar(image, ax=ax, shrink=0.65, label="Elevation (metres assumed)")
image = map_view(axes[2], residual, "Measured minus fitted elevations", "RdBu_r", (-0.5, 0.5))
fig.colorbar(image, ax=axes[2], shrink=0.65, label="Difference (m)")
save_figure(fig, "step_05_fitted_planes.png")

main_planes = [
    plane
    for plane in planes
    if plane["plan_area_m2"] >= MAIN_FACE_AREA
    and MIN_MAIN_SLOPE <= plane["slope_deg"] <= MAX_MAIN_SLOPE
]
if not main_planes:
    raise ValueError("No main roof planes meet the configured suitability assumptions.")
main_ids = [plane["face_id"] for plane in main_planes]
support = np.where(np.isin(face_seeds, main_ids), face_seeds, 0)
_, nearest = ndi.distance_transform_edt(support == 0, return_indices=True)
final_labels = np.where(roof_mask, support[tuple(nearest)], 0).astype(np.int16)
baseline = np.full(dsm.shape, np.nan)
for plane in main_planes:
    pixels = final_labels == plane["face_id"]
    baseline[pixels] = plane["a"] * x[pixels] + plane["b"] * y[pixels] + plane["c"]
blue_contrast = (blue - red) / (blue + red + 1)
glazing = (blue_contrast > BLUE_CONTRAST) & (green > red * 1.15) & (blue > 120)
brightness = rgb.mean(axis=2)
details = ndi.gaussian_filter(brightness, 1) - ndi.gaussian_filter(brightness, 6) > 15
light_regions, _ = ndi.label(details & (blue > red * 1.05) & roof_mask)
for region_id, slices in enumerate(ndi.find_objects(light_regions), 1):
    region = light_regions[slices] == region_id
    height, width = region.shape
    if (
        0.08 <= region.sum() * pixel_area <= 2.5
        and max(height, width) * pixel_size <= 2.5
        and (max(height, width) / min(height, width) <= 3)
        and (region.mean() > 0.25)
    ):
        glazing[slices] |= ndi.binary_fill_holes(region)
large_gaps = ndi.distance_transform_edt(~valid) * pixel_size > 0.5
obstructions = roof_mask & ((np.abs(dsm - baseline) > HEIGHT_DEVIATION) | glazing | large_gaps)
obstructions = ndi.binary_opening(obstructions, iterations=1)
obstructions = ndi.binary_closing(obstructions, iterations=2) & roof_mask
regions, _ = ndi.label(obstructions)
sizes = np.bincount(regions.ravel()) * pixel_area
keep_regions = np.flatnonzero((sizes >= 0.12) & (np.arange(len(sizes)) > 0))
obstructions = ndi.binary_fill_holes(np.isin(regions, keep_regions)) & roof_mask
best_score = np.full(dsm.shape, np.inf)
for plane in main_planes:
    distance = ndi.distance_transform_edt(support != plane["face_id"]) * pixel_size
    prediction = plane["a"] * x + plane["b"] * y + plane["c"]
    score = np.abs(dsm - prediction) + 0.025 * distance
    improve = roof_mask & valid & ~obstructions & (distance < 2.5) & (score < best_score)
    final_labels[improve] = plane["face_id"]
    best_score[improve] = score[improve]
local_transform = Affine(pixel_size, 0, 0, 0, -pixel_size, 0)


def mask_polygon(mask):
    return unary_union(
        [
            shape(geometry)
            for geometry, value in shapes(
                mask.astype("uint8"), mask=mask, transform=local_transform
            )
            if value
        ]
    )


roof_polygon = mask_polygon(roof_mask)
obstruction_polygon = mask_polygon(obstructions)
face_polygons = {p["face_id"]: mask_polygon(final_labels == p["face_id"]) for p in main_planes}
for plane in main_planes:
    plane["plan_area_m2"] = float((final_labels == plane["face_id"]).sum() * pixel_area)
    plane["surface_area_m2"] = plane["plan_area_m2"] / np.cos(np.radians(plane["slope_deg"]))
    plane["status"] = "automatic_main_plane"
print("Main roof faces:", len(main_planes))
print(f"Detected obstruction plan area: {obstructions.sum() * pixel_area:.2f} m²")

fig, axes = plt.subplots(1, 4, figsize=(18, 6), layout="constrained")
for ax in axes:
    map_view(ax, rgb, "")
map_view(
    axes[0],
    np.ma.masked_where(face_labels == 0, face_labels),
    "Initial roof faces",
    "tab20",
    (0.5, face_count + 0.5),
    alpha=0.65,
)
map_view(
    axes[1],
    np.ma.masked_where(final_labels == 0, final_labels),
    "Main roof faces",
    "tab20",
    (0.5, face_count + 0.5),
    alpha=0.65,
)
draw_polygon(axes[2], obstruction_polygon, color="#d83d37", alpha=0.65)
axes[2].set_title("Detected obstructions")
remaining = roof_mask & ~obstructions
map_view(
    axes[3],
    np.ma.masked_where(~remaining, np.ones_like(dsm)),
    "Remaining roof before clearances",
    "Greens",
    (0, 1),
    alpha=0.55,
)
save_figure(fig, "step_06_obstructions_faces.png")

roof_rows, roof_cols = np.where(roof_mask)
corner_rows = np.column_stack((roof_rows, roof_rows + 1, roof_rows + 1, roof_rows))
corner_cols = np.column_stack((roof_cols, roof_cols, roof_cols + 1, roof_cols + 1))
grid_ids = corner_rows * (dsm.shape[1] + 1) + corner_cols
node_ids, inverse = np.unique(grid_ids, return_inverse=True)
quads = inverse.reshape(-1, 4)
node_rows, node_cols = np.divmod(node_ids, dsm.shape[1] + 1)
node_x, node_y = (node_cols * pixel_size, -node_rows * pixel_size)
coefficient_table = np.zeros((face_count + 1, 3))
for plane in main_planes:
    coefficient_table[plane["face_id"]] = [plane["a"], plane["b"], plane["c"]]
coefficients = coefficient_table[final_labels[roof_rows, roof_cols]]
predictions = (
    coefficients[:, 0, None] * node_x[quads]
    + coefficients[:, 1, None] * node_y[quads]
    + coefficients[:, 2, None]
)
node_z = np.bincount(quads.ravel(), weights=predictions.ravel()) / np.bincount(quads.ravel())
vertices = np.column_stack((node_x, node_y, node_z))
triangles = np.vstack((quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]))
triangle_labels = np.tile(final_labels[roof_rows, roof_cols], 2)
obstruction_triangles = np.tile(obstructions[roof_rows, roof_cols], 2)
surface = np.full(dsm.shape, np.nan)
surface[roof_rows, roof_cols] = node_z[quads].mean(axis=1)
final_residual = np.where(roof_mask & valid, dsm - surface, np.nan)
observed_node_z = ndi.map_coordinates(
    filled, [node_rows - 0.5, node_cols - 0.5], order=1, mode="nearest"
)
print("Mesh vertices:", len(vertices), "Triangles:", len(triangles))
print("Missing roof pixels reconstructed:", int((roof_mask & ~valid).sum()))

fig = plt.figure(figsize=(10, 7), layout="constrained")
ax = fig.add_subplot(111, projection="3d")
ax.plot_trisurf(
    node_x, node_y, node_z, triangles=triangles, color="#8fa6ae", linewidth=0, antialiased=False
)
ax.set_box_aspect((np.ptp(node_x), np.ptp(node_y), np.ptp(node_z)))
ax.view_init(elev=45, azim=-65)
ax.set(
    xlabel="Local east (m)",
    ylabel="Local north (m)",
    zlabel="Elevation (metres assumed)",
    title="Continuous reconstructed roof mesh",
)
save_figure(fig, "step_07_roof_mesh.png")

def plane_basis(plane):
    normal = np.array([-plane["a"], -plane["b"], 1.0])
    normal /= np.linalg.norm(normal)
    u = np.array([plane["b"], -plane["a"], 0.0])
    u /= np.linalg.norm(u)
    return (u, np.cross(u, normal))


panels = []
usable_faces = {}
for plane in main_planes:
    face_id = plane["face_id"]
    a, b = (plane["a"], plane["b"])
    u, v = plane_basis(plane)
    projection = [u[0], u[1], v[0] + v[2] * a, v[1] + v[2] * b, 0, 0]
    face_uv = affine_transform(face_polygons[face_id], projection)
    obstacles_uv = affine_transform(obstruction_polygon, projection)
    usable = face_uv.buffer(-EDGE_CLEARANCE).difference(obstacles_uv.buffer(OBSTRUCTION_CLEARANCE))
    usable_faces[face_id] = usable
    if usable.is_empty:
        continue
    prepared = prep(usable)
    min_u, min_v, max_u, max_v = usable.bounds
    best_panels, best_span = ([], np.inf)
    for width, height in [(PANEL_WIDTH, PANEL_LENGTH), (PANEL_LENGTH, PANEL_WIDTH)]:
        for du in np.linspace(0, width + PANEL_GAP, 4, endpoint=False):
            for dv in np.linspace(0, height + PANEL_GAP, 4, endpoint=False):
                batch, positions = ([], [])
                for i, uu in enumerate(
                    np.arange(min_u + du, max_u - width + 1e-09, width + PANEL_GAP)
                ):
                    for j, vv in enumerate(
                        np.arange(min_v + dv, max_v - height + 1e-09, height + PANEL_GAP)
                    ):
                        rectangle = box(uu, vv, uu + width, vv + height)
                        if prepared.covers(rectangle):
                            batch.append(rectangle)
                            positions.append((i, j))
                if not batch:
                    continue
                positions = np.asarray(positions)
                occupied = np.zeros(tuple(positions.max(axis=0) + 1), dtype=bool)
                occupied[tuple(positions.T)] = True
                groups, _ = ndi.label(occupied)
                counts = np.bincount(groups.ravel())
                counts[0] = 0
                largest = counts.argmax()
                batch = [
                    rectangle
                    for rectangle, position in zip(batch, positions)
                    if groups[tuple(position)] == largest
                ]
                if len(batch) < MIN_ARRAY_PANELS:
                    continue
                span = unary_union(batch).envelope.area
                if len(batch) > len(best_panels) or (
                    len(batch) == len(best_panels) and span < best_span
                ):
                    best_panels, best_span = (batch, span)
    for rectangle in best_panels:
        uv = np.asarray(rectangle.exterior.coords)[:4]
        xyz = uv[:, 0, None] * u + uv[:, 1, None] * v + [0, 0, plane["c"]]
        panels.append(
            {
                "panel_id": len(panels) + 1,
                "face_id": face_id,
                "uv": rectangle,
                "xy": affine_transform(rectangle, [u[0], v[0], u[1], v[1], 0, 0]),
                "xyz": xyz,
            }
        )
print("Panels placed:", len(panels))
print(f"Total panel surface area: {len(panels) * PANEL_LENGTH * PANEL_WIDTH:.2f} m²")
for plane in main_planes:
    plane["panel_count"] = sum((p["face_id"] == plane["face_id"] for p in panels))

fig, ax = plt.subplots(figsize=(9, 10), layout="constrained")
map_view(
    ax, rgb, f"{len(panels)} panels; each {PANEL_LENGTH:.3f} × {PANEL_WIDTH:.3f} m on the roof"
)
draw_polygon(ax, obstruction_polygon, color="#d83d37", alpha=0.65)
for panel in panels:
    draw_polygon(ax, panel["xy"], facecolor="#c6f63c", edgecolor="#19342a", linewidth=0.6)
save_figure(fig, "step_08_panel_layout.png")

def format_map(ax):
    ax.set_xlim(
        transform.c + (roof_cols.min() - 15) * pixel_size,
        transform.c + (roof_cols.max() + 15) * pixel_size,
    )
    ax.set_ylim(
        transform.f - (roof_rows.max() + 15) * pixel_size,
        transform.f - (roof_rows.min() - 15) * pixel_size,
    )
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(axis="x", labelrotation=45)


fig, axes = plt.subplots(1, 3, figsize=(16, 6), layout="constrained")
for ax in axes:
    ax.imshow(rgb, extent=extent)
    format_map(ax)
axes[0].imshow(
    np.ma.masked_where(final_labels == 0, final_labels), extent=extent, cmap="tab20", alpha=0.55
)
axes[0].set_title("Main roof faces")
draw_polygon(axes[1], obstruction_polygon, color="#d83d37", alpha=0.65)
axes[1].set_title("Detected obstructions")
draw_polygon(axes[2], obstruction_polygon, color="#d83d37", alpha=0.65)
for panel in panels:
    draw_polygon(axes[2], panel["xy"], facecolor="#c6f63c", edgecolor="#19342a", linewidth=0.5)
axes[2].set_title(f"Automatic layout: {len(panels)} panels")
save_figure(fig, "final_layout.png", dpi=200)
fig, axes = plt.subplots(1, 3, figsize=(15, 6), layout="constrained")
for ax, values, title, cmap, limits in [
    (axes[0], slope, "Pixel slope", "viridis", (0, 65)),
    (axes[1], aspect, "Downslope aspect", "hsv", (0, 360)),
    (axes[2], final_residual, "DSM minus reconstructed roof", "RdBu_r", (-0.5, 0.5)),
]:
    image = ax.imshow(
        np.where(roof_mask, values, np.nan),
        extent=extent,
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
    )
    format_map(ax)
    ax.set_title(title)
    fig.colorbar(image, ax=ax, shrink=0.65, label="Metres assumed" if ax is axes[2] else "Degrees")
save_figure(fig, "terrain_and_residual.png", dpi=180)

scene = go.Figure(
    go.Mesh3d(
        x=node_x,
        y=node_y,
        z=node_z,
        i=triangles[:, 0],
        j=triangles[:, 1],
        k=triangles[:, 2],
        intensity=triangle_labels,
        intensitymode="cell",
        colorscale="Viridis",
        showscale=False,
        name="Reconstructed roof",
        showlegend=True,
        flatshading=True,
    )
)
object_triangles = triangles[obstruction_triangles]
scene.add_trace(
    go.Mesh3d(
        x=node_x,
        y=node_y,
        z=observed_node_z + 0.025,
        i=object_triangles[:, 0],
        j=object_triangles[:, 1],
        k=object_triangles[:, 2],
        color="#d83d37",
        name="Detected obstructions",
        showlegend=True,
    )
)
panel_vertices = (
    np.concatenate([panel["xyz"] + [0, 0, 0.05] for panel in panels])
    if panels
    else np.empty((0, 3))
)
panel_triangles = np.asarray(
    [[4 * i, 4 * i + 1, 4 * i + 2] for i in range(len(panels))]
    + [[4 * i, 4 * i + 2, 4 * i + 3] for i in range(len(panels))],
    dtype=int,
).reshape(-1, 3)
scene.add_trace(
    go.Mesh3d(
        x=panel_vertices[:, 0],
        y=panel_vertices[:, 1],
        z=panel_vertices[:, 2],
        i=panel_triangles[:, 0],
        j=panel_triangles[:, 1],
        k=panel_triangles[:, 2],
        color="#c6f63c",
        name="Panels",
        showlegend=True,
    )
)
scene.update_layout(
    title="Automatic roof reconstruction and panel layout",
    height=750,
    scene={
        "aspectmode": "data",
        "xaxis_title": "Local east (m)",
        "yaxis_title": "Local north (m)",
        "zaxis_title": "Elevation (metres assumed)",
    },
    margin={"l": 0, "r": 0, "t": 50, "b": 0},
)
scene.write_html(OUTPUT_DIR / "roof_3d.html", include_plotlyjs=True, auto_open=False)
with (OUTPUT_DIR / "roof_mesh.obj").open("w", encoding="utf-8") as stream:
    for xx, yy, zz in vertices:
        stream.write(f"v {xx:.5f} {yy:.5f} {zz:.5f}\n")
    for triangle in triangles + 1:
        stream.write("f " + " ".join(map(str, triangle)) + "\n")
print("3D model:", OUTPUT_DIR / "roof_3d.html")

def save_raster(filename, values, dtype, nodata):
    settings = profile.copy()
    settings.update(driver="GTiff", dtype=dtype, count=1, nodata=nodata, compress="deflate")
    with rasterio.open(OUTPUT_DIR / filename, "w", **settings) as target:
        target.write(values.astype(dtype), 1)


save_raster("roof_mask.tif", roof_mask, "uint8", None)
save_raster("obstructions.tif", obstructions, "uint8", None)
save_raster("roof_faces.tif", final_labels, "int16", 0)
for filename, values in [
    ("elevation.tif", dsm),
    ("slope.tif", slope),
    ("aspect.tif", aspect),
    ("reconstructed_roof.tif", surface),
    ("residual.tif", final_residual),
]:
    save_raster(
        filename, np.where(roof_mask & np.isfinite(values), values, -9999), "float32", -9999
    )


def save_geojson(filename, items):
    features = []
    for geometry, properties in items:
        projected = affine_transform(geometry, [1, 0, 0, 1, transform.c, transform.f])
        geographic = transform_geom(profile["crs"], "EPSG:4326", mapping(projected))
        features.append({"type": "Feature", "geometry": geographic, "properties": properties})
    content = {"type": "FeatureCollection", "features": features}
    (OUTPUT_DIR / filename).write_text(json.dumps(content, indent=2), encoding="utf-8")


save_geojson(
    "roof_faces.geojson", [(face_polygons[plane["face_id"]], plane) for plane in main_planes]
)
save_geojson(
    "obstructions.geojson",
    [
        (g, {"id": i + 1, "type": "obstruction_candidate"})
        for i, g in enumerate(polygon_parts(obstruction_polygon))
    ],
)
save_geojson(
    "panels.geojson",
    [
        (
            panel["xy"],
            {
                "panel_id": panel["panel_id"],
                "face_id": panel["face_id"],
                "length_m": PANEL_LENGTH,
                "width_m": PANEL_WIDTH,
            },
        )
        for panel in panels
    ],
)
with (OUTPUT_DIR / "roof_planes.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(main_planes[0]))
    writer.writeheader()
    writer.writerows(main_planes)
print("Results saved in:", OUTPUT_DIR)

all_inside = all((usable_faces[p["face_id"]].covers(p["uv"]) for p in panels))
all_clear = all((not p["xy"].intersects(obstruction_polygon) for p in panels))
overlap_area = (
    sum((p["xy"].area for p in panels)) - unary_union([panel["xy"] for panel in panels]).area
)
dimension_errors = []
minimum_gap = np.inf
for i, panel in enumerate(panels):
    lengths = np.linalg.norm(np.roll(panel["xyz"], -1, axis=0) - panel["xyz"], axis=1)
    dimension_errors.append(
        np.max(
            np.abs(
                np.sort(lengths) - np.sort([PANEL_LENGTH, PANEL_LENGTH, PANEL_WIDTH, PANEL_WIDTH])
            )
        )
    )
    for other in panels[i + 1 :]:
        if other["face_id"] == panel["face_id"]:
            minimum_gap = min(minimum_gap, panel["uv"].distance(other["uv"]))
edges = np.vstack([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]])
_, edge_counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
cross_products = np.cross(
    vertices[triangles[:, 1]] - vertices[triangles[:, 0]],
    vertices[triangles[:, 2]] - vertices[triangles[:, 0]],
)
checks = {
    "panels_present": len(panels) > 0,
    "all_panels_inside_usable_areas": bool(all_inside),
    "all_panels_clear_of_detected_obstructions": bool(all_clear),
    "no_panel_overlaps": bool(abs(overlap_area) < 1e-07),
    "panel_dimensions_correct_in_3d": bool(max(dimension_errors, default=0) < 1e-07),
    "panel_gaps_respected": bool(minimum_gap >= PANEL_GAP - 1e-07),
    "all_roof_pixels_reconstructed": bool(np.isfinite(surface[roof_mask]).all()),
    "mesh_edges_have_at_most_two_triangles": bool(edge_counts.max() <= 2),
    "mesh_triangles_face_upward": bool((cross_products[:, 2] > 0).all()),
}
evaluation = roof_mask & valid & ~obstructions
rmse = float(np.sqrt(np.mean(final_residual[evaluation] ** 2)))
summary = {
    "status": "automatic_results_pending_final_QA",
    "roof_plan_area_m2": float(roof_mask.sum() * pixel_area),
    "obstruction_candidate_plan_area_m2": float(obstructions.sum() * pixel_area),
    "main_roof_regions": len(main_planes),
    "panel_count": len(panels),
    "panel_surface_area_m2": len(panels) * PANEL_LENGTH * PANEL_WIDTH,
    "surface_to_DSM_rmse_excluding_detected_obstructions_m": rmse,
    "missing_roof_pixels_reconstructed": int((roof_mask & ~valid).sum()),
    "analysis_crs": profile["crs"].to_string(),
    "mesh_origin_easting": transform.c,
    "mesh_origin_northing": transform.f,
    "vertical_units": "metres_assumed",
    "vertical_scale": VERTICAL_SCALE,
    "panel_length_m": PANEL_LENGTH,
    "panel_width_m": PANEL_WIDTH,
    "edge_clearance_m": EDGE_CLEARANCE,
    "obstruction_clearance_m": OBSTRUCTION_CLEARANCE,
    "panel_gap_m": PANEL_GAP,
    "geometry_checks": checks,
    "minimum_same_face_panel_gap_m": float(minimum_gap) if np.isfinite(minimum_gap) else None,
    "maximum_panel_dimension_error_m": float(max(dimension_errors, default=0)),
    "parameters": {
        "ground_percentile": GROUND_PERCENTILE,
        "min_main_slope_deg": MIN_MAIN_SLOPE,
        "max_main_slope_deg": MAX_MAIN_SLOPE,
        "min_height": MIN_HEIGHT,
        "green_limit": GREEN_LIMIT,
        "min_face_area": MIN_FACE_AREA,
        "main_face_area": MAIN_FACE_AREA,
        "direction_bin": DIRECTION_BIN,
        "height_deviation": HEIGHT_DEVIATION,
        "blue_contrast": BLUE_CONTRAST,
        "edge_distance": EDGE_DISTANCE,
        "minimum_array_panels": MIN_ARRAY_PANELS,
    },
    "limitation": "Checks validate generated geometry against detected masks; final RGB QA is pending.",
}
(OUTPUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
for name, passed in checks.items():
    print("PASS" if passed else "FAIL", name)
if not all(checks.values()):
    raise ValueError("One or more geometry checks failed. Inspect summary.json.")
print(f"Surface-to-DSM RMSE outside detected obstructions: {rmse:.3f} m")
print("Open this file in a browser:", OUTPUT_DIR / "roof_3d.html")

fig, ax = plt.subplots(figsize=(11, 5), layout="constrained")
ax.axis("off")
ax.set_title(f"Geometry checks: {len(panels)} panels", pad=15)
table = ax.table(
    cellText=[
        [name.replace("_", " "), "PASS" if passed else "FAIL"] for name, passed in checks.items()
    ],
    colLabels=["Check", "Result"],
    colWidths=[0.85, 0.15],
    cellLoc="left",
    loc="center",
)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.8)
for index, passed in enumerate(checks.values(), 1):
    table[index, 1].set_facecolor("#dff0df" if passed else "#ffd6d6")
save_figure(fig, "step_10_geometry_checks.png")


