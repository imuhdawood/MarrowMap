import numpy as np
from geopandas import GeoDataFrame
from shapely import MultiPolygon, Polygon, Point, make_valid
from shapely.affinity import scale as shapely_scale


def clean_geometry(gdf : GeoDataFrame, log_fn) -> GeoDataFrame:
    cleaned = gdf.copy()

    invalid = cleaned.loc[~cleaned.geometry.is_valid].make_valid()
    log_fn(f"Megs with invalid Geometry # {len(invalid)}")
    for idx in invalid.index:
        largest_poly = max(list(invalid[idx].geoms), key=lambda geom: geom.area)
        cleaned.geometry[idx] = largest_poly

    invalid = cleaned.loc[~cleaned.geometry.is_valid].make_valid()
    log_fn(f"Megs with invalid Geometry after correction # {len(invalid)}")

    non_poly = cleaned[cleaned.geometry.type != "Polygon"]
    log_fn(f"Megs with Multipolygon geometry # {len(non_poly)}")
    for idx in non_poly.index:
        largest_poly = max(list(non_poly.geometry[idx].geoms), key=lambda geom: geom.area)
        cleaned.geometry[idx] = largest_poly

    return cleaned.drop_duplicates(subset="geometry")

def rescale_polygons(gdf: GeoDataFrame, mpp: float) -> GeoDataFrame:
    """
    Multiply polygon coordinates by 1 / microns_per_pixel.
    Assumes all geometries are simple Polygons.
    """
    scale = 1.0 / mpp
    geoms = []

    for geom in gdf.geometry:
        coords = np.asarray(geom.exterior.coords)
        scaled = coords * scale
        geoms.append(Polygon(scaled))
    return gdf.assign(geometry=geoms)