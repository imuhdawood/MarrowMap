import numpy as np
from geopandas import GeoDataFrame
from shapely import MultiPolygon, Polygon, Point, make_valid
from shapely.affinity import scale as shapely_scale

from shapely.ops import unary_union
from shapely.affinity import rotate, scale, translate


def align_sample_polygons(group):
    """
    For a given group of polygons (from one sample), fix invalid geometries,
    determine if the sample is horizontally oriented, and then rotate/translate
    the polygons accordingly.
    """
    group = group.copy()
    
    # Fix invalid geometries
    group['geometry'] = group['geometry'].apply(lambda geom: geom.buffer(0) if not geom.is_valid else geom)
    
    # Combine all geometries into one to get overall bounds
    sample_union = unary_union(list(group['geometry']))
    minx, miny, maxx, maxy = sample_union.bounds
    width, height = maxx - minx, maxy - miny

    if width > height:
        # Rotate all polygons by 90 degrees around the sample's centroid
        origin = sample_union.centroid
        group['geometry_aligned'] = group['geometry'].apply(lambda g: rotate(g, 90, origin=origin))
        new_union = unary_union(list(group['geometry_aligned']))
        new_minx, new_miny, _, _ = new_union.bounds
        group['geometry_aligned'] = group['geometry_aligned'].apply(lambda g: translate(g, xoff=-new_minx, yoff=-new_miny))
    else:
        # No rotation needed, just translate so that (minx, miny) becomes (0,0)
        group['geometry_aligned'] = group['geometry'].apply(lambda g: translate(g, xoff=-minx, yoff=-miny))
    return group


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