from spatialdata import SpatialData
from spatialdata.models import SpatialElement

from spatialdata.transformations import (
    Identity,
    set_transformation,
    get_transformation
)

def get_intrinsic_cs(
    sdata: SpatialData, element: SpatialElement | str, name: str | None = None
) -> str:
    """Gets the name of the intrinsic coordinate system of an element

    Args:
        sdata: A SpatialData object
        element: `SpatialElement`, or its key
        name: Name to provide to the intrinsic coordinate system if not existing. By default, uses the element id.

    Returns:
        Name of the intrinsic coordinate system
    """
    if name is None:
        name = f"_{element if isinstance(element, str) else id(element)}_intrinsic"

    if isinstance(element, str):
        element = sdata[element]

    for cs, transform in get_transformation(element, get_all=True).items():
        if isinstance(transform, Identity):
            return cs

    set_transformation(element, Identity(), name)
    return name



def to_intrinsic(
    sdata: SpatialData, element: SpatialElement | str, element_cs: SpatialElement | str
) -> SpatialElement:
    """Transforms a `SpatialElement` into the intrinsic coordinate system of another `SpatialElement`

    Args:
        sdata: A SpatialData object
        element: `SpatialElement` to transform, or its key
        element_cs: `SpatialElement` of the target coordinate system, or its key

    Returns:
        The `SpatialElement` after transformation in the target coordinate system
    """
    if isinstance(element, str):
        element = sdata[element]
    cs = get_intrinsic_cs(sdata, element_cs)
    return sdata.transform_element_to_coordinate_system(element, cs)