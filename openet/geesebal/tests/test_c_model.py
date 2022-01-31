import ee
import pytest

import openet.geesebal.model as model
import openet.geesebal.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


@pytest.mark.parametrize(
    'rn, ndvi, albedo, lst_dem, ndwi, expected',
    [
        [100, 0.5, 0.2, 300, 0.1, 100 * 0.5],
        [100, 0.5, 0.2, 300, -0.1, 13.3],
    ]
)
def test_soil_heat_flux(rn, ndvi, albedo, lst_dem, ndwi, expected, tol=0.1):
    output = utils.constant_image_value(model.soil_heat_flux(
        rn=ee.Image.constant(rn), ndvi=ee.Image.constant(ndvi),
        albedo=ee.Image.constant(albedo),
        lst_dem=ee.Image.constant(lst_dem),
        ndwi=ee.Image.constant(ndwi)))
    assert abs(output - expected) <= tol
