import datetime
# import pprint

import ee
import pytest

import openet.geesebal.openet_landsat as landsat
import openet.geesebal.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
SCENE_ID = 'LC08_044033_20170716'
# SCENE_ID = 'LC08_042035_20150713'
SCENE_DT = datetime.datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_TIME = utils.millis(SCENE_DT)


def landsat_image(ultra_blue=0.2, blue=0.2, green=0.2, red=0.2, nir=0.7,
                  swir1=0.2, swir2=0.2, tir=300):
    """Construct a fake Landsat 8 SR image with renamed bands"""
    return ee.Image.constant([ultra_blue, blue, green, red, nir, swir1, swir2, tir]) \
        .rename(['ultra_blue', 'blue', 'green', 'red', 'nir',
                 'swir1', 'swir2', 'tir']) \


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.2, 0.1, -0.333333],    # 0.163636
        [0.2, 9.0 / 55, -0.1],    # 0.163636
        [0.2, 0.2, 0.0],
        # [0.1, 11.0 / 90,  0.1], # 0.122222
        [0.2, 0.3, 0.2],
        # [0.1, 13.0 / 70, 0.3],  # 0.185714
        [0.3, 0.7, 0.4],
        # [0.2, 0.6, 0.5],
        [0.2, 0.8, 0.6],
        # [0.1, 17.0 / 30, 0.7],  # 0.566666
        [0.1, 0.9, 0.8],
        [0.2, 0.7, 0.55555555],
    ]
)
def test_ndvi_values(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.ndvi(
        landsat_image(red=red, nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_ndvi_band_name():
    output = landsat.ndvi(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'ndvi'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        # CGM - These tests aren't all necessary, just playing around
        # Negative NDVI
        [0.2, 0.1, 0],
        # 0 NDVI
        [0.2, 0.2, 0],
        # NDVIs from 0.0 to 0.05 should end up as 0
        [0.2, 0.205, 0.0],
        [0.2, 0.220, 0.0],
        # Positive NDVI
        # [0.2, 0.3, 0.15],
        # [0.3, 0.7, 0.35],
        [0.2, 0.7, 0.505555],
        # [0.2, 0.8, 0.55],
        # [0.1, 0.9, 0.75],
    ]
)
def test_fipar_values(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.fipar(
        landsat_image(red=red, nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_fipar_band_name():
    output = landsat.fipar(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'fipar'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        # CGM - Why are these all negative?
        [0.2, 0.1, 0],
        [0.2, 0.2, 0],
        [0.2, 0.3, -0.325038],
        [0.3, 0.7, -0.861566],
        [0.2, 0.7, -1.408641],
        [0.2, 0.8, -1.597015],
        [0.1, 0.9, -2.772589],
    ]
)
def test_lai_values(red, nir, expected, tol=0.0001):
    output = utils.constant_image_value(landsat.lai(
        landsat_image(red=red, nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_lai_band_name():
    output = landsat.lai(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'lai'


@pytest.mark.parametrize(
    'green, nir, expected',
    [
        # Water should go positive
        [0.2, 0.1, 0.333333],
        [0.2, 0.2, 0.0],
        [0.2, 0.8, -0.6],
    ]
)
def test_ndwi_values(green, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.ndwi(
        landsat_image(green=green, nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_ndwi_band_name():
    output = landsat.ndwi(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'ndwi'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        # CGM - The commented out number is the LAI
        # CGM - These tests don't make sense because of the weird LAI values
        [0.2, 0.1, 0.95],
        [0.2, 0.2, 0.95],
        # [0.2, 0.3, 0],  # -0.325038
        # [0.3, 0.7, 0],  # -0.861566
        # [0.2, 0.7, 0],  # -1.408641
        # [0.2, 0.8, 0],  # -1.597015
        # [0.1, 0.9, 0],  # -2.772589
    ]
)
def test_emissivity_values(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.emissivity(
        landsat_image(red=red, nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_emissivity_band_name():
    output = landsat.emissivity(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'emissivity'


@pytest.mark.parametrize(
    'tir, red, nir, expected',
    [
        # CGM - This was just the number returned
        # I still need to come up with inputs values to test the edge cases
        #   but I want to figure out the LAI issue first.
        [300, 0.2, 0.7, 302.48],
    ]
)
def test_lst_values(tir, red, nir, expected, tol=0.01):
    output = utils.constant_image_value(landsat.lst(
        landsat_image(tir=tir, red=red, nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_lst_band_name():
    output = landsat.lst(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'lst'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.2, 0.7, 0.535714],
        [0.1, 0.9, 0.689],
    ]
)
def test_savi_values(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.savi(
        landsat_image(red=red, nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_savi_band_name():
    output = landsat.savi(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'savi'


# def albedo_l457():
#     assert False


def test_albedo_l457_band_name():
    output = landsat.albedo_l457(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'albedo'


@pytest.mark.parametrize(
    'nir, expected',
    [
        # CGM - The expected value is not 0.2 since the coefficients don't sum to 100
        [0.2, 0.1998],
        [0.7, 0.3403],
    ]
)
def test_albedo_l8_values(nir, expected, tol=0.0001):
    """The default image is all 0.2 except NIR"""
    output = utils.constant_image_value(landsat.albedo_l8(
        landsat_image(nir=nir)).rename(['output']))['output']
    assert abs(output - expected) <= tol


def test_albedo_l8_band_name():
    output = landsat.albedo_l8(landsat_image()).getInfo()['bands'][0]['id']
    assert output == 'albedo'


# def cloud_mask_sr_l457():
#     assert False
#
#
# def cloud_mask_sr_l8():
#     assert False
#
#
# def cloud_mask_C2_l457():
#     assert False
#
#
# def cloud_mask_C2_l8():
#     assert False
