import datetime
import logging
import pprint

import ee
import pytest

import openet.geesebal as geesebal
import openet.geesebal.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
COLL_ID = 'LANDSAT/LC08/C01/T1_SR/'
SCENE_ID = 'LC08_044033_20170716'
SCENE_TIME = 1500230731090
SCENE_DT = datetime.datetime.utcfromtimestamp(SCENE_TIME / 1000.0)
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_DOY = int(SCENE_DT.strftime('%j'))
SCENE_0UTC_DT = datetime.datetime.strptime(SCENE_DATE, '%Y-%m-%d')
SUN_ELEVATION = 64.3
TEST_POINT = (-121.5265, 38.7399)
# TEST_POINT = (-120.113, 36.336)


# # Should these be test fixtures instead?
# # I'm not sure how to make them fixtures and allow input parameters
# def l8_image(ultra_blue=0.2, blue=0.2, green=0.2, red=0.2, nir=0.7,
#              swir1=0.2, swir2=0.2, bt=300):
#     """Construct a fake Landsat 8 SR image with renamed bands"""
#     return ee.Image.constant([ultra_blue, blue, green, red, nir, swir1, swir2, bt]) \
#         .rename(['ultra_blue', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst']) \
#         .set({
#             'system:index': SCENE_ID,
#             'system:time_start': SCENE_TIME,
#             'system:id': f'{COLL_ID}/{SCENE_ID}',
#             'SUN_ELEVATION': SUN_ELEVATION,
#         })


# def l7_image(blue=0.2, green=0.2, red=0.2, nir=0.7,
#              swir1=0.2, swir2=0.2, bt=300):
#     """Construct a fake Landsat 8 SR image with renamed bands"""
#     return ee.Image.constant([blue, green, red, nir, swir1, swir2, bt]) \
#         .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst']) \
#         .set({
#             'system:index': SCENE_ID,
#             'system:time_start': SCENE_TIME,
#             'system:id': COLL_ID + SCENE_ID',
#             'SUN_ELEVATION': SUN_ELEVATION,
#         })


def default_image(albedo=0.2, emissivity=0.964, lai=1.4, lst=300,
                  ndvi=0.5, ndwi=-0.1, savi=0.5):
    # First construct a fake 'prepped' input image
    # Map the values to an actual Landsat image
    return ee.Image(COLL_ID + SCENE_ID) \
        .select([0, 1, 2, 3, 4, 5, 6]) \
        .multiply(0).double() \
        .add(ee.Image.constant([albedo, emissivity, lai, lst, ndvi, ndwi, savi])) \
        .rename(['albedo', 'emissivity', 'lai', 'lst', 'ndvi', 'ndwi', 'savi']) \
        .set({
            'system:index': SCENE_ID,
            'system:time_start': SCENE_TIME,
            'system:id': COLL_ID + SCENE_ID,
            'SUN_ELEVATION': SUN_ELEVATION,
        })
    # return ee.Image.constant([albedo, emissivity, lai, lst, ndvi, ndwi, savi]) \
    #     .rename(['albedo', 'emissivity', 'lai', 'lst', 'ndvi', 'ndwi', 'savi']) \
    #     .set({
    #         'system:index': SCENE_ID,
    #         'system:time_start': SCENE_TIME,
    #         'system:id': COLL_ID + SCENE_ID,
    #         'SUN_ELEVATION': SUN_ELEVATION,
    # })


# Setting etr_source and etr_band on the default image to simplify testing
#   but these do not have defaults in the Image class init
def default_image_args(albedo=0.2, emissivity=0.964, lai=1.4, lst=300,
                       ndvi=0.5, ndwi=-0.1, savi=0.5,
                       # meteorology_source_inst='NASA/NLDAS/FORA0125_H002',
                       # meteorology_source_daily='IDAHO_EPSCOR/GRIDMET',
                       # elev_source='USGS/SRTMGL1_003',
                       ndvi_cold=5, ndvi_hot=10, lst_cold=20, lst_hot=20,
                       et_reference_source=15, et_reference_band='etr',
                       et_reference_factor=0.85,
                       et_reference_resample='nearest',
                       ):
    return {
        'image': default_image(albedo=albedo, emissivity=emissivity, lai=lai,
                               lst=lst, ndvi=ndvi, ndwi=ndwi, savi=savi),
        # 'meteorology_source_inst': meteorology_source_inst,
        # 'meteorology_source_daily': meteorology_source_daily,
        # 'elev_source': elev_source,
        'ndvi_cold': ndvi_cold, 'ndvi_hot': ndvi_hot,
        'lst_cold': lst_cold, 'lst_hot': lst_hot,
        'et_reference_source': et_reference_source,
        'et_reference_band': et_reference_band,
        'et_reference_factor': et_reference_factor,
        'et_reference_resample': et_reference_resample,
    }


def default_image_obj(albedo=0.2, emissivity=0.964, lai=1.4, lst=300,
                      ndvi=0.5, ndwi=-0.1, savi=0.5,
                      # meteorology_source_inst='NASA/NLDAS/FORA0125_H002',
                      # meteorology_source_daily='IDAHO_EPSCOR/GRIDMET',
                      # elev_source='USGS/SRTMGL1_003',
                      ndvi_cold=5, ndvi_hot=10, lst_cold=20, lst_hot=10,
                      et_reference_source=15, et_reference_band='etr',
                      et_reference_factor=0.85,
                      et_reference_resample='nearest',
                      ):
    return geesebal.Image(**default_image_args(
        albedo=albedo, emissivity=emissivity, lai=lai, lst=lst,
        ndvi=ndvi, ndwi=ndwi, savi=savi,
        # meteorology_source_inst=meteorology_source_inst,
        # meteorology_source_daily=meteorology_source_daily,
        # elev_source=elev_source,
        ndvi_cold=ndvi_cold, ndvi_hot=ndvi_hot,
        lst_cold=lst_cold, lst_hot=lst_hot,
        et_reference_source=et_reference_source,
        et_reference_band=et_reference_band,
        et_reference_factor=et_reference_factor,
        et_reference_resample=et_reference_resample,
    ))


def test_Image_init_default_parameters():
    m = geesebal.Image(default_image())
    assert m._meteorology_source_inst == 'NASA/NLDAS/FORA0125_H002'
    assert m._meteorology_source_daily == 'IDAHO_EPSCOR/GRIDMET'
    assert m._elev_source == 'USGS/SRTMGL1_003'
    assert m._ndvi_cold == 5
    assert m._ndvi_hot == 10
    assert m._lst_cold == 20
    assert m._lst_hot == 20
    # assert m._et_reference_source == None
    # assert m._et_reference_band == None
    # assert m._et_reference_factor == None
    # assert m._et_reference_resample == None


# CGM - I should write image properties and image bands fixures
def test_Image_ndvi_properties():
    """Test band name and if properties are set on the image"""
    output = utils.getinfo(default_image_obj().ndvi)
    assert output['bands'][0]['id'] == 'ndvi'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_ndvi_defaults(expected=0.1, tol=0.001):
    output = utils.point_image_value(
        ee.Image(default_image_obj(ndvi=expected).ndvi), TEST_POINT)
    assert abs(output - expected) <= tol


def test_Image_et_properties():
    """Test band name and if properties are set on the image"""
    output = utils.getinfo(default_image_obj().et)
    assert output['bands'][0]['id'] == 'et'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


# CGM - This test probably won't work since running the model with an image
#   that is constant is going to work
# def test_Image_et_defaults(expected=0, tol=0.001):
#     output = utils.point_image_value(
#         ee.Image(default_image_obj().et), TEST_POINT)
#     assert abs(output - expected) <= tol


def test_Image_et_reference_properties():
    """Test if properties are set on the reference ET image"""
    output =  utils.getinfo(default_image_obj().et_reference)
    print(default_image_obj())
    assert output['bands'][0]['id'] == 'et_reference'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.parametrize(
    'source, band, factor, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 1, TEST_POINT, 11.2],
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 0.85, TEST_POINT, 11.2 * 0.85],
        ['projects/earthengine-legacy/assets/projects/climate-engine/cimis/daily',
         'ETr_ASCE', 1, TEST_POINT, 10.125],
        [10, 'FOO', 1, TEST_POINT, 10.0],
        [10, 'FOO', 0.85, TEST_POINT, 8.5],
    ]
)
def test_Image_et_reference_sources(source, band, factor, xy, expected,
                                    tol=0.001):
    """Test getting reference ET values for a single date at a real point"""
    output = utils.point_image_value(default_image_obj(
        et_reference_source=source, et_reference_band=band,
        et_reference_factor=factor).et_reference, xy)
    assert abs(output - expected) <= tol


def test_Image_et_fraction_properties():
    """Test if properties are set on the ET fraction image"""
    output =  utils.getinfo(default_image_obj().et_fraction)
    assert output['bands'][0]['id'] == 'et_fraction'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_mask_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().mask)
    assert output['bands'][0]['id'] == 'mask'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


# def test_Image_mask_values():
#     assert utils.constant_image_value(default_image_obj().mask) == 1


def test_Image_time_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().time)
    assert output['bands'][0]['id'] == 'time'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_time_values():
    """The time band should have the 0 UTC time in it for interpolation"""
    output = utils.point_image_value(default_image_obj().time, TEST_POINT)
    assert output == utils.millis(SCENE_0UTC_DT)


# def test_Image_calculate_properties():
#     """Test if properties are set on the output image"""
#     output =  utils.getinfo(default_image_obj().calculate(['ndvi']))
#     assert output['properties']['system:index'] == SCENE_ID
#     assert output['properties']['system:time_start'] == SCENE_TIME
#     assert output['properties']['image_id'] == COLL_ID + SCENE_ID
#
#
# def test_Image_calculate_variables_default():
#     output = utils.getinfo(default_image_obj().calculate())
#     assert set([x['id'] for x in output['bands']]) == {'ndvi', 'lst', 'et', 'et_fraction'}
#     # assert set([x['id'] for x in output['bands']]) == {'et'}
#
#
# def test_Image_calculate_variables_custom():
#     variables = {'ndvi'}
#     output = utils.getinfo(default_image_obj().calculate(variables))
#     assert set([x['id'] for x in output['bands']]) == variables
#
#
# def test_Image_calculate_variables_all():
#     variables = {'et', 'et_fraction', 'lst', 'ndvi', 'mask', 'time'}
#     # variables = {'et', 'et_fraction', 'et_reference', 'mask', 'ndvi', 'time'}
#     output = utils.getinfo(default_image_obj().calculate(
#         variables=list(variables)))
#     assert set([x['id'] for x in output['bands']]) == variables
#
#
# def test_Image_calculate_values():
#     """Test if the calculate method returns values"""
#     output_img = default_image_obj().calculate(['et'])
#     # output_img = default_image_obj().calculate(['et', 'et_reference', 'et_fraction'])
#     assert utils.constant_image_value(output_img.select(['et'])) > 0
#     # assert utils.constant_image_value(output_img.select(['et_reference'])) > 0
#     # assert utils.constant_image_value(output_img.select(['et_fraction'])) > 0
#
#
# def test_Image_calculate_variables_valueerror():
#     """Test if calculate method raises a valueerror for invalid variables"""
#     with pytest.raises(ValueError):
#         utils.getinfo(default_image_obj().calculate(['FOO']))



def test_Image_from_landsat_c1_sr_default_image():
    """Test that the classmethod is returning a class object"""
    output = geesebal.Image.from_landsat_c1_sr(
        ee.Image('LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716'))
    assert type(output) == type(default_image_obj())


@pytest.mark.parametrize(
    'image_id',
    [
        # 'LANDSAT/LT04/C01/T1_SR/LT04_044033_19830812',
        'LANDSAT/LT05/C01/T1_SR/LT05_044033_20110716',
        'LANDSAT/LE07/C01/T1_SR/LE07_044033_20170708',
        'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716',
    ]
)
def test_Image_from_landsat_c1_sr_landsat_image(image_id):
    """Test instantiating the class from a real Landsat images"""
    output = utils.getinfo(geesebal.Image.from_landsat_c1_sr(
        ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


# CGM - I'm not sure why these don't raise exceptions
# def test_Image_from_landsat_c1_sr_exception():
#     """Test instantiating the class for an invalid image ID"""
#     with pytest.raises(Exception):
#         utils.getinfo(geesebal.Image.from_landsat_c1_sr(ee.Image('DEADBEEF'))._index)


def test_Image_from_landsat_c2_sr_default_image():
    """Test that the classmethod is returning a class object"""
    output = geesebal.Image.from_landsat_c2_sr(
        ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_038031_20130828'))
    assert type(output) == type(default_image_obj())


@pytest.mark.parametrize(
    'image_id',
    [
        # 'LANDSAT/LT04/C02/T1_L2/LT04_044033_19830812',
        'LANDSAT/LT05/C02/T1_L2/LT05_044033_20110716',
        'LANDSAT/LE07/C02/T1_L2/LE07_044033_20170708',
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
    ]
)
def test_Image_from_landsat_c2_sr_landsat_image(image_id):
    """Test instantiating the class from a real Landsat images"""
    output = utils.getinfo(geesebal.Image.from_landsat_c2_sr(ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


# def test_Image_from_landsat_c2_sr_exception():
#     """Test instantiating the class for an invalid image ID"""
#     with pytest.raises(Exception):
#         utils.getinfo(geesebal.Image.from_landsat_c2_sr(ee.Image('DEADBEEF'))._index)


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716',
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
    ]
)
def test_Image_from_image_id(image_id):
    """Test instantiating the class using the from_image_id method"""
    output = utils.getinfo(geesebal.Image.from_image_id(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]
    assert output['properties']['image_id'] == image_id


def test_Image_from_method_kwargs():
    """Test that the init parameters can be passed through the helper methods"""
    assert geesebal.Image.from_landsat_c1_sr(
        'LANDSAT/LC08/C01/T1_SR/LC08_042035_20150713',
        elev_source='FOO')._elev_source == 'FOO'
