import math

import ee

# from openet.geesebal import openet_landsat as landsat
from openet.geesebal import utils

DEG2RAD = math.pi / 180.0


def et(image, ndvi, ndwi, lst, albedo, emissivity, savi,
       meteo_inst_source, meteo_daily_source, elev_product,
       ndvi_cold, ndvi_hot, lst_cold, lst_hot,
       time_start, geometry_image, proj, coords
       ):

    """
    Daily Evapotranspiration [mm day-1].

    Parameters
    ----------
    image : ee.Image
        Landsat image.
    ndvi : ee.Image
        Normalized difference vegetation index.
    ndwi : ee.Image
        Normalized difference water index.
    lst : ee.Image
        Land Surface Temperature [K].
    albedo : ee.Image
        Surface albedo.
    emissivity : ee.Image
        Broad-band surface emissivity.
    savi : ee.Image
        Soil-adjusted vegetation index.
    meteo_inst_source : ee.ImageCollection
        Meteorological dataset [inst]
    meteo_daily_source : ee.ImageCollection
        Meteorological dataset [daily]
    elev_product : ee.Image
    ndvi_cold : ee.Number, int
        NDVI Percentile value to determinate cold pixel.
    ndvi_hot : ee.Number, int
        NDVI Percentile value to determinate hot pixel.
    lst_cold : ee.Number, int
        Lst Percentile value to determinate cold pixel.
    lst_hot : ee.Number, int
        Lst Percentile value to determinate hot pixel.
    time_start : str
        Image property: time start of the image.
    geometry_image : ee.Geometry
        Image geometry.
    proj : ee.Image
        Landsat image projection.
    coords : ee.Image
        Landsat image Latitude and longitude.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Laipelt2021] L. Laipelt, R. Kayser, A. Fleischmann, A. Ruhoff,
        W. Bastiaanssen, T. Erickson, F. Melton,
        Long-term monitoring of evapotranspiration using the SEBAL
        algorithm and Google Earth Engine cloud computing,
        ISPRS Journal of Photogrammetry and Remote Sensing, Vol 178,
        https://doi.org/10.1016/j.isprsjprs.2021.05.018

    """
    # Image properties
    date = ee.Date(time_start)
    year = ee.Number(date.get('year'))
    month = ee.Number(date.get('month'))
    # day = ee.Number(date.get('day'))
    hour = ee.Number(date.get('hour'))
    minutes = ee.Number(date.get('minutes'))

    # Endmembers
    p_top_NDVI= ee.Number(ndvi_cold)
    p_coldest_Ts = ee.Number(lst_cold)
    p_lowest_NDVI = ee.Number(ndvi_hot)
    p_hottest_Ts = ee.Number(lst_hot)

    # Meteorology parameters
    tmin, tmax, tair, ux, rh, rso_inst, rso24h = meteorology(
        time_start, meteo_inst_source, meteo_daily_source,
    )

    # Elevation data
    dem_product = ee.Image(elev_product)
    elev = dem_product.select('elevation')

    # Sun elevation
    sun_elevation = ee.Number(image.get('SUN_ELEVATION'))

    # Terrain cos
    cos_zn = cos_terrain(time_start, dem_product, hour, minutes, coords)
    # LL: iterative process or endmembers selections may return empty values
    # in this case, return an empty image instead of broken the code
    try:
        # Air temperature correction
        tair_dem = tair_dem_correction(tmin, tmax, elev)

        # Land surface temperature correction
        lst_dem = lst_correction(
            time_start, lst, elev, tair_dem, rh, sun_elevation,
            hour, minutes, coords)

        # Cold pixel for wet conditions repretation of the image
        d_cold_pixel = cold_pixel(
            ndvi, ndwi, lst_dem, year, month, p_top_NDVI, p_coldest_Ts,
            geometry_image, coords, proj, elev
        )

        ts_cold_number = ee.Number(d_cold_pixel.get('temp'))
        ts_cold_elev = ee.Number(d_cold_pixel.get('elev'))

        # Instantaneous net radiation using reanalysis dataset
        rad_inst = radiation_inst(
            elev, lst, emissivity, albedo, tair, rh, rso_inst,
            sun_elevation, cos_zn, ts_cold_elev
        )

        # Instantaneous soil heat flux (g)
        g_inst = soil_heat_flux(rad_inst, ndvi, albedo, lst_dem, ndwi)

        # Daily ney radiation (rn)
        rad_24h = radiation_24h(
            time_start, tmax, tmin, elev,
            sun_elevation, cos_zn, rso24h
        )

        # Hot pixel
        d_hot_pixel = fexp_hot_pixel(
            time_start, ndvi, ndwi, lst_dem, rad_inst, g_inst, year, month,
             p_lowest_NDVI, p_hottest_Ts, geometry_image, coords, proj
        )

        # Instantaneous sensible heat flux (h)
        h_inst = sensible_heat_flux(
            savi, ux, ts_cold_number, d_hot_pixel, lst_dem,
            lst, elev, geometry_image
        )

        # Daily evapotranspiration (et)
        et_24hr = daily_et(h_inst, g_inst, rad_inst, lst_dem, rad_24h)

    except Exception as e:
        # CGM - We should probably log the exception so the user knows,
        #   but this will cause problems when mapping over a collection
        print(f'Unhandled Exception: {e}')

        # Return a masked image
        et_24hr = ee.Image.constant(0).updateMask(0).rename('et')

    return et_24hr


def meteorology(time_start, meteo_inst_source, meteo_daily_source):

    """
    Parameters
    ----------
    time_start : str
        Image property: time start of the image.
    meteo_inst_source: ee.ImageCollection, str
        Instantaneous meteorological data.
    meteo_daily_source :  ee.ImageCollection, str
        Daily meteorological data.

    Returns
    -------
    ee.Image

    Notes
    -----
    Accepted collections:
    Inst : NASA/NLDAS/FORA0125_H002
    Daily : IDAHO_EPSCOR/GRIDMET

    References
    ----------
    """
    time_start = ee.Number(time_start)

    meteorology_daily = ee.ImageCollection(meteo_daily_source)\
        .filterDate(ee.Date(time_start).advance(-1, 'day'), ee.Date(time_start))
    meteorology_inst_collection = ee.ImageCollection(meteo_inst_source)

    # Linear interpolation
    previous_time = time_start.subtract(2*60*60*1000)
    next_time = time_start.add(2*60*60*1000)
    previous_image = meteorology_inst_collection\
        .filterDate(previous_time, time_start)\
        .limit(1, 'system:time_start', False).first()
    next_image = meteorology_inst_collection\
        .filterDate(time_start, next_time)\
        .limit(1, 'system:time_start', True).first()
    image_previous_time = ee.Number(previous_image.get('system:time_start'))
    image_next_time = ee.Number(next_image.get('system:time_start'))
    delta_time = time_start.subtract(image_previous_time)\
        .divide(image_next_time.subtract(image_previous_time))

    # Daily variables
    # Incoming shorwave down [W m-2]
    swdown24h = meteorology_daily.select('srad').first().rename('short_wave_down')

    tmin = meteorology_daily.select('tmmn').first().rename('tmin')
    tmax = meteorology_daily.select('tmmx').first().rename('tmax')

    # Instantaneous short wave radiation [W m-2]
    rso_inst = next_image.select('shortwave_radiation')\
        .subtract(previous_image.select('shortwave_radiation'))\
        .multiply(delta_time).add(previous_image.select('shortwave_radiation'))\
        .rename('rso_inst')

    # Specific humidity [Kg Kg-1]
    q_med = next_image.select('specific_humidity')\
        .subtract(previous_image.select('specific_humidity'))\
        .multiply(delta_time).add(previous_image.select('specific_humidity'))

    # Air temperature [K]
    tair_c = next_image.select('temperature')\
        .subtract(previous_image.select('temperature'))\
        .multiply(delta_time).add(previous_image.select('temperature'))\
        .rename('tair')

    # Wind speed u [m s-1]
    wind_u = next_image.select('wind_u')\
        .subtract(previous_image.select('wind_u'))\
        .multiply(delta_time).add(previous_image.select('wind_u'))

    # Wind speed u [m s-1]
    wind_v = next_image.select('wind_v')\
        .subtract(previous_image.select('wind_v'))\
        .multiply(delta_time).add(previous_image.select('wind_v'))

    wind_med = wind_u.expression(
        'sqrt(ux_u ** 2 + ux_v ** 2)', {'ux_u': wind_u, 'ux_v': wind_v},
    ).rename('ux')

    # Wind speed [m s-1] (FAO56 Eqn 47)
    wind_med = wind_med.expression(
        'ux * (4.87) / log(67.8 * z - 5.42)', {'ux': wind_med, 'z': 10.0})

    # Pressure [kPa]
    p_med = next_image.select('pressure')\
        .subtract(previous_image.select('pressure'))\
        .multiply(delta_time).add(previous_image.select('pressure'))\
        .divide(ee.Number(1000))

    # Actual vapor pressure [kPa] (Shuttleworth Eqn 2.10)
    ea = p_med.expression('(1 / 0.622) * Q * P', {'Q': q_med, 'P': p_med})

    # Saturated vapor pressure [kPa] (FAO56 Eqn 11)
    esat = tair_c.expression(
        '0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))', {'T_air': tair_c})

    # Relative humidity (%)  (FAO56 Eqn 10)
    rh = ea.divide(esat).multiply(100).rename('RH')

    # Resample
    tmin = tmin.subtract(273.15).resample('bilinear')
    tmax = tmax.subtract(273.15).resample('bilinear')
    rso_inst = rso_inst.resample('bilinear')
    tair_c = tair_c.resample('bilinear')
    wind_med = wind_med.resample('bilinear')
    rh = rh.resample('bilinear')
    swdown24h = swdown24h.resample('bilinear')

    return [tmin, tmax, tair_c, wind_med, rh, rso_inst, swdown24h]


def tao_sw(dem, tair, rh, sun_elevation, cos_zn):

    """
    Correct declivity and aspect effects from Land Surface Temperature.

    Parameters
    ----------
    dem : ee.Image
        Elevation product data [m].
    tair : ee.Image
        Air temperature [Celsius].
    rh : ee.Image
        Relative Humidity [%]
    sun_elevation : ee.Number, ee.Image
        Sun elevation angle.
    cos_zn : ee.Number, ee.Image
        Solar zenith angle cos.

    Returns
    -------
    ee.Image

    References
    ----------
    """
    # Atmospheric pressure [kPa] (FAO56 Eqn 7)
    pres = dem.expression(
        '101.3 * ((293 - (0.0065 * Z)) / 293) ** 5.26 ', {'Z': dem})

    # Saturated vapor pressure [kPa] (FAO56 Eqn 11)
    es = tair.expression(
        '0.6108 * exp((17.27 * tair) / (tair + 237.3))', {'tair': tair})

    # Actual vapor pressure [kPa]  (FAO56 Eqn 10)
    ea = es.multiply(rh).divide(100).rename('ea')

    # Water in the atmosphere [mm] (Garrison and Adler (1990))
    w = ea.expression(
        '(0.14 * EA * PATM) + 2.1', {'PATM': pres, 'EA': ea})

    # Solar zenith angle over a horizontal surface
    solar_zenith = ee.Number(90).subtract(sun_elevation)

    solar_zenith_radians = solar_zenith.multiply(DEG2RAD)
    # Cos only in flat areas
    #cos_theta = solar_zenith_radians.cos()

    # Broad-band atmospheric transmissivity (ASCE-EWRI (2005))
    tao_sw_img = pres.expression(
        '0.35 + 0.627 * exp(((-0.00146 * P) / (Kt * ct)) - (0.075 * (W / ct) ** 0.4))',
        {'P': pres, 'W': w, 'Kt': 1.0, 'ct': cos_zn},
    )

    return tao_sw_img.rename('tao_sw')


def cos_terrain(time_start, dem, hour, minutes, coords):

    """
    Cosine zenith angle elevation (Allen et al. (2006)).

    Parameters
    ----------
    time_start : str
        Image property: time start of the image.
    dem : ee.Image
        Elevation product data [m].
    hour : ee.Number, int
        Hour.
    minutes : ee.Number, int
        Minutes.
    coords : ee.Image
        Latitude and longitude of the image.

    Returns
    -------
    ee.Image

    References
    ----------
    """
    # Day of the year
    doy = ee.Date(time_start).getRelative('day', 'year').add(1)

    # Slope and aspect
    slope_aspect = ee.Terrain.products(dem)

    # Variables
    B = doy.subtract(81).multiply(360 / 365)
    delta = ee.Number(23.45 * DEG2RAD).sin().multiply(B.multiply(DEG2RAD).sin()).asin()
    s = slope_aspect.select('slope').multiply(DEG2RAD)
    gamma = slope_aspect.select('aspect').subtract(180).multiply(DEG2RAD)
    phi = coords.select('latitude').multiply(DEG2RAD)

    # Constants
    delta = ee.Image(delta)
    a = delta.sin().multiply(phi.cos()).multiply(s.sin()).multiply(gamma.cos())\
        .subtract(delta.sin().multiply(phi.sin().multiply(s.cos())))
    b = delta.cos().multiply(phi.cos()).multiply(s.cos())\
        .add(delta.cos().multiply(phi.sin().multiply(s.sin()).multiply(gamma.cos())))
    c = delta.cos().multiply(s.sin()).multiply(gamma.sin())

    # Centroid image
    longitude_center = coords.select('longitude')

    delta_gtm = longitude_center.divide(ee.Image(15)).int()

    # Local hour time
    lht = delta_gtm.add(hour).add(minutes.divide(60))

    w = lht.subtract(12).multiply(15).multiply(DEG2RAD)

    # Cosine  zenith angle elevation
    cos_zn = w.expression(
        '-a + b * w_cos + c * w_sin',
        {'a': a, 'b': b, 'c': c, 'w_cos': w.cos(), 'w_sin': w.sin()},
    )

    return cos_zn

def tair_dem_correction(tmin, tmax, dem):

    """
    Correct Air temperature for mountain areas

    Parameters
    ----------
    tmin : ee.Image
        Minimum temperature [Celsius].
    tmax : ee.Image
        Maximum temperature [Celsius].
    dem : ee.Image
        Elevation product data [m].

    Returns
    -------
    ee.Image

    References
    ---------
    """
    tmin_dem = tmin.expression(
        'temp - 0.0065 * (dem - alt_meteo)',
        {'temp': tmin, 'dem': dem.select('elevation'), 'alt_meteo': ee.Number(2)})
    tmax_dem = tmax.expression(
        'temp - 0.0065 * (dem - alt_meteo)',
        {'temp': tmax, 'dem': dem.select('elevation'), 'alt_meteo': ee.Number(2)})

    tair_dem = tmin_dem.add(tmax_dem).divide(2)

    return tair_dem.rename('tair_dem')


def lst_correction(time_start, lst, dem, tair, rh, sun_elevation, hour, minutes,
                   coords):

    """
    Correct declivity and aspect effects from Land Surface Temperature.

    Parameters
    ----------
    time_start : str
        Image property: time start of the image.
    lst : ee.Image
        Land surface temperature [K]
    dem : ee.Image
        Elevation product data [m].
    tair : ee.Image
        Air temperature [Celsius].
    rh : ee.Image
        Relative Humidity [%].
    sun_elevation : ee.Number, ee.Image
        Sun elevation angle.
    hour : ee.Number, int
        Hour.
    minutes : ee.Number, int
        Minutes.
    coords : ee.Image
        Landsat image Latitude and longitude.

    Returns
    -------
    ee.Image

    References
    ----------
    """
    # Solar constant [W m-2]
    gsc = ee.Number(1367)

    # Day of the year
    doy = ee.Date(time_start).getRelative('day', 'year').add(1)

    # Inverse relative distance earth-sun (FAO56 Eqn 23)
    dr = doy.multiply(2 * math.pi / 365).cos().multiply(0.033).add(1)

    # Atmospheric pressure [kPa] (FAO56 Eqn 7)
    pres = lst.expression(
        '101.3 * ((293 - (0.0065 * Z)) / 293) ** 5.26 ', {'Z': dem})

    # Solar zenith angle over a horizontal surface
    solar_zenith = ee.Number(90).subtract(sun_elevation)

    solar_zenith_radians = solar_zenith.multiply(DEG2RAD)
    cos_theta = solar_zenith_radians.cos()

    # Air density [Kg m-3]
    air_dens = lst.expression(
        '(1000 * Pair)/(1.01 * LST * 287)', {'Pair': pres, 'LST': lst})

    # Temperature lapse rate (0.0065)
    temp_lapse_rate = ee.Number(0.0065)

    # Added temperature lapse rate
    temp_corr = lst.add(dem.select('elevation').multiply(temp_lapse_rate))

    cos_zn = cos_terrain(time_start, dem, hour, minutes, coords)

    # Broad-band atmospheric transmissivity (ASCE-EWRI (2005))
    tao_sw_img = tao_sw(dem, tair, rh, sun_elevation, cos_zn)

    # CGM - Check if the following expression can be simplified as:
    # 'Temp_corr + (Gsc * dr * Transm_corr) * (cos_zn - cos_zenith_flat) / (air_dens * 1004 * 0.050)',

    # Corrected Land Surface temperature [K] (Zaafar and Farah (2020) Eqn 2)
    lst_dem = lst.expression(
        'Temp_corr + (Gsc * dr * Transm_corr * cos_zn - Gsc * dr * Transm_corr * cos_zenith_flat) / '
        '(air_dens * 1004 * 0.050)',
        {'Temp_corr': temp_corr, 'Gsc': gsc, 'dr': dr, 'Transm_corr': tao_sw_img,
         'cos_zenith_flat': cos_theta, 'cos_zn': cos_zn, 'air_dens': air_dens},
    )

    return lst_dem.rename('lst_dem')


def lc_mask(month, year, geometry_image, mask_img):

    """
    Filtering pre-candidates pixels using a Land cover mask.

    Parameters
    ----------
    month : ee.Number, int
        Month.
    year : ee.Number, int
        Year.
    geometry_image : ee.Geometry
        Landsat image geometry.
    mask_img : ee.Image

    Returns
    -------
    ee.Image

    References
    ----------
    """
    # Conditions
    cdl_year_min = 2008
    cdl_year_max = 2021

    year_condition = ee.Number(year).max(cdl_year_min).min(cdl_year_max)

    isWinter = ee.Number(month.eq(1).Or(month.eq(2)).Or(month.eq(3))
                         .Or(month.eq(11)).Or(month.eq(12)))

    start = ee.Date.fromYMD(year_condition, 1, 1)
    end = ee.Date.fromYMD(year_condition, 12, 31)

    # Select classification corresponding to the year of the imageq
    lc = ee.ImageCollection('USDA/NASS/CDL').select('cropland')\
        .filter(ee.Filter.date(start, end)).first()

    # Filter cropland classes 1
    crop1 = lc.updateMask(lc.lt(61))
    crop1 = crop1.where(crop1, 1).unmask(0)

    # Filter cropland classes 2
    crop2 = lc.updateMask(lc.gte(196))
    crop2 = crop2.where(crop2, 1).unmask(0)

    # Land cover mask - total croplands
    # CGM - Should probably rename to something different than the function name
    lc_mask = crop1.add(crop2)

    lc_mask = lc_mask.updateMask(lc_mask.eq(1))

    # Check if there are more than 3000 pixels in the land cover masks
    # otherwise land cover mask is not applied (return a full scene mask)
    count_land_cover_pixels = lc_mask.rename('land_cover_pixels')\
        .reduceRegion(reducer=ee.Reducer.count(), scale=30,
                      geometry=geometry_image, maxPixels=10e14)
    n_count_lc = ee.Number(count_land_cover_pixels.get('land_cover_pixels'))

    mask = ee.Algorithms.If(n_count_lc.gte(3000), lc_mask, mask_img)

    mask = ee.Algorithms.If(isWinter.eq(1), mask_img, mask)
    # mask = landsat_image.select(0).updateMask(1)

    return ee.Image(mask)


def homogeneous_mask(ndvi, proj):

    """
    Homogeneous mask for endmembers selection (Allen et al. (2013)).

    Parameters
    ----------

    ndvi : ee.Image
        Normalized difference vegetation index (ndvi).
    proj : ee.Dictionary
        Landsat image projection.

    Returns
    -------
    ee.Image

    References
    ----------

    """

    sd_ndvi = ndvi\
        .reduceNeighborhood(reducer=ee.Reducer.stdDev(),
                            kernel=ee.Kernel.square(radius=3, units='pixels'),
                            skipMasked=False)\
        .reproject(proj)\
        .updateMask(1)

    sd_mask = sd_ndvi.updateMask(sd_ndvi.lte(0.15))

    return ee.Image(sd_mask)


def cold_pixel(ndvi, ndwi, lst_dem, year, month, ndvi_cold, lst_cold,
               geometry_image, coords, proj, dem):

    """
    Simplified CIMEC method to select the cold pixel

    Parameters
    ----------
    ndvi : ee.Image
        Normalized difference vegetation index.
    ndwi : ee.Image
        Normalized difference water index.
    lst_dem : ee.Image
        Land surface temperature [K].
    year : ee.Number, int
        Year.
    month : ee.Number, int
        Month.
    ndvi_cold : ee.Number, int
        NDVI Percentile value to determinate cold pixel.
    lst_cold : ee.Number, int
        LST Percentile value to determinate cold pixel.
    geometry_image : ee.Geometry
        Landsat image geometry.
    coords : ee.Image
        Latitude and longitude coordinates of the image.
    proj : ee.Dictionary
        Landsat image projection.
    dem : ee.Image
        Elevation data [m]

    Returns
    -------
    ee.Dictionary

    Notes
    -----
    Based on Allen et al (2013) procedure to represent extreme conditions
    to use in METRIC (adaptable for SEBAL) using endmembers candidates from
    pre-defined percentiles of LST and NDVI.

    References
    ----------

    ..[Allen2013] Allen, R.G., Burnett, B., Kramber, W., Huntington, J.,
        Kjaersgaard, J., Kilic, A., Kelly, C., Trezza, R., (2013).
        Automated Calibration of the METRIC-Landsat Evapotranspiration Process.
        JAWRA J. Am. Water Resour. Assoc. 49, 563–576.

    """
    # Pre-filter
    pos_ndvi = ndvi.updateMask(ndvi.gt(0)).rename('post_ndvi')
    ndvi_neg = pos_ndvi.multiply(-1).rename('ndvi_neg')

    lst_neg = lst_dem.multiply(-1).rename('lst_neg').rename('lst_neg')
    lst_nw = lst_dem.updateMask(ndwi.lte(0)).rename('lst_nw')

    # Crea a 0 mask
    mask = ndvi.select(0).updateMask(1)
    # Land cover mask
    land_cover_mask = lc_mask(month, year, geometry_image, mask)

    # Creates a homogeneous ndvi mask
    stdev_ndvi = homogeneous_mask(ndvi, proj)

    images = pos_ndvi.addBands([ndvi, ndvi_neg, pos_ndvi, lst_neg, lst_nw,
                                coords, dem.toFloat()])

    d_perc_top_NDVI = images.select('ndvi_neg')\
        .updateMask(land_cover_mask)\
        .updateMask(stdev_ndvi)\
        .reduceRegion(reducer=ee.Reducer.percentile([ndvi_cold]),
                      geometry=geometry_image, scale=30, maxPixels=1e9)\
        .combine(ee.Dictionary({'ndvi_neg': 100}), overwrite=False)

    n_perc_top_NDVI = ee.Number(d_perc_top_NDVI.get('ndvi_neg'))

    i_top_NDVI = images.updateMask(land_cover_mask).updateMask(stdev_ndvi)\
        .updateMask(images.select('ndvi_neg').lte(n_perc_top_NDVI))

    d_perc_low_LST = i_top_NDVI.select('lst_nw')\
        .updateMask(land_cover_mask)\
        .updateMask(stdev_ndvi)\
        .reduceRegion(reducer=ee.Reducer.percentile([lst_cold]),
                      geometry=geometry_image, scale=30, maxPixels=1e9)\
        .combine(ee.Dictionary({'lst_nw': 350}), overwrite=False)

    n_perc_low_LST = ee.Number(d_perc_low_LST.get('lst_nw'))

    i_cold_lst = i_top_NDVI.updateMask(land_cover_mask)\
        .updateMask(i_top_NDVI.select('lst_nw').lte(n_perc_low_LST))

    # Filtes
    c_lst_cold20 = i_cold_lst.updateMask(images.select('lst_nw').gte(200))
    c_lst_cold20_int = c_lst_cold20.select('lst_nw').min(1).max(1).int().rename('int')
    c_lst_cold20 = c_lst_cold20.addBands(c_lst_cold20_int)

    sum_final_cold_pix = c_lst_cold20.select('int')\
        .reduceRegion(reducer=ee.Reducer.sum(), geometry=geometry_image,
                      scale=30, maxPixels=1e9)
    n_sum_final_cold_pix = ee.Number(sum_final_cold_pix.get('int'))

    def function_def_pixel(f):
        return f.setGeometry(ee.Geometry.Point([f.get('longitude'), f.get('latitude')]))

    # Get Cold Pixel (random)
    fc_cold_pix = c_lst_cold20.stratifiedSample(1, 'int', geometry_image, 30)\
        .map(function_def_pixel)
    n_Ts_cold = ee.Number(fc_cold_pix.aggregate_first('lst_nw'))
    n_long_cold = ee.Number(fc_cold_pix.aggregate_first('longitude'))
    n_lat_cold = ee.Number(fc_cold_pix.aggregate_first('latitude'))
    n_ndvi_cold = ee.Number(fc_cold_pix.aggregate_first('ndvi'))
    n_dem_cold = ee.Number(fc_cold_pix.aggregate_first('elevation'))

    # Dictionary
    d_cold_pixel = ee.Dictionary({
        'temp': n_Ts_cold,
        'ndvi': n_ndvi_cold,
        'x': n_long_cold,
        'y': n_lat_cold,
        'sum': n_sum_final_cold_pix,
        'elev': n_dem_cold
    }).combine(ee.Dictionary({'temp': 0, 'ndvi': 0, 'x': 0, 'y': 0, 'sum': 0, 'elev': 0}),
               overwrite=False)

    return d_cold_pixel


def radiation_inst(dem, lst, emissivity, albedo, tair, rh, swdown_inst,
                   sun_elevation, cos_terrain, ts_cold_elev):

    """
    Instantaneous Net Radiation [W m-2]

    Parameters
    ----------
    dem : ee.Image
        Digital elevation product [m].
    lst : ee.Image
        Land surface temperature [k].
    emissivity : ee.Image
        Emissivity.
    albedo : ee.Image
        Albedo.
    tair : ee.Image
        Air temperature [Celsius].
    rh : ee.Image
        Relative Humidity [%].
    swdown_inst : ee.Image
        Instantaneous Short Wave radiation [W m-2].
    sun_elevation : ee.Number, int
        Sun elevation information.
    cos_terrain : ee.Image
        Solar zenith angle cos (aspect/slope).
    ts_cold_elev : ee.Number, int
        Elevation at the cold pixel [m].

    Returns
    -------
    ee.Image

    References
    ----------

    """

    rad_long_up = lst.expression(
        'emi * 5.67e-8 * (lst ** 4)', {'emi': emissivity, 'lst': lst})

    tao_sw_img = tao_sw(dem, tair, rh, sun_elevation, cos_terrain)

    log_taosw = tao_sw_img.log()

    # CGM - This variable/line is not used
    #   Could remove ts_cold_elev function input if this line was removed
    delta_z = ee.Image((dem.select('elevation').subtract(ee.Number(ts_cold_elev)))
                       .multiply(0.0065))

    rad_long_down = lst.expression(
        '(0.85 * (- log_taosw) ** 0.09) * 5.67e-8 * (n_Ts_cold ** 4)',
        {'log_taosw': log_taosw, 'n_Ts_cold': tair.add(273.15)},
    )

    # Rso inst aspect/slope
    solar_zenith = ee.Number(90).subtract(sun_elevation)
    solar_zenith_radians = solar_zenith.multiply(DEG2RAD)
    cos_zeni = solar_zenith_radians.cos()
    swdown_inst_dem = swdown_inst.multiply(cos_terrain.divide(cos_zeni))

    rn_inst = lst.expression(
        '((1 - alfa) * Rs_down) + Rl_down - Rl_up - ((1 - e_0) * Rl_down) ',
        {'alfa': albedo, 'Rs_down': swdown_inst_dem, 'Rl_down': rad_long_down,
         'Rl_up': rad_long_up, 'e_0': emissivity},
    )

    return rn_inst.rename('rn_inst')


def soil_heat_flux(rn, ndvi, albedo, lst_dem, ndwi):

    """
    Instantaneous Soil Heat Flux [W m-2]

    Parameters
    ----------
    rn : ee.Image
        Instantaneous Net Radiation [W m-2]
    ndvi : ee.Image
        Normalized difference vegetation index.
    albedo : ee.Image
        Albedo.
    lst_dem : ee.Image
        Land surface temperature [K].
    ndwi : ee.Image
        Normalized difference water index.

    Returns
    -------
    ee.Image

    References
    ----------

    """

    g = rn.expression(
        'rn * (lst - 273.15) * (0.0038 + (0.0074 * albedo)) * '
        '(1 - 0.98 * (ndvi ** 4)) ',
        {'rn': rn, 'ndvi': ndvi, 'albedo': albedo, 'lst': lst_dem},
    )

    g = g.where(ndwi.gt(0), rn.multiply(0.5))

    return g.rename('g_inst')


def radiation_24h(time_start, tmax, tmin, elev, sun_elevation, cos_terrain, rso24h):

    """
    Daily Net radiation [W m-2] - FAO56

    Parameters
    ----------
    time_start : ee.Date
        Date information of the image.
    tmax : ee.Image
        Maximum air temperature [Celsius].
    tmin : ee.Image
        Minimum air temperature [Celsius].
    elev : ee.Image
        Digital Elevation information [m].
    sun_elevation : ee.Number, int
        Sun elevation information.
    cos_terrain : ee.Image
        Solar zenith angle cos (aspect/slope).
    rso24h : ee.Image
        Daily Short wave radiation [W m-2]

    Returns
    -------
    ee.Image

    References
    ----------
    .. [FAO56] Allen, R., Pereira, L., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration: Guidelines for computing crop water
       requirements. FAO Irrigation and Drainage Paper (Vol. 56).

    """

    # Convert to MJ m-2
    rs = rso24h.multiply(0.0864).rename('Rs')

    # Solar constant [MJ m-2]
    gsc = 0.0820

    # Day of the year
    doy = ee.Date(time_start).getRelative('day', 'year').add(1)

    # Inverse relative distance earth-sun (FAO56 Eqn 23)
    dr = tmax.expression(
        '1 + (0.033 * cos((2 * pi / 365) * doy))', {'doy': doy, 'pi': math.pi})

    # Solar declination [rad] (FAO56 Eqn 24)
    sd = tmax.expression(
        '0.40928 * sin(((2 * pi / 365) * doy) - 1.39)', {'doy': doy, 'pi': math.pi})

    # Latitude of the image
    lat = tmax.pixelLonLat().select(['latitude']).multiply(DEG2RAD)\
        .rename('latitude')

    #  Sunset hour angle [rad] (FAO56 Eqn 25)
    ws = tmax.expression('acos(-tan(Lat) * tan(Sd))', {'Lat': lat, 'Sd': sd})

    # Extraterrestrial radiation [MJ m-2 d-1] (FAO56 Eqn 21)
    rad_a = tmax.expression(
        'Ws * sin(Lat) * sin(Sd) + cos(Lat) * cos(Sd) * sin(Ws)',
        {'Ws': ws, 'Lat': lat, 'Sd': sd},
    )

    ra = tmax.expression(
        '((24 * 60) / pi) * Gsc * Dr * rad_a',
        {'pi': math.pi, 'Gsc': gsc, 'Dr': dr, 'rad_a': rad_a},
    )
    # Simplified clear sky solar formulation [MJ m-2 d-1] (FAO56 Eqn 37)
    rso = tmax.expression('(0.75 + 2E-5 * z) * Ra', {'z': elev, 'Ra': ra})

    # Net shortwave radiation [MJ m-2 d-1] (FAO56 Eqn 38)
    rns = tmax.expression('(1 - albedo) * Rs', {'Rs': rs, 'albedo': 0.23})

    # Actual vapor pressure [MJ m-2 d-1] (FAO56 Eqn 11)
    ea = tmax.expression(
        '0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))', {'T_air': tmin})

    # Rso slope/aspect
    solar_zenith = ee.Number(90).subtract(sun_elevation)
    solar_zenith_radians = solar_zenith.multiply(DEG2RAD)
    cos_zeni = solar_zenith_radians.cos()

    rso24h_dem = rso.multiply(cos_terrain.divide(cos_zeni))

    # Net longwave radiation [MJ m-2 d-1] (FAO56 Eqn 39)
    rnl = tmax.expression(
        '4.901E-9 * ((Tmax ** 4 + Tmin ** 4) / 2) * (0.34 - 0.14 * sqrt(ea)) * '
        '(1.35 * (Rs / Rso) - 0.35)',
        {'Tmax': tmax.add(273.15), 'Tmin': tmin.add(273.15), 'ea': ea,
         'Rs': rs, 'Rso': rso24h_dem},
    )

    # Net radiation [MJ m-2 d-1] (FAO56 Eqn 40)
    rn = tmax.expression('Rns - Rnl', {'Rns': rns, 'Rnl': rnl})

    # Convert to W m-2
    rn = rn.multiply(11.6)

    return rn.rename('rad_24h')


def fexp_hot_pixel(time_start, ndvi, ndwi, lst_dem, rn, g, year, month,
                   ndvi_hot, lst_hot, geometry_image, coords, proj):

    """
    Simplified CIMEC method to select the hot pixel


    Parameters
    ----------
    time_start : ee.Date
        Date information of the image.
    ndvi : ee.Image
        Normalized difference vegetation index.
    ndwi : ee.Image
        Normalized difference water index.
    lst_dem : ee.Image
        Land surface temperature [K].
    rn : ee.Image
        Instantaneous Net Radiation [W m-2]
    g : ee.Image
        Instantaneous Soil heat flux [W m-2]
    year : ee.Number, int
        Year.
    month : ee.Number, int
        Month.
    ndvi_hot : ee.Number, int
        NDVI Percentile value to determinate hot pixel.
    lst_hot : ee.Number, int
        LST Percentile value to determinate hot pixel.
    geometry_image : ee.Geometry
        Image geometry.
    coords : ee.Image
        Latitude and longitude coordinates of the image.
    proj : ee.Dictionary
        Landsat image projection.

    Returns
    -------
    ee.Dictionary

    Notes
    -----
    Based on Allen et al (2013) procedure to represent extreme conditions
    to use in METRIC (adaptable for SEBAL) using endmembers candidates from
    pre-defined percentiles of LST and NDVI.

    References
    ----------

    .. [Allen2013] Allen, R.G., Burnett, B., Kramber, W., Huntington, J.,
        Kjaersgaard, J., Kilic, A., Kelly, C., Trezza, R., (2013).
        Automated Calibration of the METRIC-Landsat Evapotranspiration Process.
        JAWRA J. Am. Water Resour. Assoc. 49, 563–576.
    ..
    """

    # Pre-filter
    pos_ndvi = ndvi.updateMask(ndvi.gt(0)).rename('post_ndvi')
    ndvi_neg = pos_ndvi.multiply(-1).rename('ndvi_neg')

    lst_neg = lst_dem.multiply(-1).rename('lst_neg')
    lst_nw = lst_dem.updateMask(ndwi.lte(0)).rename('lst_nw')

    # Create a 0 mask
    mask = ndvi.select(0).updateMask(1)

    # Land cover mask
    land_cover_mask = lc_mask(month, year, geometry_image, mask)

    # Create a homogeneous ndvi mask
    stdev_ndvi = homogeneous_mask(ndvi, proj)

    images = pos_ndvi.addBands([
        ndvi, ndvi_neg, rn, g, pos_ndvi, lst_neg, lst_nw, coords])

    d_perc_down_ndvi = images.select('post_ndvi')\
        .updateMask(land_cover_mask)\
        .updateMask(stdev_ndvi)\
        .reduceRegion(reducer=ee.Reducer.percentile([ndvi_hot]),
                      geometry=geometry_image, scale=30, maxPixels=1e9)\
        .combine(ee.Dictionary({'post_ndvi': 100}), overwrite=False)
    n_perc_low_NDVI = ee.Number(d_perc_down_ndvi.get('post_ndvi'))

    i_low_NDVI = images.updateMask(land_cover_mask)\
        .updateMask(images.select('post_ndvi').lte(n_perc_low_NDVI))

    d_perc_top_lst = i_low_NDVI.select('lst_neg')\
        .updateMask(land_cover_mask)\
        .updateMask(stdev_ndvi)\
        .reduceRegion(reducer=ee.Reducer.percentile([lst_hot]),
                      geometry=geometry_image, scale=30, maxPixels=1e9)\
        .combine(ee.Dictionary({'lst_neg': 350}), overwrite=False)

    n_perc_top_lst = ee.Number(d_perc_top_lst.get('lst_neg'))

    i_top_LST = i_low_NDVI.updateMask(land_cover_mask).updateMask(stdev_ndvi)\
        .updateMask(i_low_NDVI.select('lst_neg').lte(n_perc_top_lst))

    c_lst_hot_int = i_top_LST.select('lst_nw').min(1).max(1).int().rename('int')
    c_lst_hotpix = i_top_LST.addBands(c_lst_hot_int)

    sum_final_hot_pix = c_lst_hotpix.select('int')\
        .reduceRegion(reducer=ee.Reducer.sum(), geometry=geometry_image,
                      scale=30, maxPixels=1e9)
    n_sum_final_hot_pix = ee.Number(sum_final_hot_pix.get('int'))

    # Precipitation product
    gridmet = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')\
        .filterDate(ee.Date(time_start).advance(-60, 'days'), ee.Date(time_start))

    etr_60mm = gridmet.select('etr').sum()
    precipt_60mm = gridmet.select('pr').sum()
    ratio = precipt_60mm.divide(etr_60mm)

    # Temperature adjustement offset (Allen2013 Eqn 8)
    Tfac = etr_60mm.expression('2.6 - 13 * ratio', {'ratio': ratio})

    Tfac = ee.Image(Tfac.where(ratio.gt(0.2), 0)).rename('Tfac')

    c_lst_hotpix = c_lst_hotpix.addBands(Tfac)

    def function_def_pixel(f):
        return f.setGeometry(ee.Geometry.Point([f.get('longitude'), f.get('latitude')]))

    # Get Hot Pixel (random)
    fc_hot_pix = c_lst_hotpix.stratifiedSample(1, 'int', geometry_image, 30)\
        .map(function_def_pixel)

    n_Ts_hot = ee.Number(fc_hot_pix.aggregate_first('lst_nw'))
    n_long_hot = ee.Number(fc_hot_pix.aggregate_first('longitude'))
    n_lat_hot = ee.Number(fc_hot_pix.aggregate_first('latitude'))
    n_ndvi_hot = ee.Number(fc_hot_pix.aggregate_first('ndvi'))
    n_Rn_hot = ee.Number(fc_hot_pix.aggregate_first('rn_inst'))
    n_G_hot = ee.Number(fc_hot_pix.aggregate_first('g_inst'))
    n_Tfac = ee.Number(fc_hot_pix.aggregate_first('Tfac'))

    # Dictionary
    d_hot_pixel = ee.Dictionary({
        'temp': n_Ts_hot,
        'tfac': n_Tfac,
        'x': n_long_hot,
        'y': n_lat_hot,
        'rn': n_Rn_hot,
        'g': n_G_hot,
        'ndvi': n_ndvi_hot,
        'sum': n_sum_final_hot_pix,
    }).combine(ee.Dictionary({'temp': 0, 'tfac': 0, 'x': 0, 'y': 0,
                              'rn': 0, 'g': 0, 'ndvi': 0, 'sum': 0}),
               overwrite=False)

    return d_hot_pixel


def sensible_heat_flux(savi, ux, ts_cold_number, d_hot_pixel,
                       lst_dem, lst, dem, geometry_image):

    """
    Instantaneous Sensible Heat Flux [W m-2]

    Parameters
    ----------
    savi : ee.Image
        Soil-adjusted vegetation index.
    ux : ee.Image
        Wind speed [m s-1].
    ts_cold_number : ee.Number
        Cold pixel value [K].
    d_hot_pixel : ee.Dictionary
        Rn, G, lat, lon, Ts hot pixel.
    lst_dem : ee.Image
        Land surface temperature (aspect/slope correction) [K].
    lst : ee.Image
        Land surface temperature [K].
    dem : ee.Image
        Digital elevation product [m].
    geometry_image : ee.Geometry
        Image geometry.

    Returns
    -------
    ee.Image

    References
    ----------

    .. [Bastiaanssen1998] Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A.,
        Holtslag, A.A.M., 1998. A remote sensing surface energy balance
        algorithm for land (SEBAL): 1. Formulation. J. Hydrol. 212–213, 198–212.
    .. [Allen2002] Allen, R., Bastiaanssen., W.G.M. 2002.
        Surface Energy Balance Algorithms for Land. Idaho Implementation.
        Advanced Training and Users Manual. 2002.

    """
    # Vegetation height [m]
    n_veg_height = ee.Number(0.5)

    # Wind speed height [m]
    n_zx = ee.Number(2)

    # Blending height [m]
    n_height = ee.Number(200)

    # Air specific heat [J kg-1 K-1]
    n_Cp = ee.Number(1004)

    # Von Karman’s constant
    n_K = ee.Number(0.41)

    # Filtering low lalues of wind speed
    wind_speed_std = ux.rename('ux')\
        .reduceRegion(reducer=ee.Reducer.stdDev(), geometry=geometry_image,
                      scale=10000, maxPixels=1e9)\
        .combine(ee.Dictionary({'ux': 0}), overwrite=False)

    n_wind_speed_std = ee.Number(wind_speed_std.get('ux'))

    # LL : Values less than 1.5 m s-1 tend to generate instability in
    # the iterative process to estimate aerodynamic resistance.
    # Standard Deviation is added in this situations.

    ux = ux.where(ux.lt(1.5), ux.add(n_wind_speed_std))

    # Slope/ Aspect
    slope_aspect = ee.Terrain.products(dem)
    n_Ts_cold = ee.Number(ts_cold_number)
    n_Ts_hot = ee.Number(d_hot_pixel.get('temp')).subtract(ee.Number(d_hot_pixel.get('tfac')))
    n_G_hot = ee.Number(d_hot_pixel.get('g'))
    n_Rn_hot = ee.Number(d_hot_pixel.get('rn'))
    n_long_hot = ee.Number(d_hot_pixel.get('x'))
    n_lat_hot = ee.Number(d_hot_pixel.get('y'))
    p_hot_pix = ee.Geometry.Point([n_long_hot, n_lat_hot])

    # Momentum roughness length at the weather station. (Allen2002 Eqn 28)
    n_zom = n_veg_height.multiply(0.123)

    # Friction velocity at the weather station. (Allen2002 Eqn 37)
    i_ufric_ws = lst.expression(
        '(n_K * ux) / log(n_zx / n_zom)',
        {'n_K': n_K, 'n_zx': n_zx, 'n_zom': n_zom, 'ux': ux},
    )

    # Wind speed at blending height at the weather station.  (Allen2002 Eqn 29)
    i_u200 = lst.expression(
        'i_ufric_ws * log(n_height / n_zom) / n_K',
        {'i_ufric_ws': i_ufric_ws, 'n_height': n_height, 'n_zom': n_zom, 'n_K': n_K}
    )

    # Momentum roughness length for each pixel.
    i_zom = lst.expression(
        'exp((5.62 * SAVI) - 5.809)', {'SAVI': savi},
    )

    # Momentum roughness slope/aspect Correction.  (Allen2002  A12 Eqn9)
    i_zom = i_zom.expression(
        'zom * (1 + (slope - 5) / 20)',
        {'zom': i_zom, 'slope': slope_aspect.select('slope')},
    )

    # Friction velocity for each pixel. (Allen2002 Eqn 30)
    i_ufric = lst.expression(
        '(n_K * u200) / log(height / i_zom)',
        {'u200': i_u200, 'height': n_height, 'i_zom': n_zom, 'n_K': n_K},
    )

    # Heights [m] above the zero plane displacement.
    z1 = ee.Number(0.01)
    z2 = ee.Number(2)

    # Aerodynamic resistance to heat transport (Allen2002 Eqn 26)
    i_rah = i_ufric.expression(
        '(log(z2 / z1)) / (i_ufric * 0.41)',
        {'z2': z2, 'z1': z1, 'i_ufric': i_ufric},
    ).rename(['rah'])

    n_ro_hot = n_Ts_hot.multiply(-0.0046).add(2.5538)

    # Iterative Process
    # Sensible heat flux at the hot pixel
    n_H_hot = ee.Number(n_Rn_hot).subtract(ee.Number(n_G_hot))
    n = ee.Number(1)
    n_dif = ee.Number(1)
    n_dif_min = ee.Number(0.1)
    list_dif = ee.List([])
    list_dT_hot = ee.List([])
    list_rah_hot = ee.List([])
    list_coef_a = ee.List([])
    list_coef_b = ee.List([])

    for n in range(15):
        d_rah_hot = i_rah\
            .reduceRegion(reducer=ee.Reducer.first(), geometry=p_hot_pix,
                          scale=30, maxPixels=9000000000)\
            .combine(ee.Dictionary({'rah': 0}), overwrite=False)

        # LL : To avoid 'Max (NaN) cannot be less than min (NaN)' erros in
        # cases which iterative process not converge
        n_rah_hot = ee.Number(d_rah_hot.get('rah'))\
            .multiply(100).short().divide(100)

        # Near surface temperature difference in hot pixel (dT = Tz1 – Tz2)
        # dThot=Hhot*rah/(ρCp)
        n_dT_hot = n_H_hot.multiply(n_rah_hot).divide(n_ro_hot.multiply(n_Cp))

        # Near surface temperature difference in cold pixel (dT = Tz1 – Tz2)
        n_dT_cold = ee.Number(0)
        # dT = aTs + b

        # Angular coefficient
        n_coef_a = (n_dT_cold.subtract(n_dT_hot))\
            .divide(n_Ts_cold.subtract(n_Ts_hot))

        # Linear coefficient
        n_coef_b = n_dT_hot.subtract(n_coef_a.multiply(n_Ts_hot))

        # dT for each pixel
        i_dT_int = lst.expression(
            '(n_coef_a * i_lst_med_dem) + n_coef_b',
            {'n_coef_a': n_coef_a, 'n_coef_b': n_coef_b, 'i_lst_med_dem': lst_dem},
        )

        # LL - Just testing others approaches for mountains areas
        '''
        dt_std = i_dT_int.select('dt').reduceRegion(
            reducer=ee.Reducer.stdDev(),
            geometry=geometry_image,
            scale=30,
            maxPixels=1e12)


        #i_dT_int = i_dT_int.where(slope_aspect.select('slope').gt(10),i_dT_int.add(ee.Image.constant(dt_std.get('dt')).multiply(0.5)))
        '''

        # Air temperature (Ta) for each pixel (Ta = Ts-dT)
        i_Ta = lst.expression(
            'i_lst_med - i_dT_int', {'i_lst_med': lst, 'i_dT_int': i_dT_int})

        # ro=-0.0046.*Ta+2.5538
        i_ro = i_Ta.expression('(-0.0046 * i_Ta) + 2.5538', {'i_Ta': i_Ta})

        # Sensible heat flux (H) for each pixel - iteration
        i_H_int = i_dT_int.expression(
            '(i_ro * n_Cp * i_dT_int) / i_rah',
            {'i_ro': i_ro, 'n_Cp': n_Cp, 'i_dT_int': i_dT_int, 'i_rah': i_rah},
        )

        # Monin-Obukhov length (L) - iteration
        i_L_int = i_dT_int.expression(
            '-(i_ro * n_Cp * (i_ufric ** 3) * i_lst_med) / (0.41 * 9.81 * i_H_int)',
            {'i_ro': i_ro, 'n_Cp': n_Cp, 'i_ufric': i_ufric, 'i_lst_med': lst,
             'i_H_int': i_H_int},
        )

        # Limiting L values to avoid errors in rah.
        i_L_int = i_L_int.where(i_L_int.lt(-1000), 1000)

        # Stability corrections for stable conditions
        i_psim_200 = lst.expression(
            '-5 * (height / i_L_int)', {'height': 200.0, 'i_L_int': i_L_int},
        )
        i_psih_2 = lst.expression(
            '-5 * (height / i_L_int)', {'height': 2.0, 'i_L_int': i_L_int},
        )
        i_psih_01 = lst.expression(
            '-5 * (height / i_L_int)', {'height': 0.1, 'i_L_int': i_L_int},
        )

        # x for different height
        i_x200 = i_L_int.expression(
            '(1 - (16 * (height / i_L_int))) ** 0.25',
            {'height': 200.0, 'i_L_int': i_L_int}
        )
        i_x2 = i_L_int.expression(
            '(1 - (16 * (height / i_L_int))) ** 0.25',
            {'height': 2.0, 'i_L_int': i_L_int}
        )
        i_x01 = i_L_int.expression(
            '(1 - (16 * (height / i_L_int))) ** 0.25',
            {'height': 0.1, 'i_L_int': i_L_int}
        )

        # Stability corrections for unstable conditions
        i_psimu_200 = i_x200.expression(
            '2 * log((1 + i_x200) / 2) + log((1 + i_x200 ** 2) / 2) - '
            '2 * atan(i_x200) + 0.5 * pi',
            {'i_x200': i_x200, 'pi': math.pi},
        )
        i_psihu_2 = i_x2.expression(
            '2 * log((1 + i_x2 ** 2) / 2)', {'i_x2': i_x2})
        i_psihu_01 = i_x01.expression(
            '2 * log((1 + i_x01 ** 2) / 2)', {'i_x01': i_x01})

        i_psim_200 = i_psim_200.where(i_L_int.lt(0), i_psimu_200)
        i_psih_2 = i_psih_2.where(i_L_int.lt(0), i_psihu_2)
        i_psih_01 = i_psih_01.where(i_L_int.lt(0), i_psihu_01)
        i_psim_200 = i_psim_200.where(i_L_int.eq(0), 0)
        i_psih_2 = i_psih_2.where(i_L_int.eq(0), 0)
        i_psih_01 = i_psih_01.where(i_L_int.eq(0), 0)

        # Corrected value for the friction velocity.
        i_ufric = i_ufric.expression(
            '(u200 * 0.41) / (log(height / i_zom) - i_psim_200)',
            {'u200': i_u200, 'height': n_height,
             'i_zom': i_zom, 'i_psim_200': i_psim_200},
        )

        #Limiting minimum ufric values
        i_ufric = i_ufric.where(i_ufric.lt(0.02),0.02)

        # Corrected value for the aerodinamic resistance to the heat transport
        i_rah_unstable = i_rah.expression(
            '(log(z2 / z1) - psi_h2 + psi_h01) / (i_ufric * 0.41)',
            {'z2': z2, 'z1': z1, 'i_ufric': i_ufric,
             'psi_h2': i_psih_2, 'psi_h01': i_psih_01},
        ).rename('rah')

        i_rah = i_rah.where(i_L_int.lt(0),i_rah_unstable)

        if n == 1:
            n_dT_hot_old = n_dT_hot
            n_rah_hot_old = n_rah_hot
            n_dif = ee.Number(1)

        if n > 1:
            n_dT_hot_abs = n_dT_hot.abs()
            n_dT_hot_old_abs = n_dT_hot_old.abs()
            n_rah_hot_abs = n_rah_hot.abs()
            n_rah_hot_old_abs = n_rah_hot_old.abs()
            n_dif = (n_dT_hot_abs.subtract(n_dT_hot_old_abs)
                     .add(n_rah_hot_abs).subtract(n_rah_hot_old_abs)).abs()
            n_dT_hot_old = n_dT_hot
            n_rah_hot_old = n_rah_hot
            # insert each iteration value into a list

        list_dif = list_dif.add(n_dif)
        list_coef_a = list_coef_a.add(n_coef_a)
        list_coef_b = list_coef_b.add(n_coef_b)
        list_dT_hot = list_dT_hot.add(n_dT_hot)
        list_rah_hot = list_rah_hot.add(n_rah_hot)

    # Final aerodynamic resistance to heat transport [s m-1].
    i_rah_final = i_rah.rename('rah')

    # Final near surface temperature difference [K]
    i_dT_final = i_dT_int.rename('dT')

    # Final sensible heat flux [W m-2]
    i_H_final = i_H_int.expression(
        '(i_ro * n_Cp * i_dT_int) / i_rah',
        {'i_ro': i_ro, 'n_Cp': n_Cp, 'i_dT_int': i_dT_final,
         'i_rah': i_rah_final},
    )

    #LL - Need more analysis.
    '''
    # Evapotranspiration (advection)
    et_ad = i_dT_int.expression(
        '(8*(1+ws/100))/(log((0.67*3-n_zom)/i_zom)**2)',{
        'ws':ux,
        'n_zom':n_zom,
        'i_zom':i_zom
    }
    ).rename('et_ad')
    '''

    return i_H_final.rename('h_inst')


def daily_et(h_inst, g_inst, rn_inst, lst_dem, rad_24h):

    """
    Daily Evapotranspiration [mm day-1]

    Parameters
    ----------
    h_inst : ee.Image
        Instantaneous Sensible heat flux [W m-2].
    g_inst : ee.Image
        Instantaneous Soil heat flux [W m-2].
    rn_inst : ee.Image
        Instantaneous Net radiation [ W m-2].
    lst_dem : ee.Image
        Land surface temperature (aspect/slope correction) [K].
    rad_24h : ee.Image
        Daily Net Radiation [W m-2].

    Returns
    -------
    ee.Image

    References
    ----------

    .. [Bastiaanssen1998] Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A.,
        Holtslag, A.A.M., 1998. A remote sensing surface energy balance
        algorithm for land (SEBAL): 1. Formulation. J. Hydrol. 212–213, 198–212.

    """

    # Instantaneous Latent Heat flux [W m-2]
    le_inst = h_inst.expression(
        '(i_Rn - i_G - i_H)',
        {'i_Rn': rn_inst, 'i_G': g_inst, 'i_H': h_inst},
    )

    # Latent heat of vaporization or the heat
    # absorbed when a kilogram of water evaporates [J/kg].
    i_lambda = h_inst.expression(
        '(2.501 - 0.002361 * (Ts - 273.15))', {'Ts': lst_dem},
    )

    # Evaporative fraction
    i_FE = h_inst.expression(
        'i_lambda_ET / (i_Rn - i_G)',
        {'i_lambda_ET': le_inst, 'i_Rn': rn_inst, 'i_G': g_inst},
    )
    i_FE = i_FE.clamp(0, 1)

    i_ET24h_calc = i_FE.expression(
        '(0.0864 * i_FE * Rn24hobs) / i_lambda',
        {'i_FE': i_FE, 'i_lambda': i_lambda, 'Rn24hobs': rad_24h},
    )

    # Filtering et values
    i_ET24h_calc = i_ET24h_calc\
        .where(i_ET24h_calc.gte(-1).And(i_ET24h_calc.lt(0)), 0.01)
    i_ET24h_calc = i_ET24h_calc.updateMask(i_ET24h_calc.gte(0))
    i_ET24h_calc = i_ET24h_calc.updateMask(i_ET24h_calc.lte(9))

    return i_ET24h_calc.rename('et')


def et_fraction(time_start, et, et_reference_source, et_reference_band,
                et_reference_factor):

    """ET Fraction

    Parameters
    ----------
    time_start : ee.Image
        Instantaneous Sensible heat flux [W m-2].
    et : ee.Image
        Daily evapotranspiration (et) [mm day-1]
    et_reference_source : ee.ImageCollection, str
        ET reference collection
    et_reference_band : str
        ETr band name.
    et_reference_factor : ee.Number, int
        ETr factor.

    Returns
    -------
    ee.Image

    References
    ----------
    """
    date = ee.Date(time_start)
    start_date = ee.Date(utils.date_to_time_0utc(date))

    eto = ee.ImageCollection(et_reference_source)\
        .select(et_reference_band)\
        .filterDate(start_date, start_date.advance(1, 'day'))
    et_reference_img = ee.Image(eto.first())
    et_reference_img = et_reference_img.multiply(et_reference_factor)

    et_fraction = et.divide(et_reference_img).rename('et_fraction')

    return et_fraction
