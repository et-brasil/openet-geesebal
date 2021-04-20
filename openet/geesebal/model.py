import math

import ee

# from openet.geesebal import openet_landsat as landsat
from openet.geesebal import utils


def et(image,ndvi,ndwi,lst,albedo,lai,meteo_inst_source,meteo_daily_source,
       elev_product,ndvi_cold,ndvi_hot,lst_cold,lst_hot,emissivity,savi, time_start,geometry_image,crs,transform,coords
       ):

    #GET INFORMATIONS FROM IMAGE
    _date=ee.Date(time_start)
    _year=ee.Number(_date.get('year'))
    _month=ee.Number(_date.get('month'))
    _day=ee.Number(_date.get('day'))
    _hour=ee.Number(_date.get('hour'))
    _minuts = ee.Number(_date.get('minutes'))

    #ENDMEMBERS
    p_top_NDVI=ee.Number(ndvi_cold)
    p_coldest_Ts=ee.Number(lst_cold)
    p_lowest_NDVI=ee.Number(ndvi_hot)
    p_hottest_Ts=ee.Number(lst_hot)

    #METEOROLOGY PARAMETERS - GRIDMET AND NLDAS
    tmin,tmax,tair,ux,rh,rso_inst,rso24h=meteorology(image,time_start,meteo_inst_source,meteo_daily_source) #<----LOOK AFTER

    #SRTM DATA ELEVATION
    dem_product = ee.Image(elev_product)
    elev=dem_product.select('elevation')

    #sun elevation
    sun_elevation=ee.Number(image.get('SUN_ELEVATION'))

    #LL: iterative process or endmembers selections may return empty values
    #in this case, return an empty image instead of broken the code
    try:
        #LAND SURFACE TEMPERATURE CORRECTION
        lst_dem=lst_correction(
            image,time_start,ndwi,lst,elev,tair,rh,
            sun_elevation,_hour,_minuts,coords)

        #COLD PIXEL
        d_cold_pixel=cold_pixel(
            image,ndvi,ndwi,lst_dem,_year,p_top_NDVI,p_coldest_Ts,
            geometry_image,crs,transform,coords,time_start,elev,_hour,_minuts)

        #T COLD
        ts_cold_number = ee.Number(d_cold_pixel.get('temp'))

        #NET RADIATION
        rad_inst=radiation_inst(
            image,time_start,elev,lst,emissivity,albedo,
            tair,rh,rso_inst,sun_elevation)

        #SOIL HEAT FLUX (G)
        g_inst=soil_heat_flux(rad_inst,ndvi,albedo,lst_dem)

        #DAILY NET RADIATION (FAO56)
        rad_24h=radiation_24h(image,time_start,tmax,tmin,elev,rso24h)

        #HOT PIXEL
        d_hot_pixel=fexp_hot_pixel(
            image,time_start,ndvi,ndwi,lst_dem,rad_inst,g_inst,
            _year,p_lowest_NDVI,p_hottest_Ts,
            geometry_image,crs,transform,coords,elev,_hour,_minuts)

        #SENSIBLE HEAT FLUX (H) [W M-2]
        h_inst=sensible_heat_flux(
            image,savi,ux,rh,rad_24h,ts_cold_number,d_hot_pixel,
            lst_dem,lst,elev,geometry_image,crs,transform)

        #DAILY EVAPOTRANSPIRATION [MM DAY-1]
        et=daily_et(image,h_inst,g_inst,rad_inst,lst_dem,rad_24h)

    except:
        et=ee.Image.constant(0).updateMask(0) #return a masked image

    return et


def meteorology(image_region,time_start,meteo_inst_source,meteo_daily_source):

    meteorology_inst_collection=ee.ImageCollection(meteo_inst_source)

    #linear interpolation
    time_start=ee.Number(time_start)
    previous_time=time_start.subtract(2*60*60*1000)
    next_time=time_start.add(2*60*60*1000)

    meteorology_daily = ee.ImageCollection(meteo_daily_source)\
        .filter(ee.Filter.date(ee.Date(time_start).advance(-1,'day'),ee.Date(time_start)))

    previous_image=meteorology_inst_collection\
        .filter(ee.Filter.date(previous_time,time_start))\
        .limit(1, 'system:time_start', False).first()

    next_image=meteorology_inst_collection\
        .filter(ee.Filter.date(time_start,next_time))\
        .limit(1, 'system:time_start', True).first()

    image_previous_time= ee.Number(previous_image.get('system:time_start'))

    image_next_time=ee.Number(next_image.get('system:time_start'))

    delta_time=(time_start.subtract(image_previous_time)).divide(image_next_time.subtract(image_previous_time))

    #daily variables

    #incoming shorwave down (W m2)
    Swdown24h = (meteorology_daily.select('srad').first()).rename('short_wave_down')

    tmin=(meteorology_daily.select('tmmn').first()).rename('tmin')

    tmax=(meteorology_daily.select('tmmx').first()).rename('tmax')

    #instantaneous
    rso_inst=next_image.select('shortwave_radiation')\
        .subtract(previous_image.select('shortwave_radiation'))\
        .multiply(delta_time).add(previous_image.select('shortwave_radiation')).rename('rso_inst')

    #SPECIFIC HUMIDITY [KG KG-1]
    q_med =next_image.select('specific_humidity')\
        .subtract(previous_image.select('specific_humidity'))\
        .multiply(delta_time).add(previous_image.select('specific_humidity'))

    #AIR TEMPERATURE [K]
    tair_C = next_image.select('temperature')\
        .subtract(previous_image.select('temperature'))\
        .multiply(delta_time).add(previous_image.select('temperature')).rename('tair')

    #WIND SPEED [M S-1]
    wind_u = next_image.select('wind_u')\
        .subtract(previous_image.select('wind_u'))\
        .multiply(delta_time).add(previous_image.select('wind_u')).rename('wind_u')

    wind_v = next_image.select('wind_v')\
        .subtract(previous_image.select('wind_v'))\
        .multiply(delta_time).add(previous_image.select('wind_v')).rename('wind_v')

    wind_med=wind_u.expression(
        'sqrt((ux_u)**2 + (ux_v) ** 2)',{
            'ux_u':wind_u.select('wind_u'),
            'ux_v':wind_v.select('wind_v')
            }).rename('ux')

    wind_med=wind_med.expression(
        'ux*(4.87)/log(67.8*z-5.42)',{
            'ux': wind_med,
            'z':ee.Number(10)})

    #PRESSURE [PA] CONVERTED TO KPA
    p_med= next_image.select('pressure')\
        .subtract(previous_image.select('pressure'))\
        .multiply(delta_time).add(previous_image.select('pressure')).divide(ee.Number(1000))

    #ACTUAL VAPOR PRESSURE [KPA]
    ea=p_med.expression('(1/0.622)*Q*P',{
        'Q': q_med,
        'P':p_med}).rename('ea')

    #SATURATED VAPOR PRESSURE [KPA]
    esat=tair_C.expression('0.6108*(exp((17.27*T_air)/(T_air+237.3)))',{
        'T_air':tair_C}).rename('esat')

    #RELATIVE HUMIDITY (%)
    RH=ea.divide(esat).multiply(100).rename('RH')

    #Resample

    tmin=tmin.subtract(273.15).resample('bilinear')
    tmax=tmax.subtract(273.15).resample('bilinear')
    tair_C=tair_C.resample('bilinear')
    wind_med=wind_med.resample('bilinear')
    rh=RH.resample('bilinear')
    Swdown24h=Swdown24h.resample('bilinear')

    return [tmin,tmax,tair_C,wind_med,rh,rso_inst,Swdown24h]

def tao_sw(landsat_image,time_start,dem,tair,rh,sun_elevation):

    """Correct declivity and aspect effects from Land Surface Temperature"""

    #ATMOSPHERIC PRESSURE [KPA]
    #SHUTTLEWORTH (2012)
    pres = landsat_image.expression(
        '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
            'Z' : dem}).rename('p_atm')

    #SATURATION VAPOR PRESSURE (es) [KPA]
    es = landsat_image.expression(
        '0.6108 *(exp( (17.27 * tair) / (tair + 237.3)))', {
            'tair': tair}).rename('es')

    #ACTUAL VAPOR PRESSURE (ea) [KPA]
    ea = es.multiply(rh).divide(100).rename('ea')

    #WATER IN THE ATMOSPHERE [mm]
    #Garrison and Adler (1990)
    w = landsat_image.expression(
        '(0.14 * EA * PATM) + 2.1', {
            'PATM' : pres,
            'EA' : ea}).rename('w_atm')
    #print(sun_elevation.getInfo())
    #SOLAR ZENITH ANGLE OVER A HORZONTAL SURFACE
    solar_zenith = ee.Number(90).subtract(sun_elevation)

    #DEGREE TO RADIAN
    degree2radian = 0.01745

    solar_zenith_radians = solar_zenith.multiply(degree2radian)
    cos_theta = solar_zenith_radians.cos()

    #BROAD-BAND ATMOSPHERIC TRANSMISSIVITY (tao_sw)
    #ASCE-EWRI (2005)
    tao_sw_img = landsat_image.expression(
        '0.35 + 0.627 * exp(((-0.00146 * P)/(Kt * cos_theta)) - (0.075 * (W / cos_theta)**0.4))', {
            'P' : pres,
            'W': w,
            'Kt' : ee.Number(1),
            'cos_theta' : cos_theta}).rename('tao_sw')

    return tao_sw_img

def cos_terrain(landsat_image,time_start,dem,hour,minuts,coords):

    dateStr = ee.Date(time_start)

    #DAY OF YEAR
    doy = dateStr.getRelative('day', 'year').add(1)

    #DEGREE TO RADIAN
    degree2radian = 0.01745

    #COS ZENITH ANGLE SUN ELEVATION #ALLEN ET AL. (2006)
    slope_aspect = ee.Terrain.products(dem)
    B=(ee.Number(360).divide(ee.Number(365))).multiply(doy.subtract(ee.Number(81)))
    delta = ee.Image(ee.Number(23.45).multiply(degree2radian).sin().asin().multiply(B.multiply(degree2radian).sin()))
    s =slope_aspect.select('slope').multiply(degree2radian)
    gamma = (slope_aspect.select('aspect').subtract(180)).multiply(degree2radian)
    phi = coords.select('latitude').multiply(degree2radian)

    #CONSTANTS ALLEN ET AL. (2006)
    a = ee.Image((delta.sin().multiply(phi.cos()).multiply(s.sin()).multiply(gamma.cos())).subtract(delta.sin().multiply(phi.sin().multiply(s.cos()))))
    b = (delta.cos().multiply(phi.cos()).multiply(s.cos())).add(delta.cos().multiply(phi.sin().multiply(s.sin()).multiply(gamma.cos())))
    c= (delta.cos().multiply(s.sin()).multiply(gamma.sin()))

    #GET IMAGE CENTROID
    longitude_center=coords.select('longitude')
    #DELTA GTM
    delta_gtm =longitude_center.divide(ee.Image(15)).int()

    min_to_hour=ee.Image(minuts).divide(60)

    #LOCAL HOUR TIME
    lht = ee.Image(hour).add(delta_gtm).add(min_to_hour)

    hour_a = (lht.subtract(12)).multiply(15)

    w = hour_a.multiply(degree2radian)

    cos_zn =landsat_image.expression(
        '-a +b*w_cos +c*w_sin',{
            'a': a,
            'b': b,
            'c': c,
            'w_cos': w.cos(),
            'w_sin': w.sin()})

    return cos_zn

def lst_correction(landsat_image,time_start,ndwi,lst,dem,tair,rh,sun_elevation,hour,minuts,coords):

    """Correct declivity and aspect effects from Land Surface Temperature"""

    #SOLAR CONSTANT [W M-2]
    gsc = ee.Number(1367)
    dateStr = ee.Date(time_start)

    #DAY OF YEAR
    doy = dateStr.getRelative('day', 'year').add(1)
    pi=ee.Number(math.pi)

    #INVERSE RELATIVE  DISTANCE EARTH-SUN
    d1 =  ee.Number(2).multiply(pi).divide(ee.Number(365))
    d2 = d1.multiply(doy)
    d3 = d2.cos()
    dr = ee.Number(1).add(ee.Number(0.033).multiply(d3))

    #ATMOSPHERIC PRESSURE [KPA]
    #SHUTTLEWORTH (2012)
    pres = landsat_image.expression(
        '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
            'Z' : dem}).rename('p_atm')

    #SOLAR ZENITH ANGLE OVER A HORZONTAL SURFACE
    solar_zenith = ee.Number(90).subtract(sun_elevation)

    #DEGREE TO RADIAN
    degree2radian = 0.01745

    solar_zenith_radians = solar_zenith.multiply(degree2radian)
    cos_theta = solar_zenith_radians.cos()

    #BROAD-BAND ATMOSPHERIC TRANSMISSIVITY (tao_sw)
    #ASCE-EWRI (2005)
    tao_sw_img =tao_sw(landsat_image,time_start,dem,tair,rh,sun_elevation)

    #AIR DENSITY [KG M-3]
    air_dens = landsat_image.expression(
        '(1000* Pair)/(1.01*LST*287)',{
            'Pair': pres,
            'LST': lst})

    #TEMPERATURE LAPSE RATE (0.0065)
    Temp_lapse_rate= ee.Number(0.0065)

    #LAND SURFACE TEMPERATURE CORRECTION DEM [K]
    Temp_corr= lst.add(dem.select('elevation').multiply(Temp_lapse_rate))

    cos_zn=cos_terrain(landsat_image,time_start,dem,hour,minuts,coords)

    #LAND SURFACE TEMPERATURE WITH ASPECT/SLOPE CORRECTION [K]
    lst_dem=landsat_image.expression(
        '(Temp_corr + (Gsc * dr * Transm_corr * cos_zn -Gsc * dr * Transm_corr * cos_zenith_flat) / (air_dens * 1004 * 0.050))',{
            'Temp_corr':Temp_corr,
            'Gsc':gsc,
            'dr':dr,
            'Transm_corr':tao_sw_img,
            'cos_zenith_flat':cos_theta,
            'cos_zn':cos_zn,
            'air_dens':air_dens}).rename('lst_dem')

    return lst_dem

def lc_mask(landsat_image,year,geometry_image):

    """Filtering pre-candidates pixels using a Land cover mask"""

    #CONDITIONS
    year_condition=ee.Algorithms.If(ee.Number(year).lte(2007),2008,year)

    start = ee.Date.fromYMD(year_condition, 1, 1)
    end = ee.Date.fromYMD(year_condition, 12, 31)

    #select classification corresponding to the year of the image
    lc = ee.ImageCollection('USDA/NASS/CDL').select('cropland')\
        .filter(ee.Filter.date(start, end)).first()

    #filter cropland classes 1
    crop1 = lc.updateMask(lc.lte(62))
    crop1 = crop1.where(crop1, 1).unmask(0)

    #filter cropland classes 2
    crop2 = lc.updateMask(lc.gte(196))
    crop2 = crop2.where(crop2, 1).unmask(0)

    #land cover mask - total croplands
    lc_mask = crop1.add(crop2)

    lc_mask = lc_mask.updateMask(lc_mask.eq(1))

    count_land_cover_pixels = lc_mask.rename('land_cover_pixels').reduceRegion(
        reducer=ee.Reducer.count(),
        scale= 30,
        geometry= geometry_image,
        maxPixels= 10e14)
    n_count_lc=ee.Number(count_land_cover_pixels.get('land_cover_pixels'))

    mask=ee.Algorithms.If(n_count_lc.gte(3000),lc_mask,landsat_image.select(0).updateMask(1))

    return ee.Image(mask)

def cold_pixel(landsat_image,ndvi,ndwi,lst_dem,year,ndvi_cold,lst_cold,
               geometry_image,_crs,_transform,coords,time_start,dem,hour,minuts):

    """Simplified CIMEC method to select the cold pixel"""

    #pre-filter
    pos_ndvi = ndvi.updateMask(ndvi.gt(0)).rename('post_ndvi')
    ndvi_neg =  pos_ndvi.multiply(-1).rename('ndvi_neg')

    lst_neg = lst_dem.multiply(-1).rename('lst_neg').rename('lst_neg')
    lst_nw = lst_dem.updateMask(ndwi.lte(0)).rename('lst_nw')

    #land cover mask
    land_cover_mask=lc_mask(landsat_image,year,geometry_image)

    image=pos_ndvi.addBands([ndvi,ndvi_neg,pos_ndvi,lst_neg,lst_nw,coords]) #.reproject(crs=_crs, crsTransform=_transform)

    d_perc_top_NDVI=image.select('ndvi_neg').updateMask(land_cover_mask).reduceRegion(
        reducer=ee.Reducer.percentile([ndvi_cold]),
        geometry= geometry_image,
        #crs=_crs,
        #crsTransform=_transform,
        scale= 30,
        maxPixels=1e9)
    n_perc_top_NDVI= ee.Number(d_perc_top_NDVI.get('ndvi_neg'))

    i_top_NDVI=image.updateMask(land_cover_mask).updateMask(image.select('ndvi_neg').lte(n_perc_top_NDVI))

    d_perc_low_LST = i_top_NDVI.updateMask(land_cover_mask).select('lst_nw').reduceRegion(
        reducer= ee.Reducer.percentile([lst_cold]),
        geometry=geometry_image,
        scale= 30,
        maxPixels=1e9)
    n_perc_low_LST = ee.Number(d_perc_low_LST.get('lst_nw'))

    i_cold_lst = i_top_NDVI.updateMask(land_cover_mask).updateMask(i_top_NDVI.select('lst_nw').lte(n_perc_low_LST))
    #Filtes
    c_lst_cold20 =  i_cold_lst.updateMask(image.select('lst_nw').gte(200))
    c_lst_cold20_int=c_lst_cold20.select('lst_nw').int().rename('int')
    c_lst_cold20=c_lst_cold20.addBands(c_lst_cold20_int)

    med_lst_cold20 = c_lst_cold20.select('lst_nw').reduceRegion(
        reducer=  ee.Reducer.median(),
        geometry= geometry_image,
        scale= 30,
        maxPixels= 1e9)
    n_med_lst_cold20 = ee.Number(med_lst_cold20.get('lst_nw'))

    sum_final_cold_pix = c_lst_cold20.select('int').reduceRegion(
        reducer=  ee.Reducer.sum(),
        geometry= geometry_image,
        scale= 30,
        maxPixels=1e9)
    n_sum_final_cold_pix = ee.Number(sum_final_cold_pix.get('int'))

    def function_def_pixel(f):
        return f.setGeometry(ee.Geometry.Point([f.get('longitude'), f.get('latitude')]))
    #Define Cold Pixel (random)
    fc_cold_pix = c_lst_cold20.stratifiedSample(1, "int", geometry_image, 30).map(function_def_pixel)
    n_Ts_cold = ee.Number(fc_cold_pix.aggregate_first('lst_nw'))
    n_long_cold = ee.Number(fc_cold_pix.aggregate_first('longitude'))
    n_lat_cold = ee.Number(fc_cold_pix.aggregate_first('latitude'))
    n_ndvi_cold = ee.Number(fc_cold_pix.aggregate_first('ndvi'))

    #Dictionary
    d_cold_pixel = ee.Dictionary({
        'temp': n_Ts_cold,
        'ndvi': n_ndvi_cold,
        'x':n_long_cold,
        'y': n_lat_cold,
        'sum': n_sum_final_cold_pix
    })
    return d_cold_pixel

def radiation_inst(landsat_image,time_start,dem,lst,emissivity,albedo,
                   tair,rh,swdown_inst,sun_elevation):

    """Instantaneous Net Radiation (W m-2)"""

    rad_long_up=landsat_image.expression(
        'emi * stefBol * (LST ** 4)', {
            'emi' : emissivity,
            'stefBol': ee.Image(5.67e-8),
            'LST': lst})

    tao_sw_img=tao_sw(landsat_image,time_start,dem,tair,rh,sun_elevation)

    log_taosw = tao_sw_img.log()

    rad_long_down = landsat_image.expression(
        '(0.85 * (- log_taosw) ** 0.09) * stefBol * (n_Ts_cold ** 4)', {
            'log_taosw' : log_taosw,
            'stefBol': ee.Number(5.67e-8),
            'n_Ts_cold' :tair.add(273.15)}).rename('Rl_down')

    rn_inst = landsat_image.expression(
        '((1-alfa) * Rs_down) + Rl_down - Rl_up - ((1 - e_0) * Rl_down) ', {
            'alfa' : albedo,
            'Rs_down': swdown_inst,
            'Rl_down' : rad_long_down,
            'Rl_up' : rad_long_up,
            'e_0' : emissivity}).rename('rn_inst')

    return rn_inst

def soil_heat_flux(rn,ndvi,albedo,lst_dem):

    """Instantaneous Soil Heat Flux (W m-2)"""

    g = rn.expression(
        'rn * (lst - 273.15) * ( 0.0038 + (0.0074 * albedo)) *  (1 - 0.98 * (ndvi ** 4)) ', {
            'rn' : rn,
            'ndvi': ndvi,
            'albedo': albedo,
            'lst': lst_dem}).rename('g_inst')

    return g

def radiation_24h(image,time_start,tmax,tmin,elev,rso24h):

    """Daily Net radiation (W m-2) - FAO56 """

    #CONVERT TO MJ m^2
    rs=rso24h.multiply(0.0864).rename('Rs')
    gsc = 0.0820 #MJ m2

    dateStr = ee.Date(time_start)

    doy = dateStr.getRelative('day', 'year').add(1)
    pi=ee.Number(math.pi)
    #Relative Earth–Sun distance
    dr = image.expression(
        '1 + (0.033 * cos((2 * pi/365) * doy) )',
        {'doy': doy, 'pi': pi})

    sd = image.expression(
        "0.40928 * sin( (((2 * pi) / 365) * doy) - 1.39 )",
        {'doy': doy, 'pi': pi})

    degree2radian = 0.01745

    lat = image.pixelLonLat().select(['latitude']).multiply(degree2radian).rename("latitude")

    ws = image.expression(
        'acos(-tan(Lat) * tan(Sd))',
        {'Lat':lat, 'Sd':sd})

    rad_a = image.expression(
        'Ws * sin(Lat) * sin(Sd) + cos(Lat) * cos(Sd) * sin(Ws)',
        {'Ws':ws, 'Lat':lat, 'Sd':sd}).rename('rad_a')

    ra = image.expression(
        '((24 * 60) / pi) * Gsc * Dr * rad_a',
        {'pi':pi,'Gsc': gsc, 'Dr':dr, 'rad_a':rad_a}).rename('ra')

    rso = image.expression(
        '(0.75 + 2E-5 * z) * Ra',
        {'z':elev, 'Ra':ra}).rename('Rso')

    rns = image.expression(
        '(1 -albedo) * Rs',
        {'Rs': rs, 'albedo': ee.Number(0.23)}).rename('rns')

    ea = image.expression(
        '0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
        'T_air': tmin}).rename('ea_rad')

    rnl = image.expression(
        'o * ((Tmax ** 4 + Tmin ** 4)/2) * (0.34 - 0.14 * sqrt(ea)) * (1.35*(Rs/Rso)- 0.35)',
        {'o': ee.Number(4.901E-9), 'Tmax':tmax.add(273.15), 'Tmin':tmin.add(273.15), 'ea':ea, 'Rs': rs, 'Rso': rso} ).rename('rnl')

    Rn = image.expression(
        'Rns - Rnl',
        {'Rns': rns, 'Rnl': rnl}).rename('rad_24h')

    Rn=Rn.multiply(ee.Number(11.6))

    return Rn

def fexp_hot_pixel(landsat_image,time_start,ndvi,ndwi,lst_dem,rn,g,year,
                   ndvi_hot,lst_hot,geometry_image,crs,transform,coords,
                   dem,hour,minuts):

    """Simplified CIMEC method to select the hot pixel"""

    #pre-filter
    pos_ndvi = ndvi.updateMask(ndvi.gt(0)).rename('post_ndvi')
    ndvi_neg =  pos_ndvi.multiply(-1).rename('ndvi_neg')

    lst_neg = lst_dem.multiply(-1).rename('lst_neg')
    lst_nw = lst_dem.updateMask(ndwi.lte(0)).rename('lst_nw')

    #land cover mask
    land_cover_mask=lc_mask(landsat_image,year,geometry_image)

    image=pos_ndvi.addBands([ndvi,ndvi_neg,rn,g,pos_ndvi,lst_neg,lst_nw,coords])

    d_perc_down_ndvi = image.select('post_ndvi').updateMask(land_cover_mask).reduceRegion(
        reducer= ee.Reducer.percentile([ndvi_hot]),
        geometry= geometry_image,
        scale= 30,
        maxPixels=1e9)
    n_perc_low_NDVI = ee.Number(d_perc_down_ndvi.get('post_ndvi'))

    i_low_NDVI = image.updateMask(land_cover_mask).updateMask(image.select('post_ndvi').lte(n_perc_low_NDVI))

    d_perc_top_lst = i_low_NDVI.updateMask(land_cover_mask).select('lst_neg').reduceRegion(
        reducer=ee.Reducer.percentile([lst_hot]),
        geometry= geometry_image,
        scale= 30,
        maxPixels= 1e9)
    n_perc_top_lst = ee.Number(d_perc_top_lst.get('lst_neg'))

    i_top_LST = i_low_NDVI.updateMask(land_cover_mask).updateMask(i_low_NDVI.select('lst_neg').lte(n_perc_top_lst))

    c_lst_hot_int=i_top_LST.select('lst_nw').int().rename('int')
    c_lst_hotpix=i_top_LST.addBands(c_lst_hot_int)

    sum_final_hot_pix = c_lst_hotpix.select('int').reduceRegion(
        reducer=  ee.Reducer.sum(),
        geometry= geometry_image,
        scale= 30,
        maxPixels= 1e9,)
    n_sum_final_hot_pix = ee.Number(sum_final_hot_pix.get('int'))

    #PRECIPITATION CORRECTION
    gridmet=ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")

    date_init=ee.Date(time_start)

    date60_previous=date_init.advance(-60,'days')
    precipt=gridmet.select('pr').filterDate(date60_previous,date_init)

    etr=gridmet.select('etr').filterDate(date60_previous,date_init)

    etr_60mm=etr.sum()

    precipt_60mm=precipt.sum()

    ratio=precipt_60mm.divide(etr_60mm)

    Tfac=landsat_image.expression(
        '2.6 - 13*ratio',{
            'ratio': ratio })

    Tfac=ee.Image(Tfac.where(ratio.gt(0.2),0)).rename('Tfac')

    c_lst_hotpix=c_lst_hotpix.addBands(Tfac)

    def function_def_pixel(f):
        return f.setGeometry(ee.Geometry.Point([f.get('longitude'), f.get('latitude')]))
    #Define Hot Pixel (random)
    fc_hot_pix = c_lst_hotpix.stratifiedSample(1, "int", geometry_image, 30).map(function_def_pixel)
    n_Ts_hot = ee.Number(fc_hot_pix.aggregate_first('lst_nw')) # //transforma de objeto para numero
    n_long_hot = ee.Number(fc_hot_pix.aggregate_first('longitude'))
    n_lat_hot = ee.Number(fc_hot_pix.aggregate_first('latitude'))
    n_ndvi_hot = ee.Number(fc_hot_pix.aggregate_first('ndvi'))
    n_Rn_hot = ee.Number(fc_hot_pix.aggregate_first('rn_inst'))
    n_G_hot = ee.Number(fc_hot_pix.aggregate_first('g_inst'))
    n_Tfac= ee.Number(fc_hot_pix.aggregate_first('Tfac'))

    d_hot_pixel = ee.Dictionary({
        'temp': ee.Number(n_Ts_hot).subtract(ee.Number(n_Tfac)),
        'x': n_long_hot,
        'y': n_lat_hot,
        'rn': n_Rn_hot,
        'g': n_G_hot,
        'ndvi': n_ndvi_hot,
        'sum': n_sum_final_hot_pix,
    })

    return d_hot_pixel

def sensible_heat_flux(landsat_image,savi,ux,rh,rad_24h,
                       ts_cold_number,d_hot_pixel,lst_dem,lst,dem,
                       geometry_image,crs,transform):

    """Instantaneous Sensible Heat Flux (W m-2)"""

    #parameters
    n_veg_hight = ee.Number(2)
    n_zx = ee.Number(2)
    n_hight = ee.Number(200)
    n_Cp = ee.Number(1004)#Air specific heat (J/kg/K)
    n_K = ee.Number(0.41)  #k is von Karman’s constant
    slope_aspect = ee.Terrain.products(dem)

    n_Ts_cold = ee.Number(ts_cold_number)
    n_Ts_hot = ee.Number(d_hot_pixel.get('temp'))
    n_G_hot = ee.Number(d_hot_pixel.get('g'))
    n_Rn_hot = ee.Number(d_hot_pixel.get('rn'))
    n_long_hot = ee.Number(d_hot_pixel.get('x'))
    n_lat_hot = ee.Number(d_hot_pixel.get('y'))
    p_hot_pix =  ee.Geometry.Point([n_long_hot, n_lat_hot])

    #Momentum roughness length (zom) at the weather station:
    n_zom = n_veg_hight.multiply(0.12) # //% at the weather station

    #Friction velocity at the weather station
    i_ufric_ws = landsat_image.expression(
        '(n_K * ux)/ log(n_zx /n_zom)', {
            'n_K': n_K,
            'n_zx': n_zx,
            'n_zom': n_zom,
            'ux': ux }).rename('ufric_ws')

    #Wind speed at blending height at the weather station
    i_u200 = landsat_image.expression(
        'i_ufric_ws *  (log(n_hight/n_zom)/n_K)', {
            'i_ufric_ws' : i_ufric_ws,
            'n_hight' : n_hight,
            'n_zom' : n_zom,
            'n_K' : n_K}).rename('i_u200')

    #Momentum roughness length (zom) for each pixel:
    i_zom = landsat_image.expression(
        'exp((5.62 * (SAVI))-5.809)', {
            'SAVI' : savi,}).rename('zom')

    #correction
    i_zom=i_zom.expression(
        'zom*(1+(slope-5)/20)',{
            'zom': i_zom,
            'slope':slope_aspect.select('slope')
        })

    #Friction velocity for each pixel:
    i_ufric = landsat_image.expression(
        '(n_K *u200) /(log(hight/i_zom))', {
            'u200' : i_u200,
            'hight': n_hight,
            'i_zom':n_zom,
            'n_K': n_K }).rename('u_fr')

    # Aerodynamic resistance to heat transport (rah)
    z1= ee.Number(0.1)
    z2= ee.Number(2)  #vegetation
    i_rah = i_ufric.expression(
        '(log(z2/z1))/(i_ufric*0.41)', {
            'z2' : z2,
            'z1': z1,
            'i_ufric':i_ufric }).rename('rah')

    i_rah_first = i_rah.rename('rah_first')

    n_ro_hot= (ee.Number(-0.0046).multiply(n_Ts_hot)).add(ee.Number(2.5538))

    #Iterative Process
    #Sensible heat flux at the hot pixel
    n_H_hot = ee.Number(n_Rn_hot).subtract(ee.Number(n_G_hot))
    n= ee.Number(1)
    n_dif= ee.Number(1)
    n_dif_min = ee.Number(0.1)
    list_dif = ee.List([])
    list_dT_hot = ee.List([])
    list_rah_hot = ee.List([])
    list_coef_a = ee.List([])
    list_coef_b = ee.List([])
    for n in range(20):
        d_rah_hot = i_rah.reduceRegion(
            reducer= ee.Reducer.first(),
            geometry= p_hot_pix,
            scale= 30,
            maxPixels=9000000000)

        # LL : To avoid 'Max (NaN) cannot be less than min (NaN)' erros in cases which iterative process
        # not converge
        n_rah_hot =   ee.Number(d_rah_hot.get('rah')).multiply(100).short().divide(100)

        #Near surface temperature difference in hot pixel (dT = Tz1 – Tz2)
        # dThot= Hhot*rah/(ρCp)
        n_dT_hot = (n_H_hot.multiply(n_rah_hot)).divide(n_ro_hot.multiply(n_Cp))

        # Near surface temperature difference in cold pixel (dT = Tz1 – Tz2)
        n_dT_cold = ee.Number(0)
        # dT =  aTs + b
        #Angular coefficient
        n_coef_a = (n_dT_cold.subtract(n_dT_hot)).divide(n_Ts_cold.subtract(n_Ts_hot))

        #Linear coefficient
        n_coef_b = n_dT_hot.subtract(n_coef_a.multiply(n_Ts_hot))

        #dT for each pixel
        i_dT_int = landsat_image.expression(
            '(n_coef_a * i_lst_med_dem) + n_coef_b', {
                'n_coef_a' : n_coef_a,
                'n_coef_b': n_coef_b,
                'i_lst_med_dem':lst_dem }).rename('dT')

        #Air temperature (Ta) for each pixel (Ta = Ts-dT)
        i_Ta = landsat_image.expression(
            'i_lst_med - i_dT_int', {
                'i_lst_med' : lst,
                'i_dT_int': i_dT_int})

        #ro (ρ) - air density (kg/m3)
        i_ro = i_Ta.expression(
            '(-0.0046 * i_Ta) + 2.5538', {
                'i_Ta' : i_Ta}).rename('ro')  #// ro=-0.0046.*Ta+2.5538

        #Sensible heat flux (H) for each pixel - iteration
        i_H_int = i_dT_int.expression(
            '(i_ro*n_Cp*i_dT_int)/i_rah', {
                'i_ro' : i_ro,
                'n_Cp': n_Cp,
                'i_dT_int':i_dT_int,
                'i_rah':i_rah }).rename('H')

        #Monin-Obukhov length (L) - iteration
        i_L_int = i_dT_int.expression(
            '-(i_ro*n_Cp*(i_ufric**3)*i_lst_med)/(0.41*9.81*i_H_int)',{
                'i_ro' : i_ro,
                'n_Cp': n_Cp,
                'i_ufric':i_ufric,
                'i_lst_med':lst,
                'i_H_int':i_H_int }).rename('L')

        #Stability corrections for momentum and heat transport
        img = landsat_image

        #stability corrections for stable conditions
        i_psim_200 = img.expression(
            '-5*(hight/i_L_int)', {'hight' : ee.Number(200),'i_L_int': i_L_int}).rename('psim_200')
        i_psih_2 = img.expression(
            '-5*(hight/i_L_int)',{'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('psih_2')
        i_psih_01 = img.expression(
            '-5*(hight/i_L_int)',{'hight' : ee.Number(0.1),'i_L_int': i_L_int}).rename('psih_01')

        #x for different height
        i_x200 = i_L_int.expression(
            '(1-(16*(hight/i_L_int)))**0.25',
            {'hight' : ee.Number(200),'i_L_int': i_L_int}).rename('i_x200')
        i_x2 = i_L_int.expression(
            '(1-(16*(hight/i_L_int)))**0.25',
            {'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('i_x2')
        i_x01 = i_L_int.expression(
            '(1-(16*(hight/i_L_int)))**0.25',
            {'hight' : ee.Number(0.1),'i_L_int': i_L_int})

        #stability corrections for unstable conditions
        i_psimu_200 = i_x200.expression(
            '2*log((1+i_x200)/2)+log((1+i_x200**2)/2)-2*atan(i_x200)+0.5*pi',
            {'i_x200' : i_x200,'pi': ee.Number(3.14159265)})
        i_psihu_2 = i_x2.expression(
            '2*log((1+i_x2**2)/2)',
            {'i_x2' : i_x2})
        i_psihu_01 = i_x01.expression(
            '2*log((1+i_x01**2)/2)',
            {'i_x01' : i_x01})

        i_psim_200 = i_psim_200.where(i_L_int.lt(0), i_psimu_200)
        i_psih_2 = i_psih_2.where(i_L_int.lt(0), i_psihu_2)
        i_psih_01 = i_psih_01.where(i_L_int.lt(0), i_psihu_01)
        i_psim_200 = i_psim_200.where(i_L_int.eq(0), 0)
        i_psih_2 = i_psih_2.where(i_L_int.eq(0), 0)
        i_psih_01 = i_psih_01.where(i_L_int.eq(0), 0)

        if n==1:
            i_psim_200_exp = i_psim_200
            i_psih_2_exp = i_psih_2
            i_psih_01_exp = i_psih_01
            i_L_int_exp = i_L_int
            i_H_int_exp = i_H_int
            i_dT_int_exp = i_dT_int
            i_rah_exp = i_rah

        #Corrected value for the friction velocity (u_asterisk)
        #u_asterisk=(u200.*0.41)./(log(hight./zom_pixel)-psi_m200)
        i_ufric = i_ufric.expression(
            '(u200*0.41)/(log(hight/i_zom)-i_psim_200)',
            {'u200' : i_u200,'hight': n_hight, 'i_zom':i_zom,'i_psim_200': i_psim_200}).rename('ufric_star')

        #Corrected value for the aerodinamic resistance to the heat transport (rah)
        #rah=(log(z2/z1)-psi_h2+psi_h01)./(u_asterisk*0.41)
        i_rah = i_rah.expression(
            '(log(z2/z1)-psi_h2+psi_h01)/(i_ufric*0.41)',
            {'z2' : z2,'z1': z1, 'i_ufric':i_ufric, 'psi_h2':i_psih_2, 'psi_h01':i_psih_01}).rename('rah')

        if n==1:
            n_dT_hot_old = n_dT_hot
            n_rah_hot_old = n_rah_hot
            n_dif = ee.Number(1)

        if n > 1:
            n_dT_hot_abs = n_dT_hot.abs()
            n_dT_hot_old_abs = n_dT_hot_old.abs()
            n_rah_hot_abs = n_rah_hot.abs()
            n_rah_hot_old_abs = n_rah_hot_old.abs()
            n_dif=(n_dT_hot_abs.subtract(n_dT_hot_old_abs).add(n_rah_hot_abs).subtract(n_rah_hot_old_abs)).abs()
            n_dT_hot_old = n_dT_hot
            n_rah_hot_old = n_rah_hot
            #insert each iteration value into a list

        list_dif = list_dif.add(n_dif)
        list_coef_a = list_coef_a.add(n_coef_a)
        list_coef_b = list_coef_b.add(n_coef_b)
        list_dT_hot = list_dT_hot.add(n_dT_hot)
        list_rah_hot = list_rah_hot.add(n_rah_hot)

    i_rah_final = i_rah.rename('rah')
    i_dT_final = i_dT_int.rename('dT')
    i_H_final = i_H_int.expression(
        '(i_ro*n_Cp*i_dT_int)/i_rah',
        {'i_ro' : i_ro,'n_Cp': n_Cp, 'i_dT_int':i_dT_final, 'i_rah':i_rah_final }).rename('H')

    return i_H_final

def daily_et(landsat_image,h_inst,g_inst,rn_inst,lst_dem,rad_24h):

    """Daily Evapotranspiration (mm day-1)"""

    le_inst = h_inst.expression(
        '(i_Rn-i_G-i_H_fim)', {
            'i_Rn' : rn_inst,
            'i_G': g_inst,
            'i_H_fim':h_inst }).rename('le_inst')

    #Latent heat of vaporization or the heat
    #absorbed when a kilogram of water evaporates - lambda (J/kg).
    i_lambda = h_inst.expression(
        '(2.501-0.002361*(Ts-273.15))', {'Ts' : lst_dem })

    #Instantaneous value of ET (mm/H)
    i_ET_inst = h_inst.expression(
        '0.0036 * (i_lambda_ET/i_lambda)', {
            'i_lambda_ET' : le_inst,
            'i_lambda' : i_lambda  }).rename('et_inst')

    #Evaporative fraction
    i_FE = h_inst.expression(
        'i_lambda_ET/(i_Rn-i_G)', {
            'i_lambda_ET' : le_inst,
            'i_Rn' : rn_inst,
            'i_G' : g_inst }).rename('FE')
    i_FE=i_FE.clamp(0,1)

    i_ET24h_calc = i_FE.expression(
        '(0.0864 *i_FE * Rn24hobs)/(i_lambda)', {
            'i_FE' : i_FE,
            'i_lambda' : i_lambda,
            'Rn24hobs' : rad_24h
        }).rename('et')

    #filtering et values
    i_ET24h_calc=i_ET24h_calc.where(i_ET24h_calc.gte(-1).And(i_ET24h_calc.lt(0)),0.01)
    i_ET24h_calc=i_ET24h_calc.updateMask(i_ET24h_calc.gte(0))
    i_ET24h_calc=i_ET24h_calc.updateMask(i_ET24h_calc.lte(9))

    return i_ET24h_calc

def et_fraction(landsat_image,time_start,et,
                et_reference_source,et_reference_band,et_reference_factor):

    """ET Fraction"""

    _date=ee.Date(time_start)
    _start_date = ee.Date(utils.date_to_time_0utc(_date))
    _end_date = _start_date.advance(1, 'day')

    eto=(ee.ImageCollection(et_reference_source).select(et_reference_band).filterDate(_start_date,_end_date))
    et_reference_img = ee.Image(eto.first())
    et_reference_img=et_reference_img.multiply(et_reference_factor)

    et_fraction=et.divide(et_reference_img).rename('et_fraction')

    return et_fraction
