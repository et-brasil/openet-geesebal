import ee


def ndvi(landsat_image):
    """Normalized difference vegetation index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return ee.Image(landsat_image).normalizedDifference(['nir', 'red'])\
        .rename(['ndvi']).unmask(0)


def fipar(landsat_image):
    """Fraction of intercepted photosynthetically active radiation

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    fipar = m1*ndvi + b1 (Eqn)
    m1 =1
    b1 = -0.05

    References
    ----------
    .. [Fisher2008] J. Fisher, K. Tu, D. Baldocchi,
       Global estimates of the land-atmosphere water flux based on monthly
       AVHRR and ISLSCP-II data, validated at 16 FLUXNET sites,
       Remote Sensing of Environment,
       https://doi.org/10.1016/j.rse.2007.06.025

    """
    ndvi_clamp = ndvi(landsat_image).clamp(0.0, 1.0)


    return ndvi_clamp.multiply(1).subtract(0.05).clamp(0.0, 1.0)\
        .rename('fipar')


def lai(landsat_image):
    """Leaf area index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Ross1976] J. Ross, Radiative transfer in plant communities,
       In J. L. Monteith (Ed.),Vegetation and the atmosphere.

    """
    return ee.Image(landsat_image).expression('-log(1 - fIPAR) / (KPAR)', {
            'fIPAR': fipar(ee.Image(landsat_image)),
            'KPAR': ee.Number(0.5)
        }).rename('lai')


def ndwi(landsat_image):
    """Normalized difference water index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return ee.Image(landsat_image).normalizedDifference(['green', 'nir'])\
        .rename('ndwi').unmask(0)


def emissivity(landsat_image):
    """Broad-band surface emissivity

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return lai(landsat_image).multiply(0.01).add(0.95).min(0.98)\
        .rename('emissivity')


def lst(landsat_image):
    """Land surface temperature (lst)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    lai_img = lai(landsat_image)
    ndvi_img = ndvi(landsat_image)

    # Narrow band transmissivity
    e_NB = landsat_image.expression('0.97 + (0.0033 * LAI)', {'LAI': lai_img})
    e_NB = e_NB.where(lai_img.gt(3), 0.98).rename('e_NB')

    # Water and Snow filter
    e_NB = e_NB.where(ndvi_img.lt(0), 0.99).rename('e_NB')
    e_NB = e_NB.where(ndvi_img.lt(0), 0.99).rename('e_NB')

    lst = landsat_image.expression(
        'Tb / ( 1 + ( ( comp_onda * Tb / fator) * log_eNB))', {
            'Tb': landsat_image.select('tir'),
            'comp_onda': ee.Number(1.115e-05),
            'log_eNB': e_NB.log(),
            'fator': ee.Number(1.438e-02),
        }).rename('lst')

    return lst


def savi(landsat_image):
    """Soil adjusted vegetation index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    savi = landsat_image.expression(
        '((1 + 0.5) * (B5 - B4)) / (0.5 + (B5 + B4))', {
            'B4': landsat_image.select('red'),
            'B5': landsat_image.select('nir'),
        }).rename('savi')
    savi = savi.where(savi.gt(0.689), 0.689)

    return savi


def albedo_l457(landsat_image):
    """Albedo (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Tasumi2008] M. Tasumi, R. Allen, R; Trezza,
       At-Surface Reflectance and Albedo from Satellite for Operational
       Calculation of Land Surface Energy Balance, Journal of Hydrology,
       https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51)

    """
    albedo = landsat_image.expression(
        '(0.254 * B1) + (0.149 * B2) + (0.147 * B3) + (0.311 * B4) + '
        '(0.103 * B5) + (0.036 * B7)', {
            'B1' : landsat_image.select(['blue']),
            'B2' : landsat_image.select(['green']),
            'B3' : landsat_image.select(['red']),
            'B4' : landsat_image.select(['nir']),
            'B5' : landsat_image.select(['swir1']),
            'B7' : landsat_image.select(['swir2']),
        }).rename('albedo')

    return albedo


def albedo_l89(landsat_image):
    """Albedo (Landsat 8/9)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Ke2016] Y. Ke, J. Im, S. Park, H. Gong,
       Downscaling of MODIS One Kilometer Evapotranspiration
       Using Landsat-8 Data and Machine Learning Approaches,
       Remote Sensing, https://doi.org/10.3390/rs8030215

    """

    # CGM - These coefficients don't sum to 1.0?
    albedo = landsat_image.expression(
        '(0.130 * B1) + (0.115 * B2) + (0.143 * B3) + (0.180 * B4) + '
        '(0.281 * B5) + (0.108 * B6) + (0.042 * B7)', {
            'B1' : landsat_image.select(['ultra_blue']),
            'B2' : landsat_image.select(['blue']),
            'B3' : landsat_image.select(['green']),
            'B4' : landsat_image.select(['red']),
            'B5' : landsat_image.select(['nir']),
            'B6' : landsat_image.select(['swir1']),
            'B7' : landsat_image.select(['swir2']),
        }).rename('albedo')

    # # CGM - Just curious if the sum reducer would work
    # albedo = landsat_image\
    #     .select(['ultra_blue', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2'])\
    #     .multiply([0.130, 0.115, 0.143, 0.180, 0.281, 0.108, 0.042])\
    #     .reduce(ee.Reducer.sum())\
    #     .rename(['albedo'])

    return albedo


def cloud_mask_sr_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 66  (01000010)
    water value = 68   (01000100)

    References
    ----------

    """
    quality = landsat_image.select('pixel_qa')
    c01 = quality.eq(66)  # Clear (01000010)
    c02 = quality.eq(68)  # Water (01000100)
    mask = c01.Or(c02)

    return mask


def cloud_mask_sr_l8(landsat_image):
    """Cloud mask (Landsat 8)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 322  (00101000010)
    water value = 324   (00101000100)

    References
    ----------

    """
    quality = landsat_image.select('pixel_qa')
    c01 = quality.eq(322)
    c02 = quality.eq(324)
    c03 = quality.eq(1346) #       (10101000010)
    mask = c01.Or(c02).Or(c03)

    return mask


def cloud_mask_C2_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 5440  (0001010101000000)
    water value = 5504   (0001010110000000)

    References
    ----------

    """
    quality = landsat_image.select('QA_PIXEL')
    c01 = quality.eq(5440)
    c02 = quality.eq(5504)
    mask = c01.Or(c02)

    return mask


def cloud_mask_C2_l89(landsat_image):
    """Cloud mask (Landsat 8/9)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 21824  (0001010101000000)
    water value = 21952   (0001010110000000)

    References
    ----------

    """
    quality = landsat_image.select('QA_PIXEL')
    c01 = quality.eq(21824)
    c02 = quality.eq(21952)
    mask = c01.Or(c02)

    return mask
