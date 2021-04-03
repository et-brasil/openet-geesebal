import ee


def ndvi(landsat_image):
    """Normalized difference vegetation index"""
    return ee.Image(landsat_image).normalizedDifference(['nir', 'red'])\
        .rename(['ndvi'])


def fipar(landsat_image):
    """Fraction of intercepted photosynthetically active radiation"""
    ndvi_clamp = ndvi(landsat_image).clamp(0.0, 1.00)

    return ndvi_clamp.multiply(1).subtract(ee.Number(0.05)).clamp(0, 1)\
        .rename('fipar')


def lai(landsat_image):
    """Leaf area index"""
    return ee.Image(landsat_image).expression('log(1 - fIPAR)/(KPAR)', {
            'fIPAR': fipar(ee.Image(landsat_image)),
            'KPAR': ee.Number(0.5)
        }).rename('lai')


def ndwi(landsat_image):

    return ee.Image(landsat_image).normalizedDifference(['green', 'nir'])\
        .rename('ndwi')


def emissivity(landsat_image):
    """Broad-band surface emissivity"""
    lai_img = lai(landsat_image)
    
    return landsat_image.expression('0.95 + 0.01 * LAI', {'LAI': lai_img})\
        .where(lai_img.gt(3), 0.98).rename('emissivity')


def lst(landsat_image):
    """Land surface temperature"""
    lai_img = lai(landsat_image)
    ndvi_img = ndvi(landsat_image)

    # Narrow band transmissivity
    e_NB = landsat_image.expression('0.97 + (0.0033 * LAI)', {'LAI': lai_img})
    e_NB = e_NB.where(lai_img.gt(3), 0.98).rename('e_NB')
    
    # Water and Snow filter
    e_NB = e_NB.where(ndvi_img.lt(0), 0.99).rename('e_NB')
    e_NB = e_NB.where(ndvi_img.lt(0), 0.99).rename('e_NB')

    lst = landsat_image.expression(
        'Tb / ( 1+ ( ( comp_onda * Tb / fator) * log_eNB))', {
            'Tb': landsat_image.select('tir'),
            'comp_onda': ee.Number(1.115e-05),
            'log_eNB': e_NB.log(),
            'fator': ee.Number(1.438e-02),
        }).rename('lst')
    
    return lst


def savi(landsat_image):
    """Soil adjusted vegetation index"""
    savi = landsat_image.expression(
        '((1 + 0.5)*(B5 - B4)) / (0.5 + (B5 + B4))', {
            'B4': landsat_image.select('red'),
            'B5': landsat_image.select('nir'),
        }).rename('savi')
    savi = savi.where(savi.gt(0.689), 0.689)
    
    return savi


def albedo_l457(landsat_image):
    """Albedo (Landsat 4/5/7)"""
    albedo = landsat_image.expression(
        '(0.254*B1) + (0.149*B2) + (0.147*B3) + (0.311*B4) + (0.103*B5) + (0.036*B7)', {
            'B1' : landsat_image.select(['blue']),
            'B2' : landsat_image.select(['green']),
            'B3' : landsat_image.select(['red']),
            'B4' : landsat_image.select(['nir']),
            'B5' : landsat_image.select(['swir1']),
            'B7' : landsat_image.select(['swir2']),
        }).rename('albedo')
    
    return albedo


def albedo_l8(landsat_image):
    """Albedo (Landsat 8)"""
    albedo = landsat_image.expression(
        '(0.130*B1) + (0.115*B2) + (0.143*B3) + (0.180*B4) + (0.281*B5) + (0.108*B6) + (0.042*B7)', {
            'B1' : landsat_image.select(['ultra_blue']),
            'B2' : landsat_image.select(['blue']),
            'B3' : landsat_image.select(['green']),
            'B4' : landsat_image.select(['red']),
            'B5' : landsat_image.select(['nir']),
            'B6' : landsat_image.select(['swir1']),
            'B7' : landsat_image.select(['swir2']),
        }).rename('albedo')

    return albedo


def cloud_mask_sr_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)"""
    quality = landsat_image.select('pixel_qa')
    c01 = quality.eq(66)
    c02 = quality.eq(68)
    mask = c01.Or(c02)

    return mask


def cloud_mask_sr_l8(landsat_image):
    """Cloud mask (Landsat 8)"""
    quality = landsat_image.select('pixel_qa')
    c01 = quality.eq(322)
    c02 = quality.eq(324)
    c03 = quality.eq(1346)
    mask = c01.Or(c02).Or(c03)

    return mask

def cloud_mask_C2_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)"""
    quality = landsat_image.select('QA_PIXEL')
    c01 = quality.eq(5440)
    c02 = quality.eq(5504)
    mask = c01.Or(c02)

    return mask

def cloud_mask_C2_l8(landsat_image):
    """Cloud mask (Landsat 8)"""
    quality = landsat_image.select('QA_PIXEL')
    c01 = quality.eq(21824)
    c02 = quality.eq(21952)
    c03 = quality.eq(1346)
    mask = c01.Or(c02).Or(c03)

    return mask    
