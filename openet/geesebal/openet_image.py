import pprint

import ee

from openet.geesebal import openet_landsat as landsat
from openet.geesebal import model
from openet.geesebal import utils


def lazy_property(fn):
    """Decorator that makes a property lazy-evaluated
    https://stevenloria.com/lazy-properties/
    """
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazy_property


class Image():
    """Google Earth Engine SEBAL - GEESEBAL for Landsat image"""

    def __init__(
            self, image,
            meteorology_source_inst='NASA/NLDAS/FORA0125_H002',
            meteorology_source_daily='IDAHO_EPSCOR/GRIDMET',
            elev_source='USGS/SRTMGL1_003',
            ndvi_cold=5,
            ndvi_hot=10,
            lst_cold=20,
            lst_hot=20,
            **kwargs,
            ):

        """Construct a generic GEESEBAL Image

        Parameters
        ----------
        image : ee.Image
            A "prepped" GEESEBAL input image.
            Image must have bands:
                ndvi, lai, savi, lst, emissivity, ndwi, albedo
            Image must have properties:
                SUN_ELEVATION, system:id, system:index, system:time_start
        meteorology_source_inst : str, optional
            Instantaneous meteorology source collection ID.
            Collection supported:
                NASA/NLDAS/FORA0125_H002
            Meteorology collection must have bands:
                tair, ux, rh, rso_inst
        meteorology_source_daily : str
            Daily meteorology source collection ID.
            Collection supported:
                IDAHO_EPSCOR/GRIDMET
            Meteorology collection must have bands:
                tmin,tmax, rso24h
        elev_source : str, optional
            Elevation source image ID.
        ndvi_cold : int, ee.Number, optional
            NDVI Percentile value to determinate cold pixel
        ndvi_hot : int, ee.Number, optional
            NDVI Percentile value to determinate hot pixel
        lst_cold : int, ee.Number, optional
            LST Percentile value to determinate cold pixel
        lst_hot : int, ee.Number, optional
            LST Percentile value to determinate hot pixel
        kwargs : dict, optional
            et_reference_source : str, float
                Reference ET source (the default is None).
                Parameter is required if computing 'et_fraction' or 'et_reference'.
            et_reference_band : str
                Reference ET band name (the default is None).
                Parameter is required if computing 'et_fraction' or 'et_reference'.
            et_reference_factor : float, None
                Reference ET scaling factor.  The default is None which is
                equivalent to 1.0 (or no scaling).
            et_reference_resample : {'nearest', 'bilinear', 'bicubic', None}
                Reference ET resampling.  The default is None which is
                equivalent to nearest neighbor resampling.

        Notes
        -----
        Standard percentiles are from Allen et al. (2013)
        """

        self.image = image

        # Copy system properties
        self._id = self.image.get('system:id')
        self._index = self.image.get('system:index')
        self._time_start = self.image.get('system:time_start')
        self._properties = {
            'system:index': self._index,
            'system:time_start': self._time_start,
            'image_id': self._id,
        }
        # Build SCENE_ID from the (possibly merged) system:index
        scene_id = ee.List(ee.String(self._index).split('_')).slice(-3)
        self._scene_id = ee.String(scene_id.get(0)).cat('_')\
            .cat(ee.String(scene_id.get(1))).cat('_')\
            .cat(ee.String(scene_id.get(2)))

        # Build WRS2_TILE from the scene_id
        self._wrs2_tile = ee.String('p').cat(self._scene_id.slice(5, 8))\
            .cat('r').cat(self._scene_id.slice(8, 11))

        # Set server side date/time properties using the 'system:time_start'
        self._date = ee.Date(self._time_start)
        self._year = ee.Number(self._date.get('year'))
        self._month = ee.Number(self._date.get('month'))
        self._start_date = ee.Date(utils.date_to_time_0utc(self._date))
        self._end_date = self._start_date.advance(1, 'day')
        self._doy = ee.Number(self._date.getRelative('day', 'year')).add(1).int()

        # Model input parameters
        self._meteorology_source_inst = meteorology_source_inst
        self._meteorology_source_daily = meteorology_source_daily
        self._elev_source = elev_source
        self._ndvi_cold = ndvi_cold
        self._ndvi_hot = ndvi_hot
        self._lst_cold = lst_cold
        self._lst_hot = lst_hot

        # Reference ET parameters
        try:
            self.et_reference_source = kwargs['et_reference_source']
        except:
            self.et_reference_source = None
        try:
            self.et_reference_band = kwargs['et_reference_band']
        except:
            self.et_reference_band = None
        try:
            self.et_reference_factor = kwargs['et_reference_factor']
        except:
            self.et_reference_factor = None
        try:
            self.et_reference_resample = kwargs['et_reference_resample']
        except:
            self.et_reference_resample = None

        # Check reference ET parameters
        if (self.et_reference_factor and
                not utils.is_number(self.et_reference_factor)):
            raise ValueError('et_reference_factor must be a number')
        if self.et_reference_factor and self.et_reference_factor < 0:
            raise ValueError('et_reference_factor must be greater than zero')
        resample_methods = ['nearest', 'bilinear', 'bicubic']
        if (self.et_reference_resample and
                self.et_reference_resample.lower() not in resample_methods):
            raise ValueError('unsupported et_reference_resample method')

        self.geometry = self.image.select(0).geometry()
        self.proj = self.image.select(0).projection()
        self.latlon = ee.Image.pixelLonLat().reproject(self.proj)
        self.coords = self.latlon.select(['longitude', 'latitude' ])

        # Image projection and geotransform
        self.crs = image.projection().crs()
        self.transform = ee.List(ee.Dictionary(
            ee.Algorithms.Describe(image.projection())).get('transform'))

    @classmethod
    def from_image_id(cls, image_id, **kwargs):
        """Constructs an GEESEBAL Image instance from an image ID

        Parameters
        ----------
        image_id : str
            An earth engine image ID.
            (i.e. 'LANDSAT/LC08/C01/T1_SR/LC08_044033_20170716')
        kwargs
            Keyword arguments to pass through to model init.

        Returns
        -------
        new instance of Image class
        """

        # DEADBEEF - Should the supported image collection IDs and helper
        # function mappings be set in a property or method of the Image class?
        collection_methods = {
            'LANDSAT/LT04/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LT05/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LE07/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LC08/C01/T1_SR': 'from_landsat_c1_sr',
            'LANDSAT/LT04/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LT05/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LE07/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LC08/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LC09/C02/T1_L2': 'from_landsat_c2_sr',
        }

        try:
            method_name = collection_methods[image_id.rsplit('/', 1)[0]]
        except KeyError:
            raise ValueError(f'unsupported collection ID: {image_id}')
        except Exception as e:
            raise Exception(f'unhandled exception: {e}')

        method = getattr(Image, method_name)

        return method(ee.Image(image_id), **kwargs)

    @classmethod
    def from_landsat_c1_sr(cls, sr_image, cloudmask_args={}, **kwargs):
        """Returns a GEESEBAL Image instance from a Landsat Collection 1 SR image

        Parameters
        ----------
        sr_image : ee.Image, str
            A raw Landsat Collection 1 SR image or image ID.
        cloudmask_args : dict
            keyword arguments to pass through to cloud mask function
        kwargs : dict
            Keyword arguments to pass through to Image init function

        Returns
        -------
        Image
        """

        sr_image = ee.Image(sr_image)

        # Use the SATELLITE property identify each Landsat type
        # This should be equivalent tot he SPACECRAFT_ID property in collection 2
        spacecraft_id = ee.String(sr_image.get('SATELLITE'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_4': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'pixel_qa'],
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'pixel_qa'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'pixel_qa'],
            'LANDSAT_8': ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10',
                          'pixel_qa'],
        })
        output_bands = ee.Dictionary({
            'LANDSAT_4': ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir',
                          'pixel_qa'],
            'LANDSAT_5': ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir',
                          'pixel_qa'],
            'LANDSAT_7': ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir',
                          'pixel_qa'],
            'LANDSAT_8': ['ultra_blue', 'blue', 'green', 'red', 'nir',
                          'swir1', 'swir2', 'tir', 'pixel_qa'],
        })
        scalars = ee.Dictionary({
            'LANDSAT_4': [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.1, 1.0],
            'LANDSAT_5': [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.1, 1.0],
            'LANDSAT_7': [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.1, 1.0],
            'LANDSAT_8': [0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
                          0.0001, 0.0001, 0.1, 1.0],
        })
        prep_image = sr_image\
            .select(input_bands.get(spacecraft_id), output_bands.get(spacecraft_id))\
            .multiply(ee.Image.constant(ee.List(scalars.get(spacecraft_id))))\
            .set({'SPACECRAFT_ID': spacecraft_id})

        # TODO: Restructure these to avoid the "If" calls if possible
        albedo = ee.Algorithms.If(
            spacecraft_id.compareTo(ee.String('LANDSAT_8')),
            landsat.albedo_l457(prep_image),
            landsat.albedo_l89(prep_image))

        cloud_mask = ee.Algorithms.If(
            spacecraft_id.compareTo(ee.String('LANDSAT_8')),
            landsat.cloud_mask_sr_l457(sr_image),
            landsat.cloud_mask_sr_l8(sr_image))

        # # Default the cloudmask flags to True if they were not
        # # Eventually these will probably all default to True in openet.core
        # if 'shadow_flag' not in cloudmask_args.keys():
        #     cloudmask_args['shadow_flag'] = True
        # if 'snow_flag' not in cloudmask_args.keys():
        #     cloudmask_args['snow_flag'] = True
        # cloud_mask = openet.core.common.landsat_c1_sr_cloud_mask(
        #     sr_image, **cloudmask_args)

        # Build the input image
        input_image = ee.Image([
            landsat.ndvi(prep_image),
            landsat.lai(prep_image),
            landsat.savi(prep_image),
            landsat.lst(prep_image),
            landsat.emissivity(prep_image),
            landsat.ndwi(prep_image),
            albedo,
        ])

        # Calculate sun elevation from zenith
        sun_elevation = ee.Number(90)\
            .subtract(ee.Number(sr_image.get('SOLAR_ZENITH_ANGLE')))

        # Apply the cloud mask and add properties
        input_image = input_image.updateMask(cloud_mask).set({
            'system:index': sr_image.get('system:index'),
            'system:time_start': sr_image.get('system:time_start'),
            'system:id': sr_image.get('system:id'),
            'SUN_ELEVATION':sun_elevation,
        })

        # Instantiate the class
        return cls(input_image, **kwargs)

    @classmethod
    def from_landsat_c2_sr(cls, sr_image, cloudmask_args={}, **kwargs):
        """Returns a GEESEBAL Image instance from a Landsat Collection 2 SR image

        Parameters
        ----------
        sr_image : ee.Image, str
            A raw Landsat Collection 2 SR image or image ID.
        cloudmask_args : dict
            keyword arguments to pass through to cloud mask function
        kwargs : dict
            Keyword arguments to pass through to Image init function

        Returns
        -------
        Image
        """

        sr_image = ee.Image(sr_image)

        # Use the SPACECRAFT_ID property identify each Landsat type
        spacecraft_id = ee.String(sr_image.get('SPACECRAFT_ID'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_4': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
                          'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_5': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
                          'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_7': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
                          'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_8': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5',
                          'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
            'LANDSAT_9': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5',
                          'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
        })
        output_bands = ee.Dictionary({
            'LANDSAT_4': ['blue', 'green', 'red', 'nir',
                          'swir1', 'swir2', 'lst', 'QA_PIXEL'],
            'LANDSAT_5': ['blue', 'green', 'red', 'nir',
                          'swir1', 'swir2', 'lst', 'QA_PIXEL'],
            'LANDSAT_7': ['blue', 'green', 'red', 'nir',
                          'swir1', 'swir2', 'lst', 'QA_PIXEL'],
            'LANDSAT_8': ['ultra_blue', 'blue', 'green', 'red', 'nir',
                          'swir1', 'swir2', 'lst', 'QA_PIXEL'],
            'LANDSAT_9': ['ultra_blue', 'blue', 'green', 'red', 'nir',
                          'swir1', 'swir2', 'lst', 'QA_PIXEL'],
        })
        scalars = ee.Dictionary({
            'LANDSAT_4': [0.0000275, 0.0000275, 0.0000275, 0.0000275,
                          0.0000275, 0.0000275, 0.00341802, 1],
            'LANDSAT_5': [0.0000275, 0.0000275, 0.0000275, 0.0000275,
                          0.0000275, 0.0000275, 0.00341802, 1],
            'LANDSAT_7': [0.0000275, 0.0000275, 0.0000275, 0.0000275,
                          0.0000275, 0.0000275, 0.00341802, 1],
            'LANDSAT_8': [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275,
                          0.0000275, 0.0000275, 0.00341802, 1],
            'LANDSAT_9': [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275,
                          0.0000275, 0.0000275, 0.00341802, 1],
        })
        offsets = ee.Dictionary({
            'LANDSAT_4': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 1],
            'LANDSAT_5': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 1],
            'LANDSAT_7': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 1],
            'LANDSAT_8': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 1],
            'LANDSAT_9': [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 1],
        })

        prep_image = sr_image \
            .select(input_bands.get(spacecraft_id), output_bands.get(spacecraft_id)) \
            .multiply(ee.Image.constant(ee.List(scalars.get(spacecraft_id)))) \
            .add(ee.Image.constant(ee.List(offsets.get(spacecraft_id))))

        # CGM - Need to come up with a more robust approach,
        #   but this seems to work for now
        albedo = ee.Algorithms.If(
            ee.List(['LANDSAT_8', 'LANDSAT_9']).contains(spacecraft_id),
            landsat.albedo_l89(prep_image),
            landsat.albedo_l457(prep_image),
        )
        cloud_mask = ee.Algorithms.If(
            ee.List(['LANDSAT_8', 'LANDSAT_9']).contains(spacecraft_id),
            landsat.cloud_mask_C2_l89(sr_image),
            landsat.cloud_mask_C2_l457(sr_image),
        )
        # albedo = ee.Algorithms.If(
        #     spacecraft_id.compareTo(ee.String('LANDSAT_8')),
        #     landsat.albedo_l457(prep_image),
        #     landsat.albedo_l89(prep_image))
        # cloud_mask = ee.Algorithms.If(
        #     spacecraft_id.compareTo(ee.String('LANDSAT_8')),
        #     landsat.cloud_mask_C2_l457(sr_image),
        #     landsat.cloud_mask_C2_l89(sr_image))

        # # Default the cloudmask flags to True if they were not
        # # Eventually these will probably all default to True in openet.core
        # if 'cirrus_flag' not in cloudmask_args.keys():
        #     cloudmask_args['cirrus_flag'] = True
        # if 'dilate_flag' not in cloudmask_args.keys():
        #     cloudmask_args['dilate_flag'] = True
        # if 'shadow_flag' not in cloudmask_args.keys():
        #     cloudmask_args['shadow_flag'] = True
        # if 'snow_flag' not in cloudmask_args.keys():
        #     cloudmask_args['snow_flag'] = True
        # cloud_mask = openet.core.common.landsat_c2_sr_cloud_mask(
        #     sr_image, **cloudmask_args)

        # Build the input image
        # Don't compute LST since it is being provided
        input_image = ee.Image([
            prep_image.select(['lst']),
            landsat.ndvi(prep_image),
            landsat.lai(prep_image),
            landsat.savi(prep_image),
            landsat.emissivity(prep_image),
            landsat.ndwi(prep_image),
            albedo,
        ])

        # Apply the cloud mask and add properties
        input_image = input_image.updateMask(cloud_mask)\
            .set({'system:index': sr_image.get('system:index'),
                  'system:time_start': sr_image.get('system:time_start'),
                  'system:id': sr_image.get('system:id'),
                  'SUN_ELEVATION': sr_image.get('SUN_ELEVATION'),
            })

        # Instantiate the class
        return cls(input_image, **kwargs)

    def calculate(self, variables=['ndvi', 'lst', 'et', 'et_fraction']):
        """Return a multiband image of calculated variables

        Parameters
        ----------
        variables : list

        Returns
        -------
        ee.Image
        """

        output_images = []
        for v in variables:
            if v.lower() == 'et':
                output_images.append(self.et.float())
            elif v.lower() == 'et_fraction':
                output_images.append(self.et_fraction.float())
            elif v.lower() == 'lst':
                output_images.append(self.lst.float())
            elif v.lower() == 'ndvi':
                output_images.append(self.ndvi.float())
            # CGM - Calculate method must support the mask and time bands
            elif v.lower() == 'mask':
                output_images.append(self.mask)
            elif v.lower() == 'time':
                output_images.append(self.time)
            else:
                raise ValueError(f'unsupported variable: {v}')

        return ee.Image(output_images).set(self._properties)

    @lazy_property
    def ndvi(self,):

        return self.image.select(['ndvi']).set(self._properties)

    @lazy_property
    def emissivity(self,):

        return self.image.select(['emissivity']).set(self._properties)

    @lazy_property
    def ndwi(self,):

        return self.image.select(['ndwi']).set(self._properties)

    @lazy_property
    def lai(self,):

        return self.image.select(['lai']).set(self._properties)

    @lazy_property
    def albedo(self,):

        return self.image.select(['albedo']).set(self._properties)

    @lazy_property
    def savi(self,):

        return self.image.select(['savi']).set(self._properties)

    @lazy_property
    def lst(self,):

        return self.image.select(['lst']).set(self._properties)

    @lazy_property
    def et(self,):

        et = model.et(image=self.image,
                        ndvi=self.ndvi,
                        ndwi=self.ndwi,
                        lst=self.lst,
                        albedo=self.albedo,
                        emissivity=self.emissivity,
                        savi=self.savi,
                        # lai=self.lai,
                        meteo_inst_source=self._meteorology_source_inst,
                        meteo_daily_source=self._meteorology_source_daily,
                        elev_product=self._elev_source,
                        ndvi_cold=self._ndvi_cold,
                        ndvi_hot=self._ndvi_hot,
                        lst_cold=self._lst_cold,
                        lst_hot=self._lst_hot,
                        time_start=self._time_start,
                        geometry_image=self.geometry,
                        proj=self.proj,
                        coords=self.coords,
                        )

        return et.set(self._properties)

    # @lazy_property
    # def et_fraction(self):
    #
    #     et_fr = model.et_fraction(
    #         self.image,
    #         self._time_start,
    #         self.et,
    #         self.et_reference_source, self.et_reference_band,
    #         self.et_reference_factor,
    #     )
    #
    #     return et_fr.set(self._properties)

    @lazy_property
    def et_reference(self):
        """Reference ET for the image date"""
        if utils.is_number(self.et_reference_source):
            # Interpret numbers as constant images
            # CGM - Should we use the ee_types here instead?
            #   i.e. ee.ee_types.isNumber(self.et_reference_source)
            et_reference_img = ee.Image.constant(self.et_reference_source)
        elif type(self.et_reference_source) is str:
            # Assume a string source is an image collection ID (not an image ID)
            et_reference_coll = ee.ImageCollection(self.et_reference_source)\
                .filterDate(self._start_date, self._end_date)\
                .select([self.et_reference_band])
            et_reference_img = ee.Image(et_reference_coll.first())
            if self.et_reference_resample in ['bilinear', 'bicubic']:
                et_reference_img = et_reference_img\
                    .resample(self.et_reference_resample)
        else:
            raise ValueError('unsupported et_reference_source: {}'.format(
                self.et_reference_source))

        if self.et_reference_factor:
            et_reference_img = et_reference_img.multiply(self.et_reference_factor)

        # Map ETr values directly to the input (i.e. Landsat) image pixels
        # The benefit of this is the ETr image is now in the same crs as the
        #   input image.  Not all models may want this though.
        # Note, doing this will cause the reference ET to be cloud masked.
        # CGM - Should the output band name match the input ETr band name?
        return self.ndvi.multiply(0).add(et_reference_img)\
            .rename(['et_reference']).set(self._properties)

    @lazy_property
    def et_fraction(self):
        """Fraction of reference ET (equivalent to the Kc)"""
        return self.et.divide(self.et_reference) \
            .rename(['et_fraction']).set(self._properties)

    # CGM - The mask band is currently needed for the time band
    # If the model does not do any additional masking we might be able to
    #   build the mask from the NDVI or QA band instead.
    @lazy_property
    def mask(self,):
        """Mask of all active pixels (based on the final et)"""

        mask = self.et.multiply(0).add(1).updateMask(1).uint8().rename(['mask'])

        return  mask.set(self._properties)

    # CGM - The image class must have a "time" method for the interpolation
    # I'm not sure if it needs to be built from the active pixels mask
    #   or could be built from the NDVI or QA band instead.
    @lazy_property
    def time(self,):
        """Return an image of the 0 UTC time (in milliseconds)"""

        time = self.mask \
            .double().multiply(0).add(utils.date_to_time_0utc(self._date)) \
            .rename(['time'])

        return time.set(self._properties)
