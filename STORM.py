"""
STORM - STOchastic Rainfall Model STORM produces realistic regional or
watershed rainfall under various climate scenarios based on empirical-
stochastic selection of historical rainfall characteristics.

Based on Singer, M. B., and Michaelides, K. (2017), Deciphering the expression
of climate change within the Lower Colorado River basin by stochastic
simulation of convective rainfall.

version name: STORM_v1_17_py

Author: Michael Singer 2017. Python transliteration by Daniel Hobley 26/03/18
following Matlab version of STORM dated 03/23/2018
"""

import os
import time
import numpy as np
import datetime as date
from matplotlib.path import Path
from scipy.stats import genextreme, fisk
from six.moves import range
try:
    import shapefile as shp
except ImportError:
    print('Install the pyshp package to run STORM! ' +
          'e.g., conda install pyshp or pip install pyshp')
    raise ImportError

from matplotlib.pyplot import imshow, plot, show, scatter, axis, figure

PLOT_EACH_STORM = False
# Use this to visualise every storm over the catchment. Warning: will be slow
# and tedious for long runs or many storms


def storm(mode, numsims, numsimyrs, seasons, ptot_scenario,
          storminess_scenario, ptot_scenario2, storminess_scenario2,
          ET_scenario,
          storminess_scaling_factor=0.05,
          storm_stepchange=0.25,
          storminess_scaling_factor2=0.05,
          storm_stepchange2=0.25,
          ptot_scaling_factor=0.05,
          ptot_scaling_factor2=0.05,
          ET_scaling_factor=0.25,
          numcurves=11,
          orographic_limits=(1350., 1500.)):
    """
    USAGE:

    >>> storm(mode, numsims, numsimyrs, seasons, ptot_scenario,
    ...       storminess_scenario, ptot_scenario2, storminess_scenario2,
    ...       ET_scenario, storminess_scaling_factor=0.05,
    ...       storm_stepchange=0.25, storminess_scaling_factor2=0.05,
    ...       storm_stepchange2=0.25, ptot_scaling_factor=0.05,
    ...       ptot_scaling_factor2=0.05, ET_scaling_factor=0.25)

    *mode* is a string in single quotations containing one of these two
    options:
    'Validation' for validating STORM's output against historical observations
    'Simulation' for simulating stochastic rainfall under various climate
    scenarios
    *mode* affects where STORM generates output--either for a set of
    gauge locations for comparison with observational data (Validation) or on
    an aribitrary watershed grid (Simulation).

    *numsims* is the integer number of x-year (*numsimyrs*) simulations to be
    run as a batch. Each simulation will produce its own output folder that is
    timestamped with relevant simulation information

    *numsimyrs* is the integer number of years in each simulation. Recommended
    value >=30

    *seasons* is number specifying whether there are 1 or 2 seasons. If
    there are two seasons, pdfs will be required specifying rainfall
    characteristics for the 2nd season (Ptot, Duration, etc). Also, ptot and
    storminess scenarios will be required arguments in the function (see
    below). Currently, season one runs from May 1 to Sept 30 and season two
    runs from Oct 1 to Apr 30, but these dates can be changed.

    *ptot_scenario* is a string in single quotations specifying the type of
    Ptot simulation to be run. This variable characterizes the degree of
    overall wetness as measured by annual precipitation totals. There are
    currently 5 options here:
      'ptotC': for control (stationary) conditions akin to observations of
               precipitation totals
      'ptot+':  for step change of increased overall wetness from observed
      'ptot-':  for step change of decreased overall wetness from observed
      'ptotT+': for progressive increasing trend in storminess from observed
      'ptotT-':  for progressive decreasing trend in storminess from observed

    *storminess_scenario* is a string in single quotations specifying the type
    of Storminess simulation to be run. This variable characterizes the degree
    of overall storminess as measured by storm intensity for a particular
    duration. There are currently 5 options here:
      'stormsC': for no change in storminess -- use intensity-duration
                 relationships from historical records
      'storms+': for step change of increased storminess from observed
      'storms-': for step change of decreased storminess from observed
      'stormsT+': for progressive increasing trend in storminess from observed
      'stormsT-': for progressive decreasing trend in storminess from observed.

    *ptot_scenario2* is a string in single quotations specifying the type of
    Ptot simulation to be run. This variable characterizes the degree of
    overall wetness as measured by annual precipitation totals. There are
    currently 5 options here:
      'ptot2C': for control (stationary) conditions akin to observations of
                precipitation totals
      'ptot2+': for step change of increased overall wetness from observed
      'ptot2-': for step change of decreased overall wetness from observed
      'ptot2T+': for progressive increasing trend in storminess from observed
      'ptot2T-': for progressive decreasing trend in storminess from observed
    For seasons = 1, use None for this argument.

    *storminess_scenario2* is a string in single quotations specifying the type
    of Storminess simulation to be run. This variable characterizes the degree
    of overall storminess as measured by storm intensity for a particular
    duration. There are currently 5 options here:
      'storms2C': for no change in storminess -- use intensity-duration
                  relationships from historical records
      'storms2+': for step change of increased storminess from observed
      'storms2-': for step change of decreased storminess from observed
      'storms2T+:' for progressive increasing trend in storminess from observed
      'storms2T-': for progressive decreasing trend in storminess from
                   observed.
    For seasons == 1, use None for this argument.

    *ET_scenario* is a string in single quotations specifying the type of
    Storminess simulation to be run. This variable characterizes the degree of
    evaporative demand and its reflection of climate change. There are
    currently 5 options here:
      'ETC' for no change in ET - use historical records
      'ET+' for step change of increased ET from observed
      'ET-' for step change of decreased ET from observed
    Note: TRENDS and SEASONAL CHANGES IN ET HAVE NOT YET BEEN CODED!!
      'ET_T+' for progressive increasing trend in ET from observed
      'ET_T-' for progressive decreasing trend in ET from observed.

    If used, the following input arguments control the rates of change and
    sizes of steps specified by the preceeding arguments:

    *storminess_scaling_factor* is a float that specifies the fractional change
        in intensity per year when storm_trend is applied in
        storminess_scenario.
    *storm_stepchange* is a float that specifies the value of fractional step
        change in intensity when storms+ or storms- are applied in
        storminess_scenario.
    *storminess_scaling_factor2* is a float that specifies the fractional
        change in intensity per year when storm_trend is applied in
        storminess_scenario2.
    *storm_stepchange2* is a float that specifies the value of fractional step
        change in intensity when storms+ or storms- are applied in
        storminess_scenario2.
    *ptot_scaling_factor* is a float that specifies the fractional change in
        wetness per year when storm_trend is applied in ptot_scenario.
    *ptot_scaling_factor2* is a float that specifies the fractional change in
        wetness per year when storm_trend is applied in ptot_scenario2.
    *ET_scaling_factor* is a float that specifies the value of fractional step
        change in ET when ET+ or ET- are applied in ET_scenario.

    Note that the orographic settings of STORM are at present hard-wired to
    those of Walnut Gulch. Modify the code directly if you wish to change the
    orographic rainfall dynamics.

    The code requires the following variables to be created as csv files and
    placed in a folder called 'model_input':

    Pdfs listed below should be generated after exploring the appropriate
    distribution type. If these are not of the distribution types detailed in
    Singer et al., 2018, the code below will have to be altered. These should
    be specified in CSV files, following the Matlab conventions of [shape (if
    needed), mu, sigma, lower_truncation_thresh (if needed),
    upper_truncation_thresh (if needed)]

    'Ptot_pdf' - A pdf fitted to all available station precip totals data (mm).
    'Ptot_pdf2' - A pdf fitted to all available station precip totals data for
        season 2 (mm).
    'Ptot_pdf_cc_up' - Historical pdf shifted toward wetter conditions
    'Ptot_pdf_cc_up2' - Historical pdf for season 2 shifted toward wetter
        conditions
    'Ptot_pdf_cc_down' - Historical pdf shifted toward drier conditions
    'Ptot_pdf_cc_down2' - Historical pdf for season 2 shifted toward drier
        conditions
    'Duration_pdf' - A pdf fitted to all available storm duration data (min)
    'Duration_pdf2' - A pdf fitted to all available storm duration data of
        season 2 (min)
    'Area_pdf' - A pdf fitted to all available storm area data (m^2)
    'Area_pdf2' - A pdf fitted to all available storm area data for season 2
        (m^2)
    'Int_arr_pdf' - A pdf fitted to all available interarrival times data (hrs)
    'Int_arr_pdf2' - A pdf fitted to all available interarrival times data for
        season 2 (hrs)
    'Recess_pdf' - A pdf of recession coefficients of intensity with distance
        from storm center (mm/hr/km).

    'Easting' Actual rain gauge location. Used when mode = 'Validation'
    'Northing' Actual rain gauge location. Used when mode = 'Validation'
    'gauges' Numbers of actual rain gauges. Used when mode = 'Validation'
    'gauge_elev' Elevations of actual rain gauges. Elevation value is extracted
        from field called 'RASTERVALU'. Used when  mode = 'Validation'
    'boundary/boundary.shp' Shapefile of the watershed boundary
    'point_elevations.shp' Shapefile of the elevations at each 'gauging' point
        on the grid within the watershed boundary
    'X' Longitudinal data for each grid point. Used when mode = 'Simulation'
    'Y' Latitudinal data for each grid point. Used when mode = 'Simulation'.
        X & Y will be sampled below to determine storm center location.
    'Storm_depth_data' Storm depth data (mm). Used when mode = 'Validation'.
    'Intensity_data' Intensity data (mm/hr). Used when mode = 'Validation'.
    'Duration_data' Duration data (min) for use when mode = 'Validation'
    'fuzz' A vector of fuzzy tolerance values for intensity selection, assumed
        based on differences between PI-PD curves(mm/hr).
    'bndry_buffer.shp' A shapefile of the watershed boundary buffer (5000 m).
        Used for storm center determination allowing for partial storm
        coverage.
    'ET_monthly_day' A matrix of daytime ET values organized with one month per
        column, drawn from data or calculations (mm).
    'ET_monthly_night' A matrix of nighttime ET values organized with one month
        per column, drawn from data or calculations (mm).
    """

    # check the inputs all make sense:
    assert mode in ('Validation', 'Simulation'), 'mode is invalid'
    assert seasons in (1, 2), 'seasons must be 1 or 2'
    assert ptot_scenario in ('ptotC', 'ptot+', 'ptot-', 'ptotT+', 'ptotT-')
    assert storminess_scenario in ('stormsC', 'storms+', 'storms-',
                                   'stormsT+', 'stormsT-')
    if seasons == 2:
        assert ptot_scenario2 in ('ptot2C', 'ptot2+', 'ptot2-', 'ptot2T+',
                                  'ptot2T-')
        assert storminess_scenario2 in ('storms2C', 'storms2+', 'storms2-',
                                        'storms2T+', 'storms2T-')
    else:
        assert (ptot_scenario2 is None) and (storminess_scenario2 is None)
    if ET_scenario in ('ET_T+', 'ET_T-'):
        print('Trend options for ET are not yet available')
        raise NameError
    else:
        assert ET_scenario in ('ETC', 'ET+', 'ET-')

    # Restatement of arguments in function
    XXX = ('Thank you for using STORM. Your variable settings are: ')
    print(XXX)
    XXX = ['Type of Model Run = ', mode]
    print(XXX)
    XXX = ['Number of Simulations = ' + str(numsims)]
    print(XXX)
    XXX = ['Number of Simulation Years = ' + str(numsimyrs)]
    print(XXX)
    if seasons == 2:
        XXX = ('Number of Seasons = Two')
    elif seasons == 1:
        XXX = ('Number of Seasons = One')
    print(XXX)

    if ptot_scenario == 'ptotC':
        XXX = ('Scenario of Total Precipitation (Season 1) = Control climate')
    elif ptot_scenario == 'ptot+':
        XXX = ('Scenario of Total Precipitation (Season 1) = Step change ' +
               'increase in seasonal precipitation')
    elif ptot_scenario == 'ptot-':
        XXX = ('Scenario of Total Precipitation (Season 1) = Step change ' +
               'decrease in seasonal precipitation')
    elif ptot_scenario == 'ptotT+':
        XXX = ('Scenario of Total Precipitation (Season 1) = Postive trend ' +
               'in seasonal precipitation')
    elif ptot_scenario == 'ptotT-':
        XXX = ('Scenario of Total Precipitation (Season 1) = Negative trend ' +
               'in seasonal precipitation')
    print(XXX)

    if storminess_scenario == 'stormsC':
        XXX = ('Scenario of Precipitation Intensity (Season 1) = Control ' +
               'climate')
    elif storminess_scenario == 'storms+':
        XXX = ('Scenario of Precipitation Intensity (Season 1) = Step ' +
               'change increase in seasonal storminess')
    elif storminess_scenario == 'storms-':
        XXX = ('Scenario of Precipitation Intensity (Season 1) = Step ' +
               'change decrease in seasonal storminess')
    elif storminess_scenario == 'stormsT+':
        XXX = ('Scenario of Precipitation Intensity (Season 1) = Postive ' +
               'trend in seasonal storminess')
    elif storminess_scenario == 'stormsT-':
        XXX = ('Scenario of Precipitation Intensity (Season 1) = Negative ' +
               'trend in seasonal storminess')
    print(XXX)

    if ptot_scenario2 == 'ptot2C':
        XXX = ('Scenario of Total Precipitation (Season 2) = Control climate')
    elif ptot_scenario2 == 'ptot2+':
        XXX = ('Scenario of Total Precipitation (Season 2) = Step change ' +
               'increase in seasonal precipitation')
    elif ptot_scenario2 == 'ptot2-':
        XXX = ('Scenario of Total Precipitation (Season 2) = Step change ' +
               'decrease in seasonal precipitation')
    elif ptot_scenario2 == 'ptot2T+':
        XXX = ('Scenario of Total Precipitation (Season 2) = Postive trend ' +
               'in seasonal precipitation')
    elif ptot_scenario2 == 'ptot2T-':
        XXX = ('Scenario of Total Precipitation (Season 2) = Negative trend ' +
               'in seasonal precipitation')
    print(XXX)

    if storminess_scenario2 == 'storms2C':
        XXX = ('Scenario of Precipitation Intensity (Season 2) = Control ' +
               'climate')
    elif storminess_scenario2 == 'storms2+':
        XXX = ('Scenario of Precipitation Intensity (Season 2) = Step ' +
               'change increase in seasonal storminess')
    elif storminess_scenario2 == 'storms2-':
        XXX = ('Scenario of Precipitation Intensity (Season 2) = Step ' +
               'change decrease in seasonal storminess')
    elif storminess_scenario2 == 'storms2T+':
        XXX = ('Scenario of Precipitation Intensity (Season 2) = Postive ' +
               'trend in seasonal storminess')
    elif storminess_scenario2 == 'storms2T-':
        XXX = ('Scenario of Precipitation Intensity (Season 2) = Negative ' +
               'trend in seasonal storminess')
    print(XXX)

    if ET_scenario == 'ETC':
        XXX = ('Scenario of ET = Control climate')
    elif ET_scenario == 'ET+':
        XXX = ('Scenario of ET = Step change increase in PET')
    elif ET_scenario == 'ET-':
        XXX = ('Scenario of ET = Step change decrease in PET')
    print(XXX)

    # load the necessary driving distributions and data
    cwd = os.getcwd()
    datapath = os.path.join(cwd, 'model_input')

    # These ptot_pdf's are simply mu, sigma
    if ptot_scenario == 'ptotC':
        Ptot_pdf = np.loadtxt(os.path.join(datapath, 'Ptot_pdf.csv'))
        # This is the pdf fitted to all available station precip data
        # (normal dist). It will be sampled below.
        mu = Ptot_pdf[0]
        sigma = Ptot_pdf[1]
    elif ptot_scenario == 'ptot+':
        Ptot_pdf_cc_up = np.loadtxt(os.path.join(
            datapath, 'Ptot_pdf_cc_up.csv'))
        # This is the pdf fitted to all available station
        # precip data (normal dist). It will be sampled below.
        mu = Ptot_pdf_cc_up[0]
        sigma = Ptot_pdf_cc_up[1]
    elif ptot_scenario == 'ptot-':
        Ptot_pdf_cc_down = np.loadtxt(os.path.join(
            datapath, 'Ptot_pdf_cc_down.csv'))
        # This is the pdf fitted to all available station
        # precip data (normal dist). It will be sampled below.
        mu = Ptot_pdf_cc_down[0]
        sigma = Ptot_pdf_cc_down[1]
    elif ptot_scenario == 'ptotT+':
        Ptot_pdf = np.loadtxt(os.path.join(datapath, 'Ptot_pdf.csv'))
        # This is the pdf fitted to all available station precip data
        # (normal dist). It will be modified as a trend and sampled below.
        mu = Ptot_pdf[0]
        sigma = Ptot_pdf[1]
    elif ptot_scenario == 'ptotT-':
        Ptot_pdf = np.loadtxt(os.path.join(datapath, 'Ptot_pdf.csv'))
        # This is the pdf fitted to all available station precip data
        # (normal dist). It will be modified as a trend and sampled below.
        mu = Ptot_pdf[0]
        sigma = Ptot_pdf[1]

    # #### matlab's GEV is (shape_param, scale(sigma), pos(mu))
    # note that in Scipy, we must add a minus to the shape param for a GEV
    # to match Matlab's implementation
    Duration_pdf = np.loadtxt(os.path.join(datapath, 'Duration_pdf.csv'))
    # formatting will be 1D: shape_param, mu, sigma, low_trunc, hi_trunc
    # note this is different from the .mat file
    Duration_pdf_GEV = {'shape': -Duration_pdf[0], 'sigma': Duration_pdf[2],
                        'mu': Duration_pdf[1]}
    try:
        Duration_pdf_GEV['trunc_interval'] = [Duration_pdf[3],
                                              Duration_pdf[4]]
    except IndexError:
        Duration_pdf_GEV['trunc_interval'] = [0., np.inf]
    # This is the pdf fitted to all available station duration data (GEV
    # dist). It will be sampled below.
    Area_pdf = np.loadtxt(os.path.join(datapath, 'Area_pdf.csv'))
    # formatting will be 1D: mu, sigma, low_trunc, hi_trunc
    Area_pdf_EV = {'shape': 0., 'sigma': Area_pdf[1],
                   'mu': Area_pdf[0]}
    try:
        Area_pdf_EV['trunc_interval'] = [Area_pdf[2],
                                         Area_pdf[3]]
    except IndexError:
        Area_pdf_EV['trunc_interval'] = [0., np.inf]
    # This is the pdf fitted to all available station area data (EV dist). It
    # will be sampled below.
    Int_arr_pdf = np.loadtxt(os.path.join(datapath, 'Int_arr_pdf.csv'))
    # formatting will be 1D: shape_param, mu, sigma, low_trunc, hi_trunc
    # note this is different from the .mat file
    Int_arr_pdf_GEV = {'shape': -Int_arr_pdf[0], 'sigma': Int_arr_pdf[2],
                       'mu': Int_arr_pdf[1]}
    try:
        Int_arr_pdf_GEV['trunc_interval'] = [Int_arr_pdf[3],
                                             Int_arr_pdf[4]]
    except IndexError:
        Int_arr_pdf_GEV['trunc_interval'] = [0., np.inf]
    # This is the pdf fitted to all available station interarrival time data
    # (GEV dist). It will be sampled below.

    # Load pdfs for season two
    if seasons == 2:
        Duration_pdf2 = np.loadtxt(os.path.join(datapath, 'Duration_pdf2.csv'))
        # This is the pdf for season2 fitted to all available station duration
        # data (log-logistic dist). It will be sampled below.
        # This is a log-logistic, or Fisk, distribution. This differs Matlab to
        # Python. Scipy treats the dist with a c-param, where:
        # c = 1/matlab_sigma
        # scale = exp(matlab_mean)
        beta_fisk = 1./Duration_pdf2[1]  # following WP's notation
        c_fisk = beta_fisk
        scale_fisk = np.exp(Duration_pdf2[0])
        Duration_pdf_Fisk2 = {'c': c_fisk, 'scale': scale_fisk}
        try:
            Duration_pdf_Fisk2['trunc_interval'] = [Duration_pdf2[2],
                                                    Duration_pdf2[3]]
        except IndexError:
            Duration_pdf_Fisk2['trunc_interval'] = [0., np.inf]
        Area_pdf2 = np.loadtxt(os.path.join(datapath, 'Area_pdf2.csv'))
        # This is the pdf fitted for season2 to all available station area
        # data (EV dist). It will be sampled below.
        Area_pdf_EV2 = {'shape': 0., 'sigma': Area_pdf2[1],
                        'mu': Area_pdf2[0]}
        try:
            Area_pdf_EV2['trunc_interval'] = [Area_pdf2[2],
                                                  Area_pdf2[3]]
        except IndexError:
            Area_pdf_EV2['trunc_interval'] = [0., np.inf]
        # This is the pdf fitted to all available station area data (EV dist).
        # It will be sampled below.
        Int_arr_pdf2 = np.loadtxt(os.path.join(datapath, 'Int_arr_pdf2.csv'))
        # This is the pdf for season2 fitted to all available station
        # interarrival time data (GEV dist). It will be sampled below.
        Int_arr_pdf_GEV2 = {'shape': -Int_arr_pdf2[0],
                            'sigma': Int_arr_pdf2[2],
                            'mu': Int_arr_pdf2[1]}
        try:
            Int_arr_pdf_GEV2['trunc_interval'] = [Int_arr_pdf2[3],
                                                  Int_arr_pdf2[4]]
        except IndexError:
            Int_arr_pdf_GEV2['trunc_interval'] = [0, np.inf]
        # This is the pdf fitted to all available station interarrival time
        # data (GEV dist). It will be sampled below.
        if ptot_scenario2 == 'ptot2C':
            Ptot_pdf2 = np.loadtxt(os.path.join(datapath, 'Ptot_pdf2.csv'))
            # This is the pdf fitted to all available station precip data
            # (gamma dist). It will be sampled below.
        elif ptot_scenario2 == 'ptot2+':
            Ptot_pdf_cc_up2 = np.loadtxt(os.path.join(
                datapath, 'Ptot_pdf_cc_up2.csv'))
            # This is the pdf for season2 fitted to all available station
            # precip data (gamma dist). It will be sampled below.
        elif ptot_scenario2 == 'ptot2-':
            Ptot_pdf_cc_down2 = np.loadtxt(os.path.join(
                datapath, 'Ptot_pdf_cc_down2.csv'))
            # This is the pdf for season2 fitted to all available station
            # precip data (gamma dist). It will be sampled below.
        elif ptot_scenario2 == 'ptot2T+':
            Ptot_pdf2 = np.loadtxt(os.path.join(datapath, 'Ptot_pdf2.csv'))
            # This is the pdf fitted to all available station precip data
            # (gamma dist). It will be modified as a trend and sampled below.
            gamma_a2 = 1.65466
            # this is actually the shape parameter 'a' for a gamma distribution
            gamma_b2 = 52.4768
            # this is actually the shape parameter 'b' for a gamma distribution
        elif ptot_scenario2 == 'ptot2T-':
            Ptot_pdf2 = np.loadtxt(os.path.join(datapath, 'Ptot_pdf2.csv'))
            # This is the pdf fitted to all available station precip data
            # (gamma dist). It will be modified as a trend and sampled below.
            gamma_a2 = 1.65466
            gamma_b2 = 52.4768

    Recess_pdf = np.loadtxt(os.path.join(datapath, 'Recess_pdf.csv'))
    # This is the pdf of storm gradient recession coefficiencts from Morin et
    # al, 2005 (normal dist). It will be sampled below.
    # i.e., mu, sigma
    if mode == 'Validation':
        Easting = np.loadtxt(os.path.join(datapath, 'Easting.csv'))
        # This is the Longitudinal data for each gauge.
        Northing = np.loadtxt(os.path.join(datapath, 'Northing.csv'))
        # This is the Latitudinal data for each gauge. It will be sampled
        # below.
        gauges = np.loadtxt(os.path.join(datapath, 'gauges.csv'))
        # This is the list of gauge numbers. It will be sampled below.
        gauge_elev = np.loadtxt(os.path.join(datapath, 'gauge_elev.csv'))
        # This is the list of gauge numbers. It will be sampled below.
        Zz = gauge_elev
        numgauges = len(gauges)
    elif mode == 'Simulation':
        a = shp.Reader(os.path.join(datapath, 'boundary', 'boundary.shp'))
        # This is shapefile of the watershed boundary
        b = shp.Reader(os.path.join(
            datapath, 'point_elevations', 'point_elevations',
            'point_elevations.shp'))
        # This is shapefile of the elevations at each 'gauging' point on the
        # grid within the watershed boundary
        coords = a.shapeRecords()[0].shape.__geo_interface__['coordinates'][0]
        Yy = [Y for (X, Y) in coords]
        Xx = [X for (X, Y) in coords]
        polypts_xy = [[X, Y] for (X, Y) in coords]
        Zz = [recs.record[4] for recs in b.shapeRecords()]
        (X1, Y1) = np.meshgrid(
            np.linspace(min(Xx), max(Xx),
                        int(round((max(Xx)-min(Xx))/1000.))),
            np.linspace(min(Yy), max(Yy),
                        int(round((max(Yy)-min(Yy))/1000.))))
        # creates a mesh with 1000m spacings
        isin = Path(polypts_xy).contains_points(
            [(zx, zy) for (zx, zy) in zip(X1.flat, Y1.flat)]).reshape(X1.shape)
        # isin=inpolygon(X1(:),Y1(:),Xx,Yy)
        Yin = Y1[isin]
        Xin = X1[isin]
        # creates a list of points, Xin, Yin on the grid that are within the
        # watershed boundary and which will be used as 'gauging' output
        # This list is NOT grid structured.
        np.savetxt('Xin.txt', Xin)
        np.savetxt('Yin.txt', Yin)
        numgauges = len(Xin)
        # number of rain gauges in the basin.
        # NOTE: In this version this produces output on a grid, rather than at
        # real gauge locations.
    X = np.loadtxt(os.path.join(datapath, 'X.csv'))
    # This is the Longitudinal data for each grid point. It will be sampled
    # below to determine storm center location.
    Y = np.loadtxt(os.path.join(datapath, 'Y.csv'))
    # This is the Latitudinal data for each grid point. It will be sampled
    # below to determine storm center location.

    fuzz = np.loadtxt(os.path.join(datapath, 'fuzz.csv'))
    # This is a vector of fuzzy tolerace values for intensity selection.
    ET_monthly_day = np.loadtxt(os.path.join(datapath,
                                             'ET_monthly_day.csv'))
    # This is a matrix of averaged daytime values of ET (mm) grouped as one
    # column per month.
    ET_monthly_night = np.loadtxt(os.path.join(datapath,
                                               'ET_monthly_night.csv'))
    # This is a matrix of averaged nighttime values of ET (mm) grouped as one
    # column per month.

    Gauges = np.arange(numgauges, dtype=int)
    print(os.path.join(datapath, 'bndry_buffer', 'bndry_buffer.shp'))
    c = shp.Reader(os.path.join(datapath, 'bndry_buffer', 'bndry_buffer.shp'))
    # This is a shapefile of the watershed boundary buffer (5000 m)
    coordsc = c.shapeRecords()[0].shape.__geo_interface__['coordinates'][0]
    Yyy = [Yi for (Xi, Yi) in coordsc]
    Xxx = [Xi for (Xi, Yi) in coordsc]
    polypts_xyc = [[Xi, Yi] for (Xi, Yi) in coordsc]
    (Xx1, Yy1) = np.meshgrid(
        np.linspace(min(Xxx), max(Xxx), int(round((max(Xxx)-min(Xxx))/500.))),
        np.linspace(min(Yyy), max(Yyy), int(round((max(Yyy)-min(Yyy))/500.))))
    # creates a mesh with 500m spacings
    isinc = Path(polypts_xyc).contains_points(
        [(zx, zy) for (zx, zy) in zip(Xx1.flat, Yy1.flat)]).reshape(Xx1.shape)
    Yyin = Yy1[isinc]
    Xxin = Xx1[isinc]
    # creates a list of points, Xxin, Yyin on the grid that are within the
    # watershed boundary buffer and which will be used as storm center
    # locations

    # These are elevation ranges for the 3 orographic groups in Walnut Gulch,
    # AZ. They will need to be adapted for other basins (e.g., via hypsometric
    # curve).
    # The below is an abortive effort to generalise the orography away from WG:
    # Orogrp_list = []
    # for i in range(len(orographic_limits) + 1):
    #     if i == 0:
    #         oromin = int(round(min(Zz)))
    #         oromax = orographic_limits[0]
    #     elif i == len(orographic_limits):
    #         oromin = orographic_limits[-1]
    #         oromax = int(round(max(Zz)))
    #     else:
    #         oromin = orographic_limits[i-1]
    #         oromax = orographic_limits[i]
    #     Orogrp_list.append(np.arange(oromin, oromax, 1))
    OroGrp1 = np.arange(int(np.floor(min(Zz))), 1350+1, 1)
    OroGrp2 = np.arange(1351, 1500+1, 1)
    OroGrp3 = np.arange(1501, int(np.ceil(max(Zz))), 1)

    # lambda, kappa, and C are parameters of the intensity-duration curves
    # from WGEW of the form:
    # intensity1 = lambda*exp(delta*duration)+kappa*exp(gamma*duration)+C
    delta = -0.508
    gamma = -0.008
    # formerly lambda in the Matlab code
    lambda1 = [642.2, 578.0, 513.8, 449.5, 385.3, 321.1, 256.9, 192.7, 128.4,
               64.1, 21.]
    kappa = [93.1, 83.8, 74.5, 65.2, 55.9, 46.6, 37.2, 27.9, 18.6, 9.3, 0.9]
    C = [4.5, 4., 3.5, 3., 2.5, 2., 1.5, 1., 0.5, 0.25, 0.05]
    # constant value of intensity for each recession curve at high duration
    # values
    # season2
    # NOTE: C=0 for season2
    delta2 = -0.04
    gamma2 = -0.001038
    fracs = (1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
    lambda2 = [72.71*frac for frac in fracs]
    kappa2 = [5.432*frac for frac in fracs]
    # NOTE: These curves were obtained by finding the maximum intensity
    # value
    # for each duration, followed by LOWESS smoothing with a span of 0.1.
    # Next we fit a double
    # exponential form to this LOWESS fit and obtained the curve
    # parameters,lambda, kappa, and C (C=0 for season 2).
    # Finally, we obtain the other curves as percentiles of the maximum
    # fitted curve.

    max_numstorms = 10000
    # This is for initializing matrices. Trailing zeros are deleted from
    # the matrix at the end of the code.

    ####################################################################
    # Everything above this line is model input. Model operation begins below
    # this line

    # Create named and timestamped model output folder for the simulation set
    t1 = date.datetime.now()
    month = t1.strftime("%B")[:3]
    t2 = str(t1.day)+month+str(t1.year)+'_'+str(t1.hour)+str(t1.minute)
    # for naming output directories and files by current date/time
    tx = (mode[:3] + str(numsims) + '_' + str(numsimyrs) + 'Y_' +
          str(seasons) + 'S_' + ptot_scenario + '_' + storminess_scenario +
          '_' + str(ptot_scenario2) + '_' + str(storminess_scenario2) + '_' +
          ET_scenario + '_')
    tx0 = os.path.join(cwd, 'model_output', tx) + t2
    namecounter = 1
    try:
        os.mkdir(tx0)
    except OSError:
        tx0 += '_' + str(t1.second)  # adds differentiation if needed.
        os.mkdir(tx0)

    for simulations in range(numsims):
        # number of simulations to run. A separate timestamped output folder
        # is created for each
        t0 = time.time()
        t1a = date.datetime.now()
        montha = t1a.strftime("%B")[:3]
        t2a = str(t1a.day)+montha+str(t1a.year)+'_'+str(t1a.hour)+str(
            t1a.minute)
        tx2 = os.path.join(tx0, t2a)
        try:
            os.mkdir(tx2)
        except OSError:
            tx2 += '_' + str(t1a.second)
            os.mkdir(tx2)
            # pass  # dir already exists

        # next section of code simulates ET on a 12-hour basis (average day
        # and night values for each month of each simulation year) based on
        # sampling a pdf composed from regionally measured data

        jan_index = np.logical_not(np.isnan(ET_monthly_day[:, 0]))
        feb_index = np.logical_not(np.isnan(ET_monthly_day[:, 1]))
        mar_index = np.logical_not(np.isnan(ET_monthly_day[:, 2]))
        apr_index = np.logical_not(np.isnan(ET_monthly_day[:, 3]))
        may_index = np.logical_not(np.isnan(ET_monthly_day[:, 4]))
        jun_index = np.logical_not(np.isnan(ET_monthly_day[:, 5]))
        jul_index = np.logical_not(np.isnan(ET_monthly_day[:, 6]))
        aug_index = np.logical_not(np.isnan(ET_monthly_day[:, 7]))
        sep_index = np.logical_not(np.isnan(ET_monthly_day[:, 8]))
        oct_index = np.logical_not(np.isnan(ET_monthly_day[:, 9]))
        nov_index = np.logical_not(np.isnan(ET_monthly_day[:, 10]))
        dec_index = np.logical_not(np.isnan(ET_monthly_day[:, 11]))
        ET_final = [np.array([np.nan, ])]
        for yooo in range(numsimyrs):
            yoo = str(yooo)
            ET_day_jan = np.random.choice(ET_monthly_day[jan_index, 0], 31)
            ET_night_jan = np.random.choice(ET_monthly_night[jan_index, 0], 31)
            ET_day_feb = np.random.choice(ET_monthly_day[feb_index, 1], 28)
            ET_night_feb = np.random.choice(ET_monthly_night[feb_index, 1], 28)
            ET_day_mar = np.random.choice(ET_monthly_day[mar_index, 2], 31)
            ET_night_mar = np.random.choice(ET_monthly_night[mar_index, 2], 31)
            ET_day_apr = np.random.choice(ET_monthly_day[apr_index, 3], 30)
            ET_night_apr = np.random.choice(ET_monthly_night[apr_index, 3], 30)
            ET_day_may = np.random.choice(ET_monthly_day[may_index, 4], 31)
            ET_night_may = np.random.choice(ET_monthly_night[may_index, 4], 31)
            ET_day_jun = np.random.choice(ET_monthly_day[jun_index, 5], 30)
            ET_night_jun = np.random.choice(ET_monthly_night[jun_index, 5], 30)
            ET_day_jul = np.random.choice(ET_monthly_day[jul_index, 6], 31)
            ET_night_jul = np.random.choice(ET_monthly_night[jul_index, 6], 31)
            ET_day_aug = np.random.choice(ET_monthly_day[aug_index, 7], 31)
            ET_night_aug = np.random.choice(ET_monthly_night[aug_index, 7], 31)
            ET_day_sep = np.random.choice(ET_monthly_day[sep_index, 8], 30)
            ET_night_sep = np.random.choice(ET_monthly_night[sep_index, 8], 30)
            ET_day_oct = np.random.choice(ET_monthly_day[oct_index, 9], 31)
            ET_night_oct = np.random.choice(ET_monthly_night[oct_index, 9], 31)
            ET_day_nov = np.random.choice(ET_monthly_day[nov_index, 10], 30)
            ET_night_nov = np.random.choice(ET_monthly_night[nov_index, 10],
                                            30)
            ET_day_dec = np.random.choice(ET_monthly_day[dec_index, 11], 31)
            ET_night_dec = np.random.choice(ET_monthly_night[dec_index, 11],
                                            31)
            jan_interleave = np.empty((2*ET_day_jan.size,),
                                      dtype=ET_day_jan.dtype)
            jan_interleave[0::2] = ET_day_jan
            jan_interleave[1::2] = ET_night_jan
            feb_interleave = np.empty((2*ET_day_feb.size,),
                                      dtype=ET_day_feb.dtype)
            feb_interleave[0::2] = ET_day_feb
            feb_interleave[1::2] = ET_night_feb
            mar_interleave = np.empty((2*ET_day_mar.size,),
                                      dtype=ET_day_mar.dtype)
            mar_interleave[0::2] = ET_day_mar
            mar_interleave[1::2] = ET_night_mar
            apr_interleave = np.empty((2*ET_day_apr.size,),
                                      dtype=ET_day_apr.dtype)
            apr_interleave[0::2] = ET_day_apr
            apr_interleave[1::2] = ET_night_apr
            may_interleave = np.empty((2*ET_day_may.size,),
                                      dtype=ET_day_may.dtype)
            may_interleave[0::2] = ET_day_may
            may_interleave[1::2] = ET_night_may
            jun_interleave = np.empty((2*ET_day_jun.size,),
                                      dtype=ET_day_jun.dtype)
            jun_interleave[0::2] = ET_day_jun
            jun_interleave[1::2] = ET_night_jun
            jul_interleave = np.empty((2*ET_day_jul.size,),
                                      dtype=ET_day_jul.dtype)
            jul_interleave[0::2] = ET_day_jul
            jul_interleave[1::2] = ET_night_jul
            aug_interleave = np.empty((2*ET_day_aug.size,),
                                      dtype=ET_day_aug.dtype)
            aug_interleave[0::2] = ET_day_aug
            aug_interleave[1::2] = ET_night_aug
            sep_interleave = np.empty((2*ET_day_sep.size,),
                                      dtype=ET_day_sep.dtype)
            sep_interleave[0::2] = ET_day_sep
            sep_interleave[1::2] = ET_night_sep
            oct_interleave = np.empty((2*ET_day_oct.size,),
                                      dtype=ET_day_oct.dtype)
            oct_interleave[0::2] = ET_day_oct
            oct_interleave[1::2] = ET_night_oct
            nov_interleave = np.empty((2*ET_day_nov.size,),
                                      dtype=ET_day_nov.dtype)
            nov_interleave[0::2] = ET_day_nov
            nov_interleave[1::2] = ET_night_nov
            dec_interleave = np.empty((2*ET_day_dec.size,),
                                      dtype=ET_day_dec.dtype)
            dec_interleave[0::2] = ET_day_dec
            dec_interleave[1::2] = ET_night_dec

            cat_chain = np.concatenate(
                [may_interleave, jun_interleave, jul_interleave,
                 aug_interleave, sep_interleave, oct_interleave,
                 nov_interleave, dec_interleave, jan_interleave,
                 feb_interleave, mar_interleave, apr_interleave])
            ET_final.append(cat_chain.copy())

        ET_final = np.concatenate(ET_final)
        ET_matrix = ET_final[1:]
        # This matrix contains a complete time series of day/night ET values
        # for the length of the simulation.
        if ET_scenario == 'ET+':
            ET_matrix += ET_matrix*ET_scaling_factor
        if ET_scenario == 'ET-':
            ET_matrix -= ET_matrix*ET_scaling_factor

        # Initialize output matrices
        # Storm matrix = [Storm# StormArea StormDuration Int_DurCurve#
        #                 StormIntensity #GaugesHit RecessionVal StormTotal
        #                 UTM_Longitude UTM_Latitude Year SimTime]
        # The values in Storm_matrix have the following corresponding units:
        # [# m^2 min # mm/hr mm/hr/km mm m m y hr]

        Storm_matrix = np.zeros((max_numstorms*numsimyrs, 12))
        # ^based on presumed maximum number of storms per year (see above)
        Ptot_ann_global = np.zeros(numsimyrs)
        Ptot_ann_global2 = np.zeros(numsimyrs)

        # Gauge matrix = [Year Storm# StormIntensity StormDuration StormTotal
        #                 Ann_Cum_PTotal Int_Arr_Time SimTime]
        # The values in Gauge_matrix have the following corresponding units:
        # [y # mm/hr min mm mm hr hr]

        Gauge_matrix = []
        for ii in range(numgauges):
            # init matrices of storm output for each gauge
            # in the Python iteration, we're going to store this as a list
            Gauge_matrix.append(np.zeros((max_numstorms*numsimyrs, 8)))

        # initialize all variables (concatenated matrices of generated output)
        Intensity_local_all = [0., ]
        Storm_totals_all = [0., ]
        Duration_local_all = [0., ]
        Gauges_hit_all = [0., ]
        Ptotal_local = []
        Storm_total_local_year_across_yrs = []
        storm_count = 0
        master_storm_count = 0
        storm_trend = 0
        storm_trend2 = 0
        calendar_time = 0.  # tracks simulation time per year in hours

        for syear in range(numsimyrs):
            calendar_year_time = 0.  # tracks simulation time per year in hours
            leftovers = 0.
            Storm_total_local_year_across_yrs.append([])
            # add a new storage slot
            storm_trend += storminess_scaling_factor
            Ptotal = 0.  # ok
            if ptot_scenario == 'ptotT+':
                mu += mu*ptot_scaling_factor
                # scale the mean of the PTotal distribution each year by the
                # scaling factor defined above and add it to the prior mean
                # value.
            elif ptot_scenario == 'ptotT-':
                mu -= mu*ptot_scaling_factor
                # scale the mean of the PTotal distribution each year by the
                # scaling factor defined above and add it to the prior mean
                # value.
            Ptot_pdf_norm = {'sigma': sigma, 'mu': mu}
            Ptot_ann_global[syear] = np.random.normal(
                loc=Ptot_pdf_norm['mu'], scale=Ptot_pdf_norm['sigma'])
            # samples from the appropriate normal distribution and saves global
            # value of Ptot (that must be equalled or exceeded) for each year

            Storm_total_local_year = []
            ann_cum_Ptot_gauge = np.zeros(numgauges)
            for storm in range(max_numstorms):
                int_arr_val = genextreme.rvs(c=Int_arr_pdf_GEV['shape'],
                                             loc=Int_arr_pdf_GEV['mu'],
                                             scale=Int_arr_pdf_GEV['sigma'])
                try:
                    int_arr_val = int_arr_val.clip(
                        Int_arr_pdf_GEV['trunc_interval'][0],
                        Int_arr_pdf_GEV['trunc_interval'][1])
                except KeyError:
                    pass
                # Samples from distribution of interarrival times (hr). This
                # can be used to develop STORM output for use in rainfall-
                # runoff models or any water balance application.
                rain_int_gauge = np.zeros(numgauges)
                center_val_X = np.random.choice(Xxin)
                # sample uniformly from storm center matrix from grid with
                # 500 m spacings within a 5000 m buffer around the basin.
                center_val_Y = np.random.choice(Yyin)
                North = center_val_Y
                East = center_val_X
                area_val = genextreme.rvs(c=Area_pdf_EV['shape'],
                                          loc=Area_pdf_EV['mu'],
                                          scale=Area_pdf_EV['sigma'])
                try:
                    area_val = area_val.clip(
                        Area_pdf_EV['trunc_interval'][0],
                        Area_pdf_EV['trunc_interval'][1])
                except KeyError:
                    pass
                # Samples from distribution of storm areas in m2
                cx = East
                # value of coord should be set to storm center selected (same
                # below)
                cy = North
                r = np.sqrt(area_val/np.pi)
                # value here should be selected based on area above in meters
                # to match the UTM values in North and East vectors.
                if mode == 'Validation':
                    mask_name = (
                        (Easting[:]-cx)**2 + (Northing[:]-cy)**2) <= r**2
                    # determine which gauges are hit by Euclidean geometry
                elif mode == 'Simulation':
                    mask_name = ((Xin[:]-cx)**2 + (Yin[:]-cy)**2) <= r**2
                    # determine which grid locations are hit by Euclidean
                    # geometry
                if np.sum(mask_name) == 0:
                    # this short circuits the storm loop in the case that the
                    # storm does not affect any 'gauging' location
                    Storm_total_local_year.append(np.zeros(numgauges))
                    Storm_total_local_year_across_yrs[-1].append(
                        Storm_total_local_year[-1].copy())
                    continue
                storm_count += 1
                master_storm_count += 1
                gauges_hit = Gauges[mask_name]
                num_gauges_hit = len(gauges_hit)
                Storm_matrix[master_storm_count, 0] = master_storm_count
                Storm_matrix[master_storm_count, 1] = area_val
                Storm_matrix[master_storm_count, 8] = cx
                Storm_matrix[master_storm_count, 9] = cy
                Storm_matrix[master_storm_count, 10] = syear
                # this routine below allows for orography in precip by first
                # determining the closest gauge and then determining its
                # orographic grouping
                if mode == 'Validation':
                    gdist = (Easting[:]-cx)**2 + (Northing[:]-cy)**2
                elif mode == 'Simulation':
                    gdist = (Xin[:]-cx)**2 + (Yin[:]-cy)**2
                closest_gauge = round(Zz[gdist.argmin()])
                # this will be compared against orographic gauge groupings to
                # determine the appropriate set of intensity-duration curves
                Storm_matrix[master_storm_count, 5] = num_gauges_hit
                Gauges_hit_all.extend(gauges_hit)

                # this routine below determines to which orographic group the
                # closest gauge to the storm center belongs to, and censors the
                # number of curves accordingly missing top curve in GR1, top
                # and bottom curves for GR2, and bottom curve for GR3
                if closest_gauge in OroGrp1:
                    # new version of orography compares local 'gauge' elevation
                    # to elevation bands called OroGrp, defined above
                    baa = 'a'
                elif closest_gauge in OroGrp2:
                    baa = 'b'
                elif closest_gauge in OroGrp3:
                    baa = 'c'
                int_dur_curve_num = np.arange(numcurves, dtype=int)
                duration_val = genextreme.rvs(c=Duration_pdf_GEV['shape'],
                                              loc=Duration_pdf_GEV['mu'],
                                              scale=Duration_pdf_GEV['sigma'])
                try:
                    duration_val = duration_val.clip(
                        Duration_pdf_GEV['trunc_interval'][0],
                        Duration_pdf_GEV['trunc_interval'][1])
                except KeyError:
                    pass
                duration_val = round(duration_val)
                # round to nearest minute for consistency with measured data
                Storm_matrix[master_storm_count, 2] = duration_val
                calendar_year_time += int_arr_val + duration_val/60.  # hours
                calendar_time += int_arr_val + duration_val/60.  # hours
                Storm_matrix[master_storm_count, 11] = calendar_time
                # stored cumulative time of simulation

                # original curve# probs for 30%-20%-10%:
                # [0.0636 0.0727 0.0819 0.0909 0.0909 0.0909 0.0909 0.0909
                # 0.1001 0.1090 0.1182]
                # original curve# probs are modified as below
                if baa == 'a':
                    weights = [0.0318, 0.0759, 0.0851, 0.0941, 0.0941, 0.0941,
                               0.0941, 0.0941, 0.1033, 0.1121, 0.1213]
                    # weights to reflect probabilities that favor selection of
                    # lower curves.
                elif baa == 'b':
                    weights = [0.0478, 0.0778, 0.0869, 0.0959, 0.0959, 0.0959,
                               0.0959, 0.0959, 0.1051, 0.1141, 0.0888]
                elif baa == 'c':
                    weights = [0.0696, 0.0786, 0.0878, 0.0968, 0.0968, 0.0968,
                               0.0968, 0.0968, 0.1060, 0.1149, 0.0591]
                int_dur_curve_val = np.random.choice(
                    int_dur_curve_num, p=weights)
                # each of these is adapted based on its orographic grouping

                Storm_matrix[master_storm_count, 3] = int_dur_curve_val
                intensity_val = (
                    lambda1[int_dur_curve_val]*np.exp(delta*duration_val) +
                    kappa[int_dur_curve_val]*np.exp(gamma*duration_val) +
                    C[int_dur_curve_val])
                # these curves are based on empirical data from WG.
                fuzz_int_val = np.random.choice(fuzz)
                intensity_val2 = intensity_val + fuzz_int_val
                # allowing for +/-5 mm/hr fuzzy tolerance around selected
                # intensity
                if intensity_val2 < 0.5:  # cannot have zero or -ve intensity
                    intensity_val = 0
                elif intensity_val2 < 1. and intensity_val >= 0.5:
                    intensity_val = 1.
                else:
                    intensity_val = intensity_val2
                if storminess_scenario == 'stormsT+':
                    intensity_val = round(
                        intensity_val + intensity_val*storm_trend)
                    # storminess trend is applied and its effect rises each yr
                    # of simulation
                elif storminess_scenario == 'stormsT-':
                    intensity_val = round(
                        intensity_val - intensity_val*storm_trend)
                elif storminess_scenario == 'storms+':
                    intensity_val = round(
                        intensity_val + intensity_val*storm_stepchange)
                    # storminess change is applied as a step change uniformly
                    # over the simulation
                elif storminess_scenario == 'storms-':
                    intensity_val = round(
                        intensity_val - intensity_val*storm_stepchange)
                Storm_matrix[master_storm_count, 4] = intensity_val

                # area to determine which gauges are hit
                recess_val = np.random.normal(
                    loc=Recess_pdf[0], scale=Recess_pdf[1])
                # this pdf of recession coefficients determines how intensity
                # declines with distance from storm center (see below)
                Storm_matrix[master_storm_count, 6] = recess_val
                for j in range(numgauges):
                    # determine cartesian distances to all hit gauges and
                    # associated intensity values at each gauge hit by storm
                    if j in gauges_hit:
                        if mode == 'Validation':
                            gauge_dist = np.sqrt(((Easting[j]-cx)**2 +
                                                  (Northing[j]-cy)**2))
                        elif mode == 'Simulation':
                            gauge_dist = np.sqrt(((Xin[j]-cx)**2 +
                                                  (Yin[j]-cy)**2))
                        gauge_dist_km = gauge_dist/1000.  # change to km
                        rain_int_gauge[j] = intensity_val*np.exp(
                            -2.*(recess_val**2)*gauge_dist_km**2)
                        # Rodriguez-Iturbe et al., 1986; Morin et al., 2005 but
                        # sampled from a distribution
                    else:
                        rain_int_gauge[j] = 0.
                        # give zero intensity to gauges not hit
                    ann_cum_Ptot_gauge[j] += rain_int_gauge[j]*duration_val/60.

                    Gauge_matrix[j][master_storm_count, 0] = syear  # year
                    Gauge_matrix[j][master_storm_count, 1] = master_storm_count
                    # ^storm number
                    Gauge_matrix[j][master_storm_count, 2] = rain_int_gauge[j]
                    Gauge_matrix[j][master_storm_count, 3] = duration_val
                    Gauge_matrix[j][master_storm_count, 4] = (
                        rain_int_gauge[j] * duration_val/60.)
                    # ^storm total
                    Gauge_matrix[j][master_storm_count, 5] = (
                        ann_cum_Ptot_gauge[j])
                    # ^annual cumulative total (Ptot)
                    Gauge_matrix[j][master_storm_count, 6] = int_arr_val
                    # ^interarrival time in hrs
                    Gauge_matrix[j][master_storm_count, 7] = calendar_time
                    # ^total sim time in hours

                Intensity_local_all.extend(rain_int_gauge)
                # ^collect into vector of all simulated intensities (at all
                # gauges)
                dur_step = np.full(numgauges, duration_val)
                Duration_local_all.extend(dur_step)
                Storm_total_local_year.append(rain_int_gauge*duration_val/60.)
                # collect storm total data for all gauges into rows by storm
                Storm_totals_all.extend(rain_int_gauge*duration_val/60.)
                Storm_matrix[master_storm_count, 7] = (
                    intensity_val*duration_val/60.)
                Ptotal = np.zeros(numgauges)

                # This section permits optional visualisation of each
                # storm...
                if PLOT_EACH_STORM:
                    if ~np.allclose(rain_int_gauge, 0.):
                        if mode == 'Simulation':
                            plot(Xin, Yin, 'r.', alpha=0.3)
                            scatter(Xin, Yin, s=rain_int_gauge)
                        else:
                            plot(Easting, Northing, 'r.', alpha=0.3)
                            scatter(Easting, Northing, s=rain_int_gauge)
                        plot(East, North, 'ko')
                        axis('equal')
                        show()

                for stormtot in Storm_total_local_year:
                    Ptotal[:] += stormtot
                    # sum the annual storm total at each gauge
                if np.nanmedian(Ptotal) > Ptot_ann_global[syear]:
                    # once the median of all gauges exceeds the selected annual
                    # storm total, a new simulation year begins
                    Ptotal_local.append(Ptotal.copy())
                    break  # end storm loop and start a new simulation year
                Storm_total_local_year_across_yrs[-1].append(
                    Storm_total_local_year[-1].copy())
                # collect all local annual storm totals for each gauge.
                # end of the storm loop
            stormcount = storm  # + 1
            ##################################################################
            # SEASON 2
            if seasons == 2:
                storm_trend2 += storminess_scaling_factor2
                if ptot_scenario2 == 'ptot2C':
                    Ptot_ann_global2[syear] = np.random.gamma(
                        Ptot_pdf2[0], scale=Ptot_pdf2[1])
                    # samples from the control normal distribution and saves
                    # global value of Ptot (that must be equalled or exceeded)
                    # for each year
                if ptot_scenario2 == 'ptot2+':
                    Ptot_ann_global2[syear] = np.random.gamma(
                        Ptot_pdf_cc_up2[0], scale=Ptot_pdf_cc_up2[1])
                    # samples from the increased wetness normal distribution
                    # and saves global value of Ptot (that must be equalled or
                    # exceeded) for each year
                if ptot_scenario2 == 'ptot2-':
                    Ptot_ann_global2[syear] = np.random.gamma(
                        Ptot_pdf_cc_down2[0], scale=Ptot_pdf_cc_down2[1])
                    # samples from the decreased wetness normal distribution
                    # and saves global value of Ptot (that must be equalled or
                    # exceeded) for each year
                if ptot_scenario2 == 'ptot2T+':
                    gamma_a2 += ptot_scaling_factor2
                    # create new normal pdf with the updated mean scaled by the
                    # climate change trend
                    Ptot_ann_global2[syear] = np.random.gamma(
                        gamma_a2, scale=gamma_b2)
                    # samples from normal distribution and saves global value
                    # of Ptot (that must be equalled or exceeded) for each year
                elif ptot_scenario == 'ptot2T-':
                    gamma_a2 -= ptot_scaling_factor2
                    # create new normal pdf with the updated mean scaled by the
                    # climate change trend
                    Ptot_ann_global2[syear] = np.random.gamma(
                        gamma_a2, scale=gamma_b2)
                    # samples from normal distribution and saves global value
                    # of Ptot (that must be equalled or exceeded) for each year

                # add a new storage slot
                for storm in range(stormcount, max_numstorms):
                    int_arr_val = genextreme.rvs(
                        c=Int_arr_pdf_GEV2['shape'],
                        loc=Int_arr_pdf_GEV2['mu'],
                        scale=Int_arr_pdf_GEV2['sigma'])
                    try:
                        int_arr_val = int_arr_val.clip(
                            Int_arr_pdf_GEV2['trunc_interval'][0],
                            Int_arr_pdf_GEV2['trunc_interval'][1])
                    except KeyError:
                        pass
                    # Samples from distribution of interarrival times (hr).
                    # This can be used to develop STORM output for use in
                    # rainfall-runoff models or any water balance application.
                    rain_int_gauge = np.zeros(numgauges)
                    center_val_X = np.random.choice(Xxin)
                    # sample uniformly from storm center matrix from grid with
                    # 500 m spacings within a 5000 m buffer around the basin.
                    center_val_Y = np.random.choice(Yyin)
                    North = center_val_Y
                    East = center_val_X
                    area_val = genextreme.rvs(c=Area_pdf_EV2['shape'],
                                              loc=Area_pdf_EV2['mu'],
                                              scale=Area_pdf_EV2['sigma'])
                    try:
                        area_val = area_val.clip(
                            Area_pdf_EV2['trunc_interval'][0],
                            Area_pdf_EV2['trunc_interval'][1])
                    except KeyError:
                        pass
                    # Samples from distribution of storm areas in m2
                    cx = East
                    # value of coord should be set to storm center selected
                    # (same below)
                    cy = North
                    r = np.sqrt(area_val/np.pi)
                    # value here should be selected based on area above in
                    # meters to match the UTM values in North and East vectors.
                    if mode == 'Validation':
                        mask_name = (
                            (Easting[:]-cx)**2 + (Northing[:]-cy)**2) <= r**2
                        # determine which gauges are hit by Euclidean geometry
                    elif mode == 'Simulation':
                        mask_name = ((Xin[:]-cx)**2 + (Yin[:]-cy)**2) <= r**2
                        # determine which grid locations are hit by Euclidean
                        # geometry
                    if np.sum(mask_name) == 0:
                        # this short circuits the storm loop in the case that
                        # the storm does not affect any 'gauging' location
                        Storm_total_local_year.append(np.zeros(numgauges))
                        Storm_total_local_year_across_yrs[-1].append(
                            Storm_total_local_year[-1].copy())
                        continue
                    storm_count += 1
                    master_storm_count += 1
                    gauges_hit = Gauges[mask_name]
                    num_gauges_hit = len(gauges_hit)
                    Storm_matrix[master_storm_count, 0] = master_storm_count
                    Storm_matrix[master_storm_count, 1] = area_val
                    Storm_matrix[master_storm_count, 8] = cx
                    Storm_matrix[master_storm_count, 9] = cy
                    Storm_matrix[master_storm_count, 10] = syear
                    # this routine below allows for orography in precip by
                    # first determining the closest gauge and then determining
                    # its orographic grouping
                    if mode == 'Validation':
                        gdist = (Easting[:]-cx)**2 + (Northing[:]-cy)**2
                    elif mode == 'Simulation':
                        gdist = (Xin[:]-cx)**2 + (Yin[:]-cy)**2
                    closest_gauge = int(round(Zz[gdist.argmin()]))
                    # this will be compared against orographic gauge groupings
                    # to determine the appropriate set of intensity-duration
                    # curves
                    Storm_matrix[master_storm_count, 5] = num_gauges_hit
                    Gauges_hit_all.extend(gauges_hit)

                    # this routine below determines to which orographic group
                    # the closest gauge to the storm center belongs to, and
                    # censors the number of curves accordingly missing top
                    # curve in GR1, top and bottom curves for GR2, and bottom
                    # curve for GR3
                    if closest_gauge in OroGrp1:
                        # new version of orography compares local 'gauge'
                        # elevation to elevation bands called OroGrp, defined
                        # above
                        baa = 'a'
                    elif closest_gauge in OroGrp2:
                        baa = 'b'
                    elif closest_gauge in OroGrp3:
                        baa = 'c'
                    int_dur_curve_num = np.arange(numcurves, dtype=int)
                    duration_val = fisk.rvs(
                        c=Duration_pdf_Fisk2['c'],
                        scale=Duration_pdf_Fisk2['scale'])
                    try:
                        duration_val = duration_val.clip(
                            Duration_pdf_Fisk2['trunc_interval'][0],
                            Duration_pdf_Fisk2['trunc_interval'][1])
                    except KeyError:
                        pass
                    duration_val = round(duration_val)
                    # round to nearest minute for consistency with measured
                    # data
                    Storm_matrix[master_storm_count, 2] = duration_val
                    calendar_year_time += int_arr_val + duration_val/60.  # hrs
                    calendar_time += int_arr_val + duration_val/60.  # hours
                    Storm_matrix[master_storm_count, 11] = calendar_time
                    # stored cumulative time of simulation

                    # original curve# probs for 30%-20%-10%:
                    # [0.0636 0.0727 0.0819 0.0909 0.0909 0.0909 0.0909 0.0909
                    # 0.1001 0.1090 0.1182]
                    # original curve# probs are modified as below
                    if baa == 'a':
                        weights = [0.0318, 0.0759, 0.0851, 0.0941, 0.0941,
                                   0.0941, 0.0941, 0.0941, 0.1033, 0.1121,
                                   0.1213]
                        # weights to reflect probabilities that favor selection
                        # of lower curves.
                    elif baa == 'b':
                        weights = [0.0478, 0.0778, 0.0869, 0.0959, 0.0959,
                                   0.0959, 0.0959, 0.0959, 0.1051, 0.1141,
                                   0.0888]
                    elif baa == 'c':
                        weights = [0.0696, 0.0786, 0.0878, 0.0968, 0.0968,
                                   0.0968, 0.0968, 0.0968, 0.1060, 0.1149,
                                   0.0591]
                    int_dur_curve_val = np.random.choice(int_dur_curve_num,
                                                         p=weights)
                    # each of these is adapted based on its orographic grouping

                    Storm_matrix[master_storm_count, 3] = int_dur_curve_val
                    intensity_val = (
                        lambda2[int_dur_curve_val]*np.exp(delta2 *
                                                          duration_val) +
                        kappa2[int_dur_curve_val]*np.exp(gamma2 *
                                                         duration_val) +
                        C[int_dur_curve_val])
                    # No C this time
                    # these curves are based on empirical data from WG.
                    fuzz_int_val = np.random.choice(fuzz)
                    intensity_val2 = intensity_val + fuzz_int_val
                    # allowing for +/-5 mm/hr fuzzy tolerance around selected
                    # intensity
                    if intensity_val2 < 0.5:
                        # cannot have zero or negative intensity
                        intensity_val = 0
                    elif intensity_val2 < 1. and intensity_val >= 0.5:
                        intensity_val = 1.
                    else:
                        intensity_val = intensity_val2
                    if storminess_scenario == 'storms2T+':
                        intensity_val = round(
                            intensity_val + intensity_val*storm_trend2)
                        # storminess trend is applied and its effect rises each
                        # year of simulation
                    elif storminess_scenario == 'storms2T-':
                        intensity_val = round(
                            intensity_val - intensity_val*storm_trend2)
                    elif storminess_scenario == 'storms2+':
                        intensity_val = round(
                            intensity_val + intensity_val*storm_stepchange2)
                        # storminess change is applied as a step change
                        # uniformly over the simulation
                    elif storminess_scenario == 'storms2-':
                        intensity_val = round(
                            intensity_val - intensity_val*storm_stepchange2)
                    Storm_matrix[master_storm_count, 4] = intensity_val

                    # area to determine which gauges are hit
                    recess_val = np.random.normal(
                        loc=Recess_pdf[0], scale=Recess_pdf[1])
                    # this pdf of recession coefficients determines how
                    # intensity declines with distance from storm center (see
                    # below)
                    Storm_matrix[master_storm_count, 6] = recess_val
                    for j in range(numgauges):
                        # determine cartesian distances to all hit gauges and
                        # associated intensity values at each gauge hit by the
                        # storm
                        if j in gauges_hit:
                            if mode == 'Validation':
                                gauge_dist = np.sqrt(((Easting[j]-cx)**2 +
                                                      (Northing[j]-cy)**2))
                            elif mode == 'Simulation':
                                gauge_dist = np.sqrt(((Xin[j]-cx)**2 +
                                                      (Yin[j]-cy)**2))
                            gauge_dist_km = gauge_dist/1000.  # change to km
                            rain_int_gauge[j] = intensity_val*np.exp(
                                -2.*(recess_val**2)*gauge_dist_km**2)
                            # Rodriguez-Iturbe et al., 1986; Morin et al., 2005
                            # but sampled from a distribution
                        else:
                            rain_int_gauge[j] = 0.
                            # give zero intensity to gauges not hit
                        ann_cum_Ptot_gauge[j] += (
                            rain_int_gauge[j]*duration_val/60.)

                        Gauge_matrix[j][master_storm_count, 0] = syear  # year
                        Gauge_matrix[j][
                            master_storm_count, 1] = master_storm_count
                        # ^storm number
                        Gauge_matrix[j][
                            master_storm_count, 2] = rain_int_gauge[j]
                        Gauge_matrix[j][master_storm_count, 3] = duration_val
                        Gauge_matrix[j][master_storm_count, 4] = (
                            rain_int_gauge[j] * duration_val/60.)
                        # ^storm total
                        Gauge_matrix[j][
                            master_storm_count, 5] = ann_cum_Ptot_gauge[j]
                        # ^annual cumulative total (Ptot)
                        Gauge_matrix[j][master_storm_count, 6] = int_arr_val
                        # ^interarrival time in hrs
                        Gauge_matrix[j][master_storm_count, 7] = calendar_time
                        # ^total sim time in hours

                    Intensity_local_all.extend(rain_int_gauge)
                    # ^collect into vector of all simulated intensities (at all
                    # gauges)
                    dur_step = np.full(numgauges, duration_val)
                    Duration_local_all.extend(dur_step)
                    Storm_total_local_year.append(
                        rain_int_gauge*duration_val/60.)
                    # collect storm totals for all gauges into rows by storm
                    Storm_totals_all.extend(rain_int_gauge*duration_val/60.)
                    Storm_matrix[master_storm_count, 7] = (
                        intensity_val*duration_val/60.)
                    Ptotal = np.zeros(numgauges)

                    # This section permits optional visualisation of each
                    # storm...
                    if PLOT_EACH_STORM:
                        if ~np.allclose(rain_int_gauge, 0.):
                            if mode == 'Simulation':
                                plot(Xin, Yin, 'r.', alpha=0.3)
                                scatter(Xin, Yin, s=rain_int_gauge)
                            else:
                                plot(Easting, Northing, 'r.', alpha=0.3)
                                scatter(Easting, Northing, s=rain_int_gauge)
                            plot(East, North, 'ko')
                            axis('equal')
                            show()

                    for stormtot in Storm_total_local_year:
                        Ptotal[:] += stormtot
                        # sum the annual storm total at each gauge

                    if np.nanmedian(Ptotal) > Ptot_ann_global2[syear]:
                        # once the median of all gauges exceeds the selected
                        # annual storm total, a new simulation year begins
                        Ptotal_local.append(Ptotal.copy())
                        break  # end storm loop and start a new simulation year
                    Storm_total_local_year_across_yrs[-1].append(
                        Storm_total_local_year[-1].copy())
                    # collect all local annual storm totals for each gauge.
                    # end of the storm loop
                    # ...& end of the second season loop
                    # now tab back to the simulations loop, forever ago

        trunc_at = np.argmax(np.isclose(Storm_matrix[1:, 8], 0.)) + 1
        Storm_matrix_out = Storm_matrix[:trunc_at, :]
        # ^gets rid of trailing zeros from initialized matrix
        # Turn Gauge_matrix from a list of arrays to a 3D array.
        # Dims are numgauges, master_storm_ID, property_ID
        Gauge_matrix_out = np.empty((len(Gauge_matrix),
                                     Gauge_matrix[0].shape[0],
                                     Gauge_matrix[0].shape[1]), dtype=float)
        rw = 0
        for GM in Gauge_matrix:
            Gauge_matrix_out[rw, :, :] = GM
            rw += 1
        # use the master_storm_count property row (1) to trim
        AB = np.argmax(np.isclose(Gauge_matrix_out[0, 1:, 1], 0.)) + 1
        Gauge_matrix_out = Gauge_matrix_out[:, :AB, :]

        leftstuff2 = np.amax(Gauge_matrix_out[0, :, 7])
        # ^simulated time in hours
        nn = Gauge_matrix_out.shape[1]
        leftovers = 365.25 * numsimyrs * 24. - leftstuff2
        # remaining annual time in hours (not occupied by storms or interstorm
        # periods)
        leftovers22 = leftovers / nn
        # adjustment to add to each interarrival time in hours
        Gauge_matrix_1back = Gauge_matrix_out[0, :, :]
        Gauge_matrix_out[0, :, 6] += leftovers22
        cal_time2 = 0.
        for jjj in range(nn):
            Gauge_matrix_out[0, jjj, 7] = (
                cal_time2 + Gauge_matrix_out[0, jjj, 6] +
                Gauge_matrix_out[0, jjj, 3]/60.)
            cal_time2 = Gauge_matrix_out[0, jjj, 7]
        Gauge_matrix_out[:, :, 6] = Gauge_matrix_out[0, :, 6]
        Gauge_matrix_out[:, :, 7] = Gauge_matrix_out[0, :, 7]

        Storm_matrix_out[:, 11] = Gauge_matrix_out[0, :, 7]
        Gauges_hit_all[Gauges_hit_all == 0] = np.nan
        GZ = np.logical_not(np.isclose(Intensity_local_all, 0.))
        Intensity_all = np.array(Intensity_local_all)[GZ]
        Duration_all = np.array(Duration_local_all)[GZ]
        Storm_totals_all_out = np.array(Storm_totals_all)[GZ]
        Ptot_ann_global_out = Ptot_ann_global[1:]
        if seasons == 2:
            Ptot_ann_global_out2 = Ptot_ann_global2[1:]
        Gauges_hit_all_out = np.array(Gauges_hit_all)[1:]

    # save all the output data:
        savename = [
            'Ptot_ann', 'Storm_matrix', 'Gauges_hit', 'Storm_totals',
            'Intensity', 'Duration', 'ET_matrix']
        savemod = [
            'global_', '', 'all_', 'all_', 'selected_', 'selected_',
            '']
        savedata = [
            Ptot_ann_global_out, Storm_matrix_out,
            Gauges_hit_all_out, Storm_totals_all_out, Intensity_all,
            Duration_all, ET_matrix]
        if seasons == 2:
            savename.append('Ptot_ann2')
            savemod.append('global2_')
            savedata.append(Ptot_ann_global_out2)
        for nm, mod, dat in zip(savename, savemod, savedata):
            np.savetxt(os.path.join(tx2, (
                nm + '_' + str(numsims) + 'sims_' + str(numsimyrs) + 'y_' +
                mod + t2)), dat)

        for kk in range(numgauges):
            np.savetxt(os.path.join(tx2, (
                'Gauge_matrix' + str(kk) + '_' + str(numsims) + 'sims_' +
                str(numsimyrs) + 'y_' + t2)), Gauge_matrix_out[kk, :, :])

        boo = time.time() - t0
        runtime_seconds = boo
        runtime_minutes = boo/60.

########
    return Storm_matrix_out, Gauge_matrix, Duration_all

if __name__ == "__main__":
    # storm(mode='Validation', numsims=1, numsimyrs=2, seasons=1,
    #       ptot_scenario='ptotC', storminess_scenario='stormsC',
    #       ptot_scenario2=None, storminess_scenario2=None, ET_scenario='ETC')

    # # a slightly longer, more complex exercise of the code:
    # storm(mode='Simulation', numsims=2, numsimyrs=20, seasons=2,
    #       ptot_scenario='ptot+', storminess_scenario='stormsT-',
    #       ptot_scenario2='ptot2T-', storminess_scenario2='storms2+',
    #       ET_scenario='ET+')

    SM, GM, Dur = storm(mode='Validation', numsims=1, numsimyrs=1, seasons=2,
          ptot_scenario='ptotC', storminess_scenario='stormsC',
          # ptot_scenario2=None, storminess_scenario2=None,
          ptot_scenario2='ptot2C', storminess_scenario2='storms2C',
          ET_scenario='ETC')
    print('Total number of storms:', SM.shape[0])
