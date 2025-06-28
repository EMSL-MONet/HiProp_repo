# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 20:47:03 2025

@author: chak803
"""
import numpy as np
from rosetta import Rosetta
rose33 = Rosetta(rosetta_version=3, model_code=3)


def soil_ptf(sand_pct, clay_pct, bulkd):
    """
    Predicts soil properties using the Rosetta Soil texture triangle.

    Parameters:
    - sand_pct (float): Percentage of sand in the soil.
    - clay_pct (float): Percentage of clay in the soil.
    - bulkd (float): Bulk density of the soil in kg/m^3.

    Returns:
    - swilt (float): Wilting point (residual volumetric water content) in saturation unit
    - phi_max (float): Saturated volumetric water content
    - kstex (float): Saturated hydraulic conductivity in m/hour.
    - m (float): Fraction of soil water content corresponding to matrix potential.

    References:
    The function uses the Rosetta Soil model for soil texture prediction.
    For more information, refer to: https://github.com/usda-ars-ussl/rosetta-soil?tab=readme-ov-file
    """
    silt_pct = 100 - (sand_pct+clay_pct)
    texture_data = np.array([[sand_pct, silt_pct, clay_pct, bulkd/1000]])
    means, stdev = rose33.predict(texture_data)

    theta_r, theta_s, log10_a, log10_n, log10_ksat = means.T

    n = (10**log10_n[0])
    m = 1 - 1/n
    swilt = theta_r[0]/theta_s[0]
    phi_max = theta_s[0]
    kstex = 10**log10_ksat[0] * 0.01 / 24  # m/hour
    # assuming field capacity is mid point between wilting and saturation moisture conditions
    sfc = ((theta_r[0]+theta_s[0])/2)/theta_s[0]
    return swilt, phi_max, sfc, kstex, m


def simulate_seasonal_rainfall(total_time, mean_rate, mean_rainfall):
    """
    Simulate rainfall over time. The scale parameter in the exponential distribution represents the
    average amount of rainfall per event and average interval between consecutive rainfall events.

    Parameters:
        mean_rate (float): Event rate (events per unit time).
        total_time (int): Total time period.
        mean_rainfall (float): Mean rainfall amount per event.

    Returns:
        rainfall (ndarray): Simulated rainfall over time.
    """
    # Define seasonal patterns for event rate (events per unit time) and mean rainfall amount per event
    spring_rate, summer_rate, autumn_rate, winter_rate = mean_rate
    spring_mean_rainfall, summer_mean_rainfall, autumn_mean_rainfall, winter_mean_rainfall = mean_rainfall

    # Initialize rainfall array
    rainfall = np.zeros(total_time)

    # Simulate rainfall events over time
    for t in range(total_time):
        # Determine the current season based on time
        month = (t // (24 * 30)) % 12
        if 3 <= month < 6:  # Spring (March, April, May)
            rate = spring_rate
            mean_rainfall = spring_mean_rainfall
        elif 6 <= month < 9:  # Summer (June, July, August)
            rate = summer_rate
            mean_rainfall = summer_mean_rainfall
        elif 9 <= month < 12:  # Autumn (September, October, November)
            rate = autumn_rate
            mean_rainfall = autumn_mean_rainfall
        else:  # Winter (December, January, February)
            rate = winter_rate
            mean_rainfall = winter_mean_rainfall

        # Generate rainfall events according to the current event rate and mean rainfall
        if np.random.rand() < rate:
            # Sample the time until the next event from an exponential distribution
            rng = np.random.default_rng(seed=42)  # Use a fixed seed for reproducibility
            time_until_next_event = rng.exponential(scale=1/rate)

            # Ensure the event occurs within the total time period
            if t + time_until_next_event < total_time:
                # Sample the amount of rain for the event from an exponential distribution
                rainfall[t + int(time_until_next_event)] += rng.exponential(scale=mean_rainfall)

    return rainfall


# @jit(nopython=True)
def soil_wat_bal(theta, dt,  rain, soil_dict):
    """
    Computes the soil water balance for one time step based on the specified parameters.
  
    Parameters
    ----------
    s : float
        Current soil moisture level.
    dt : float
        Time step for the computation (in hours).
    rain : float
        Rainfall rate (in meters per hour).
    AGG : float
        Aggregation factor affecting soil permeability.
    soil_dict : dict containing soil parameters in the following order:
            - param_bulkd: bulk density of the soil (kg/m^3)
            - sand_pct: percentage of sand content in the soil
            - silt_pct: percentage of silt content in the soil
            - clay_pct: percentage of clay content in the soil
            - phi_max: maximum soil moisture retention capacity (m^3/m^3)
            - porosity: soil moisture at the permanent wilting point (m^3/m^3)
            - alpha : shape parameter for Aggregate-hydraulic conductivity relationship
            - beta : shape parameter for Aggregate-hydraulic conductivity relationship
            - sh: soil saturation at hygroscopic point (m^3/m^3)
            - sstar: soil saturation where ET decreases with soil moisture (m^3/m^3)
            - et_max: maximum evapotranspiration rate [m/day]
            - ew: evaporation rate [m/day]
            - zr: rooting depth (m)
            - rstar: fraction of rain intercepted by canopy
            - kstex: textural saturated hydraulic conductivity of the soil (m/s)
            - m_prime: parameter affecting soil water retention
            - m: parameter affecting soil water retention
            - pH: soil pH
            - swilt: soil saturation at wilting point (m^3/m^3)

    Returns
    -------
    dict
        A dictionary containing the computed values for the following parameters:
            - 'Rain [m/day]': Rainfall rate
            - 'Runoff [m/day]': Rate of runoff
            - 'saturation': Current soil moisture level after computation
            - 'ET [m/day]': Evapotranspiration rate
    
    """
    # Computes evapotranspiration
    def evapotranspiration(s, soil_dict):
        """
        Computes evapotranspiration
        """
        if s <= soil_dict['swilt']:
            et = 0
        elif soil_dict['swilt'] < s <= soil_dict['phitex']:
            et = soil_dict['ew'] + (soil_dict['et_max'] - soil_dict['ew']) * (s - soil_dict['swilt']) / \
                (soil_dict['phitex'] - soil_dict['swilt'])
        else:
            et = soil_dict['et_max']
        return et

    # Soil water storage capacity
    swc = soil_dict['zr']  # m

    # Canopy interception
    tf = rain - rain * soil_dict['rstar']  # Througfall   # m/day

    # Add througfall
    theta = theta + tf / swc   # Dimensionless
    # Verify if there is runoff
    if theta > soil_dict['phitex']:
        q = (theta - soil_dict['phitex']) * swc
        theta = soil_dict['phitex']
        Kh = soil_dict['kstex']
    else:
        q = 0
        # Calculate Effective Saturation (Se)
        theta_r = soil_dict['swilt']
        theta_s = soil_dict['phitex']
        Se = (theta - theta_r) / (theta_s - theta_r)
        Se = np.clip(Se, 0, 1)  # Ensure Se is within [0,1]

        # Calculate Hydraulic Conductivity (Kh) using Mualem-van Genuchten
        # Calculate Hydraulic Conductivity (Kh) using Mualem-van Genuchten
        # Avoid division by zero when Se=0
        m = soil_dict['m']
        if Se > 0:
            Kh = soil_dict['kstex'] * (Se ** 0.5) * (1 - (1 - Se ** (1 / m)) ** m) ** 2  # m/day
        else:
            Kh = 0
    et = evapotranspiration(theta, soil_dict)
    theta = theta - (et+Kh) * dt / swc
    theta = max(theta, soil_dict['swilt'])  # Ensure s doesn't go below sh

    # s = np.arange(0, 1.01, 0.01)
    # ET = np.array([evapotranspiration(s, soil_dict) for s in s])
    # plt.figure()
    # plt.plot(s, ET)

    output = np.array([theta, rain, q, et, Kh])
    return output
