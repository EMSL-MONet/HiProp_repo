import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
# import numpy as np
# import seaborn as sns

def generate_daily_soil_temperature(MAT, Trange, days=365, northern_hemisphere=True):
    # Generate sequence of values from 0 to pi
    PIseq = np.linspace(0, np.pi, days)
    
    # Adjust for hemisphere: positive for northern, negative for southern
    coef = 1.5 if northern_hemisphere else -1.5
    
    # Calculate soil temperature for each day
    soil_temps = (Trange / 2) * np.sin((2 * PIseq) - coef) + MAT
    
    return soil_temps


def calculate_annual_to_daily_ct(seqDAY, peakDAY, sdNPP, annNPP):
    """
    Calculate daily carbon inputs distribution using a normal distribution.

    Parameters:
    seqDAY (ndarray): Array of julian day integers (1 to 365).
    peakDAY (int): Julian day when carbon inputs peak.
    sdNPP (float): Standard deviation of the normal distribution around peakDAY.
    annNPP (float): Annual Net Primary Productivity in gC m^-2 yr^-1.

    Returns:
    ndarray: Array of daily total carbon inputs.
    """
    daily_prob = norm.pdf(seqDAY, peakDAY, sdNPP)
    daily_prob /= daily_prob.sum()
    daily_ct = daily_prob * annNPP
    return daily_ct

def calculate_aboveground_ct(daily_ct, RtoSi, aHARVj):
    """
    Calculate daily aboveground carbon input.

    Parameters:
    daily_ct (ndarray): Daily total carbon inputs.
    RtoSi (float): Root-to-shoot ratio of the material.
    aHARVj (float): Fraction of aboveground material harvested.

    Returns:
    ndarray: Array of daily aboveground carbon inputs.
    """
    return daily_ct * (1 / (RtoSi + 1)) * (1 - aHARVj)

def calculate_belowground_ct(daily_ct, RtoSi, bHARVj, depth, Rdep50, Rdepmax):
    """
    Calculate daily belowground carbon input.

    Parameters:
    daily_ct (ndarray): Daily total carbon inputs.
    RtoSi (float): Root-to-shoot ratio of the material.
    bHARVj (float): Fraction of belowground material harvested.
    depth (float): Depth of soil layer simulated (cm).
    Rdep50 (float): Depth at which 50% of root biomass is proportioned (cm).
    Rdepmax (float): Maximum rooting depth (cm).

    Returns:
    ndarray: Array of daily belowground carbon inputs.
    """
    bCT = daily_ct * (RtoSi / (RtoSi + 1)) * (1 - bHARVj)
    proportion_factor = (depth * (Rdep50 + Rdepmax)) / (Rdepmax * (Rdep50 + depth))
    return bCT * proportion_factor

def add_noise(carbon_input, noise_level=0.05):
    """
    Add noise to simulate inter-day variation.

    Parameters:
    carbon_input (ndarray): Daily carbon input data.
    noise_level (float): Maximum amplitude of the noise as a fraction of carbon_input.

    Returns:
    ndarray: Noisy carbon input data.
    """
    noise = np.random.normal(0, noise_level * carbon_input, size=carbon_input.shape)
    return carbon_input + noise



def generate_forcings(generate_plot=False):
    daily_mean_rainfall = pd.read_csv('daily_mean_rainfall.csv')
    daily_mean_rainfall=daily_mean_rainfall[daily_mean_rainfall['year'] == 2020].reset_index()
    daily_rainfall = daily_mean_rainfall['secPrecipBulk'].values*0.001  # m/d
# =============================================================================
#     daily_soil_temperatures
# =============================================================================

    # Example parameters
    num_pt = len(daily_rainfall)
    mean_annual_temperature = 10  # mean annual temperature in degrees Celsius
    temperature_range = 18        # temperature range in degrees Celsius

    # Generate soil temperature
    daily_soil_temperatures = generate_daily_soil_temperature(mean_annual_temperature, temperature_range, days=num_pt)

    # Optionally, plot the results
    if generate_plot:
        fig, ax1 = plt.subplots(figsize=(7, 4))

        # Plot daily soil temperatures
        ax1.plot(daily_soil_temperatures, label="Daily Soil Temperature", color='tab:blue')
        ax1.set_xlabel("Day of the Year")
        ax1.set_ylabel("Soil Temperature (°C)", color='tab:blue')
        ax1.tick_params(axis='y', labelcolor='tab:blue')

        # Create a twin y-axis to plot rainfall
        ax2 = ax1.twinx()
        ax2.bar(range(len(daily_rainfall)), daily_rainfall, label="Daily Rainfall", color='tab:gray', alpha=0.5)
        ax2.set_ylabel("Rainfall (m/day)", color='tab:gray')
        ax2.tick_params(axis='y', labelcolor='tab:gray')

        # Add grid and legend
        fig.tight_layout()
        ax1.legend(loc="upper left")
        ax2.legend(loc="upper right")
        
        # plt.plot(daily_temp_h05_v503['soilTempMean'], label="Obs Daily Soil Temperature")

        # plt.title("Simulated Daily Soil Temperature Over a Year")
        # plt.legend()
        plt.grid(True)
        plt.savefig('figs/daily_soil_temperatures.png', dpi=300, bbox_inches='tight')
        plt.show()
# =============================================================================
#   Calculate daily carbon inputs
# =============================================================================
    # Example parameters for simulation
    annNPP = 450  # Annual NPP in gC m^-2 yr^-1
    peakNPP_DAY = 182  # Julian day of peak NPP (mid-year)
    peakLitter_DAY = 200  # Julian day of peak litterfall; offset to simulate phase difference
    sdNPP = 30     # Standard deviation for the NPP distribution spread
    RtoSi = 0.5    # Ratio of root biomass to shoot biomass
    aHARV = 0.1    # Proportion of aboveground biomass harvested
    bHARV = 0.1    # Proportion of belowground biomass harvested
    depth = 15     # Soil depth for carbon input calculation in cm
    Rdep50 = 30    # Depth at which 50% of roots are distributed in cm
    Rdepmax = 100  # Maximum rooting depth in cm

    # Generate a sequence of Julian days (1 through 365)
    seqDAY = np.arange(1, num_pt+1)

    # Calculate daily carbon inputs
    daily_ct_npp = calculate_annual_to_daily_ct(seqDAY, peakNPP_DAY, sdNPP, annNPP)
    daily_ct_litter = calculate_annual_to_daily_ct(seqDAY, peakLitter_DAY, sdNPP, annNPP)
    jaCT = calculate_aboveground_ct(daily_ct_litter, RtoSi, aHARV)  # Use litterfall timing for aboveground carbon
    jbCT = calculate_belowground_ct(daily_ct_npp, RtoSi, bHARV, depth, Rdep50, Rdepmax)  # Use NPP timing for belowground

    # Add noise to simulate inter-day variation
    jaCT_noisy = add_noise(jaCT, noise_level=0.05)
    jbCT_noisy = add_noise(jbCT, noise_level=0.05)
    total_plant_input= jaCT_noisy+jbCT_noisy

    total_plant_input= jaCT_noisy+jbCT_noisy
    if generate_plot:
        # # Plot NPP alongside noisy carbon input results over the year
        plt.figure(figsize=(10, 5))
        # plt.plot(seqDAY, daily_ct_npp, label="Daily NPP", color='tab:blue')
        # plt.plot(seqDAY, jaCT_noisy, label="Noisy Aboveground Carbon Input (jaCT)", color='tab:green')
        # plt.plot(seqDAY, jbCT_noisy, label="Noisy Belowground Carbon Input (jbCT)", linestyle='dashed', color='tab:red')
        plt.plot(seqDAY, total_plant_input, label="Total Plant Carbon Input (jbCT)", linestyle='dashed', color='black')
        plt.xlabel("Day of the Year")
        plt.ylabel(r"Plant carbon Input (gC m$^{-2}$ day$^{-1}$)")
        # plt.title("Simulated Daily Carbon Inputs with Noise and Phase Offset Over a Year")
        # plt.legend()
        plt.grid(True)
    
    return daily_rainfall, daily_soil_temperatures, total_plant_input