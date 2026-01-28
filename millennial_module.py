# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 21:19:32 2025

@author: chak803
"""
from soil_hydrology_module import  soil_wat_bal

from scipy.integrate import solve_ivp
import pandas as pd
import numpy as np
from scipy.interpolate import pchip_interpolate
import numdifftools as nd
from scipy import optimize
import matplotlib.pyplot as plt
import MEMES_functions as ms
import importlib
importlib.reload(ms)

def load_par():
    soil_texture = pd.read_excel('soil_texture.xlsx')
    soil_texture = soil_texture[soil_texture['source'].isin(['Intact', 'Repacked','HiHydroSoil v2.0'])].copy()

    # Read parameters
    parameters_file = pd.read_table("soilpara_in_fit.txt", delimiter=' ', header=None)
    parameters_file.columns = ["parameter", "value"]
    parameters = dict(zip(parameters_file.iloc[:, 0], parameters_file.iloc[:, 1]))
    
    
    parameters['pH'] = 7  # pH
    parameters['bulkd'] = 1350  # bulk density in kg soil m-3
    parameters['param_pc'] = 0.86  # slope of mineral C - clay relationship
    parameters['depth'] = 0.3  # core sample depth in meter
    
    
    kinetic_dict = {
        'param_pb': parameters['param_pb'],
        'param_pa': parameters['param_pa'],
        'param_pi': parameters['param_pi'],
        'alpha_pl': parameters['alpha_pl'],
        'eact_pl': parameters['eact_pl'],
        'kaff_pl': parameters['kaff_pl'],
        'kaff_des': parameters['kaff_des'],
        'rate_pa': parameters['rate_pa'],
        'rate_break': parameters['rate_break'],
        'rate_leach': parameters['rate_leach'],
        'rate_ma': parameters['rate_ma'],
        'rate_bd': parameters['rate_bd'],
        'alpha_lb': parameters['alpha_lb'],
        'eact_lb': parameters['eact_lb'],
        'kaff_lb': parameters['kaff_lb'],
        'cue_t': parameters['cue_t'],
        'tae_ref': parameters['tae_ref'],
        'param_p1': parameters['param_p1'],
        'param_p2': parameters['param_p2'],
        'cue_ref': parameters['cue_ref']
    }
    
    
    soil_dict = {
        'bulkd': parameters['bulkd'],  # kg soil/ m3
        'sand_pct': 50,
        'silt_pct': 30,
        'clay_pct': 20,
        'param_pc': 0.86,
        'lambda_val': parameters['lambda'],
        'matpot': parameters['matpot'],
        'kamin': parameters['kamin'],
        'depth': parameters['depth'],
        'phitex': 0.6,
        'et_max': 0.005,  # m/d
        'ew': 0.0005,  # m/d
        'zr': parameters['depth'],  # rooting depth
        'rstar': 0.1,
        'kstex': 1,  # m/d
        'm': 0.286,
        'pH': parameters['pH'],
        'swilt': 0.08,
        'coloc_factor':2.8,
        'rel_opt_theta':0.65,
        'b': 0.75
    }

    soil_dict.keys()
    return soil_dict, kinetic_dict, soil_texture

def Yan2018(theta,phi, clay,coloc_factor, rel_opt_theta,b):
    ktheta = 0.1
    # def compute_a(clay, coloc_factor):
    #     if clay < 0.046 * 100 / coloc_factor:
    #         return 0
    #     elif clay <= 104.6/2.8:
    #         return (coloc_factor * clay / 100) - 0.046
    #     else:
    #         return 1
    def compute_a(clay, coloc_factor):
        if clay <= (1+0.046)*100/coloc_factor:
            return (coloc_factor * clay / 100) - 0.046
        else:
            return 1
    theta_opt =rel_opt_theta*phi 
    ns=2
    
    if theta<theta_opt:
        fm= ((ktheta+theta_opt)/(ktheta+theta))*(theta/theta_opt)**(1+compute_a(clay,coloc_factor)*ns)
    else:
        fm = ((phi - theta)/(phi-theta_opt))**b
    return fm


def cal_derivatives(t, state, npp, st, sw, kinetic_dict, soil_dict, SS, default_moisture=True):
    """
    Function to compute the derivatives of state variables in a soil carbon model.
    Parameters:
    - t: Current time.
    - state: Array containing the state variables (POM, LMWC, AGG, MIC, MAOM) or (POM, LMWC, AGG, MIC, MAOM, CO2)
      depending on SS.
    - npp: Net Primary Production data.
    - st: Soil temperature data.
    - rain: rain in m per hour
    - efficiency: Calculates the carbon use efficiency
    - kinetic_dict: Array containing kinetic parameters for decomposition.
    - soil_dict: Array containing soil parameters.
    - SS: Boolean indicating whether the derivatives are returned for steady-state (True) or not (False).
    - dynamic_phi: Boolean indicating whether to use dynamic phi (True) or not (False).
    Returns:
    - Array containing the derivatives of state variables.
    """
    if SS:
        POM, LMWC, AGG, MIC, MAOM = state
    else:
        POM, LMWC, AGG, MIC, MAOM, CO2 = state
    # Extracting forcing_var data
    temperature = np.interp(t, np.arange(0, len(st)), st)
    plant_input = np.interp(t, np.arange(0, len(npp)), npp)
    theta = np.interp(t, np.arange(0, len(sw)), sw)

    porosity = soil_dict['phitex']

    param_claysilt = soil_dict['silt_pct']+soil_dict['clay_pct']

    # Soil type properties
    kaff_lm = np.exp(-kinetic_dict['param_p1'] * soil_dict['pH'] - kinetic_dict['param_p2']) * kinetic_dict['kaff_des']
    param_qmax = soil_dict['depth'] * soil_dict['bulkd'] * param_claysilt * soil_dict['param_pc']

    if default_moisture:
        scalar_wd = (theta / porosity)**(0.5)
        spot = np.exp(soil_dict['lambda_val'] * -soil_dict['matpot'])
        so2 = (soil_dict['kamin'] + (1 - soil_dict['kamin']) * ((porosity - theta) / porosity)**0.5)
        scalar_wb = spot * so2 * scalar_wd
    else:
        coloc_a = soil_dict['coloc_factor']*soil_dict['clay_pct']/100 - 0.046 if soil_dict['clay_pct'] < 100/soil_dict['coloc_factor'] else 1  # from fig 6 of Yan 2018 Nat Comm
        scalar_wd = (theta / porosity)**(0.5*coloc_a)
        scalar_wb = Yan2018(theta, porosity, soil_dict['clay_pct'], soil_dict['coloc_factor'],
                            soil_dict['rel_opt_theta'],soil_dict['b'])
    

    # Decomposition rates
    gas_const = 8.31446
    vmax_pl = kinetic_dict['alpha_pl'] * np.exp(-kinetic_dict['eact_pl'] / (gas_const * (temperature + 273.15)))
    f_PO_LM = vmax_pl * scalar_wd * POM * MIC / (kinetic_dict['kaff_pl'] + MIC)
    f_PO_AG = kinetic_dict['rate_pa'] * scalar_wd * POM
    f_AG_break = kinetic_dict['rate_break'] * scalar_wd * AGG

    f_LM_leach = kinetic_dict['rate_leach'] * scalar_wd * LMWC

    f_LM_MA = scalar_wd * kaff_lm * LMWC * (1 - MAOM / param_qmax)
    f_MA_LM = kinetic_dict['kaff_des'] * MAOM / param_qmax
    vmax_lb = kinetic_dict['alpha_lb'] * np.exp(-kinetic_dict['eact_lb'] / (gas_const * (temperature + 273.15)))
    f_LM_MB = vmax_lb* scalar_wb * MIC * LMWC / (kinetic_dict['kaff_lb'] + LMWC)
    f_MB_turn = kinetic_dict['rate_bd'] * MIC**2.0
    f_MA_AG = kinetic_dict['rate_ma'] * scalar_wd * MAOM
    f_MB_atm = f_LM_MB * (1 - (kinetic_dict['cue_ref'] - kinetic_dict['cue_t'] * (temperature - kinetic_dict['tae_ref'])))

    # mass balance
    dPOM = plant_input * kinetic_dict['param_pi'] + f_AG_break * kinetic_dict['param_pa'] - f_PO_AG - f_PO_LM
    dLMWC = plant_input * (1. - kinetic_dict['param_pi']) - f_LM_leach + f_PO_LM - f_LM_MA - \
        f_LM_MB + f_MB_turn * (1. - kinetic_dict['param_pb']) + f_MA_LM
    dAGG = f_MA_AG + f_PO_AG - f_AG_break
    dMIC = f_LM_MB - f_MB_turn - f_MB_atm
    dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * kinetic_dict['param_pb'] - f_MA_AG + f_AG_break * (1. - kinetic_dict['param_pa'])
    fluxes = {
        'f_AG_break': f_AG_break,
        'f_POM_to_AGG': f_PO_AG,
        'f_POM_to_DOM': f_PO_LM,
        'f_MAOM_to_DOM': f_MA_LM,
        'f_DOM_to_MAOM': f_LM_MA,
        'f_DOM_leach': f_LM_leach,
        'f_DOM_to_MIC': f_LM_MB,
        'f_MB_turn':f_MB_turn,
        'f_MB_atm':f_MB_atm
    }
    if SS:
        return np.array([dPOM, dLMWC, dAGG, dMIC, dMAOM])
    else:
        return np.array([dPOM, dLMWC, dAGG, dMIC, dMAOM, f_MB_atm]), fluxes


def cal_steady_state(forcing_var, kinetic_dict, soil_dict, calc_eigen,default_moisture=True):
    """
    Calculate the steady state of a system using Newton-Krylov method and perform stability analysis if requested.

    Parameters:
    - derivs (callable): Function that computes the derivatives of the system.
    - forcing_var (dict): Dictionary containing forcing_var functions for temperature, moisture, and NPP.
    - totalK (callable):Compute the total hydraulic conductivity.
    - kinetic_dict (dict): Dictionary containing kinetic parameters.
    - soil_dict (dict): Dictionary containing soil parameters.
    - SS (bool, optional): Flag indicating whether to perform steady-state calculation. Defaults to True.
    - dynamic_phi (bool, optional): Flag indicating whether to use dynamic phi. Defaults to False.
    - calc_eigen (bool, optional): Flag indicating whether to perform stability analysis. Defaults to True.
    Returns:
    - sol (array): Steady state solution of the system.
    - eig (array): Eigenvalues of the Jacobian matrix if stability analysis is performed.
    Notes:
    - The function internally defines another function, `SS_derivatives`, which computes the derivatives for
      steady-state calculation.
    - It uses the Newton-Krylov method to find the steady state.
    - If stability analysis is requested (`calc_eigen=True`), it computes the Jacobian matrix and its eigenvalues.
    """

   
    num_years = 1
    forcing_var_time = np.arange(1, num_years * 365 + 1)
    st = np.interp(forcing_var_time, forcing_var_time, np.tile(np.mean(forcing_var['forc_st']), len(forcing_var_time)))
    sw =  np.interp(forcing_var_time, forcing_var_time, np.tile(np.mean(forcing_var['forc_sw']), len(forcing_var_time)))
    npp = np.interp(forcing_var_time, forcing_var_time, np.tile(np.mean(forcing_var['forc_npp']), len(forcing_var_time)))
    
    t = 10
    def fun(state):
        dCdt = cal_derivatives(t, state,npp, st, sw, kinetic_dict, soil_dict,  SS=True,default_moisture=default_moisture)
        return dCdt
    
    # def fun(y): return maas_bal_fun(t, y, npp, st, sw, kinetic_dict, soil_dict,SS=True)

    sol = optimize.newton_krylov(fun, np.array([500, 5, 500, 10, 1000]), maxiter=2000)

    # stability analysis
    if calc_eigen:
        jac_fun = nd.Jacobian(fun)
        jac = jac_fun(sol)
        eig, vec = np.linalg.eig(jac)
        return sol, eig
    else:
        return sol

def dynamic_simulation(forcing_var, kinetic_dict, soil_dict,
              state={'POM': 1, 'LMWC': 1, 'AGG': 1, 'MIC': 1, 'MAOM': 1, 'CO2': 0},
              SS=False, default_moisture=True):
    """
    Function to run a soil carbon model simulation.
    Parameters:
    - forcing_var: Dictionary containing forcing_var data (temperature, soil water content, net primary production).
    - derivs: Function that computes the derivatives of state variables.
    - kinetic_dict: Array containing kinetic parameters for decomposition.
    - soil_dict: Array containing soil parameters.
    - num_years: Number of years to run the simulation for (default is 200).
    - state: Initial state of the model (default is equal initial values for POM, LMWC, AGG, MIC, MAOM, and zero CO2).
    - use_method: Method to use for solving the differential equations ('solve_ivp' or 'integration_step').
    - SS: Boolean indicating whether the derivatives are returned for steady-state (True) or not (False).
    - dynamic_phi: Boolean indicating whether to use dynamic phi (True) or not (False).
    Returns:
    - Output of the model simulation.
    """

    init_Cpool = np.array(list(state.values()))
    
    num_steps = len(forcing_var)


    forcing_var_time = np.arange(0, num_steps, 1)

    st_dt = pchip_interpolate(np.arange(len(forcing_var['forc_st'])), forcing_var['forc_st'], forcing_var_time)

    npp_dt = pchip_interpolate(np.arange(len(forcing_var['forc_npp'])), forcing_var['forc_npp'], forcing_var_time)

    sw_dt = pchip_interpolate(np.arange(len(forcing_var['forc_sw'])), forcing_var['forc_sw'], forcing_var_time)
    
    def maas_bal_fun(t, state,npp_dt, st_dt, sw_dt, kinetic_dict, soil_dict, SS, default_moisture):
        dCdt, _ = cal_derivatives(t, state,npp_dt, st_dt, sw_dt, kinetic_dict, soil_dict, SS, default_moisture)
        return dCdt
        
    t_span = (0, num_steps) # simulation period
    sol = solve_ivp(
                    maas_bal_fun, 
                    t_span, 
                    init_Cpool, 
                    args=(npp_dt, st_dt, sw_dt, kinetic_dict, soil_dict, SS,default_moisture), 
                    t_eval=np.arange(0, num_steps)
                )

    y = np.hstack((sol.t.reshape(-1, 1), sol.y.T))

    millennialC = pd.DataFrame(y, columns=["time [day]","POM [gC m-2]", "LMWC [gC m-2]",
                               "AGG [gC m-2]", "MIC [gC m-2]", "MAOM [gC m-2]", "CO2 [gC m-2]"])


    t = millennialC['time [day]']
    fluxes_all = {}

    for i in range(len(t)):

        state = millennialC.loc[i,['POM [gC m-2]', 'LMWC [gC m-2]', 'AGG [gC m-2]',
               'MIC [gC m-2]', 'MAOM [gC m-2]', 'CO2 [gC m-2]']].values

        # dCdt, fluxes = cal_derivatives(t[i], state, forcing_var, millennial_params_val, SS=False)
        # Compute fluxes using the derivative function
        _, fluxes = cal_derivatives(t[i], state,npp_dt, st_dt, sw_dt, kinetic_dict, soil_dict, SS=False, default_moisture=default_moisture)

        # Store fluxes in a structured format
        for key, value in fluxes.items():
            if key not in fluxes_all:
                fluxes_all[key] = []
            fluxes_all[key].append(value)

    for key in fluxes_all:
        fluxes_all[key] = np.array(fluxes_all[key])
    fluxes_df= pd.DataFrame(fluxes_all)
    millennialC=pd.concat((millennialC,fluxes_df), axis=1)    
    return millennialC, fluxes_df





def plot_pools(millennialC):
    # plot C pools--------------
    pools = ['POM [gC m-2]', 'LMWC [gC m-2]', 'AGG [gC m-2]',
       'MIC [gC m-2]', 'MAOM [gC m-2]']
    millennialC['time [year]']=millennialC['time [day]']/365
    fig1 = plt.figure(figsize=(10, 8))
    plt.subplot(321)
    
    for col in pools:
        plt.plot(millennialC['time [year]'], millennialC[col], label=col, linewidth=2)
    plt.ylabel(r'C stocks (gC  m$^{-2}$)')
    # plt.yscale('log')
    plt.xlabel("Time [years]") 
    plt.legend(loc='best', fontsize=8)
    # bbox_to_anchor=(1, 1))
    # plt.ylim(0, max(output[:, 1:].max(axis=0)))

    plt.subplot(322)
    plt.plot(millennialC['time [year]'], millennialC['CO2 [gC m-2]'], linewidth=2)
    plt.ylabel(r"CO$_2$ [gC $m^{-2}$]")
    plt.xlabel("Time [years]") 


    def f(x):
        return np.interp(x, millennialC['time [day]'], millennialC['CO2 [gC m-2]'])
    df = nd.Derivative(f, order=2)
    resp = df(millennialC['time [day]'].values)
    plt.subplot(323)
    plt.plot(millennialC['time [year]'], resp, '-', linewidth=1)
    plt.ylabel(r"Respiration rate [gC m$^{-2}$ d$^{-1}$]")
    plt.xlabel("Time [years]") 

    plt.tight_layout()
    
    # fluxes= [ 'f_AG_break','f_POM_to_AGG', 'f_POM_to_DOM', 'f_MAOM_to_DOM', 
    #          'f_DOM_to_MAOM',  'f_DOM_leach', 'f_DOM_to_MIC','f_MB_turn']
    # fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(10, 5), sharex=False)
    # axes=axes.flatten()
    # for ax, col in zip(axes, fluxes):
    #     ax.plot(millennialC['time [day]'], millennialC[col], alpha=0.8)
    #     ax.set_title(f"{col}\n[gC m$^{{-2}}$ d$^{{-1}}$]")
    #     ax.grid(True)

    # plt.tight_layout()

    return fig1

# %%


def steady_state_bar_chart(new_coloc=2.8,default_moisture=True, make_plot=False):
    soil_dict, kinetic_dict, soil_texture=load_par()
    rainfall, daily_soil_temperatures, total_plant_input = ms.generate_forcings()

    # Initialize an empty list to collect results
    results = []
    for _, row in soil_texture.iterrows():
        # Extract soil parameters for the current row
        if  row['source']=='Repacked':
            coloc_factor= new_coloc# Intact soils would have more isolation, higher coloc_factor more isolation
        else:
            coloc_factor=soil_dict['coloc_factor']  #Repacked soils would have less isolation
        soil_dict_temp = soil_dict.copy()

        soil_dict_temp.update({
            'sand_pct': row['Sand'],
            'silt_pct': row['Silt'],
            'clay_pct': row['Clay'],
            'phitex': row['theta_s'],  # Assuming porosity corresponds to theta_s
            'swilt': row['theta_r'],   # Assuming wilting point corresponds to theta_r
            'kstex': row['KSAT_cm_day'] * 1e-2,  # Convert from cm/day to m/day if needed
            'm': row['m'],
            'coloc_factor': coloc_factor,
            'rel_opt_theta':row['Relative_opt_water_content'],
            'b': row['b']
        })


        SWB = np.zeros((len(rainfall), 5))
        s=0.2
        for i in range(len(rainfall)):
            result = soil_wat_bal(s, 1,  rainfall[i], soil_dict_temp)
            s = result[0]
            SWB[i] = result

        SWB_df = pd.DataFrame(SWB)
        SWB_df.columns = ['theta', 'Rain (m/d)', 'Runoff (m/d)', 'ET (m/d)', 'Kh (m/d)']
        SWB_df['time [day]'] = np.arange(len(rainfall))

        forcing_var = pd.DataFrame({'forc_st': daily_soil_temperatures,
                                    'forc_sw': SWB_df['theta'],
                                    'forc_npp': total_plant_input
                                    })

        # num_years = 200  # Number of years to repeat the data
        # forcing_var_repeated = pd.concat([forcing_var] * num_years, ignore_index=True)    
        
        # # Run the model for this soil type
        # df_sol = dynamic_simulation(forcing_var_repeated, kinetic_dict, soil_dict_temp,
        #                             state={'POM': 1, 'LMWC': 1, 'AGG': 1, 'MIC': 1, 'MAOM': 1, 'CO2': 0},
        #                             SS=False, default_moisture=default_moisture)
        # print(output.sum(axis=1))
        # fig1 = plot_pools(df_sol)
        # Extract the final pool values for the current soil type and texture
        # output = pd.DataFrame(df_sol[['POM [gC m-2]', 'LMWC [gC m-2]', 'AGG [gC m-2]',
        #                               'MIC [gC m-2]', 'MAOM [gC m-2]']].iloc[-1]).T
        
        sol = cal_steady_state(forcing_var, kinetic_dict, soil_dict_temp, 
                                    calc_eigen=False, default_moisture=default_moisture)
        
        output = pd.DataFrame([{
            'POM [gC m-2]': sol[0],
            'LMWC [gC m-2]': sol[1],
            'AGG [gC m-2]': sol[2],
            'MIC [gC m-2]': sol[3],
            'MAOM [gC m-2]': sol[4]
        }])
        

        # Add metadata columns for soil type, texture, and source
        output['Soil Type'] = row['soil']
        output['Soil texture'] = row['texture']
        output['Source'] = row['source']

        # Append to results list
        results.append(output)

    # Combine results into a single DataFrame
    ss_millennial = pd.concat(results, ignore_index=True)
    # ss_millennial.to_excel("millennial_ss.xlsx", index=False)

    # Plot
    if make_plot:
        source_order = ['Intact', 'Repacked', 'HiHydroSoil v2.0']
        texture_order = ['Loam', 'Sandy Loam', 'Silt Loam', 'Silty Clay']

        # Convert to ordered categorical types
        pivot_df = ss_millennial.copy()
        pivot_df = pivot_df.rename(columns={
            'POM [gC m-2]': r'POM',
            'LMWC [gC m-2]': r'LMWC',
            'AGG [gC m-2]': r'AGG',
            'MIC [gC m-2]': r'MIC',
            'MAOM [gC m-2]': r'MAOM'
        })
        pool_order =['POM', 'LMWC', 'AGG',  'MIC', 'MAOM']

        pivot_df['Source'] = pd.Categorical(pivot_df['Source'], categories=source_order, ordered=True)
        pivot_df['Soil texture'] = pd.Categorical(pivot_df['Soil texture'], categories=texture_order, ordered=True)
        pivot_df_sort = pivot_df.sort_values(['Soil texture','Source'])

        # Rebuild index for stacked plot
        pivot_df_sort.set_index(['Soil texture', 'Source'], inplace=True)

        x_labels = [' | '.join(label) for label in pivot_df_sort.index]
        x = np.arange(len(x_labels))
        pool_order =['POM', 'LMWC', 'AGG',  'MIC', 'MAOM']

        # Plotting
        fig, ax = plt.subplots(figsize=(10, 7))

        bottom = np.zeros(len(pivot_df_sort))
        colors = plt.cm.Set2(range(len(pool_order)))
        for i, pool in enumerate(pool_order):
            values = pivot_df_sort[pool].values
            ax.bar(x, values, bottom=bottom, label=pool, color=colors[i], edgecolor='black', alpha=1, width=0.5)
            bottom += values

        # Aesthetics
        ax.set_xticks(x)
        ax.set_yticklabels(ax.get_yticks(), fontsize=16)
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=18)
        ax.set_ylabel("Soil C Stock (gC m$^{-2}$)", fontsize=20)
        # ax.set_xlabel("Soil texture | Source", fontsize=20)
        ax.legend(title="", loc='upper left',  ncol=5, fontsize=16)
        ax.set_ylim(0,9000)
        # plt.tight_layout()
    return ss_millennial