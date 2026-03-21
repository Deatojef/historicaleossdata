#!/usr/bin/python3

from datetime import datetime, timedelta
import pytz
from operator import itemgetter, attrgetter
import re
import pandas as pd
import glob
#from openpyxl import load_workbook
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import savgol_filter
from scipy.stats import zscore
from numpy.polynomial.polynomial import Polynomial as Poly
import multiprocessing as mp
import math
import json
import copy
import psycopg2 as pg
import simplekml


##########################################
# get the various spreadsheet tab names in an excel spreadsheet
##########################################
#def getsheetnames(filepath:str = None)->list:
#    wb = load_workbook(filepath, read_only=True, keep_links=False)
#    return wb.sheetnames


#####################################
# Function to use to determine distance between two points
def distance(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))

    #r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    r = 3956 # Radius of earth in kilometers. Use 3956 for miles

    return c * r


################################
# Function for connecting to the database
################################
def connectToDatabase(dbstring):
    try:

        # If not already connected to the database, then try to connect
        con = pg.connect(dbstring)
        con.set_session(autocommit=True)

        return con

    except pg.DatabaseError as error:

        # If there was a connection error, then close these, just in case they're open
        print(f"Database error: {error}")
        return None



##########################################
# read in a JSON file and return the resulting dataframe
##########################################
def readJSONData(filename:str = None)->(pd.DataFrame, pd.DataFrame, pd.DataFrame, list, list):

    json_data = None
    with open(filename, "r") as file:
        json_data = json.loads(file.read())  
    df = pd.json_normalize(json_data[0]["packets"]);

#
#        packets: [
#             {
#                "info": "/144844h3906.12N/10306.60WO179/017/A=005366 -47T8384P EOSS BALLOON",
#                "receivetime": "2025-03-23T08:48:48.29",
#                "packettime": "08:48:44",
#                "callsign": "AE0SS-4",
#                "raw": "AE0SS-4>APZEOS,EOSS,WIDE2-1,qAO,K0SCC-8:/144844h3906.12N/10306.60WO179/017/A=005366 -47T8384P EOSS BALLOON",
#                "bearing": 179,
#                "speed_mph": 20,
#                "altitude_ft": 5366,
#                "lat_deg": 39.102,
#                "lon_deg": -103.11,
#                "vert_rate_fts": -17.23,
#                "vert_rate_ftmin": -1034,
#                "lat_rate_deg": -0.0001166667,
#                "lon_rate_deg": -1.66667e-05,
#                "temperature_k": 268.45,
#                "temperature_c": -4.7,
#                "temperature_f": 23.54,
#                "pressure_pa": 83840,
#                "pressure_atm": 0.8274,
#                "timedelta_s": 30
#            },
#            ...
#            ...

    # sort it by packettime
    df = df.sort_values(by=['packettime'], ascending=True)

    return df


##########################################
# process a data frame, returning a trimmed dataframe with additional calculated columns
##########################################
def process_df(flightname:str, df:pd.DataFrame)->(pd.DataFrame, pd.DataFrame, pd.DataFrame, None, None):


    # make a copy of this df
    #df = df.copy()
    #df = df.reset_index()

    # drop empty rows
    df = df.dropna(how='all')

    # add a flightid column
    df["flightid"] = flightname

    # Add a column that identifies those rows that contain position information (or not)
    df["position_packet"] = (df["altitude_ft"] > 0) & (df["altitude_ft"].notnull()) & (df["latitude"].notnull()) & (df["longitude"].notnull())

    # look for a packet from the beacons that report the detected burst altitude
    burst_packets = df[(df['raw'].str.contains(r'DETECTED.*BURST', case=False, na=False)) | (df['raw'].str.contains(r'DETECTED COMMANDED RELEASE', case=False, na=False))]
    detected_burst = None

    # if there are matching burst packets, then try and parse the first occurance of that and save the altitude to 'detected_burst'
    detected_altitudes = []
    for _, row in burst_packets.iterrows():
        match = re.search(r'(\d+)\s*ft', row['raw'], re.IGNORECASE)
        if match:
            detected_altitudes.append(int(match.group(1)))

    detected_burst = max(detected_altitudes) if detected_altitudes else None
        
    # callsign 
    callsigns = df["callsign"].unique().tolist()

    # new empty dataframe
    consolidated = pd.DataFrame()
    status_df = df[df["position_packet"] == False]

    # loop through those packets for each beacon creating calculated columns (ex. timedelta_s, velocity, acceleration, etc.)
    for beacon in callsigns:

        # this beacon's rows that have some sort of location information
        beacon_df = df[(df["callsign"] == beacon) & (df["position_packet"] == True)]

        # Sort this based on packettime
        beacon_df = beacon_df.sort_values(by=['packettime'], ascending=True)

        # check for duplicate packets.  We keep the first one encountered.
        beacon_df = beacon_df.drop_duplicates(subset='packettime', keep='first')

        # the time delta (in secs) between data points.
        beacon_df["timedelta_s"] = beacon_df["packettime"].diff(1).dt.total_seconds()
        #ascent['elapsed_secs'] = pd.to_timedelta(ascent['packettime'] - ascent_start).dt.total_seconds()

        # velocity values
        beacon_df["velocity_x_degs"] = (beacon_df["longitude"].diff(1) / beacon_df["timedelta_s"]).round(8)
        beacon_df["velocity_y_degs"] = (beacon_df["latitude"].diff(1) / beacon_df["timedelta_s"]).round(8)
        beacon_df["velocity_z_fts"] = (beacon_df["altitude_ft"].diff(1) / beacon_df["timedelta_s"]).round(8)
        beacon_df["velocity_z_ms"] = (beacon_df["altitude_m"].diff(1) / beacon_df["timedelta_s"]).round(8)
        beacon_df["vert_rate_ftmin"] = (beacon_df["velocity_z_fts"] * 60.0).round(1)

        # air density
        beacon_df["airdensity_slugs"] = ((beacon_df["pressure_pa"] / (287.05 * beacon_df["temperature_k"])) / 515.3788199999872).round(8)
        beacon_df["airdensity_kgm3"] = (beacon_df["pressure_pa"] / (287.05 * beacon_df["temperature_k"])).round(8)



        # add these position rows to the consolidated dataframe with all beacons
        if beacon_df.shape[0] > 0:
            consolidated = pd.concat([consolidated, beacon_df])

    # if there are rows to process...
    if consolidated.shape[0] > 0:

        # Sort this based on packettime
        consolidated = consolidated.sort_values(by=['packettime'], ascending=True, ignore_index=True).reset_index()

        # default value for Airflow
        consolidated["airflow"] = "n/a"

        # Determine the max altitude, then use that to create new column to indiate if the flight is ascending (True) or descending (False)
        altmax_idx = consolidated["altitude_ft"].idxmax()
        max_altitude_ft = consolidated.loc[altmax_idx, "altitude_ft"]
        consolidated["ascending"] = False
        consolidated.loc[0:altmax_idx, "ascending"] = True

        # find those packets that were transmitted much later, after the flight had already landed.
        stragglers = consolidated[(consolidated["timedelta_s"] > 500) & (consolidated["altitude_ft"] < 8000) & (consolidated["ascending"] == False)].index
        if stragglers.shape[0] > 0:
            #print(f"{stragglers.shape[0]=}, {stragglers[0]}")
            #for i, r in consolidated[stragglers[0]-5:].tail(10).iterrows():
            #    print(f"{i:<5} {r['packettime']}   {r['callsign']:<10} {r['altitude_ft']:<8} {r['velocity_z_fts']:<13} {r['timedelta_s']}")
            consolidated = consolidated[0:stragglers[0]+1]

        # eject obvious outliers.  Ascent rates should always be > -5 and < 50.  Descent rates should always be, < 5.
        consolidated = consolidated[   ((consolidated["velocity_z_fts"] < 50) & (consolidated["velocity_z_fts"] > -5) & (consolidated["ascending"] == True))     |      ((consolidated["velocity_z_fts"] < 5) & (consolidated["ascending"] == False))   ]

        # sort this again based on packettime and reset the index having just removed outliers and stragglers
        consolidated = consolidated.sort_values(by=['packettime'], ascending=True, ignore_index=True).reset_index()

        # create forward moving averages 
        consolidated["forward_avg"] = (consolidated[::-1]["velocity_z_fts"].rolling(3).mean().round(2))[::-1]
        consolidated["long_forward_avg"] = (consolidated[::-1]["velocity_z_fts"].rolling(10).mean().round(2))[::-1]

        # Trim off the initial part of the telemetry prior to actual launch.  We do this by looking at two different length (short & long), forward looking, moving averages.  Where both short
        # and long period MA's show an upward velocity > 5ft/s, then we consider the flight as having been "launched" and trim off any packets prior to that point.
        startingindex_locations = consolidated.loc[(consolidated["velocity_z_fts"] > 0) & (consolidated["forward_avg"] > 5) & (consolidated["long_forward_avg"] > 5) & (consolidated["ascending"] == True) & (consolidated["altitude_ft"] < 8000)].index

        if startingindex_locations.shape[0] == 0:
            starting_idx = 0
            #print(f"no starting locations:  {startingindex_locations=}")
            #for idx, r in consolidated.iterrows():
            #    print(f"{r['index']}, {r['packettime']}, {r['callsign']}, {r['timedelta_s']}, {r['altitude_ft']}, {r['ascending']}, {r['velocity_z_fts']}")
            #return (pd.DataFrame(), pd.DataFrame(), None)

        elif startingindex_locations[0] < 1:
            starting_idx = 0
        else:
            #starting_idx = startingindex_locations[0] - 6
            starting_idx = startingindex_locations[0] - 1


        #ending_points = consolidated.loc[(consolidated["velocity_z_fts"] > -5) & (consolidated["forward_avg"] > -5) & (consolidated["long_forward_avg"] > -5) & (consolidated["velocity_z_fts"] < 0) & (consolidated["ascending"] == False)].index
        ending_points = consolidated.loc[(consolidated["velocity_z_fts"] > -5) & (consolidated["forward_avg"] > -5) & (consolidated["long_forward_avg"] > -5) & (consolidated["ascending"] == False) & (consolidated["altitude_ft"] < 8000)].index
        ending_idx = consolidated.shape[0] - 1
        for i in ending_points:
            if i > starting_idx:
                ending_idx = i + 1 
                break

        # check the limits on the start and end points
        if starting_idx < 0:
            starting_idx = 0
        if ending_idx > consolidated.shape[0]-1:
            ending_idx = consolidated.shape[0]-1

        #print(f"################ {flightname}: starting data points #####################")
        #print("pos, idx, packettime, callsign, altitude_ft, velocity_z_fts, forward_avg, long_forward_avg")
        #for i, r in consolidated.iterrows():
        #    print(f"{i:<5} {r['packettime']}   {r['callsign']:<10} {r['altitude_ft']:<8} {r['velocity_z_fts']:<13} {r['forward_avg']:<8} {r['long_forward_avg']:<8}")
        #    #print(f"{i:<5} {r['packettime']}   {r['callsign']:<10} {r['altitude_ft']:<8} {r['velocity_z_fts']:<13} {r['raw']}")
    
        #print(f"################ {flightname}: starting data points #####################")
        #print("pos, idx, packettime, callsign, altitude_ft, velocity_z_fts, forward_avg, long_forward_avg")
        #for i, r in consolidated[starting_idx-5:].head(10).iterrows():
        #    w = "pre" if i < starting_idx else "post"
        #    if i == starting_idx: 
        #        w = "idx"
        #    print(f"{w:<5} {i:<5} {r['packettime']}   {r['callsign']:<10} {r['altitude_ft']:<8} {r['velocity_z_fts']:<13} {r['forward_avg']:<8} {r['long_forward_avg']:<8}")

        #print(f"################ {flightname}: Ending data points #####################")
        #print("pos, idx, packettime, callsign, altitude_ft, velocity_z_fts, forward_avg, long_forward_avg, timedelta_s")
        #for i, r in consolidated[:ending_idx+5].tail(10).iterrows():
        #    w = "pre" if i < ending_idx else "post"
        #    if i == ending_idx: 
        #        w = "idx"
        #    print(f"{w:<5} {i:<5} {r['packettime']}   {r['callsign']:<10} {r['altitude_ft']:<8} {r['velocity_z_fts']:<13} {r['forward_avg']:<8} {r['long_forward_avg']:<8} {r['timedelta_s']}")

        #print("\n")

        #print(f"num: {num_rows}, starting_idx: {starting_idx}, ending_idx: {ending_idx}")
        #print(consolidated.iloc[starting_idx-5:starting_idx+5][["packettime", "callsign", "altitude_ft", "velocity_z_fts"]])
        #print(consolidated.iloc[ending_idx-10:ending_idx+10][["packettime", "callsign", "altitude_ft", "velocity_z_fts", "forward_avg", "long_forward_avg"]])
        #print(consolidated.iloc[ending_idx:][["packettime", "callsign", "altitude_ft", "velocity_z_fts", "forward_avg", "long_forward_avg"]])

        #for idx, r in consolidated.iterrows():
        #    print(f"{r['index']}, {r['packettime']}, {r['callsign']}, {r['timedelta_s']}, {r['altitude_ft']}, {r['ascending']}, {r['velocity_z_fts']}")

        #print(f"0, {starting_idx=}, {ending_idx=}, {consolidated.shape[0]}")

        # Roughly trim off rows from the beginning and end that we don't need.  We keep an few rows before and after our start/ending points.
        #consolidated = consolidated.iloc[starting_idx:ending_idx]
        consolidated = consolidated[starting_idx:ending_idx]

        #print(f"1, {starting_idx=}, {ending_idx=}, {consolidated.shape[0]}")

        # now remove any rows from the beginning with a velocity_z_fts < 1
        #consolidated.drop(consolidated[(consolidated["velocity_z_fts"] < 1.0) & (consolidated["ascending"] == True) & (consolidated["altitude_ft"] < 8000)].index, inplace=True)

        #print(f"2, {starting_idx=}, {ending_idx=}, {consolidated.shape[0]}")
        
        # now remove any rows from the end with a velocity_z_fts > -1
        #consolidated.drop(consolidated[(consolidated["velocity_z_fts"] > -.1) & (consolidated["ascending"] == False) & (consolidated["altitude_ft"] < 8000)].index, inplace=True)

        #print(f"3, {starting_idx=}, {ending_idx=}, {consolidated.shape[0]}")
        #for idx, r in consolidated.iterrows():
        #    print(f"{r['index']}, {r['packettime']}, {r['callsign']}, {r['timedelta_s']}, {r['altitude_ft']}, {r['ascending']}, {r['velocity_z_fts']}")

        # calculate the distance from the launch location to each data point
        firstpoint = consolidated.iloc[0]
        consolidated["distance_from_launch_mi"] = distance(firstpoint["latitude"], firstpoint["longitude"], consolidated["latitude"], consolidated["longitude"]).round(2)
        consolidated["distance_from_launch_km"] = (distance(firstpoint["latitude"], firstpoint["longitude"], consolidated["latitude"], consolidated["longitude"]) * 1.609344).round(2)

        # number of rows
        num_rows = consolidated.shape[0]

        # We're now ready to split the data frame into ascending and descending portions of the flight.
        ascent = consolidated.loc[consolidated["ascending"] == True].copy()
        descent = consolidated.loc[consolidated["ascending"] == False].copy()

        # lambda function to calculate the perpendicular distance from time, altitude data point during to that line
        # negative numbers indicate a point above the line, positive values indicate a point below the line
        #dist = lambda x, y, m, b : abs(m * x - y + b) / math.sqrt(m**2 + 1) 
        dist = lambda x, y, m, b : (m * x - y + b) / math.sqrt(m**2 + 1) 

        if ascent.shape[0] > 2:


            #print(outliers[["packettime", "callsign", "altitude_ft", "velocity_z_fts", "ascending"]])
            #print(ascent.iloc[0:10][["packettime", "callsign", "altitude_ft", "velocity_z_fts", "ascending"]])

            # starting time values for each phase of the flight
            ascent_start = ascent.iloc[0]["packettime"]

            # elapsed secs for each phase of the flight
            ascent['elapsed_secs'] = pd.to_timedelta(ascent['packettime'] - ascent_start).dt.total_seconds()

            # calculate the equation for a straight line (y=mx+b) between the beginning and ending points for each phase of the flight
            m_ascent = (ascent.iloc[-1]["altitude_ft"] - ascent.iloc[0]["altitude_ft"]) / (ascent.iloc[-1]["elapsed_secs"] - ascent.iloc[0]["elapsed_secs"])
            b_ascent = ascent.iloc[0]["altitude_ft"] - m_ascent * ascent.iloc[0]["elapsed_secs"]

            # add this distance to the ascent phase of the flight
            ascent["distance_to_line"] = dist(ascent["elapsed_secs"], ascent["altitude_ft"], m_ascent, b_ascent).round(6)

            # calculate the acceleration
            ascent["acceleration_fts2"] = (ascent["velocity_z_fts"].diff() / ascent["timedelta_s"]).round(6)
            ascent["acceleration_ms2"] = (ascent["velocity_z_ms"].diff() / ascent["timedelta_s"]).round(6)

            # Drop the first two rows from the ascent phase as the velocity and acceleration values are NaN
            #ascent = ascent.drop(ascent.iloc[0:2].index)

            # Calculate the mean for each phase of the flight (note:  Pandas calculates the unbiased variance/std)
            # These are "running" mean values.  That is at each point, calculate the mean for all points prior.
            ascent["velocity_mean_fts"] = ascent["velocity_z_fts"].expanding().mean().round(6)
            ascent["acceleration_mean_fts2"] = ascent["acceleration_fts2"].expanding().mean().round(6)
            ascent["velocity_mean_ms"] = ascent["velocity_z_ms"].expanding().mean().round(6)
            ascent["acceleration_mean_ms2"] = ascent["acceleration_ms2"].expanding().mean().round(6)

            # Calculate standard deviations for each phase of the flight (note:  Pandas calculates the unbiased variance/std)
            # These are "running" std values.  That is at each point, calculate the std for all points prior.
            ascent["velocity_std_fts"] = ascent["velocity_z_fts"].expanding().std().round(6)
            ascent["acceleration_std_fts2"] = ascent["acceleration_fts2"].expanding().std().round(6)
            ascent["velocity_std_ms"] = ascent["velocity_z_ms"].expanding().std().round(6)
            ascent["acceleration_std_ms2"] = ascent["acceleration_ms2"].expanding().std().round(6)

            # Normalize the acceleration and velocity for each phase of the flight for each point, use the mean & std for all points 
            # prior to calculate the "standard score". These figures are scale invariant (i.e. no units). 
            ascent["velocity_norm_fts"] = ((ascent["velocity_z_fts"] - ascent["velocity_mean_fts"]) / ascent["velocity_std_fts"]).round(6)
            ascent["acceleration_norm_fts2"] = ((ascent["acceleration_fts2"] - ascent["acceleration_mean_fts2"]) / ascent["acceleration_std_fts2"]).round(6)
            ascent["velocity_norm_ms"] = ((ascent["velocity_z_ms"] - ascent["velocity_mean_ms"]) / ascent["velocity_std_ms"]).round(6)
            ascent["acceleration_norm_ms2"] = ((ascent["acceleration_ms2"] - ascent["acceleration_mean_ms2"]) / ascent["acceleration_std_ms2"]).round(6)

            # Determine a good degree for a polynomial curve fit.  This is based off of the VMR so that low dispersion uses a higher degree polynomial
            # and a high dispersion uses a lower degree polynomial.
            var_a = ascent["velocity_z_fts"].var()
            var_a_m = ascent["velocity_z_ms"].var()
            VMR = var_a / ascent["velocity_z_fts"].mean()
            VMR_m = var_a_m / ascent["velocity_z_ms"].mean()
            max_degree = 13
            if VMR > 1.75:
                v_variance_a = int((max_degree + VMR) / VMR)
            elif 1.75 >= VMR >= 1.0:
                v_variance_a = int(max_degree / VMR)
            else:
                v_variance_a = max_degree 

            if VMR_m > 1.75:
                v_variance_a_m = int((max_degree + VMR_m) / VMR_m)
            elif 1.75 >= VMR_m >= 1.0:
                v_variance_a_m = int(max_degree / VMR_m)
            else:
                v_variance_a_m = max_degree 

            # Find a line of best fit (using least squares) for the velocity vs altitude plot
            # numpy's polynomial fit
            v_coef_a, stats_a = Poly.fit(ascent["altitude_ft"], ascent["velocity_z_fts"], v_variance_a, full=True)
            v_coef_a_m, stats_a_m = Poly.fit(ascent["altitude_m"], ascent["velocity_z_ms"], v_variance_a_m, full=True)
            
            # acceleration data from the polynomial curve fit...add this to each phase of the flight
            ascent["velocity_curvefit_fts"] = v_coef_a(ascent["altitude_ft"]).round(6)
            ascent["velocity_curvefit_ms"] = v_coef_a_m(ascent["altitude_m"]).round(6)

            # Now determine the Reynold's transition areas for the ascent phase of the flight
            # We do this by finding the points (during the ascent) that are furthest away from the average ascent velocity.  That is,
            # we want to find those two points that are the largest positive and negative numbers for "distance_to_line".  With those
            # two points we can defined a maximum of three airflow environments.
            max_dist_idx = 0
            min_dist_idx = 0
            max_dist = 0
            min_dist = 0

            minimum_distance = 30
            altitude_min = 10000

            # assuming that any Re transition would occur at much higher MSL.  So putting a altitude floor to filter out early flight jitters.
            if ascent[(ascent["distance_to_line"] > minimum_distance) & (ascent['altitude_ft'] > altitude_min)].shape[0] > 0:
                # largest positive distance from the average ascent velocity
                max_dist_idx = ascent[(ascent["distance_to_line"] > minimum_distance) & (ascent['altitude_ft'] > altitude_min)]["distance_to_line"].idxmax()
                max_dist = ascent.loc[max_dist_idx, "distance_to_line"]

            # assuming that any Re transition would occur at much higher MSL.  So putting a altitude floor to filter out early flight jitters.
            if ascent[(ascent["distance_to_line"] < -minimum_distance) & (ascent['altitude_ft'] > altitude_min)].shape[0] > 0:
                # smallest negative distance from the average ascent velocity
                min_dist_idx = ascent[(ascent["distance_to_line"] < -minimum_distance) & (ascent['altitude_ft'] > altitude_min)]["distance_to_line"].idxmin()
                min_dist = ascent.loc[min_dist_idx, "distance_to_line"]

            #print(f"{minimum_distance=}, {min_dist_idx=}, {max_dist_idx=}, {min_dist=}, {max_dist=}, delta: {max_dist-min_dist}")

            if max_dist > 0 and min_dist < 0:
                # two points of Re transition
                if max_dist_idx > min_dist_idx:
                    ascent.loc[0:min_dist_idx, "airflow"] = "high Re"
                    ascent.loc[min_dist_idx:max_dist_idx, "airflow"] = "low Re"
                    ascent.loc[max_dist_idx:, "airflow"] = "high Re"
                else:
                    ascent.loc[0:max_dist_idx, "airflow"] = "low Re"
                    ascent.loc[max_dist_idx:min_dist_idx, "airflow"] = "high Re"
                    ascent.loc[min_dist_idx:, "airflow"] = "low Re"

            elif max_dist > 0 and min_dist == 0:
                # single point of Re transition:  from laminar to turbulent
                ascent.loc[0:max_dist_idx, "airflow"] = "low Re"
                ascent.loc[max_dist_idx:, "airflow"] = "high Re"

            elif max_dist == 0 and min_dist < 0:
                # single point of Re transition:  from turbulent to laminar.  This is the most common case.
                ascent.loc[0:min_dist_idx, "airflow"] = "high Re"
                ascent.loc[min_dist_idx:, "airflow"] = "low Re"

            else:  # they're both zero?  ...implies no deviations from the average velocity?  ...should never get here, unless all transitions were filtered out.  ;)
                pass


            rbp_idx = ascent["distance_to_line"].abs().idxmax()
            rbp_alt = ascent.loc[rbp_idx, "altitude_ft"]

            # now label data points as either turbulent or laminar for the ascent phase
            #if rbp_idx is not None:
            #    ascent.loc[0:rbp_idx, "airflow"] = "turbulent"
            #    ascent.loc[rbp_idx:, "airflow"] = "laminar"

        else:
            print(f"ERROR:  {flightname.lower()} {callsigns} ascent shape was: {ascent.shape[0]}")
            return (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), None, None)

        if descent.shape[0] > 2:

            #print(descent.iloc[-10:][["packettime", "callsign", "altitude_ft", "velocity_z_fts", "forward_avg", "long_forward_avg"]])

            # starting time values for each phase of the flight
            descent_start = descent.iloc[0]["packettime"]

            # elapsed secs for each phase of the flight
            descent['elapsed_secs'] = pd.to_timedelta(descent['packettime'] - descent_start).dt.total_seconds()

            # calculate the equation for a straight line (y=mx+b) between the beginning and ending points for each phase of the flight
            m_descent = (descent.iloc[-1]["altitude_ft"] - descent.iloc[0]["altitude_ft"]) / (descent.iloc[-1]["elapsed_secs"] - descent.iloc[0]["elapsed_secs"])
            b_descent = descent.iloc[0]["altitude_ft"] - m_descent * descent.iloc[0]["elapsed_secs"]

            # add this distance to each phase of the flight
            descent["distance_to_line"] = dist(descent["elapsed_secs"], descent["altitude_ft"], m_descent, b_descent).round(6)

            # calculate the acceleration
            descent["acceleration_fts2"] = (descent["velocity_z_fts"].diff() / descent["timedelta_s"]).round(6)
            descent["acceleration_ms2"] = (descent["velocity_z_ms"].diff() / descent["timedelta_s"]).round(6)

            # Calculate the mean for each phase of the flight (note:  Pandas calculates the unbiased variance/std)
            # These are "running" mean values.  That is at each point, calculate the mean for all points prior.
            descent["velocity_mean_fts"] = descent["velocity_z_fts"].expanding().mean().round(6)
            descent["acceleration_mean_fts2"] = descent["acceleration_fts2"].expanding().mean().round(6)
            descent["velocity_mean_ms"] = descent["velocity_z_ms"].expanding().mean().round(6)
            descent["acceleration_mean_ms2"] = descent["acceleration_ms2"].expanding().mean().round(6)

            # Calculate standard deviations for each phase of the flight (note:  Pandas calculates the unbiased variance/std)
            # These are "running" std values.  That is at each point, calculate the std for all points prior.
            descent["velocity_std_fts"] = descent["velocity_z_fts"].expanding().std().round(6)
            descent["acceleration_std_fts2"] = descent["acceleration_fts2"].expanding().std().round(6)
            descent["velocity_std_ms"] = descent["velocity_z_ms"].expanding().std().round(6)
            descent["acceleration_std_ms2"] = descent["acceleration_ms2"].expanding().std().round(6)

            # Normalize the acceleration and velocity for each phase of the flight for each point, use the mean & std for all points 
            # prior to calculate the "standard score". These figures are scale invariant (i.e. no units). 
            descent["velocity_norm_fts"] = ((descent["velocity_z_fts"] - descent["velocity_mean_fts"]) / descent["velocity_std_fts"]).round(6)
            descent["acceleration_norm_fts2"] = ((descent["acceleration_fts2"] - descent["acceleration_mean_fts2"]) / descent["acceleration_std_fts2"]).round(6)
            descent["velocity_norm_ms"] = ((descent["velocity_z_ms"] - descent["velocity_mean_ms"]) / descent["velocity_std_ms"]).round(6)
            descent["acceleration_norm_ms2"] = ((descent["acceleration_ms2"] - descent["acceleration_mean_ms2"]) / descent["acceleration_std_ms2"]).round(6)

            # Determine a good degree for a polynomial curve fit.  This is based off of the VMR so that low dispersion uses a higher degree polynomial
            # and a high dispersion uses a lower degree polynomial.
            var_d = descent["velocity_z_fts"].var()
            var_d_m = descent["velocity_z_ms"].var()
            VMR = var_d / descent["velocity_z_fts"].mean()
            VMR_m = var_d_m / descent["velocity_z_ms"].mean()
            max_degree = 13
            if VMR > 1.75:
                v_variance_d = int((max_degree + VMR) / VMR)
            elif 1.75 >= VMR >= 1.0:
                v_variance_d = int(max_degree / VMR)
            else:
                v_variance_d = max_degree 

            if VMR_m > 1.75:
                v_variance_d_m = int((max_degree + VMR_m) / VMR_m)
            elif 1.75 >= VMR_m >= 1.0:
                v_variance_d_m = int(max_degree / VMR_m)
            else:
                v_variance_d_m = max_degree 

            # Find a line of best fit (using least squares) for the velocity vs altitude plot
            # numpy's polynomial fit
            v_coef_d, stats_d = Poly.fit(descent["altitude_ft"], descent["velocity_z_fts"], v_variance_d, full=True)
            v_coef_d_m, stats_d_m = Poly.fit(descent["altitude_m"], descent["velocity_z_ms"], v_variance_d_m, full=True)
            
            # acceleration data from the polynomial curve fit...add this to each phase of the flight
            descent["velocity_curvefit_fts"] = v_coef_d(descent["altitude_ft"]).round(6)
            descent["velocity_curvefit_ms"] = v_coef_d_m(descent["altitude_m"]).round(6)

        else:
            print(f"ERROR: {flightname.lower()} {callsigns} descent shape was: {descent.shape[0]}")

        #for idx, r in ascent.iterrows():
        #    print(f"{r['index']}, {r['callsign']}, {r['altitude_ft']}, {r['velocity_z_fts']}, {r['airflow']}, {r['distance_to_line']}")

        # drop the temp columns added during the analysis
        ascent = ascent.drop(['distance_to_line', 'forward_avg', 'long_forward_avg', 'timedelta_s'], axis=1)
        descent = descent.drop(['distance_to_line', 'forward_avg', 'long_forward_avg', 'timedelta_s'], axis=1)

        # print out some stats
        print(f"{flightname.lower()} {callsigns}   final rows: {num_rows}, trim: {starting_idx}:{ending_idx}, ascent: {ascent.shape[0]}, descent: {descent.shape[0]}, max_alt: {max_altitude_ft}, rbp_idx: {rbp_idx}, rbp_alt: {rbp_alt}, curvefit_deg: {v_variance_a}, VMR: {VMR:.2f}, curvefit_deg_m: {v_variance_a_m}, VMR_m: {VMR_m:.2f}")


        # return the two pandas dataframes and the degree of the polynomial used to fit a curve to the ascent velocity
        return (ascent, descent, status_df, v_variance_a, detected_burst)

    else:
        print(f"No data to return")
        return (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), None, None)



##########################################
# createPlot 
# this takes the ascent and descent pandas dataframes and saves a plot of the flight to the 'filename' provided as a PNG.
##########################################
def createPlot(filename:str, ascent: pd.DataFrame, descent: pd.DataFrame, degree: int)->None:

    if ascent.shape[0] < 2 or descent.shape[0] < 2:
        return

    # callsign 
    callsigns = [x.lower() for x in pd.concat([ascent["callsign"], descent["callsign"]]).unique().tolist()]

    # find the Reynolds transition points
    Re_points = ascent[ascent["airflow"].shift() != ascent["airflow"]].index

    # if there's more than one transition point, then lop off the first one as that's likely just the first data point.
    if Re_points.shape[0] > 1:
        Re_points = Re_points[1:]

    # burst point location (for the Altitude plot)
    burst_idx = ascent["altitude_ft"].idxmax()
    burst_alt_ft = ascent.loc[burst_idx, "altitude_ft"]
    burst_packettime = ascent.iloc[-1]["packettime"]

    # Altitude limits for the Velocity and Altitude plots
    ylims = [0, ascent.loc[burst_idx, "altitude_ft"] + 10000]

    # title of the plot
    title = filename.split(".")[0]
    title = title.split("/")
    title = title[len(title) - 1]

    ## Figure
    fig = plt.figure(1, figsize=(17, 15))
    fig.suptitle(f"30-Second Flight Data: {title}, {callsigns}", fontsize=14)
    gs = GridSpec(2, 3, figure=fig)
    ax =  fig.add_subplot(gs[0, 2])     # the velocity plot
    ax2 = fig.add_subplot(gs[0, :-1])    # the altitude plot
    ax3 = fig.add_subplot(gs[1, :-1])   # the acceleration plot
    ax4 = fig.add_subplot(gs[1, 2])   # ACF


    #####################
    ## Velocity plot
    #####################
    ax.set_xlabel(f'Velocity ($ft/min$)', fontsize=8)
    ax.set_ylabel(f'Altitude ($ft$)', fontsize=8)
    ax.set_title(f'Ascent Velocity vs. Altitude', fontsize=12)
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.minorticks_on()
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.set_ylim(ylims)
    ax.yaxis.set_label_position("right")
    ax.tick_params('y', labelleft=False, labelright=True)
    ax.grid(visible=True, which='minor', color='k', linestyle='--', alpha=0.1)
    ax.grid(visible=True, which='major')

    # plot the Velocity vs. Altitude
    v_ascentplot, = ax.plot(ascent["velocity_z_fts"]*60, ascent["altitude_ft"], '.', color="tab:blue", label=f"Velocity", alpha=.3)

    # polynomial curve fit for the velocity
    v_polyfit, = ax.plot(ascent["velocity_curvefit_fts"]*60, ascent["altitude_ft"], '-.', color="black", label=f"Curve Fit ({degree}deg poly)", linewidth=2)

    # Mark the Reynolds transition areas
    for r in Re_points:
        bp_alt = ascent.loc[r, "altitude_ft"]
        v_hline = ax.hlines(y=bp_alt, xmin=ascent.loc[ascent["velocity_z_fts"].idxmin(), "velocity_z_fts"]*60, xmax=ascent.loc[ascent["velocity_z_fts"].idxmax(), "velocity_z_fts"]*60, color='tab:red', linestyles="dotted", linewidth=2, label=fr"Reynolds Transition (${round(bp_alt):,}ft$)")
        v_tb = ax.text(x=ascent.loc[ascent["velocity_z_fts"].idxmax(), "velocity_z_fts"]*60, y=bp_alt-5000, s=f"{round(bp_alt):,}ft", ha="right", va="center", size=12, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

    # Legend location
    ax.legend(fontsize=8, loc='upper right')


    #####################
    ## altitude plot
    #####################
    ax2.set_ylabel(f'Altitude ($ft$)', fontsize=8)
    ax2.set_title(f'Altitude vs. Time', fontsize=12)
    ax2.set_xlabel('Time', fontsize=8)
    ax2.ticklabel_format(useOffset=False, style="plain")
    ax2.minorticks_on()
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)
    ax2.set_ylim(ylims)
    #ax2.tick_params('x', labelbottom=False)
    ax2.grid(visible=True, which='minor', color='k', linestyle='--', alpha=0.1)
    ax2.grid(visible=True, which='major')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

    start = [ascent.iloc[0]["packettime"], ascent.iloc[0]["altitude_ft"]]
    end = [ascent.iloc[-1]["packettime"], ascent.iloc[-1]["altitude_ft"]]

    # Plot the altitude vs. time
    ascentplot, = ax2.plot(ascent["packettime"], ascent["altitude_ft"], ".", color="tab:blue", label=f"Ascent")
    if descent.shape[0] > 1:
        descentplot, = ax2.plot(descent["packettime"], descent["altitude_ft"], ".", color="tab:orange", label=f"Descent")
    straight = ax2.plot((start[0], end[0]), (start[-1], end[-1]), '-.', color='black', label=f"Average", alpha=.8)

    # mark the Reynolds transition areas with a vertical line
    for r in Re_points:
        bp_alt = ascent.loc[r, "altitude_ft"]
        bp_idx = r
        bp_t = ascent.loc[r, "packettime"]
        bp_tminus = bp_t - pd.Timedelta(minutes=15)
        bp_tplus = bp_t + pd.Timedelta(minutes=15)
        alt_vline = ax2.vlines(x=ascent.loc[bp_idx, "packettime"], ymin=ascent.loc[bp_idx, "altitude_ft"] - 8000, ymax=ascent.loc[bp_idx, "altitude_ft"] + 8000, color='tab:red', linestyles="dotted", linewidth=2, label=fr"Reynolds Transition (${round(bp_alt):,}ft$)")
        alt_hline = ax2.hlines(y=bp_alt, xmin=bp_tminus, xmax=bp_tplus, color='tab:red', linewidth=2, linestyles="dotted")
        tb = ax2.text(x=ascent.loc[bp_idx, "packettime"], y=bp_alt-8000, s=f"{round(bp_alt):,}ft", ha="left", va="center", size=12, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

    # Add some text as the burst point
    burst = ax2.text(x=burst_packettime, y=burst_alt_ft+2000, s=f"{round(burst_alt_ft):,}ft", ha="center", va="bottom", size=12, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

    # Legend location
    ax2.legend(fontsize=8, loc='upper left')


    #####################
    # Acceleration Plot
    #####################
    ax3.set_ylabel(f'Normalized Acceleration', fontsize=8)
    ax3.set_title(f'Normalized Acceleration of Smoothed Velocity', fontsize=12)
    ax3.set_xlabel('Time', fontsize=8)
    ax3.ticklabel_format(useOffset=False, style="plain")
    ax3.minorticks_on()
    ax3.xaxis.set_tick_params(labelsize=8)
    ax3.yaxis.set_tick_params(labelsize=8)
    ax3.grid(visible=True, which='minor', color='k', linestyle='--', alpha=0.1)
    ax3.grid(visible=True, which='major')
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax3_nodata = ax3.text(x=.5, y=.5, s="No Data", ha="center", va="center", size=12, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

    # Normalized Acceleration
    #a_ascent, = ax3.plot(ascent["packettime"], ascent["accel_smoothed_norm"], ".", color="slateblue", label=r"Normalized Acceleration", alpha=.3)

    # mark the Reynolds break point
    #a_vline = ax3.vlines(x=bp_t, ymin=ascent.loc[ascent["accel_smoothed_norm"].idxmin(), "accel_smoothed_norm"], ymax=ascent.loc[ascent["accel_smoothed_norm"].idxmax(), "accel_smoothed_norm"], color='tab:red', linestyles="dotted", linewidth=2, label=fr"Reynolds BP (${round(bp_alt):,}ft$)")
    #a_hline = ax3.hlines(y=0, xmin=ascent.loc[ascent["packettime"].idxmin(), "packettime"], xmax=descent.loc[descent["packettime"].idxmax(), "packettime"], color='tab:gray', linewidth=2)

    # Legend location
    #ax3.legend(fontsize=8, loc='upper left')


    #####################
    # The ACF of the pre-launch data
    #####################
    ax4.set_ylabel(f'Autocorrelation', fontsize=8)
    ax4.set_title(f'Autocorrelation Function (pre-launch data only)', fontsize=12)
    ax4.set_xlabel('Lag', fontsize=8)
    ax4.ticklabel_format(useOffset=False, style="plain")
    ax4.minorticks_on()
    ax4.tick_params('y', labelleft=False, labelright=True)
    ax4.yaxis.set_label_position("right")
    ax4.xaxis.set_tick_params(labelsize=8)
    ax4.yaxis.set_tick_params(labelsize=8)
    #ax4.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    nodata = ax4.text(x=.5, y=.5, s="No Data", ha="center", va="center", size=12, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))


    # mark the Reynolds break point
    # Adjust space between the plots
    plt.subplots_adjust(hspace=.15, wspace=.15, top=.9, bottom=.05, right=.95, left=.05)
    #plt.subplots_adjust(hspace=0, wspace=0, top=.9, bottom=.05, right=.95, left=.05)

    # save this to a PNG
    plt.savefig(filename, format="png", bbox_inches='tight', dpi=100)

    # close
    plt.close()


############################################
# function to determine if a given date is under daylight saving time, returns true or false
# arguments:
# $somedate - string in a typical format, Y-m-d
def is_dst(date_str: str)->bool:

    # create a date object
    dt = datetime.strptime(date_str + " 03:00:00", '%Y-%m-%d %H:%M:%S')

    # the mountain time zone.  Probably should fix the upstream data so that the "timezone" for the location of the launch is computed here instead.
    mountain_tz = pytz.timezone('America/Denver')
    local_dt = mountain_tz.localize(dt)

    # if the dst offset is not zero, then this date_str is within DST
    is_dst = True if local_dt.dst() != timedelta(0) else False

    return is_dst


##########################################
# query the database for packets from a specific flight
##########################################
def queryDatabase(flight: pd.DataFrame = None)->pd.DataFrame:

    query_sql = """
        select distinct on (1)
            substring(a.raw from position(':' in a.raw)+1) as info,
            date_trunc('milliseconds', a.tm)::timestamp without time zone as receivetime,
            case
                when a.raw similar to '%%[0-9]{6}h%%' then
                    date_trunc('milliseconds', ((to_timestamp(a.tm::date || ' ' || substring(a.raw from position('h' in a.raw) - 6 for 6), 'YYYY-MM-DD HH24MISS')::timestamp at time zone 'UTC') at time zone %(timezone)s)::timestamp)::timestamp without time zone
                else
                    date_trunc('milliseconds', a.tm)::timestamp without time zone
            end as packettime,
            a.callsign,
            a.raw,
            round(a.bearing::numeric, 1) as bearing,
            round(a.speed_mph::numeric, 1) as speed_mph,
            round(a.speed_mph::numeric * 1.609344, 1) as speed_kph,
            round(a.altitude::numeric, 1) as altitude_ft,
            round(a.altitude::numeric * 0.3048, 2) as altitude_m,
            cast(st_y(a.location2d) as numeric(14, 10)) as latitude,
            cast(st_x(a.location2d) as numeric(14, 10)) as longitude,
            case when a.raw similar to '%% [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P%%' then
                round(32.0 + 1.8 * cast(substring(substring(substring(a.raw from ' [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P') from ' [-]{0,1}[0-9]{1,6}T') from ' [-]{0,1}[0-9]{1,6}') as decimal) / 10.0, 2)
            else
                NULL
            end as temperature_f,
            case when a.raw similar to '%% [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P%%' then
                round(cast(substring(substring(substring(a.raw from ' [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P') from ' [-]{0,1}[0-9]{1,6}T') from ' [-]{0,1}[0-9]{1,6}') as decimal) / 10.0, 2)
            else
                NULL
            end as temperature_c,
            case when a.raw similar to '%% [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P%%' then
                round(273.15 + cast(substring(substring(substring(a.raw from ' [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P') from ' [-]{0,1}[0-9]{1,6}T') from ' [-]{0,1}[0-9]{1,6}') as decimal) / 10.0, 2)
            else
                NULL
            end as temperature_k,
            case
                when a.raw similar to '%% [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P%%' then
                    round(cast(substring(substring(a.raw from '[0-9]{1,6}P') from '[0-9]{1,6}') as decimal) / 10132.5, 2)
                else
                    NULL
            end as pressure_atm, 
            case
                when a.raw similar to '%% [-]{0,1}[0-9]{1,6}T[-]{0,1}[0-9]{1,6}P%%' then
                    round(cast(substring(substring(a.raw from '[0-9]{1,6}P') from '[0-9]{1,6}') as decimal) * 10.0, 2)
                else
                    NULL
            end as pressure_pa


        from
            packets a

        where
            --a.location2d is not null
            --and a.altitude != 0
            a.tm > %(starttime)s and a.tm < %(endtime)s
            and a.callsign in %(callsignlist)s
            and a.raw not like '%%WA0GEH-10%%'

        order by
            1, 2
    ;"""


    try:

        # Connect to the database
        dbconn = pg.connect("dbname=jeff")

        if dbconn:

            # create database cursor
            cur = dbconn.cursor()

            # parameters for the SQL query
            launchdate = flight["day"]
            timezone = "MDT" if is_dst(launchdate) else "MST"
            startdate = launchdate + " 03:00:00"
            enddate = launchdate + " 23:59:59"
            beacons = flight["beacons"]

            # Execute the query
            cur.execute(query_sql, { "timezone": timezone, "starttime": startdate, "endtime": enddate, "callsignlist": tuple(beacons)})

            # fetch the returned rows
            rows = cur.fetchall()

            if rows is not None:

                # the list of columns names returned from the query
                column_names = [desc[0] for desc in cur.description]

                # list of numeric columns
                numeric_columns = ['bearing', 'speed_mph', 'speed_kph', 'altitude_ft', 'altitude_m', 'latitude', 'longitude', 'temperature_f', 'temperature_c', 'temperature_k', 'pressure_pa', 'pressure_atm']

                # create a new data frame
                df = pd.DataFrame(rows, columns=column_names)

                # Convert the specified columns to float64
                df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric).astype('float64')

            else:
                return None

            cur.close()
            dbconn.close()

            return df

        else:
            return None

    except pg.DatabaseError as error:
        print(f"Database error: {error}")


    return None



##########################################
# thread process
##########################################
def processThread(flight: dict = None)->None:

    # query data from the database
    df = queryDatabase(flight)

    # the flightname
    flightname = flight["flight"]

    # process this data frame
    a, d, s, deg, detected_burst = process_df(flightname, df)

    #print(f"################ {flightname}: ascent data #####################")
    #print("pos, idx, packettime, callsign, altitude_ft, velocity_z_fts")
    #for i, r in a.head(10).iterrows():
    #    print(f"{i:<5} {r['packettime']}   {r['callsign']:<10} {r['altitude_ft']:<8} {r['velocity_z_fts']:<13}")

    # order of the columns that we want data saved as
    new_order = [
            'flightid', 
            'callsign', 
            'receivetime', 
            'packettime', 
            'altitude_ft', 
            'altitude_m', 
            'vert_rate_ftmin', 
            'elapsed_secs', 
            'flight_phase', 
            'position_packet',
            'info', 
            'raw', 
            'bearing',
            'speed_mph', 
            'speed_kph', 
            'latitude', 
            'longitude', 
            'distance_from_launch_mi', 
            'distance_from_launch_km', 
            'temperature_f',
            'temperature_c',
            'temperature_k',
            'pressure_pa', 
            'pressure_atm', 
            'airdensity_slugs',
            'airdensity_kgm3',
            'velocity_x_degs', 
            'velocity_y_degs',
            'velocity_z_fts', 
            'velocity_z_ms', 
            'airflow', 
            'acceleration_fts2',
            'velocity_mean_fts', 
            'acceleration_mean_fts2', 
            'velocity_std_fts',
            'acceleration_std_fts2', 
            'velocity_norm_fts', 
            'acceleration_norm_fts2',
            'velocity_curvefit_fts',
            'acceleration_ms2',
            'velocity_mean_ms', 
            'acceleration_mean_ms2', 
            'velocity_std_ms',
            'acceleration_std_ms2', 
            'velocity_norm_ms', 
            'acceleration_norm_ms2',
            'velocity_curvefit_ms'
         ]
    
    # add this to the flight JSON
    jsondata = {}
    jsondata = copy.deepcopy(flight)
    jsondata["reynolds_transitions"] = []

    # add the flight_phase column, remove some unneeded ones, and reorder the columns
    if not a.empty:

        a.reset_index(drop=True, inplace=True)
        a.loc[a["ascending"] == True, "flight_phase"] = "ascending"
        a = a.drop(['index', 'ascending'], axis=1)
        a = a[new_order]

        # find the Reynolds transition points
        Re_points = a[a["airflow"].shift() != a["airflow"]].index

        # if there's more than one transition point, then lop off the first one as that's likely just the first data point.
        if Re_points.shape[0] > 1:
            Re_points = Re_points[1:]

        # add this column to the ascent dataframe
        a["reynolds_transition"] = ""

        for i in Re_points:
            prior_Re = a.loc[i-1, "airflow"].split(' ')[0]
            next_Re = a.loc[i+1, "airflow"].split(' ')[0]
            reynoldstext = prior_Re + "_to_" + next_Re

            a.loc[i, "reynolds_transition"] = reynoldstext
            if next_Re and prior_Re:
                jsondata["reynolds_transitions"].append({"transition": reynoldstext, "altitude_ft": a.loc[i, "altitude_ft"], "altitude_m": a.loc[i, "altitude_m"] })

    if not d.empty:
        d.reset_index(drop=True, inplace=True)
        d.loc[d["ascending"] == False, "flight_phase"] = "descending"
        d = d.drop(['index', 'ascending'], axis=1)
        d = d[new_order]

        # add this column to the descent dataframe
        d["reynolds_transition"] = ""

    # prior to consolidating status packets, trim off those that occurred before and after the flight.
    min_time = a['packettime'].min()
    max_time = d['packettime'].max()
    status_packets = s[(s['packettime'] >= min_time) & (s['packettime'] <= max_time)]

    # consolidate ascent and descent into a single data frame for csv, json, pickle, and kml output types
    #consolidated = pd.concat([a, d, status_packets])
    consolidated = pd.concat([a, d])
    consolidated = consolidated.sort_values(by=['packettime'])
    
    # max altitude
    max_altitude_ft = consolidated["altitude_ft"].max()
    max_altitude_m = consolidated["altitude_m"].max()
    numpoints = consolidated.shape[0]
    flighttime_seconds = (consolidated["packettime"].max() - consolidated["packettime"].min()) / pd.Timedelta(seconds=1)
    flighttime_hrs = int(flighttime_seconds / 3600)
    flighttime_mins = int((flighttime_seconds - flighttime_hrs*3600) / 60)
    flighttime_secs = int(flighttime_seconds - flighttime_hrs*3600 - flighttime_mins*60)
    flighttime = f"{flighttime_hrs}hrs {flighttime_mins}mins {flighttime_secs}secs"

    # add maximum altitude elements to the flight JSON structure
    flight["maxaltitude_ft"] = max_altitude_ft
    jsondata["maxaltitude_ft"] = max_altitude_ft
    flight["maxaltitude_m"] = max_altitude_m
    jsondata["maxaltitude_m"] = max_altitude_m

    # the detected burst keys
    flight["detected_burst"] = {}
    jsondata["detected_burst"] = {}

    # if a burst was detected by the payload then add that.
    flight["detected_burst"]["detected"] = True if detected_burst else False
    flight["detected_burst"]["burst_ft"] = detected_burst if detected_burst else 0
    flight["detected_burst"]["burst_m"] = round(detected_burst * 0.3048, 2) if detected_burst else 0
    jsondata["detected_burst"]["detected"] = True if detected_burst else False
    jsondata["detected_burst"]["burst_ft"] = detected_burst if detected_burst else 0
    jsondata["detected_burst"]["burst_m"] = round(detected_burst * 0.3048, 2) if detected_burst else 0

    flight["numpoints"] = consolidated.shape[0]
    jsondata["numpoints"] = consolidated.shape[0]
    flight["flighttime"] = flighttime
    jsondata["flighttime"] = flighttime
    flight["flighttime_secs"] = flighttime_seconds
    jsondata["flighttime_secs"] = flighttime_seconds

    # first and last rows
    firstrow = consolidated.iloc[0]
    lastrow = consolidated.iloc[-1]

    # down range distance in miles (i.e. distance between launch and landing lat/lon points)
    flight["range_distance_traveled_mi"] = lastrow["distance_from_launch_mi"]
    jsondata["range_distance_traveled_mi"] = lastrow["distance_from_launch_mi"]
    flight["range_distance_traveled_km"] = lastrow["distance_from_launch_km"]
    jsondata["range_distance_traveled_km"] = lastrow["distance_from_launch_km"]

    # launch location
    flight["launch_location"] = { "latitude": firstrow["latitude"], "longitude": firstrow["longitude"], "altitude_ft": firstrow["altitude_ft"], "altitude_m": firstrow["altitude_m"] }
    jsondata["launch_location"] = { "latitude": firstrow["latitude"], "longitude": firstrow["longitude"], "altitude_ft": firstrow["altitude_ft"], "altitude_m": firstrow["altitude_m"] }

    # landing location
    flight["landing_location"] = { "latitude": lastrow["latitude"], "longitude": lastrow["longitude"], "distance_from_launch_mi": lastrow["distance_from_launch_mi"], "distance_from_launch_km": lastrow["distance_from_launch_km"], "altitude_ft": lastrow["altitude_ft"], "altitude_m": lastrow["altitude_m"] }
    jsondata["landing_location"] = { "latitude": lastrow["latitude"], "longitude": lastrow["longitude"], "distance_from_launch_mi": lastrow["distance_from_launch_mi"], "distance_from_launch_km": lastrow["distance_from_launch_km"], "altitude_ft": lastrow["altitude_ft"], "altitude_m": lastrow["altitude_m"] }

    # save the pandas data itself to a "pickle" file in the web directory
    consolidated.to_pickle("output/pkl/" + flightname.lower() + ".pkl")

    # save the data to a csv file
    consolidated.to_csv("output/csv/" + flightname.lower() + ".csv", index=False)

    # flight metadata
    metadata = pd.json_normalize(flight);

    # save the data to a json file
    jsondata["packets"] = json.loads(consolidated.to_json(orient='records', date_unit='ms', date_format='iso', index=False))
    jsonfile = "output/json/" + flightname.lower() + ".json"
    with open(jsonfile, "w") as f:
        f.write(json.dumps(jsondata))


    # output to KML
    createKML(flightname.upper(), a, d)


    # create an individual plot and save the png file
    if not a.empty and not d.empty and deg:
        createPlot("output/png/" + flightname.lower() + ".png", a, d, deg)

    # save the metadata to a csv file
    metadata.to_csv("output/csv/" + flightname.lower() + "_metadata.csv", index=False)

    # save the output to an excel spreadsheet with seperate tabs for ascent and descent
    if not a.empty or not d.empty:
        with pd.ExcelWriter("output/xlsx/" + flightname.lower() + ".xlsx", engine='xlsxwriter', datetime_format='mm/dd/yy hh:mm:ss.000') as writer:
            metadata.to_excel(writer, index=False, sheet_name='Metadata')
            if not a.empty:
                a.to_excel(writer, index=False, sheet_name='Ascent')
            if not d.empty:
                d.to_excel(writer, index=False, sheet_name='Descent')


##########################################
# create KML files for the provided ascent & decent portions of this flight
##########################################
def createKML(flightname: str, ascent: pd.DataFrame, descent: pd.DataFrame) -> None:
    """create a KML file that can be used with Google Earth to show a 3D plot of the flight's path along with interesting data points along the way"""

    # the launch and landing points
    launch = ascent.iloc[0]
    landing = descent.iloc[-1]

    # Get the last datapoint during the ascent, that should be the approximate burst point
    burst = ascent.iloc[-1]

    # flight duration
    flighttime_seconds = (descent["packettime"].max() - ascent["packettime"].min()) / pd.Timedelta(seconds=1)
    flighttime_hrs = int(flighttime_seconds / 3600)
    flighttime_mins = int((flighttime_seconds - flighttime_hrs*3600) / 60)
    flighttime_secs = int(flighttime_seconds - flighttime_hrs*3600 - flighttime_mins*60)
    flighttime = f"{flighttime_hrs}hrs {flighttime_mins}mins {flighttime_secs}secs"


    #create a new kml object
    kml = simplekml.Kml()

    # create a new document object
    document = kml.newdocument(name=f"{flightname}")
    desc = f"""
    <h1>{flightname}</h1>
    <table style="width: 100%;" cellpadding=0 cellspacing=0 border=0>
    <tr><td style="border: solid 1px black;border-bottom: 0; padding: 10px;"><strong>Launch Local Date/Time</strong></td><td style=" padding: 10px;border: solid 1px black;border-bottom: 0; border-left: 0;">{launch['packettime']}</td></tr>
    <tr><td style="border: solid 1px black;border-bottom: 0; padding: 10px;"><strong>Flight Duration</strong></td><td style=" padding: 10px;border: solid 1px black;border-bottom: 0; border-left: 0; ">{flighttime}</td></tr>
    <tr><td style="border: solid 1px black;border-bottom: 0; padding: 10px;"><strong>Down Range Distance</strong></td><td style=" padding: 10px;border: solid 1px black;border-bottom: 0; border-left: 0; ">{landing['distance_from_launch_mi'].round(1):,}mi ({landing['distance_from_launch_km'].round(1):,}km)</td></tr>
    <tr><td style="border: solid 1px black; padding: 10px;"><strong>Approx. Burst Altitude</strong></td><td style=" padding: 10px;border: solid 1px black; border-left: 0; ">{burst['altitude_ft'].round(0):,.0f}ft ({burst['altitude_m'].round(1):,}m)</td></tr>
    </table>
    """
    document.description = f"<![CDATA[{desc}]]>"

    # folders for the various features
    paths = document.newfolder(name=f"Paths")
    waypoints = document.newfolder(name=f"Waypoints")
    poi = document.newfolder(name=f"Points of Interest")


    # create a list of coordinates from a dataframe
    coords = lambda df: list(zip(df['longitude'], df['latitude'], df['altitude_m']))



    #--------- start: ascent linestring -----------
    # create the ascent linestring
    linestring_ascent = paths.newlinestring(name=f"{flightname} Ascent")
    linestring_ascent.coords = coords(ascent)

    # ...or relativetoground
    linestring_ascent.altitudemode = simplekml.AltitudeMode.absolute

    # this will not extend the line down to the ground.
    linestring_ascent.extrude = 0

    # line color and width
    linestring_ascent.style.linestyle.color =  simplekml.Color.red
    linestring_ascent.style.linestyle.width = 3
    #--------- end: ascent linestring -----------
    


    #--------- start: descent linestring -----------
    # create the descent linestring
    linestring_descent = paths.newlinestring(name=f"{flightname} Descent")
    points = coords(descent)
    points.insert(0, (burst['longitude'], burst['latitude'], burst['altitude_m']))
    linestring_descent.coords = points

    # ...or relativetoground
    linestring_descent.altitudemode = simplekml.AltitudeMode.absolute

    # this will not extend the line down to the ground.
    linestring_descent.extrude = 0

    # line color and width
    linestring_descent.style.linestyle.color =  simplekml.Color.blue
    linestring_descent.style.linestyle.width = 3
    #--------- end: descent linestring -----------
    


    #--------- start: create launch/landing waypoints -----------

    # the launch waypoint
    launch_pnt = poi.newpoint(name=f"Launch")
    launch_pnt.coords = [(launch['longitude'], launch['latitude'], launch['altitude_m'])]
    launch_pnt.altitudemode = simplekml.AltitudeMode.absolute
    desc= f"""
    <h1>Launch</h1>
    <table style="width: 100%;" cellpadding=0 cellspacing=0 border=0>
    <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Local Time</strong></td><td style="border: solid 1px black;border-bottom: 0; border-left: 0; padding: 10px;">{launch['packettime'].strftime('%H:%M:%S')}</td></tr>
    <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Altitude</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{launch['altitude_ft'].round(0):,.0f}ft ({launch['altitude_m'].round(1):,}m)</td></tr>
    <tr><td style="padding: 10px; border: solid 1px black;"><strong>Coordinates</strong></td><td style="border: solid 1px black; border-left: 0; padding: 10px;">{launch['latitude'].round(8)}, {launch['longitude'].round(8)}</td></tr>
    </table>
    """
    launch_pnt.description = f"<![CDATA[{desc}]]>"
    launch_pnt.style.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
    launch_pnt.style.iconstyle.color = simplekml.Color.darkorange

    # the landing waypoint
    landing_pnt = poi.newpoint(name=f"Landing")
    landing_pnt.coords = [(landing['longitude'], landing['latitude'], landing['altitude_m'])]
    landing_pnt.altitudemode = simplekml.AltitudeMode.absolute
    desc = f"""
    <h1>Landing</h1>
    <table style="width: 100%;" cellpadding=0 cellspacing=0 border=0>
    <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Local Time</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{landing['packettime'].strftime('%H:%M:%S')}</td></tr>
    <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Altitude</strong></td><td style="border: solid 1px black;border-bottom: 0; border-left: 0; padding: 10px;">{landing['altitude_ft'].round(0):,.0f}ft ({landing['altitude_m'].round(1):,}m)</td></tr>
    <tr><td style="padding: 10px; border: solid 1px black;"><strong>Coordinates</strong></td><td style="border: solid 1px black; border-left: 0; padding: 10px;">{landing['latitude'].round(8)}, {landing['longitude'].round(8)}</td></tr>
    </table>
    """
    landing_pnt.description = f"<![CDATA[{desc}]]>"
    landing_pnt.style.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
    landing_pnt.style.iconstyle.color = simplekml.Color.darkorange
    #--------- end: create launch/landing waypoints -----------



    #--------- start: create burst waypoint -----------

    # the burst waypoint
    burst_pnt = poi.newpoint(name=f"Burst {burst['altitude_ft'].round(0):,.0f}ft")
    burst_pnt.coords = [(burst['longitude'], burst['latitude'], burst['altitude_m'])]
    burst_pnt.altitudemode = simplekml.AltitudeMode.absolute
    desc = f"""
    <h1>Burst: {burst['altitude_ft'].round(0):,.0f}ft</h1>
    <table style="width: 100%;" cellpadding=0 cellspacing=0 border=0>
    <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Local Time</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{burst['packettime'].strftime('%H:%M:%S')}</td></tr>
    <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Altitude</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{burst['altitude_ft'].round(0):,.0f}ft ({burst['altitude_m'].round(1):,}m)</td></tr>
    <tr><td style="padding: 10px; border: solid 1px black;"><strong>Coordinates</strong></td><td style="border: solid 1px black; border-left: 0; padding: 10px;">{burst['latitude'].round(8)}, {burst['longitude'].round(8)}</td></tr>
    </table>
    """
    burst_pnt.description = f"<![CDATA[{desc}]]>"
    burst_pnt.style.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
    burst_pnt.style.iconstyle.color = simplekml.Color.yellow


    #--------- end: create burst waypoint -----------



    #--------- start: reynolds transistion points -----------

    # get the number of Re points
    repoints = ascent[(ascent["reynolds_transition"].notnull()) & (ascent["reynolds_transition"] != "")]

    # for each reynolds transition point create a new placemarker
    for i, point in repoints.iterrows():

        # the Re waypoint
        re_name = "Turbulent-to-Laminar" if point['reynolds_transition'] == 'high_to_low' else "Laminar-to-Turbulent"
        re_pnt = poi.newpoint(name=f"{round(point['altitude_ft']/1000):,.0f}k, {re_name}")
        re_pnt.coords = [(point['longitude'], point['latitude'], point['altitude_m'])]
        re_pnt.altitudemode = simplekml.AltitudeMode.absolute
        desc = f"""
        <h1>Airflow Transitioning from {re_name}</h1>
        <table style="width: 100%;" cellpadding=0 cellspacing=0 border=0>
        <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Local Time</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{point['packettime'].strftime('%H:%M:%S')}</td></tr>
        <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Altitude</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{round(point['altitude_ft']):,.0f}ft ({round(point['altitude_m'], 1):,}m)</td></tr>
        <tr><td style="padding: 10px; border: solid 1px black;"><strong>Coordinates</strong></td><td style="border: solid 1px black; border-left: 0; padding: 10px;">{round(point['latitude'], 8)}, {round(point['longitude'], 8)}</td></tr>
        </table>
        """
        re_pnt.description = f"<![CDATA[{desc}]]>"
        re_pnt.style.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
        re_pnt.style.iconstyle.color = simplekml.Color.darkgreen
        re_pnt.style.labelstyle.color = simplekml.Color.darkgreen

    #--------- end: reynolds transistion points -----------




    # get a list of those altitudes we should target for use as waypoints
    target_altitudes = lambda df: np.arange(np.ceil(df['altitude_ft'].min() / 10000) * 10000, df['altitude_ft'].max(), 10000)
    
    # indices of those points within the ascent path that we'll use as waypoints
    ascent_waypoints_idx = []
    descent_waypoints_idx = []

    # find the nearest points to the target_altitudes within the ascent
    for target in target_altitudes(ascent):

        # find a point with an altitude_ft closest to the target point
        idx = (ascent['altitude_ft'] - target).abs().idxmin()
        ascent_waypoints_idx.append(idx)

    # find the nearest points to the target_altitudes within the descent
    for target in target_altitudes(descent):

        # find a point with an altitude_ft closest to the target point
        idx = (descent['altitude_ft'] - target).abs().idxmin()
        descent_waypoints_idx.append(idx)

    # those list of points within the ascent that we'll use as waypoints
    ascent_waypoints = ascent.loc[ascent_waypoints_idx]

    # those list of points within the ascent that we'll use as waypoints
    descent_waypoints = descent.loc[descent_waypoints_idx]

    # now add those points to the KML for the ascent
    for i, r in ascent_waypoints.iterrows():
        pnt = waypoints.newpoint(name=f"{round(r['altitude_ft'] / 1000)}k")
        pnt.coords = [(r['longitude'], r['latitude'], r['altitude_m'])]
        pnt.altitudemode = simplekml.AltitudeMode.absolute
        desc = f"""
        <h1>Waypoint</h1>
        <table style="width: 100%;" cellpadding=0 cellspacing=0 border=0>
        <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Local Time</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{r['packettime'].strftime('%H:%M:%S')}</td></tr>
        <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Altitude</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{round(r['altitude_ft']):,.0f}ft ({round(r['altitude_m'],1):,}m)</td></tr>
        <tr><td style="padding: 10px; border: solid 1px black;"><strong>Coordinates</strong></td><td style="border: solid 1px black; border-left: 0; padding: 10px;">{round(r['latitude'], 8)}, {round(r['longitude'], 8)}</td></tr>
        </table>
        """
        pnt.description = f"<![CDATA[{desc}]]>"
        pnt.style.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
        pnt.style.iconstyle.color = simplekml.Color.red


    # now add those points to the KML for the descent
    for i, r in descent_waypoints.iterrows():
        pnt = waypoints.newpoint(name=f"{round(r['altitude_ft'] / 1000)}k")
        pnt.coords = [(r['longitude'], r['latitude'], r['altitude_m'])]
        pnt.altitudemode = simplekml.AltitudeMode.absolute
        desc = f"""
        <h1>Waypoint</h1>
        <table style="width: 100%;" cellpadding=0 cellspacing=0 border=0>
        <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Local Time</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{r['packettime'].strftime('%H:%M:%S')}</td></tr>
        <tr><td style="padding: 10px; border: solid 1px black;border-bottom: 0;"><strong>Altitude</strong></td><td style="border: solid 1px black;border-left: 0; border-bottom: 0; padding: 10px;">{round(r['altitude_ft']):,.0f}ft ({round(r['altitude_m'],1):,}m)</td></tr>
        <tr><td style="padding: 10px; border: solid 1px black;"><strong>Coordinates</strong></td><td style="border: solid 1px black; border-left: 0; padding: 10px;">{round(r['latitude'], 8)}, {round(r['longitude'], 8)}</td></tr>
        </table>
        """
        pnt.description = f"<![CDATA[{desc}]]>"
        pnt.style.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
        pnt.style.iconstyle.color = simplekml.Color.blue

    #--------- end: create breadcrumb waypoints -----------


    kml.save(f"output/kml/{flightname.lower()}.kml")
     


##########################################
# read the flight list and return a list of dictionaries
##########################################
def readFlightList(filename: str)->dict:

    with open(filename, "r") as file:
        json_data = json.loads(file.read())  

    return json_data


##########################################
# main 
##########################################
def main():

    # thread list
    proclist = []

    # read in the flight list
    flightlist = readFlightList("./flightlist.json")

    try: 

        # loop the list of flights
        for flight in flightlist:

            # the flightname
            flightname = flight["flight"]

            # create a new parachute property with imperial and metric units
            chute_orig = flight["parachute"]
            chute_new = {}
            chute_new["description"] = chute_orig["description"]
            chute_new["size_ft"] = chute_orig["size"]
            chute_new["size_m"] = round(float(chute_orig["size"]) * 0.3048, 2)

            # replace the original parachute key
            flight["parachute"] = chute_new


            # the original weights property
            weights_orig = flight["weights"]
            weights_new = {}

            # loop through all of the weights key/value pairs, creating a new "weights" that contains key/value pairs for lbs and kgs
            for key, w in weights_orig.items():
                key_lb = key + "_lb"
                key_kg = key + "_kg"
                value_lb = float(w)
                value_kg = round(value_lb * 0.4535924, 2)
                weights_new[key_lb] = value_lb
                weights_new[key_kg] = value_kg

                # add the parachute weights to the parachute property as  well
                if key == "parachute":
                    flight["parachute"]["weight_lb"] = value_lb
                    flight["parachute"]["weight_kg"] = value_kg

            # replace the original weights property with the newly constructed one
            flight["weights"] = weights_new

            #if flightname == 'EOSS-361' or flightname == 'EOSS-300' or flightname == 'EOSS-374':
            #if flightname == 'EOSS-328':
            #if flightname == 'EOSS-376' or flightname == 'EOSS-321' or flightname == 'EOSS-328':
            # start a thread to process this flight's data
            p = mp.Process(name=flightname, target=processThread, args=(flight,))
            proclist.append(p)
            p.start()


    except (KeyboardInterrupt) as e:
        print("Waiting for remaining processes to finish...")
        pass

    finally:

        # join all the threads waiting for them to finish
        for p in proclist:
            p.join()


    ####
    # Now read in all the metadata info created and construct a single JSON file with that for every flight.
    ####
    print("Consolidating flight metadata...");

    # List where we'll store all the flight json
    flightjson = []

    # loop through each flight reading in the JSON just created from the threads above
    i = 0;
    for flight in flightlist:
        filename = "output/json/" + flight["flight"].lower() + ".json"
        with open(filename, "r") as file:

            # read json from file
            json_data = json.loads(file.read())

            # delete the packets key as we don't need those
            del json_data["packets"]

            # append to our list
            flightjson.append(json_data)
            i += 1;

    # write the results to a single metadata json file
    filename = "output/json/flights_metadata.json"
    with open(filename, "w") as file:
        file.write(json.dumps(flightjson))

    # create a CSV version of the flightlist metadata
    metadata = pd.json_normalize(flightjson);

    # save the metadata to a csv file
    metadata.to_csv("output/csv/flights_metadata.csv", index=False)

    # save the flight list to excel
    with pd.ExcelWriter("output/xlsx/flights_metadata.xlsx", engine='xlsxwriter', datetime_format='mm/dd/yy hh:mm:ss.000') as writer:
        metadata.to_excel(writer, index=False, sheet_name='Flightlist')

    print(f"Flights processed: {i}")

        
if __name__ == "__main__":
    main()
