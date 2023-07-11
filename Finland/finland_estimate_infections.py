# Purpose: To read and process the wastewater data from Finland to estimate the number of newly infected individuals
import pandas as pd
import numpy as np

# Read fecal shedding model
shedding = pd.read_csv('FecalSheddingModel.csv', index_col=0)

# Read Finnish data
ww_fi = pd.read_csv('https://raw.githubusercontent.com/necsi/WHN-Wastewater-Data/main/data/Finland/fi_wastewater_data.csv')

# Select columns of interest
ww = ww_fi[['Date of sample', 'Location of treatment plant', 'Normalized RNA count using the RNA standard']]
ww.columns = ['Date', 'Location', 'million_gc/1000p/day']

# Convert date to datetime
ww['Date'] = pd.to_datetime(ww['Date'])

# Using only dates after 2023-01-01
ww = ww[ww['Date'] > '2023-01-01']

# Transform dataframe such that dates are unique and every location has one column
ww = ww.pivot(index='Date', columns='Location', values='million_gc/1000p/day')

# Select only those locations for which we have weekly data
ww = ww[['Espoo', 'Helsinki', 'Joensuu', 'Jyväskylä', 'Kuopio', 'Oulu', 'Tampere', 'Turku', 'Vaasa']]

# Number of residents in each region served by the sewershed
p_Espoo = 390000
p_Helsinki = 860000
p_Joensuu = 98000
p_Jyväskylä = 154600
p_Kuopio = 90697
p_Oulu = 200000
p_Tampere = 200000
p_Turku = 300000
p_Vaasa = 69500
#p_EntireCountry = p_Espoo + p_Helsinki + p_Joensuu + p_Jyväskylä + p_Kuopio + p_Oulu + p_Tampere + p_Turku + p_Vaasa

# Bring values to billion gene copies per day
ww['Espoo'] = ww['Espoo'] * p_Espoo / 1000000
ww['Helsinki'] = ww['Helsinki'] * p_Helsinki / 1000000
ww['Joensuu'] = ww['Joensuu'] * p_Joensuu / 1000000
ww['Jyväskylä'] = ww['Jyväskylä'] * p_Jyväskylä / 1000000
ww['Kuopio'] = ww['Kuopio'] * p_Kuopio / 1000000
ww['Oulu'] = ww['Oulu'] * p_Oulu / 1000000
ww['Tampere'] = ww['Tampere'] * p_Tampere / 1000000
ww['Turku'] = ww['Turku'] * p_Turku / 1000000
ww['Vaasa'] = ww['Vaasa'] * p_Vaasa / 1000000
#ww['Entire country'] = ww['Entire country'] * p_EntireCountry / 1000000

# Add missing dates and interpolate missing values using linear interpolation
ww = ww.resample('D').asfreq()
ww = ww.interpolate(method='linear', axis = 0)

# Transform into million gc / capita / day
ww['Espoo_mil_gc/cap'] = ww['Espoo'] *1000 / p_Espoo
ww['Helsinki_mil_gc/cap'] = ww['Helsinki'] *1000 / p_Helsinki
ww['Joensuu_mil_gc/cap'] = ww['Joensuu'] *1000 / p_Joensuu
ww['Jyväskylä_mil_gc/cap'] = ww['Jyväskylä'] *1000 / p_Jyväskylä
ww['Kuopio_mil_gc/cap'] = ww['Kuopio'] *1000 / p_Kuopio
ww['Oulu_mil_gc/cap'] = ww['Oulu'] *1000 / p_Oulu
ww['Tampere_mil_gc/cap'] = ww['Tampere'] *1000 / p_Tampere
ww['Turku_mil_gc/cap'] = ww['Turku'] *1000 / p_Turku
ww['Vaasa_mil_gc/cap'] = ww['Vaasa'] *1000 / p_Vaasa

# Make new columns with 7 day rolling average
ww['Espoo_7day'] = ww['Espoo_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Helsinki_7day'] = ww['Helsinki_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Joensuu_7day'] = ww['Joensuu_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Jyväskylä_7day'] = ww['Jyväskylä_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Kuopio_7day'] = ww['Kuopio_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Oulu_7day'] = ww['Oulu_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Tampere_7day'] = ww['Tampere_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Turku_7day'] = ww['Turku_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Vaasa_7day'] = ww['Vaasa_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()

# Estimate new infections
# Make list of 14 first entries in gc in billions column of shedding
shedding_list = shedding['gc in billions'].tolist()[:14]

# Initialize new infections columns
ww['Espoo_new_inf_total'] = 0
ww['Helsinki_new_inf_total'] = 0
ww['Joensuu_new_inf_total'] = 0
ww['Jyväskylä_new_inf_total'] = 0
ww['Kuopio_new_inf_total'] = 0
ww['Oulu_new_inf_total'] = 0
ww['Tampere_new_inf_total'] = 0
ww['Turku_new_inf_total'] = 0
ww['Vaasa_new_inf_total'] = 0

# Fill in the rest of the new infections columns
for i in range(15, len(ww)):
    ww.iloc[i,ww.columns.get_loc('Espoo_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Espoo')]- ((ww.iloc[i-1,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Espoo_new_inf_total')]*shedding_list[13])))/shedding_list[0]
for i in range(13, len(ww)):
    ww.iloc[i,ww.columns.get_loc('Helsinki_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Helsinki')]- ((ww.iloc[i-1,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Helsinki_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('Joensuu_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Joensuu')]- ((ww.iloc[i-1,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Joensuu_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('Jyväskylä_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Jyväskylä')]- ((ww.iloc[i-1,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Jyväskylä_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('Kuopio_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Kuopio')]- ((ww.iloc[i-1,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Kuopio_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('Oulu_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Oulu')]- ((ww.iloc[i-1,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Oulu_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('Tampere_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Tampere')]- ((ww.iloc[i-1,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Tampere_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('Turku_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Turku')]- ((ww.iloc[i-1,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Turku_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('Vaasa_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('Vaasa')]- ((ww.iloc[i-1,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('Vaasa_new_inf_total')]*shedding_list[13])))/shedding_list[0]

# Make new column with 3-day average of new infections by taking the last day, the present day and the next day
ww['Espoo_new_inf_3day'] = (ww['Espoo_new_inf_total'].shift(1) + ww['Espoo_new_inf_total'] + ww['Espoo_new_inf_total'].shift(-1)) / 3
ww['Helsinki_new_inf_3day'] = (ww['Helsinki_new_inf_total'].shift(1) + ww['Helsinki_new_inf_total'] + ww['Helsinki_new_inf_total'].shift(-1)) / 3
ww['Joensuu_new_inf_3day'] = (ww['Joensuu_new_inf_total'].shift(1) + ww['Joensuu_new_inf_total'] + ww['Joensuu_new_inf_total'].shift(-1)) / 3
ww['Jyväskylä_new_inf_3day'] = (ww['Jyväskylä_new_inf_total'].shift(1) + ww['Jyväskylä_new_inf_total'] + ww['Jyväskylä_new_inf_total'].shift(-1)) / 3
ww['Kuopio_new_inf_3day'] = (ww['Kuopio_new_inf_total'].shift(1) + ww['Kuopio_new_inf_total'] + ww['Kuopio_new_inf_total'].shift(-1)) / 3
ww['Oulu_new_inf_3day'] = (ww['Oulu_new_inf_total'].shift(1) + ww['Oulu_new_inf_total'] + ww['Oulu_new_inf_total'].shift(-1)) / 3
ww['Tampere_new_inf_3day'] = (ww['Tampere_new_inf_total'].shift(1) + ww['Tampere_new_inf_total'] + ww['Tampere_new_inf_total'].shift(-1)) / 3
ww['Turku_new_inf_3day'] = (ww['Turku_new_inf_total'].shift(1) + ww['Turku_new_inf_total'] + ww['Turku_new_inf_total'].shift(-1)) / 3
ww['Vaasa_new_inf_3day'] = (ww['Vaasa_new_inf_total'].shift(1) + ww['Vaasa_new_inf_total'] + ww['Vaasa_new_inf_total'].shift(-1)) / 3


# Compare our estimated cases with officially estimated cases
estim_cases = pd.read_csv('https://raw.githubusercontent.com/necsi/WHN-Wastewater-Data/main/data/Finland/FinlandEstimCasesSewersheds.csv')

# Select columns of interest
estim_cases = estim_cases[['Treatment plant', 'Date of reporting', 'COVID-19 cases, estimate']]

# Rename columns
estim_cases.columns = ['Treatment_plant','Date', 'Estim_cases']

# Convert date to datetime
estim_cases['Date'] = pd.to_datetime(estim_cases['Date'])

# Using only dates after 2023-01-01
estim_cases = estim_cases[estim_cases['Date'] > '2023-01-01']

# Transform dataframe such that dates are unique and every location has one column
estim_cases = estim_cases.pivot(index='Date', columns='Treatment_plant', values='Estim_cases')

# Fill everything that is ".." with 0
estim_cases = estim_cases.replace('..', 0)

# Add missing dates
estim_cases = estim_cases.resample('D').asfreq()

# Set data entries to float
estim_cases = estim_cases.astype(float)

# Divide by 7 to get weekly average
estim_cases = estim_cases / 7

# Run trough all rows of all columns from the bottom, check if the entry is nan, and if it is, fill it with the last non-nan value
for i in range(len(estim_cases.columns)):
    for j in range(len(estim_cases)-1, 0, -1):
        if np.isnan(estim_cases.iloc[j,i]):
            estim_cases.iloc[j,i] = estim_cases.iloc[j+1,i]


# Restructure for danfo.js
df = ww.copy()

# Get new index
df = df.reset_index()

# Select the columns we need
df = df[['Date', 'Espoo_7day', 'Helsinki_7day', 'Joensuu_7day',
       'Jyväskylä_7day', 'Kuopio_7day', 'Oulu_7day', 'Tampere_7day',
       'Turku_7day', 'Vaasa_7day', 'Espoo_new_inf_3day', 'Helsinki_new_inf_3day', 'Joensuu_new_inf_3day',
       'Jyväskylä_new_inf_3day', 'Kuopio_new_inf_3day', 'Oulu_new_inf_3day',
       'Tampere_new_inf_3day', 'Turku_new_inf_3day', 'Vaasa_new_inf_3day']]

# Rename the columns
df.columns = ['Date', 'Espoo_wastewater', 'Helsinki_wastewater', 'Joensuu_wastewater', 'Jyväskylä_wastewater', 'Kuopio_wastewater', 'Oulu_wastewater', 'Tampere_wastewater', 'Turku_wastewater', 'Vaasa_wastewater', 'Espoo_inf', 'Helsinki_inf', 'Joensuu_inf', 'Jyväskylä_inf', 'Kuopio_inf', 'Oulu_inf', 'Tampere_inf', 'Turku_inf', 'Vaasa_inf']

# "Melt" the data so each row is a unique date-region combination
df_melted = df.melt(id_vars='Date', var_name='Region_and_Measure', value_name='Value')

# Separate the Region_and_Measure column into separate 'Region' and 'Measure' columns
df_melted[['Region', 'Measure']] = df_melted['Region_and_Measure'].str.split('_', expand=True, n=1)

# Drop the now redundant 'Region_and_Measure' column
df_melted = df_melted.drop(columns='Region_and_Measure')

# To make data uniform across different files, add a 'Country' column
df_melted['Country'] = 'Finland'

# Reorder the columns
df_melted = df_melted[['Country', 'Region', 'Date', 'Measure', 'Value']]


df2 = estim_cases.copy()

# Get new index
df2 = df2.reset_index()

# Select the columns we need
df2 = df2[['Date', 'Viikinmäki', 'Kuhasalo',
       'Nenäinniemi', 'Lehtoniemi', 'Taskila', 'Viinikanlahti',
       'Kakolanmäki', 'Pått']]

# Rename the columns
df2.columns = ['Date', 'Helsinki_official', 'Joensuu_official', 'Jyväskylä_official', 'Kuopio_official', 'Oulu_official', 'Tampere_official', 'Turku_official', 'Vaasa_official']

# "Melt" the data so each row is a unique date-region combination
df2_melted = df2.melt(id_vars='Date', var_name='Region_and_Measure', value_name='Value')

# Separate the Region_and_Measure column into separate 'Region' and 'Measure' columns
df2_melted[['Region', 'Measure']] = df2_melted['Region_and_Measure'].str.split('_', expand=True, n=1)

# Drop the now redundant 'Region_and_Measure' column
df2_melted = df2_melted.drop(columns='Region_and_Measure')

# To make data uniform across different files, add a 'Country' column
df2_melted['Country'] = 'Finland'

# Reorder the columns
df2_melted = df2_melted[['Country', 'Region', 'Date', 'Measure', 'Value']]


cleaned = pd.concat([df_melted, df2_melted], axis=0)

cleaned.to_csv('Finland_cleaned.csv', index=False)

# Json
cleaned.to_json('Finland_cleaned.json', orient='records')
