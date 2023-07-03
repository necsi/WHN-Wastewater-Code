import pandas as pd
import numpy as np

# Read fecal shedding model
shedding = pd.read_csv('https://raw.githubusercontent.com/necsi/WHN-Wastewater-Data/main/FecalSheddingModel.csv', index_col=0)

# Read prepared wastewater data
ww = pd.read_csv('https://raw.githubusercontent.com/necsi/WHN-Wastewater-Data/main/data/Canada/ww_BC_Canada.csv', index_col=0)


# Population served
# Metro Vancouver
p_Annacis = 1278019
p_Nw_Langley = 43976
p_Iona = 704759
p_Lions_Gate = 192404
p_Lulu = 201642
p_Metro_Vancouver = p_Annacis + p_Nw_Langley + p_Iona + p_Lions_Gate + p_Lulu
p_VCH = p_Iona + p_Lions_Gate + p_Lulu
p_FH = p_Annacis + p_Nw_Langley

# Island Health
p_Victoria = 299890
p_Nanaimo = 104273
p_Comox_Valley = 45260
p_Island_Health = p_Victoria + p_Nanaimo + p_Comox_Valley

# Interior Health
p_Kamloops = 99500
p_Kelowna = 130000
p_Penticton = 37000
p_Interior_Health = p_Kamloops + p_Kelowna + p_Penticton


# Transform into million gc / capita / day
ww['Annacis_mil_gc/cap'] = ww['Annacis'] *1000 / p_Annacis
ww['Nw. Langley_mil_gc/cap'] = ww['Nw. Langley'] *1000 / p_Nw_Langley
ww['Iona_mil_gc/cap'] = ww['Iona'] *1000 / p_Iona
ww['Lion\'s Gate_mil_gc/cap'] = ww['Lion\'s Gate'] *1000 / p_Lions_Gate
ww['Lulu_mil_gc/cap'] = ww['Lulu'] *1000 / p_Lulu

ww['Victoria_mil_gc/cap'] = ww['Victoria'] *1000 / p_Victoria
ww['Nanaimo_mil_gc/cap'] = ww['Nanaimo'] *1000 / p_Nanaimo
ww['Comox Valley_mil_gc/cap'] = ww['Comox Valley'] *1000 / p_Comox_Valley

ww['Kamloops_mil_gc/cap'] = ww['Kamloops'] *1000 / p_Kamloops
ww['Kelowna_mil_gc/cap'] = ww['Kelowna'] *1000 / p_Kelowna
ww['Penticton_mil_gc/cap'] = ww['Penticton'] *1000 / p_Penticton


# Make new columns with 7 day rolling average
ww['Annacis_7day'] = ww['Annacis_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Nw. Langley_7day'] = ww['Nw. Langley_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Iona_7day'] = ww['Iona_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Lion\'s Gate_7day'] = ww['Lion\'s Gate_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Lulu_7day'] = ww['Lulu_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()

ww['Victoria_7day'] = ww['Victoria_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Nanaimo_7day'] = ww['Nanaimo_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Comox Valley_7day'] = ww['Comox Valley_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()

ww['Kamloops_7day'] = ww['Kamloops_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Kelowna_7day'] = ww['Kelowna_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()
ww['Penticton_7day'] = ww['Penticton_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()


# Make new columns with weighted average of health authorities
ww['Metro_Vancouver_7day'] = (ww['Annacis_7day'] * p_Annacis + ww['Nw. Langley_7day'] * p_Nw_Langley + ww['Iona_7day'] * p_Iona + ww['Lion\'s Gate_7day'] * p_Lions_Gate + ww['Lulu_7day'] * p_Lulu) / p_Metro_Vancouver
ww['Island_Health_7day'] = (ww['Victoria_7day'] * p_Victoria + ww['Nanaimo_7day'] * p_Nanaimo + ww['Comox Valley_7day'] * p_Comox_Valley) / p_Island_Health
ww['Interior_Health_7day'] = (ww['Kamloops_7day'] * p_Kamloops + ww['Kelowna_7day'] * p_Kelowna + ww['Penticton_7day'] * p_Penticton) / p_Interior_Health


# Estimate new cases
# Add billion gc for all health authorities
ww['VCH'] = ww['Iona'] + ww['Lion\'s Gate'] + ww['Lulu']
ww['FH'] = ww['Annacis'] + ww['Nw. Langley']

ww['VIHA'] = ww['Victoria'] + ww['Nanaimo'] + ww['Comox Valley']

ww['IH'] = ww['Kamloops'] + ww['Kelowna'] + ww['Penticton']


# Make list of 14 first entries in gc in billions column of shedding
shedding_list = shedding['gc in billions'].tolist()[:14]

# Initialize new infections columns
ww['VCH_new_inf_total'] = 0
ww['FH_new_inf_total'] = 0
ww['VIHA_new_inf_total'] = 0
ww['IH_new_inf_total'] = 0


# Fill in the rest of the new infections columns
for i in range(13, len(ww)):
    ww.iloc[i,ww.columns.get_loc('VCH_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('VCH')]- ((ww.iloc[i-1,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('VCH_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    ww.iloc[i,ww.columns.get_loc('FH_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('FH')]- ((ww.iloc[i-1,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('FH_new_inf_total')]*shedding_list[13])))/shedding_list[0]

for i in range(231, len(ww)):
    ww.iloc[i,ww.columns.get_loc('VIHA_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('VIHA')]- ((ww.iloc[i-1,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('VIHA_new_inf_total')]*shedding_list[13])))/shedding_list[0]
for i in range(204, len(ww)):
    ww.iloc[i,ww.columns.get_loc('IH_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc('IH')]- ((ww.iloc[i-1,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc('IH_new_inf_total')]*shedding_list[13])))/shedding_list[0]
    
    
# Correct for whole Health Authority population
p_covered_VCH = 0.88141
p_covered_FH = 0.67549
p_covered_VIHA = 0.51834
p_covered_IH = 0.32467

ww['VCH_new_inf'] = ww['VCH_new_inf_total'] / p_covered_VCH
ww['FH_new_inf'] = ww['FH_new_inf_total'] / p_covered_FH
ww['VIHA_new_inf'] = ww['VIHA_new_inf_total'] / p_covered_VIHA
ww['IH_new_inf'] = ww['IH_new_inf_total'] / p_covered_IH


# Make new column with 3-day average of new infections by taking the last day, the present day and the next day
ww['VCH_new_inf_3day'] = (ww['VCH_new_inf'].shift(1) + ww['VCH_new_inf'] + ww['VCH_new_inf'].shift(-1)) / 3
ww['FH_new_inf_3day'] = (ww['FH_new_inf'].shift(1) + ww['FH_new_inf'] + ww['FH_new_inf'].shift(-1)) / 3
ww['VIHA_new_inf_3day'] = (ww['VIHA_new_inf'].shift(1) + ww['VIHA_new_inf'] + ww['VIHA_new_inf'].shift(-1)) / 3
ww['IH_new_inf_3day'] = (ww['IH_new_inf'].shift(1) + ww['IH_new_inf'] + ww['IH_new_inf'].shift(-1)) / 3


# Save 
ww.to_csv('Canada/BC_estimate_infections.csv')

# Restructure for danfo.js
df = ww.copy()

df['MetroVancouver_inf'] = df['VCH_new_inf_3day'] + df['FH_new_inf_3day']

# Get new index
df = df.reset_index()

# Select the columns we need
df = df[['Date', 'Metro_Vancouver_7day', 'Island_Health_7day', 'Interior_Health_7day', 'MetroVancouver_inf', 'VIHA_new_inf_3day', 'IH_new_inf_3day']]

# Rename the columns
df.columns = ['Date', 'MetroVancouver_wastewater', 'IslandHealth_wastewater', 'InteriorHealth_wastewater', 'MetroVancouver_inf', 'IslandHealth_inf', 'InteriorHealth_inf']

# "Melt" the data so each row is a unique date-region combination
df_melted = df.melt(id_vars='Date', var_name='Region_and_Measure', value_name='Value')

# Separate the Region_and_Measure column into separate 'Region' and 'Measure' columns
df_melted[['Region', 'Measure']] = df_melted['Region_and_Measure'].str.split('_', expand=True, n=1)

# Drop the now redundant 'Region_and_Measure' column
df_melted = df_melted.drop(columns='Region_and_Measure')

# To make data uniform across different files, add a 'Country' column
df_melted['Country'] = 'Canada'

# Reorder the columns
df_melted = df_melted[['Country', 'Region', 'Date', 'Measure', 'Value']]

# Remove all entries after today
import datetime
today = pd.Timestamp(datetime.date.today())
df_melted['Date'] = pd.to_datetime(df_melted['Date'])
df_melted = df_melted[df_melted['Date'] <= today]
df_melted.to_csv('Canada_cleaned.csv', index=False)
