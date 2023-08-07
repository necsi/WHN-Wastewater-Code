# Purpose: To read and process the wastewater data from the Netherlands to estimate the number of newly infected individuals
import pandas as pd
import numpy as np

# Read fecal shedding model
shedding = pd.read_csv('FecalSheddingModel.csv', index_col=0)

# Read Netherlands data from Github, separating columns by ;
ww_ne = pd.read_csv('https://data.rivm.nl/covid-19/COVID-19_rioolwaterdata.csv', sep=';')
#ww_ne = pd.read_csv('https://raw.githubusercontent.com/necsi/WHN-Wastewater-Data/main/data/Netherlands/nl_wastewater_data_test.csv')

# Select columns of interest
ww = ww_ne[['Date_measurement', 'RWZI_AWZI_code', 'RWZI_AWZI_name', 'RNA_flow_per_100000']]
ww.columns = ['Date', 'Code', 'Location', 'RNA_flow_per_100000']

# Convert date to datetime
ww['Date'] = pd.to_datetime(ww['Date'])

# Using only dates after 2022-01-01
ww = ww[ww['Date'] > '2022-01-01']

# Get all unique Locations
locations = ww['Location'].unique()

# Load the population data and assume string values in columns
pop = pd.read_csv('https://raw.githubusercontent.com/necsi/WHN-Wastewater-Code/main/Netherlands/Netherlands_people_served_by_sewershed.csv')

# Rename the columns
pop.columns = ['Code', 'Inhabitants']

# Add the population data to the wastewater data in an "inhabitants" column
ww = ww.merge(pop, how='left', left_on='Code', right_on='Code')

# Bring values to billion gene copies per day
ww['mil_gc_per_capita_per_day'] = ww['RNA_flow_per_100000'] / 100000000000
ww['bil_gc_per_day'] = ww['mil_gc_per_capita_per_day'] * ww['Inhabitants']
ww['bil_gc_per_day'] = ww['bil_gc_per_day'] / 1000

# Transform dataframe such that dates are unique and every location has one column
ww2 = ww.copy()
ww = ww.pivot(index='Date', columns='Location', values='bil_gc_per_day')
ww2 = ww2.pivot(index='Date', columns='Location', values='mil_gc_per_capita_per_day')

# Make new column for each location with the name plus "_mil_gc/cap"
for col in ww2.columns:
    ww2[col + '_mil_gc/cap'] = ww2[col]

# Delete the old columns
ww2 = ww2.drop(columns=locations)

# Add the new columns to the original dataframe
ww = ww.merge(ww2, how='left', left_on='Date', right_on='Date')

# Add missing dates and interpolate missing values using linear interpolation
ww = ww.resample('D').asfreq()
ww = ww.interpolate(method='linear', axis = 0)

# Make new columns with 7 day rolling average for each location
#for loc in locations:
#    ww[loc + '_7day_avg'] = ww[loc + '_mil_gc/cap'].rolling(window = 7,min_periods=1).mean()

# Make new columns with 3-day average by taking the last day, the present day and the next day
for loc in locations:
    ww[loc + '_3day_avg'] = (ww[loc + '_mil_gc/cap'].shift(1) + ww[loc + '_mil_gc/cap'] + ww[loc + '_mil_gc/cap'].shift(-1)) / 3


# Estimate new infections
# Make list of 14 first entries in gc in billions column of shedding
shedding_list = shedding['gc in billions'].tolist()[:14]

# Initialize new infections columns
for loc in locations:
    ww[loc + '_new_inf_total'] = 0


# Fill in the rest of the new infections columns
for loc in locations:
    for i in range(13, len(ww)):
        ww.iloc[i,ww.columns.get_loc(loc + '_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc(loc)]- ((ww.iloc[i-1,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[13])))/shedding_list[0]

# Make new column with 3-day average of new infections by taking the last day, the present day and the next day
for loc in locations:
    ww[loc + '_new_inf_3day'] = (ww[loc + '_new_inf_total'].shift(1) + ww[loc + '_new_inf_total'] + ww[loc + '_new_inf_total'].shift(-1)) / 3

# Restructure for Tableau
df_new = ww.copy()

# Delete last row to avoid null entries
df_new = df_new[:-1]

# Get new index
df_new = df_new.reset_index()

# Select the columns we need, that is Date, all ending in '_7day' and all ending in '_new_inf_3day'
df = df_new[['Date']]
for col in df_new.columns:
    if col.endswith('_3day_avg') or col.endswith('_new_inf_3day'):
        df[col] = df_new[col]

# Rename the columns by sustituting '_3day_avg' with '_wastewater' and '_new_inf_3day' with '_inf'
df.columns = df.columns.str.replace('_3day_avg', '_wastewater')
df.columns = df.columns.str.replace('_new_inf_3day', '_inf')

# "Melt" the data so each row is a unique date-region combination
df_melted = df.melt(id_vars='Date', var_name='Region_and_Measure', value_name='Value')

# Separate the Region_and_Measure column into separate 'Region' and 'Measure' columns
df_melted[['Region', 'Measure']] = df_melted['Region_and_Measure'].str.split('_', expand=True, n=1)

# Drop the now redundant 'Region_and_Measure' column
df_melted = df_melted.drop(columns='Region_and_Measure')

# To make data uniform across different files, add a 'Country' column
df_melted['Country'] = 'Netherlands'

# Reorder the columns
df_melted = df_melted[['Country', 'Region', 'Date', 'Measure', 'Value']]

# Change any negative values to 0
df_melted['Value'] = df_melted['Value'].clip(lower=0)

# Load the sewershed fraction data
sewershed_data = pd.read_csv("Sewershed area fraction covered by province.csv")

# Filter the first dataframe to only include rows where "Measure" is "inf"
infection_data = df_melted[df_melted["Measure"] == "inf"]

# Join the two dataframes on the "Region" column from the first dataframe and the "rwzi_name" column from the second dataframe
combined_data = infection_data.merge(sewershed_data, left_on="Region", right_on="rwzi_name", how="inner")

# Group by "Date" and "province" and calculate the sum of "Value" multiplied by "fraction" for each group
grouped_data = combined_data.groupby(["Date", "province"]).apply(lambda df: (df["Value"] * df["fraction"]).sum()).reset_index()

# Rename the columns to match the original dataframe
grouped_data.columns = ["Date", "Region", "Value"]

# Add the "Country" and "Measure" columns
grouped_data["Country"] = "Netherlands"
grouped_data["Measure"] = "inf"

# Reorder the columns to match the original dataframe
grouped_data = grouped_data[["Country", "Region", "Date", "Measure", "Value"]]

# Calculate the population served by the sewershed in each province
sewershed_data["pop_served"] = sewershed_data["population"] * sewershed_data["fraction"]

# Filter the first dataframe to only include rows where "Measure" is "wastewater"
wastewater_data = df_melted[df_melted["Measure"] == "wastewater"]

# Join the two dataframes on the "Region" column from the first dataframe and the "rwzi_name" column from the second dataframe
combined_data2 = wastewater_data.merge(sewershed_data, left_on="Region", right_on="rwzi_name", how="inner")

# For each date and province, calculate the sum of wastewater values weighted by the sewershed's population serving the province
grouped_waste = combined_data2.groupby(["Date", "province"]).apply(lambda df: (df["Value"] * df["pop_served"]).sum()).reset_index()

# For each date and province, calculate the sum of populations served by the sewersheds
grouped_pop = combined_data2.groupby(["Date", "province"])["pop_served"].sum().reset_index()

# Join the two grouped dataframes on "Date" and "province"
grouped_data2 = grouped_waste.merge(grouped_pop, on=["Date", "province"], how="inner")

# Divide the total wastewater value by the total population for each province to get the per capita wastewater value
grouped_data2["Value"] = grouped_data2[0] / grouped_data2["pop_served"]

# Drop the unnecessary columns
grouped_data2 = grouped_data2.drop(columns=[0, "pop_served"])

# Rename the columns to match the original dataframe
grouped_data2.columns = ["Date", "Region", "Value"]

# Add the "Country" and "Measure" columns
grouped_data2["Country"] = "Netherlands"
grouped_data2["Measure"] = "wastewater"

# Reorder the columns to match the original dataframe
grouped_data2 = grouped_data2[["Country", "Region", "Date", "Measure", "Value"]]

# Append the new dataframe to the inf dataframe
grouped_data3 = pd.concat([grouped_data, grouped_data2], ignore_index=True)

# Append the new dataframe to the original dataframe
updated_netherlands_data = pd.concat([grouped_data3, df_melted], ignore_index=True)

# CSV
updated_netherlands_data.to_csv('Netherlands_cleaned.csv', index=False)

# Json
updated_netherlands_data.to_json('Netherlands_cleaned.json', orient='records')
