# Purpose: To read and process the wastewater data from the US to estimate the number of newly infected individuals
import pandas as pd
import numpy as np
import pandas as pd
from sodapy import Socrata

# Read fecal shedding model
shedding = pd.read_csv('FecalSheddingModel.csv', index_col=0)

# Get US data

client = Socrata("data.cdc.gov", None)
results_data = client.get("g653-rqe2", limit=200000000)

# Convert to pandas DataFrame
ww_us = pd.DataFrame.from_records(results_data)
#ww_us = pd.read_csv('NWSS_Public_SARS-CoV-2_Concentration_in_Wastewater_Data.csv')

# Get all rows for which "normalization" column has the value "flow-population"
ww_us = ww_us[ww_us['normalization'] == 'flow-population']

# Select columns of interest
ww = ww_us.copy()
ww = ww[['date', 'key_plot_id', 'pcr_conc_smoothed']]
ww.columns = ['Date', 'key_plot_id', 'gc/capita/day']

# Convert date to datetime
ww['Date'] = pd.to_datetime(ww['Date'])

# Using only dates after 2021-06-01
ww = ww[ww['Date'] > '2021-06-01']

# Sort dates from oldest to newest for each location
ww = ww.sort_values(by=['key_plot_id', 'Date'])

# Get all unique Locations
locations = ww['key_plot_id'].unique()

# Load the population data
results_pop = client.get("2ew6-ywp6", limit=200000000)

# Convert to pandas DataFrame
pop = pd.DataFrame.from_records(results_pop)
#pop = pd.read_csv('NWSS_Public_SARS-CoV-2_Wastewater_Metric_Data.csv')

# Select only the columns of interest
pop = pop[['key_plot_id', 'population_served', 'date_end']]

# Rename the columns
pop.columns = ['key_plot_id', 'Inhabitants', 'Date']

# Convert date to datetime
pop['Date'] = pd.to_datetime(pop['Date'])

# Using only dates after 2021-06-01
pop = pop[pop['Date'] > '2021-06-01']

# Remove duplicates
pop = pop.drop_duplicates()

# Add the population data to the wastewater data in an "inhabitants" column
ww = ww.merge(pop, how='left', on=['key_plot_id', 'Date'])

# Make sure that the values are numeric
ww['gc/capita/day'] = pd.to_numeric(ww['gc/capita/day'], errors='coerce')
ww['Inhabitants'] = pd.to_numeric(ww['Inhabitants'], errors='coerce')

# Convert negative values to 0
ww['gc/capita/day'] = ww['gc/capita/day'].clip(lower=0)

# Bring values to billion gene copies per day
ww['mil_gc_per_capita_per_day'] = ww['gc/capita/day'] / 1000000
ww['bil_gc_per_day'] = ww['mil_gc_per_capita_per_day'] * ww['Inhabitants']
ww['bil_gc_per_day'] = ww['bil_gc_per_day'] / 1000

# Transform dataframe such that dates are unique and every location has one column
ww2 = ww.copy()
ww = ww.pivot(index='Date', columns='key_plot_id', values='bil_gc_per_day')
ww2 = ww2.pivot(index='Date', columns='key_plot_id', values='mil_gc_per_capita_per_day')

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
    ww[loc + '_new_inf_total'] = [0] * len(ww)

# Fill all empty cells with 0
ww=ww.fillna(0)

# Fill in the rest of the new infections columns
for loc in locations:
    for i in range(13, len(ww)):
        ww.iloc[i,ww.columns.get_loc(loc + '_new_inf_total')] = (ww.iloc[i,ww.columns.get_loc(loc)]- ((ww.iloc[i-1,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[1])+(ww.iloc[i-2,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[2])+(ww.iloc[i-3,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[3])+(ww.iloc[i-4,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[4])+(ww.iloc[i-5,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[5])+(ww.iloc[i-6,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[6])+(ww.iloc[i-7,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[7])+(ww.iloc[i-8,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[8])+(ww.iloc[i-9,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[9])+(ww.iloc[i-10,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[10])+(ww.iloc[i-11,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[11])+(ww.iloc[i-12,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[12])+(ww.iloc[i-13,ww.columns.get_loc(loc + '_new_inf_total')]*shedding_list[13])))/shedding_list[0]
#print(ww.tail())
# Make new column with 3-day average of new infections by taking the last day, the present day and the next day
for loc in locations:
    ww[loc + '_new_inf_3day'] = (ww[loc + '_new_inf_total'].shift(1) + ww[loc + '_new_inf_total'] + ww[loc + '_new_inf_total'].shift(-1)) / 3

# Restructure for Tableau
df_new = ww.copy()

# Delete last row to avoid null entries
df_new = df_new[:-1]

# Get new index
df_new = df_new.reset_index()

# Select the columns we need, that is Date, all ending in '_3day_avg' and all ending in '_new_inf_3day'
df = df_new[['Date']]
for col in df_new.columns:
    if col.endswith('_3day_avg') or col.endswith('_new_inf_3day'):
        df[col] = df_new[col]


# Rename the columns by sustituting '_3day_avg' with '_wastewater' and '_new_inf_3day' with '_inf'
df.columns = df.columns.str.replace('_3day_avg', '-wastewater')
df.columns = df.columns.str.replace('_new_inf_3day', '-inf')

# "Melt" the data so each row is a unique date-region combination
df_melted = df.melt(id_vars='Date', var_name='Region_and_Measure', value_name='Value')

# Separate the Region_and_Measure column into separate 'Region' and 'Measure' columns
df_melted[['Region', 'Measure']] = df_melted['Region_and_Measure'].str.split('-', expand=True, n=1)

# Drop the now redundant 'Region_and_Measure' column
df_melted = df_melted.drop(columns='Region_and_Measure')

# To make data uniform across different files, add a 'Country' column
df_melted['Country'] = 'United_States'

# Reorder the columns
df_melted = df_melted[['Country', 'Region', 'Date', 'Measure', 'Value']]

# Change any negative values to 0
df_melted['Value'] = df_melted['Value'].clip(lower=0)

# Load the sewershed jurisdiction data
ww_us2 = pd.DataFrame.from_records(results_pop)

sewershed_data = ww_us2[['wwtp_jurisdiction', 'key_plot_id', 'population_served']].drop_duplicates()
sewershed_data['population_served'] = pd.to_numeric(sewershed_data['population_served'], errors='coerce')

# Filter the first dataframe to only include rows where "Measure" is "inf"
infection_data = df_melted[df_melted["Measure"] == "inf"]

# Join the two dataframes on the "Region" column from the first dataframe and the "rwzi_name" column from the second dataframe
combined_data = infection_data.merge(sewershed_data, left_on="Region", right_on="key_plot_id", how="inner")

# Group by "Date" and "province" and calculate the sum of "Value" multiplied by "fraction" here assumed 1 for each group
grouped_data = combined_data.groupby(["Date", "wwtp_jurisdiction"]).apply(lambda df: (df["Value"]).sum()).reset_index()

# Rename the columns to match the original dataframe
grouped_data.columns = ["Date", "Region", "Value"]

# Add the "Country" and "Measure" columns
grouped_data["Country"] = "United_States"
grouped_data["Measure"] = "inf"

# Reorder the columns to match the original dataframe
grouped_data = grouped_data[["Country", "Region", "Date", "Measure", "Value"]]

# Filter the first dataframe to only include rows where "Measure" is "wastewater"
wastewater_data = df_melted[df_melted["Measure"] == "wastewater"]

# Join the two dataframes on the "Region" column from the first dataframe and the "rwzi_name" column from the second dataframe
combined_data2 = wastewater_data.merge(sewershed_data, left_on="Region", right_on="key_plot_id", how="inner")

# For each date and province, calculate the sum of wastewater values weighted by the sewershed's population serving the province
grouped_waste = combined_data2.groupby(["Date", "wwtp_jurisdiction"]).apply(lambda df: (df["Value"] * df["population_served"]).sum()).reset_index()

# Add a column of 0s and 1s to combined_data2 to indicate whether the value is >0 or not
combined_data2["Value>0"] = np.where(combined_data2["Value"] > 0, 1, 0)

# For each date and state, calculate the sum of populations served by the sewersheds
grouped_pop = combined_data2.groupby(["Date", "wwtp_jurisdiction"]).apply(lambda df: (df["population_served"] * df["Value>0"]).sum()).reset_index()
#grouped_pop = combined_data2.groupby(["Date", "wwtp_jurisdiction"])["population_served"].sum().reset_index()

# Join the two grouped dataframes on "Date" and state
grouped_data2 = grouped_waste.merge(grouped_pop, on=["Date", "wwtp_jurisdiction"], how="inner")

# Divide the total wastewater value by the total population for each province to get the per capita wastewater value
grouped_data2["Value"] = grouped_data2["0_x"] / grouped_data2["0_y"]

# Drop the unnecessary columns
grouped_data2 = grouped_data2.drop(columns=["0_x", "0_y"])

# Rename the columns to match the original dataframe
grouped_data2.columns = ["Date", "Region", "Value"]

# Add the "Country" and "Measure" columns
grouped_data2["Country"] = "United_States"
grouped_data2["Measure"] = "wastewater"

# Reorder the columns to match the original dataframe
grouped_data2 = grouped_data2[["Country", "Region", "Date", "Measure", "Value"]]

# Append the new dataframe to the inf dataframe
#us_data = pd.concat([grouped_data, grouped_data2], ignore_index=True)

# Append the new dataframe to the original dataframe
#updated_us_data = pd.concat([grouped_data3, df_melted], ignore_index=True)

# CSV
#updated_us_data.to_csv('United_States_cleaned.csv', index=False)

# Json
#updated_us_data.to_json('United_States_cleaned.json', orient='records')

# Define the census_population_2022 data
census_population_2022 = {
    "Alabama": 5074296, "Alaska": 733583, "Arizona": 7359197, "Arkansas": 3045637, 
    "California": 39029342, "Colorado": 5839926, "Connecticut": 3626205, "Delaware": 1018396, 
    "District of Columbia": 671803, "Florida": 22244823, "Georgia": 10912876, "Hawaii": 1440196, 
    "Idaho": 1939033, "Illinois": 12582032, "Indiana": 6833037, "Iowa": 3200517, "Kansas": 2937150, 
    "Kentucky": 4512310, "Louisiana": 4590241, "Maine": 1385340, "Maryland": 6164660, "Massachusetts": 6981974, 
    "Michigan": 10034113, "Minnesota": 5717184, "Mississippi": 2940057, "Missouri": 6177957, "Montana": 1122867, 
    "Nebraska": 1967923, "Nevada": 3177772, "New Hampshire": 1395231, "New Jersey": 9261699, "New Mexico": 2113344, 
    "New York": 19677151, "North Carolina": 10698973, "North Dakota": 779261, "Ohio": 11756058, "Oklahoma": 4019800, 
    "Oregon": 4240137, "Pennsylvania": 12972008, "Rhode Island": 1093734, "South Carolina": 5282634, 
    "South Dakota": 909824, "Tennessee": 7051339, "Texas": 30029572, "Utah": 3380800, "Vermont": 647064, 
    "Virginia": 8683619, "Washington": 7785786, "West Virginia": 1775156, "Wisconsin": 5892539, "Wyoming": 581381,
    "New York City": 8335897 
}

# Calculate coverage ratio for each state on each date
grouped_pop['coverage_ratio'] = grouped_pop.apply(lambda row: row[0] / census_population_2022.get(row['wwtp_jurisdiction'], 1), axis=1)


# Merge the grouped_data_df with the coverage ratios from grouped_pop_df
adjusted_data_df = pd.merge(grouped_data, grouped_pop[['Date', 'wwtp_jurisdiction', 'coverage_ratio']], 
                            left_on=['Date', 'Region'], right_on=['Date', 'wwtp_jurisdiction'], how='left')

# Adjust the "Value" column by dividing by the "coverage_ratio"
adjusted_data_df['Adjusted Value'] = adjusted_data_df['Value'] / adjusted_data_df['coverage_ratio']

# Drop the extra columns
adjusted_data_df.drop(columns=['wwtp_jurisdiction', 'coverage_ratio', 'Value'], inplace=True)

# Rename the "Adjusted Value" column to "Value"
adjusted_data_df.rename(columns={'Adjusted Value': 'Value'}, inplace=True)

# Append the dataframes for inf and wastewater
us_data = pd.concat([adjusted_data_df, grouped_data2], ignore_index=True)

# Get official data for each state
confirmed_cases_data = pd.read_csv('covid_confirmed_usafacts.csv')

statewise_new_cases = confirmed_cases_data.copy()
date_columns = statewise_new_cases.columns[4:]  # Exclude non-date columns

# Transpose the data to have dates as rows and the combination of counties and states as columns
statewise_new_cases = statewise_new_cases.melt(id_vars=['countyFIPS', 'County Name', 'State', 'StateFIPS'], 
                                               value_vars=date_columns, 
                                               var_name='Date', 
                                               value_name='New Cases')

# Convert the date column to datetime format
statewise_new_cases['Date'] = pd.to_datetime(statewise_new_cases['Date'], format='%Y-%m-%d')

# Sort the data by state, county, and date
statewise_new_cases = statewise_new_cases.sort_values(by=['State', 'County Name', 'Date'])

# Reset the index
statewise_new_cases = statewise_new_cases.reset_index(drop=True)

# Aggregate the new cases by state
statewise_new_cases = statewise_new_cases.groupby(['State','Date'])['New Cases'].sum().reset_index()

# For each state transform cumulative to new cases in new column
statewise_new_cases['New Cases'] = statewise_new_cases.groupby('State')['New Cases'].diff().fillna(0)

# Don't allow negative values
statewise_new_cases['New Cases'] = statewise_new_cases['New Cases'].clip(lower=0)

# Apply 7-day rolling average to new cases
statewise_new_cases['New Cases MA'] = statewise_new_cases.groupby('State')['New Cases'].rolling(7).mean().reset_index(0,drop=True)

# Drop the New Cases column and rename the New Cases MA column to "Value"
statewise_new_cases = statewise_new_cases.drop(columns=['New Cases']).rename(columns={'New Cases MA': 'Value'})

# Add the "Measure" column set to "official"
statewise_new_cases['Measure'] = 'official'

# Rename the 'State' column to 'Region' to match the US_data_w_cov_ratio.csv structure
statewise_new_cases = statewise_new_cases.rename(columns={'State': 'Region'})

# Add the "Country" column set to "United_States"
statewise_new_cases['Country'] = 'United_States'

# Reorder the columns to match the US_data_w_cov_ratio.csv structure
statewise_new_cases = statewise_new_cases[['Country', 'Region', 'Date', 'Measure', 'Value']]

# Mapping of full state names to their abbreviations
state_name_to_abbreviation = {
    'Alabama': 'AL', 'Alaska': 'AK', 'Arizona': 'AZ', 'Arkansas': 'AR', 'California': 'CA', 'Colorado': 'CO',
    'Connecticut': 'CT', 'Delaware': 'DE', 'District of Columbia': 'DC', 'Florida': 'FL', 'Georgia': 'GA', 'Hawaii': 'HI', 'Idaho': 'ID',
    'Illinois': 'IL', 'Indiana': 'IN', 'Iowa': 'IA', 'Kansas': 'KS', 'Kentucky': 'KY', 'Louisiana': 'LA',
    'Maine': 'ME', 'Maryland': 'MD', 'Massachusetts': 'MA', 'Michigan': 'MI', 'Minnesota': 'MN', 'Mississippi': 'MS',
    'Missouri': 'MO', 'Montana': 'MT', 'Nebraska': 'NE', 'Nevada': 'NV', 'New Hampshire': 'NH', 'New Jersey': 'NJ',
    'New Mexico': 'NM', 'New York': 'NY', 'North Carolina': 'NC', 'North Dakota': 'ND', 'Ohio': 'OH', 'Oklahoma': 'OK',
    'Oregon': 'OR', 'Pennsylvania': 'PA', 'Rhode Island': 'RI', 'South Carolina': 'SC', 'South Dakota': 'SD',
    'Tennessee': 'TN', 'Texas': 'TX', 'Utah': 'UT', 'Vermont': 'VT', 'Virginia': 'VA', 'Washington': 'WA',
    'West Virginia': 'WV', 'Wisconsin': 'WI', 'Wyoming': 'WY'
}

# Reverse the key-value pairs in the state_name_to_abbreviation dictionary
abbreviation_to_state_name = {v: k for k, v in state_name_to_abbreviation.items()}

# Update the state abbreviations in the dataframe to use full state names
statewise_new_cases['Region'] = statewise_new_cases['Region'].replace(abbreviation_to_state_name)

# Use only dates after 2021-06-01
statewise_new_cases = statewise_new_cases[statewise_new_cases['Date'] > '2021-06-01']

# Append the official_new_cases_data to us_data_with_ratio
us_data_combined = pd.concat([us_data, statewise_new_cases], ignore_index=True)

# Save to csv
us_data_combined.to_csv('United_States_cleaned.csv', index=False)

# Save to json
us_data_combined.to_json('United_States_cleaned.json', orient='records')
