#Get and plot wastewater data from the Netherlands

import pandas as pd

df = pd.read_csv('https://raw.githubusercontent.com/necsi/WHN-Wastewater-Data/main/data/Netherlands/nl_wastewater_data.csv')

print(df.head())

#plot the data
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

#convert date to datetime   
df['Date'] = pd.to_datetime(df['Date_measurement'], format='%Y-%m-%d')

#sort by date
df = df.sort_values(by='Date')

#convert to numeric
df['RNA_flow_100k'] = pd.to_numeric(df['RNA_flow_per_100000'], errors='coerce')


#plot for all locations
for location in df['RWZI_AWZI_name'].unique():
    if location in ['Arnhem','Utrecht','Amsterdam West']:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(df[df['RWZI_AWZI_name'] == location]['Date'], df[df['RWZI_AWZI_name'] == location]['RNA_flow_100k'])
        ax.set_xlabel('Date')
        ax.set_ylabel('RNA normalized')
        ax.set_title('Virus in Wastewater in ' + location)
        plt.show()
