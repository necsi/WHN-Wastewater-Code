import requests
from datetime import datetime, timedelta
# https://d1t7q96h7r5kqm.cloudfront.net/2023-12-04_automated_csvs/cases_by_census_region_nationwide.csv
def download_csv(base_url, start_date, days_back, file_extension):
    for i in range(days_back):
        # Calculate the date for the file
        date_for_file = start_date - timedelta(days=i)
        date_str = date_for_file.strftime("%Y-%m-%d")  # Adjust format if needed

        # Construct the file URL
        file_url = f"{base_url}{date_str}{file_extension}"

        # Attempt to download the file
        response = requests.get(file_url)
        if response.status_code == 200:
            # Save the file if found
            file_path = f'US_Biobot_data_{date_str}.csv'
            with open(file_path, 'wb') as file:
                file.write(response.content)
            return f"CSV file for {date_str} downloaded successfully at {file_path}.", file_path

    return "No CSV file found in the specified date range.", None

# Base URL and other parameters
base_url = "https://d1t7q96h7r5kqm.cloudfront.net/"  # Replace with actual base URL
start_date = datetime.now()
days_back = 14  # Number of days to go back from today
file_extension = "_automated_csvs/wastewater_by_census_region_nationwide.csv"  # File extension

# Call the function to download the CSV
result, file_path = download_csv(base_url, start_date, days_back, file_extension)

import pandas as pd

# Read the CSV into a pandas DataFrame, making sure to parse the first column as dates
# Specify the correct date format if pandas does not recognize it automatically
df = pd.read_csv(file_path, skiprows=2, parse_dates=['Date'], usecols=[0, 1, 3], 
                 names=['Date', 'Location', 'Concentration'], 
                 date_parser=lambda x: pd.to_datetime(x, format='%Y-%m-%d'))

# Clean non-numeric characters from the concentration column if necessary
# For example, if there are commas in the numbers or there are strings like '<1'
df['Concentration'] = pd.to_numeric(df['Concentration'].replace('[^0-9.]', '', regex=True), errors='coerce')

# Filter the DataFrame for the 'Nationwide' location
nationwide_df = df[df['Location'] == 'Nationwide'].sort_values('Date')

nationwide_df.to_csv('US_Biobot.csv', index=False)
