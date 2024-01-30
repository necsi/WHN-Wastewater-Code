import requests
from datetime import datetime, timedelta
import pandas as pd

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
file_extension = "_automated_csvs/wastewater_by_county.csv"  # File extension

# Call the function to download the CSV
result, file_path = download_csv(base_url, start_date, days_back, file_extension)


# Read the CSV into a pandas DataFrame, making sure to parse the dates column
# Specify the correct date format if pandas does not recognize it automatically
df = pd.read_csv(file_path, skiprows=1, parse_dates=['date'], 
                 date_parser=lambda x: pd.to_datetime(x, format='%Y-%m-%d'))

# Clean non-numeric characters from the concentration column if necessary
# For example, if there are commas in the numbers or there are strings like '<1'
df['eff_conc_sarscov2_weekly'] = pd.to_numeric(df['eff_conc_sarscov2_weekly'].replace('[^0-9.]', '', regex=True), errors='coerce')
df['eff_conc_sarscov2_weekly_rolling'] = pd.to_numeric(df['eff_conc_sarscov2_weekly_rolling'].replace('[^0-9.]', '', regex=True), errors='coerce')

# Filter the DataFrame for the 'Nationwide' location
df.to_csv('US_Biobot_counties.csv', index=False)
