
mv total_cases.csv csv_old/
mv total_deaths.csv csv_old/
mv casedistribution.csv csv_old/
mv coronavirus-cases_latest.csv csv_old/


## see https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/ for England deaths
##

## coronavirus-cases_latest.csv has England hospital cases ... consistent across the time period.

	

wget -nH https://covid.ourworldindata.org/data/ecdc/total_deaths.csv
wget -nH https://covid.ourworldindata.org/data/ecdc/total_cases.csv
wget -nH https://c19downloads.azureedge.net/downloads/csv/coronavirus-cases_latest.csv
wget -nH -O casedistribution.csv https://opendata.ecdc.europa.eu/covid19/casedistribution/csv
