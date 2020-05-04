
mv total_cases.csv csv_old/
mv total_deaths.csv csv_old/
mv casedistribution.csv csv_old/
mv coronavirus-cases_latest.csv csv_old/
	

wget -nH https://covid.ourworldindata.org/data/ecdc/total_deaths.csv
wget -nH https://covid.ourworldindata.org/data/ecdc/total_cases.csv
wget -nH https://c19downloads.azureedge.net/downloads/csv/coronavirus-cases_latest.csv
wget -nH -O casedistribution.csv https://opendata.ecdc.europa.eu/covid19/casedistribution/csv
