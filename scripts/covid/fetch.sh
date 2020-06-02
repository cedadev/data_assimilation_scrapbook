
## https://data.humdata.org/dataset/acaps-covid19-government-measures-dataset
## wget not yet set ...
## https://data.humdata.org/dataset/e1a91ae0-292d-4434-bc75-bf863d4608ba/resource/e53c0034-a320-450f-aaeb-0e476f462464/download/acaps-_covid19_government_measures_dataset.xlsx 

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



JH_COVID=/home/mjuckes/Repositories/git/import/john_hopkins/COVID-19

(cd $JH_COVID ; git pull)
cp ${JH_COVID}/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv .	
