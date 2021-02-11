temp <- read.csv2("data/SIG/10.18751-Climate-Timeseries-CHTM-1.0-swiss.txt")


library("jsonlite")
json_data <- fromJSON(file="data/SIG/ch.meteoschweiz.messnetz-klima_fr.json")
test <- fromJSON("data/SIG/ch.meteoschweiz.messnetz-klima_fr.json", flatten=TRUE)
head(json_data)

library("RNetCDF")
nc <- open.nc("data/SIG/RnormM8110_ch02.lonlat_000001010000_000012010000.nc")
