# Simple coordinate conversion between geographic and projected coordinates, for R
Download R package [here](https://github.com/miguel-porto/R-utmlatlong/releases/download/v1.0/utmlatlong_1.0.tar.gz) (linux)

Simple numeric coordinate conversion between latitude longitude and UTM coordinates (and vice-versa) with no dependencies.
It is not meant to be a replacement for the rgdal package functions, but only provide the simplest way to convert between
geographic and projected coordinates without any other dependencies, if all you want is just plain numbers.
	
## Using it
Install the package from local file and then:
```
library(utmlatlong)

# make some sample lat long coordinates
coo=cbind(
  c(37.3454, 36.4455, 14.566)
  ,c(-8.5467, 88.547, 34.6576)
)

# convert lat long to UTM WGS84
utm=ll2utmwgs84(coo)

# ...and back to lat long again!
utm2llwgs84(utm, utmXZone=attr(utm,"utmXZone"), utmYZone=attr(utm,"utmYZone"))
```

