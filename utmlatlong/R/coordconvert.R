######################################################################################################
# This code was adapted from the great C++ code taken originally from
# http://www.koders.com/cpp/fid56D52408FAC344874E65BF9A1C54F3731C96A39B.aspx
# http://www.mlab.cz/WebSVN/filedetails.php?repname=MLAB&path=%2FDesigns%2Fskrysohledac2%2FSW%2Futm.c
######################################################################################################

# Original headers in the C++ code
##################################
# File:			Utm.cpp
# RCS:			$Header: /cvsroot/stelvio/stelvio/NavStar/Utm.cpp,v 1.1 2001/03/18 20:07:03 steve_l Exp $
# Author:		Steve Loughran
# Created:		2001
# Language:		C++
# Package:		
# Status:		Experimental
# @doc

#	This is code to do UTM conversion.
#  I took this code from Jason Bevins' GPS thing which blagged the VB algorithms
#    from the Mapping Datum Transformation Software (MADTRAN) program,
#    written in PowerBasic.  To get the source code for MADTRAN, go to:
#
#      http://164.214.2.59/publications/guides/MADTRAN/index.html
#    
#	this version retains the core algorithms as static functions

cArray=c("C","D","E","F","G","H","J","K","L","M","N","P","Q","R","S","T","U","V","W","X")
deg2rad=pi/180
rad2deg=180/pi
fe=500000
ok=0.9996

ll2utm <- function(lat, lon=NULL, a,f,utmXZone=NULL) {
	xy=xy.coords(lat, lon)
	.LatLonToUtm(a,f,xy$x,xy$y,utmXZone)
}

# Converts a latlong coordinate to UTM coordinate in WGS84.
ll2utmwgs84 <- function(lat,lon=NULL,utmXZone=NULL) {
	return(ll2utm(lat, lon, 6378137.0, 1 / 298.257223563, utmXZone))
}

utm2ll <- function(easting, northing=NULL, a,f,utmXZone,utmYZone) {
	xy=xy.coords(easting, northing)
	.UtmToLatLon(a,f,utmXZone,utmYZone,xy$x,xy$y)
}

utm2llwgs84 <- function(easting,northing=NULL, utmXZone,utmYZone) {
	return(utm2ll(easting,northing, 6378137.0, 1 / 298.257223563, utmXZone, utmYZone))
}

	
CalculateESquared<-function(a,b) {
	return(((a * a) - (b * b)) / (a * a))
}

CalculateE2Squared <- function(a,b) {
	return(((a * a) - (b * b)) / (b * b))
}
	
denom <- function(es,sphi) {
	sinSphi = sin(sphi)
	return(sqrt(1.0 - es * (sinSphi * sinSphi)))
}

sphsr <- function (a,es,sphi) {
	dn = denom (es, sphi)
	return(a * (1.0 - es) / (dn * dn * dn))
}

sphsn <- function(a,es,sphi) {
	sinSphi = sin(sphi)
	return(a / sqrt(1.0 - es * (sinSphi * sinSphi)))
}
	
sphtmd <- function(ap,bp,cp,dp,ep,sphi) {
	return((ap * sphi) - (bp * sin(2.0 * sphi)) + (cp * sin(4.0 * sphi)) - (dp * sin(6.0 * sphi)) + (ep * sin(8.0 * sphi)))
}

# pass utmXZone if you want to force conversion within a given zone, leave null for correct conversion.
.LatLonToUtm <- function(a,f,lat,lon,utmXZone=NULL) {
	if(is.null(utmXZone)) {
		utmXZone=30 + floor(lon / 6.0)+1
	}

	if (lat < 84.0 && lat >= 72.0) {
		# Special case: zone X is 12 degrees from north to south, not 8.
		utmYZone = cArray[19+1]
	} else {
		utmYZone =cArray[floor((lat + 80.0) / 8.0)+1]
	}
	if (lat >= 84.0 || lat < -80.0) {
		# Invalid coordinate; the vertical zone is set to the invalid
		utmYZone = NA
	}
	latRad = lat * deg2rad
	lonRad = lon * deg2rad
	recf = 1 / f
	b = a * (recf - 1.0) / recf
	eSquared = CalculateESquared (a, b)
	e2Squared = CalculateE2Squared (a, b)
	tn = (a - b) / (a + b)
	tn3=tn*tn*tn
	tn4=tn3 * tn
	tn5=tn4 * tn
	ap = a * (1.0 - tn + 5.0 * ((tn * tn) - tn3) / 4.0 + 81.0 * (tn4 - tn5) / 64.0)
	bp = 3.0 * a * (tn - (tn * tn) + 7.0 * (tn3 - tn4) / 8.0 + 55.0 * tn5 / 64.0) / 2.0
	cp = 15.0 * a * ((tn * tn) - tn3 + 3.0 * (tn4 - tn5) / 4.0) / 16.0
	dp = 35.0 * a * (tn3 - tn4 + 11.0 * tn5 / 16.0) / 48.0
	ep = 315.0 * a * (tn4 - tn5) / 512.0
	olam = (utmXZone * 6 - 183) * deg2rad
	dlam = lonRad - olam
	s = sin (latRad)
	c = cos (latRad)
	t = s / c
	eta = e2Squared * (c * c)
	sn = sphsn(a, eSquared, latRad)
	tmd = sphtmd(ap, bp, cp, dp, ep, latRad)
	t1 = tmd * ok
	t2 = sn * s * c * ok / 2.0
	t3 = sn * s * (c * c * c) * ok * (5.0 - (t * t) + 9.0 * eta + 4.0 * (eta * eta)) / 24.0
	nfn=ifelse(latRad < 0.0, 10000000, 0)
	northing = (nfn + t1 + (dlam * dlam) * t2 + (dlam * dlam * dlam * dlam) * t3 + (dlam * dlam * dlam * dlam * dlam * dlam) + 0.5)
	t6 = sn * c * ok
	t7 = sn * (c * c * c) * (1.0 - (t * t) + eta) / 6.0
	easting = (fe + dlam * t6 + (dlam * dlam * dlam) * t7 + 0.5)
	northing=ifelse(northing >= 9999999, 9999999, northing)
	out=cbind(easting,northing)
	attr(out,"utmXZone")=utmXZone
	attr(out,"utmYZone")=utmYZone
	return(out)
}

.UtmToLatLon <- function(a,f,utmXZone,utmYZone,easting,northing) {
	recf = 1.0 / f
	b = a * (recf - 1) / recf
	eSquared = CalculateESquared(a, b)
	e2Squared = CalculateE2Squared(a, b)
	tn = (a - b) / (a + b)
	tn3=tn*tn*tn
	tn4=tn3*tn
	tn5=tn4*tn
	ap = a * (1.0 - tn + 5.0 * ((tn * tn) - tn3) / 4.0 + 81.0 * (tn4 - tn5) / 64.0)
	bp = 3.0 * a * (tn - (tn * tn) + 7.0 * (tn3 - tn4) / 8.0 + 55.0 * tn5 / 64.0) / 2.0
	cp = 15.0 * a * ((tn * tn) - tn3 + 3.0 * (tn4 - tn5) / 4.0) / 16.0
	dp = 35.0 * a * (tn3 - tn4 + 11.0 * tn5 / 16.0) / 48.0
	ep = 315.0 * a * (tn4 - tn5) / 512.0

	nfn=ifelse(toupper(utmYZone)<="M" & toupper(utmYZone)>="C", 10000000, 0)
	
	tmd = (northing - nfn) / ok
	sr = sphsr(a, eSquared, 0.0)
	ftphi = tmd / sr
	for(i in 0:4) {
		t10 = sphtmd(ap, bp, cp, dp, ep, ftphi)
		sr = sphsr(a, eSquared, ftphi)
		ftphi = ftphi + (tmd - t10) / sr
	}
	sr = sphsr(a, eSquared, ftphi)
	sn = sphsn(a, eSquared, ftphi)
	s = sin(ftphi)
	c = cos(ftphi)
	t = s / c
	eta = e2Squared * (c * c)
	de = easting - fe
	t10 = t / (2.0 * sr * sn * (ok * ok))
	t11 = t * (5.0 + 3.0 * (t * t) + eta - 4.0 * (eta * eta) - 9.0 * (t * t) * eta) / (24.0 * sr * (sn * sn * sn) * (ok * ok * ok * ok))
	lat = ftphi - (de * de) * t10 + (de * de * de * de) * t11
	t14 = 1.0 / (sn * c * ok)
	t15 = (1.0 + 2.0 * (t * t) + eta) / (6 * (sn * sn * sn) * c * (ok * ok * ok))
	dlam = de * t14 - (de * de * de) * t15
	olam = (utmXZone * 6 - 183.0) * deg2rad
	lon = olam + dlam
	lon = lon * rad2deg
	lat = lat * rad2deg
	return(cbind(lat,lon))
}

