// The earth is approximately spherical, but is a bit wider at the equator and shorter pole-to-pole.
// It deviates from a sphere by only about 0.33%

//////////////////////////////////////////
// Incorporating suncalc.js, removing Moon

/*
 (c) 2011-2015, Vladimir Agafonkin
 SunCalc is a JavaScript library for calculating sun/moon position and light phases.
 https://github.com/mourner/suncalc
*/

(function () { 'use strict';

// shortcuts for easier to read formulas

var PI   = Math.PI,
    sin  = Math.sin,
    cos  = Math.cos,
    tan  = Math.tan,
    asin = Math.asin,
    atan = Math.atan2,
    acos = Math.acos,
    rad  = PI / 180;

// sun calculations are based on http://aa.quae.nl/en/reken/zonpositie.html formulas


// date/time constants and conversions

var dayMs = 1000 * 60 * 60 * 24,
    J1970 = 2440588,
    J2000 = 2451545;

function toJulian(date) { return date.valueOf() / dayMs - 0.5 + J1970; }
function fromJulian(j)  { return new Date((j + 0.5 - J1970) * dayMs); }
function toDays(date)   { return toJulian(date) - J2000; }


// general calculations for position

var e = rad * 23.4397; // obliquity of the Earth

function rightAscension(l, b) { return atan(sin(l) * cos(e) - tan(b) * sin(e), cos(l)); }
function declination(l, b)    { return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l)); }

function azimuth(H, phi, dec)  { return atan(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi)); }
function altitude(H, phi, dec) { return asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H)); }

function siderealTime(d, lw) { return rad * (280.16 + 360.9856235 * d) - lw; }


// general sun calculations

function solarMeanAnomaly(d) { return rad * (357.5291 + 0.98560028 * d); }

function eclipticLongitude(M) {

    var C = rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M)), // equation of center
        P = rad * 102.9372; // perihelion of the Earth

    return M + C + P + PI;
}

function sunCoords(d) {

    var M = solarMeanAnomaly(d),
        L = eclipticLongitude(M);

    return {
        dec: declination(L, 0),
        ra: rightAscension(L, 0)
    };
}


var SunCalc = {};


// calculates sun position for a given date and latitude/longitude

SunCalc.getPosition = function (date, lat, lng) {

    var lw  = rad * -lng,
        phi = rad * lat,
        d   = toDays(date),

        c  = sunCoords(d),
        H  = siderealTime(d, lw) - c.ra;

    return {
        azimuth: azimuth(H, phi, c.dec),
        altitude: altitude(H, phi, c.dec)
    };
};


// sun times configuration (angle, morning name, evening name)

var times = SunCalc.times = [
    [-0.833, 'sunrise',       'sunset'      ]
];

// adds a custom time to the times config

SunCalc.addTime = function (angle, riseName, setName) {
    times.push([angle, riseName, setName]);
};


// calculations for sun times

var J0 = 0.0009;

function julianCycle(d, lw) { return Math.round(d - J0 - lw / (2 * PI)); }

function approxTransit(Ht, lw, n) { return J0 + (Ht + lw) / (2 * PI) + n; }
function solarTransitJ(ds, M, L)  { return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L); }

function hourAngle(h, phi, d) { return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d))); }

// returns set time for the given sun altitude
function getSetJ(h, lw, phi, dec, n, M, L) {

    var w = hourAngle(h, phi, dec),
        a = approxTransit(w, lw, n);
    return solarTransitJ(a, M, L);
}


// calculates sun times for a given date and latitude/longitude

SunCalc.getTimes = function (date, lat, lng) {

    var lw = rad * -lng,
        phi = rad * lat,

        d = toDays(date),
        n = julianCycle(d, lw),
        ds = approxTransit(0, lw, n),

        M = solarMeanAnomaly(ds),
        L = eclipticLongitude(M),
        dec = declination(L, 0),

        Jnoon = solarTransitJ(ds, M, L),

        i, len, time, Jset, Jrise;


    var result = {
        solarNoon: fromJulian(Jnoon),
        nadir: fromJulian(Jnoon - 0.5)
    };

    for (i = 0, len = times.length; i < len; i += 1) {
        time = times[i];

        Jset = getSetJ(time[0] * rad, lw, phi, dec, n, M, L);
        Jrise = Jnoon - (Jset - Jnoon);

        result[time[1]] = fromJulian(Jrise);
        result[time[2]] = fromJulian(Jset);
    }

    return result;
};

// export as AMD module / Node module / browser variable
if (typeof define === 'function' && define.amd) define(SunCalc);
else if (typeof module !== 'undefined') module.exports = SunCalc;
else window.SunCalc = SunCalc;

}());
      
////////////////////////

function getSunRange(date, lat, lon) {
    if (date.getUTCHours() != 0 || date.getUTCMinutes() != 0) {
        alert('illegal date arg to getSunRange');
        return;
    }
    var times = SunCalc.getTimes(date, lat, lon);
    
    if (isNaN(times.sunrise.getTime())) {
        // No sunrise/sunset;  we must be close to a pole and always sunny or always dark
        // Report alwaysSunny or alwaysDark based on whether the sun is above or below the horizon
        var altitude = SunCalc.getPosition(times.solarNoon, lat, lon).altitude;
        times[altitude > 0 ? 'alwaysSunny' : 'alwaysDark'] = true;
    } else {
        times.sunrise = (times.sunrise.getTime() - date.getTime())/3600000;
        times.sunset = (times.sunset.getTime() - date.getTime())/3600000;
    }
    return times;
}

////////////////////////////////
var HimawariProj = function() {
};

HimawariProj.prototype.a = 42164; // altitude of satellite, in kilometers from center of earth
HimawariProj.prototype.earthWidth = 12756; // kilometers diameters, measured from pole to pole
HimawariProj.prototype.earthHeight = 12714; // kilometers diameter at equator
HimawariProj.prototype.r = HimawariProj.prototype.earthWidth / 2; // use width for radius
HimawariProj.prototype.lon = 140.5231; // from http://www.satellite-calculations.com/Satellite/Catalog/catalogID.php?40267
HimawariProj.prototype.xCenter = 550;
HimawariProj.prototype.yCenter = 550;
HimawariProj.prototype.focalLengthPixels = 3605;


// Angle from center of earth to edge of earth at equator, as viewed from camera
HimawariProj.prototype.computeEarthHorizontalAngle = function() {
    return Math.asin(this.r / this.a);
}

HimawariProj.prototype.vecToPixel = function(vec) {
    var scale = this.focalLengthPixels / vec[2];
    return [vec[0] * scale + this.xCenter,
            vec[1] * scale * this.earthHeight / this.earthWidth + this.yCenter];
}
    
HimawariProj.prototype.pixelToVec = function(pixel) {
    var vec = [pixel[0] - this.xCenter,
               (pixel[1] - this.yCenter) * this.earthWidth / this.earthHeight,
               this.focalLengthPixels];
    return vec;
}

// returns pixel location of edge of earth (horizon) in the direction of angle (radians)
HimawariProj.prototype.pixelAtAngle = function(angle) {
    var earthAngle = this.computeEarthHorizontalAngle();
    var pixelVec = [
        Math.cos(angle) * Math.sin(earthAngle),
        Math.sin(angle) * Math.sin(earthAngle),
        1
    ];
    return this.vecToPixel(pixelVec);
}

HimawariProj.prototype.degToRad = function(deg) {
    return deg * Math.PI / 180;
}

HimawariProj.prototype.radToDeg = function(rad) {
    return rad * 180 / Math.PI;
}

HimawariProj.prototype.latLonToVec = function(latLon) {
    var lat = this.degToRad(latLon[0]);
    var lon = this.degToRad(latLon[1] - this.lon);
    var ecef = [this.r * Math.sin(lon) * Math.cos(lat),
                -this.r * Math.sin(lat),
                this.r * Math.cos(lon) * Math.cos(lat)];
    return [ecef[0], ecef[1], this.a - ecef[2]];
}

HimawariProj.prototype.testLatLon = function(latLon) {
    var pixel = this.latLonToPixel(latLon);
    var newLatLon = this.pixelToLatLon(pixel);
}

function square(x) { return x*x; }

HimawariProj.prototype.vecToLatLon = function(vec) {
    // Figure out distance from camera to sphere
    var mag = Math.sqrt(square(vec[0]) + square(vec[1]) + square(vec[2]));
    var L = [vec[0] / mag, vec[1] / mag, vec[2] / mag];
    
    var d = (L[2] * this.a) - Math.sqrt(square(L[2] * this.a) - square(this.a) + square(this.r));
    if (isNaN(d)) {
        return false;
    }
        
    // Find location on sphere
    var ecef = [L[0] * d, L[1] * d, this.a - L[2] * d];

    var lat = -Math.atan2(ecef[1], Math.sqrt(square(ecef[0]) + square(ecef[2])));
    var lon = Math.atan2(ecef[0], ecef[2]);
    return [this.radToDeg(lat), this.radToDeg(lon) + this.lon];
}

HimawariProj.prototype.latLonToPixel = function(latLon) {
    return this.vecToPixel(this.latLonToVec(latLon));
}

HimawariProj.prototype.pixelToLatLon = function(pixel) {
    return this.vecToLatLon(this.pixelToVec(pixel));
}

HimawariProj.prototype.pixelToSunTimes = function(date, pixel) {
    var latLon = this.pixelToLatLon(pixel);
    if (!latLon) return false;
    return getSunRange(date, latLon[0], latLon[1]);
}

// Returns false if should show all 24 hours (e.g. always sunny, or completely outside sphere)
// Otherwise, returns {
HimawariProj.prototype.rectToSunTimes = function(date, topLeft, bottomRight) {
    var samplePixels = [topLeft, [topLeft[0], bottomRight[1]], [bottomRight[0], topLeft[1]], bottomRight];
    var ret = false;
    for (var i = 0; i < samplePixels.length; i++) {
        var times = this.pixelToSunTimes(date, samplePixels[i]);
        if (times) {
            if (times.alwaysSunny) {
                return false;
            }
            if (!times.alwaysDark) {
                if (!ret) {
                    ret = times;
                } else {
                    ret.sunrise = Math.min(ret.sunrise, times.sunrise);
                    ret.sunset = Math.max(ret.sunset, times.sunset);
                }
            }
        }
    }
    return ret;
}
