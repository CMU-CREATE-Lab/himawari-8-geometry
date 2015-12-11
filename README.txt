// TODO:
// pixel to lat, lon
// off-earth pixel to lat, lon
How to go from pixel to lat, lon

shape of the earth
location of the camera
pixel in the camera converts to a ray

intersect the ray with the shape of the earth


if c is 0

d = -(L dot o) +- sqrt(   (L dot o)^2 - ||o||^2 + r^2)

if r is 1

//////////////////////////

Intersecting camera pixel vector (line) with a spherical earth

Starting from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

d: distance along line from starting point
L: direction of line (unit vector)
o: origin of line
c: center of earth
r: radius of earth

d = -(L dot (o-c)) +- sqrt((L dot (o-c))^2 - ||o-c||^2 + r^2)

Substituting c (center of earth) [x=0, y=0, z=0]

d = -(L dot o) +- sqrt((L dot o)^2 - ||o||^2 + r^2)

a: altitude of satellite above center of earth
Substituting o (origin of line) [x=0, y=0, z=a]

d = -(L.z * a) +- sqrt((L.z * a)^2 - a^2 + r^2)

Taking first intersection (closest to camera)

d = -(L.z * a) - sqrt((L.z * a)^2 - a^2 + r^2)


/////////////
