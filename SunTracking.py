from datetime import datetime, timezone
from astropy.time import Time
from astropy.coordinates import get_sun, EarthLocation, AltAz
import astropy.units as u

telescopeLongitude = 22.9638
telescopeLatitude = 40.6401
elevation = 50 * u.m

location = EarthLocation(lat=telescopeLatitude * u.deg, lon=telescopeLongitude * u.deg, height=elevation)

def automatic_tracking_better_astropy():
    time_now = Time(datetime.now(timezone.utc))  # Correctly define time

    sun_coord = get_sun(time_now)
    altaz_frame = AltAz(obstime=time_now, location=location)
    sun_altaz = sun_coord.transform_to(altaz_frame)

    print(f"Sun's Right Ascension: {sun_coord.ra.to_string(unit='hourangle', precision=5)} ({sun_coord.ra.deg:.6f}°)")
    print(f"Sun's Declination: {sun_coord.dec.deg:.6f}°")
    print(f"Height (Altitude): {sun_altaz.alt:.6f}")
    print(f"Azimuth: {sun_altaz.az:.6f}")

# Run in a loop
while True:
    automatic_tracking_better_astropy()
