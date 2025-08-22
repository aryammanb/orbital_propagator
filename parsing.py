from skyfield.api import load, EarthSatellite
from datetime import timedelta
import numpy as np
import time
import serial
wait=time.sleep
import csv
from datetime import datetime

def convert_to_gpgga(lat, lon, alt, year, month, day, hour, minute, second):
    # Convert latitude and longitude to NMEA format
    lat_direction = 'N' if lat >= 0 else 'S'
    lon_direction = 'E' if lon >= 0 else 'W'
    
    # Convert latitude and longitude to degrees and minutes
    lat_deg = abs(int(lat))
    lat_min = (abs(lat) - lat_deg) * 60
    lon_deg = abs(int(lon))
    lon_min = (abs(lon) - lon_deg) * 60

    # Convert second to float and format it with two decimal places
    second = float(second)
    utc_time = f"{hour:02d}{minute:02d}{second:05.2f}"  # Keeps seconds to two decimal places
    
    # Format the GPGGA string
    gpgga_string = (f"$GPGGA,{utc_time},{lat_deg:02d}{lat_min:07.4f},{lat_direction},"
                    f"{lon_deg:03d}{lon_min:07.4f},{lon_direction},1,10,1.0,{alt:.1f},M,"
                    f"-17.00,M,,\n")
    
    return gpgga_string
    # Calculate checksum

def process_csv(input_file):
    gpgga_strings = []
    
    with open(input_file, mode='r', newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=['Year', 'Day', 'Month', 'Hour', 'Minute', 'Second', 'Latitude', 'Longitude', 'Altitude'])
        next(reader)  # Skip the header row

        for row in reader:
            year = int(row['Year'])
            month = int(row['Month'])
            day = int(row['Day'])
            hour = int(row['Hour'])
            minute = int(row['Minute'])
            # Convert second to float and format to two decimal places
            second = float(row['Second'])
            second = f"{second:.2f}"  # Ensure two decimal places
            latitude = float(row['Latitude'])
            longitude = float(row['Longitude'])
            altitude = float(row['Altitude'])
            
            gpgga_string = convert_to_gpgga(latitude, longitude, altitude, year, month, day, hour, minute, second)
            gpgga_strings.append(gpgga_string)
    
    return gpgga_strings
serial_port = 'COM15'
baud_rate = 115200
def main():
    ser= serial.Serial(serial_port, baud_rate, timeout=1)
    input_file = 'gpgga7.csv'  # Replace with your CSV file path
    gpgga_strings = process_csv(input_file)
    i=0
    for gpgga in gpgga_strings:
        ser.write((gpgga.encode()))
        wait(1)
        print(gpgga)

    ser.close()

if __name__ == "__main__":
    main()
