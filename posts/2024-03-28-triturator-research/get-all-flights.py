#!/usr/bin/env python3
import pandas as pd
import subprocess
import os
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
from datetime import datetime, date
from zoneinfo import ZoneInfo
import re


def get_time_zones():
    time_zone_dict = {}
    with open("time_zones.csv", mode="r", encoding="utf-8") as file:
        for line in file:
            location, time_zone = line.split(",")
            time_zone_dict[location] = time_zone.strip()
    return time_zone_dict


def get_arrivals(day, month, year):
    if not os.path.exists("flight_data"):
        os.makedirs("flight_data")
    flight_data_path = f"data/flight_data/{year}-{month:02d}-{day:02d}.csv"
    if not os.path.exists(flight_data_path):
        subprocess.check_call(
            [
                "aws",
                "s3",
                "cp",
                f"s3://nao-bostraffic/Data/Arrivals/{year}-{month:02d}-{day:02d}_BOS_Arrivals.csv",
                flight_data_path,
            ]
        )
    return flight_data_path


def get_state_code_dict():
    state_code_dict = {}
    with open("data/state_code_to_name.tsv", mode="r", encoding="utf-8") as file:
        # source: https://docs.google.com/spreadsheets/d/1wU-Ibw9lOplcBMbCbfhgz3GeDx10yZ7iCQY1uHcMu88/edit#gid=0
        # skip first line which contains source
        next(file)
        csv_reader = csv.DictReader(file, delimiter="\t")
        for row in csv_reader:
            state_code_dict[row["state_code"]] = row["state_name"]
    return state_code_dict

def get_minor_airport_codes():
    minor_airport_codes = {}
    with open("data/minor_airports.tsv", mode="r", encoding="utf-8") as file:
        for line in file:
            airport, location = line.split("\t")
            minor_airport_codes[airport] = location.strip()
    return minor_airport_codes

def get_airport_codes():
    non_us_codes = defaultdict(tuple)
    us_codes = defaultdict(tuple)
    with open(
        # source: https://docs.google.com/spreadsheets/u/1/d/1eepIWOHicQsLyZsb0mSXGPTXDp3vlql-aGuy1AWJED0/htmlview#
        # TODO: Use official IATA data, this sheet has a couple of mistakes
        # that I have to account for below
        "data/Airport Codes by Country - Airport Codes List .tsv",
        mode="r",
        encoding="utf-8",
    ) as file:
        csv_reader = csv.DictReader(file, delimiter="\t")
        for row in csv_reader:
            fine_location, location, airport_code = (
                row["City"],
                row["Country "],  # note the space at the end of the key
                row["Code"],
            )
            if location == "USA":
                try:
                    if airport_code == "DCA":
                        city = "Washington"
                        state = "DC"
                    elif airport_code == "SFO":
                        city = "San Francisco"
                        state = "CA"
                    elif airport_code == "IAD":
                        city = "Washington"
                        state = "VA"
                    elif airport_code == "BWI":
                        city = "Baltimore"
                        state = "MD"
                    elif airport_code == "ATL":
                        city = "Atlanta"
                        state = "GA"
                    elif airport_code == "BUF":
                        city = "Buffalo"
                        state = "NY"
                    elif airport_code == "SJU":
                        city = "San Juan"
                        state = "PR"
                    elif airport_code == "IAG":
                        city = "Niagara Falls"
                        state = "NY"
                    elif airport_code == "TRI":
                        city = "Blountville"
                        state = "TN"
                    else:
                        bits = fine_location.split(", ")
                        city = bits[0]
                        state = bits[-1]
                        state = state.split(" ")[0]
                except:
                    continue
                us_codes[airport_code] = state
                continue
            else:
                if airport_code == "SJU":
                    city = "San Juan"
                    state = "PR"
                    us_codes[airport_code] = state
                    continue
                if "," in location:
                    country = location.split(", ")[1]
                else:
                    country = location
                non_us_codes[airport_code] = country

    return us_codes, non_us_codes


def create_all_flights_tsv():
    us_codes, non_us_codes = get_airport_codes()
    minor_airports = get_minor_airport_codes()
    state_code_dict = get_state_code_dict()
    time_zone_dict = get_time_zones()

    month_range = range(1, 13)
    day_range = range(1, 32)
    years = [2023, 2024]
    headers = [
        "Origin",
        "Origin Code",
        "Date",
        "Terminal",
        "Equipment",
        "Flight",
        "Airline",
        "Nation",
        "State",
        "Flight Time",
    ]
    missing_airport_codes = defaultdict(int)
    flight_times = defaultdict(int)
    flight_exclusions = defaultdict(int)
    included_flights = 0

    with open("all_flights.tsv", "w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        for year in years:
            for month in month_range:
                for day in day_range:
                    if month == 2 and day > 28:
                        continue
                    if month in [4, 6, 9, 11] and day == 31:
                        continue
                    if date(year, month, day) < date(2023, 4, 17):
                        continue
                    today = date.today()
                    if today < date(year, month, day):
                        break
                    try:
                        flight_data_path = get_arrivals(day, month, year)
                    except:
                        print(f"no data for {year}-{month}-{day}")
                        continue

                    with open(flight_data_path, newline="") as csvfile:
                        reader = csv.DictReader(csvfile)

                        for row in reader:
                            origin = row["Origin"]
                            airport_code = row["Origin Code"]
                            dep_time = row["Departure Time"]
                            arr_time = row["Arrival Time"]
                            dep_date = row["Departure Date"]
                            arr_date = row["Arrival Date"]
                            scheduled_arr_time = row["Scheduled Arrival Time"]
                            scheduled_arr_date = row["Scheduled Arrival Date"]
                            terminal = row["Terminal"]
                            equipment = row["Equipment"]
                            flight = row["Flight"]
                            airline = row["Airline"]
                            status = row["Status"]
                            if "DIVERTED" in status:
                                flight_exclusions["Diverted"] += 1
                                continue
                            elif "Canceled" in status:
                                flight_exclusions["Canceled"] += 1
                                continue
                            elif "En Route" in status:
                                flight_exclusions["En Route"] += 1
                                continue
                            elif "Unknown" in status:
                                if (
                                    arr_time == scheduled_arr_time
                                    and arr_date == scheduled_arr_date
                                ):
                                    pass # flight seems fine
                                else:
                                    flight_exclusions["Unknown status, irregular arrival time"] += 1
                                    continue # Requires further investigation #FIXME #BUG
                            if airport_code in us_codes:
                                location = us_codes[airport_code]
                                if location == "La":
                                    location = "LA"
                                try:
                                    state = state_code_dict[location]
                                except:
                                    state = None
                                country = "United States"
                            elif airport_code in non_us_codes:
                                location = non_us_codes[airport_code]
                                country = location
                                state = None
                            elif (
                                airport_code not in us_codes
                                and airport_code not in non_us_codes
                            ):
                                try:
                                    location = minor_airports[
                                        airport_code
                                    ]
                                    if location in state_code_dict:
                                        state = state_code_dict[location]
                                        country = "United States"
                                    else:
                                        country = location
                                        state = None
                                except:
                                    print(f"Missing airport code: {airport_code}")
                                    missing_airport_codes[airport_code] += 1
                                    flight_exclusions[
                                        "Missing Airport Code"
                                    ] += 1
                                    continue
                            #try:
                            #    arr_date = datetime.strptime(arr_date, "%Y-%m-%d")
                            #except:
                            #    arr_date = datetime.strptime(arr_date, "%B %d, %Y")
                            if not arr_time:
                                flight_exclusions[
                                    "No Arrival Time provided"
                                ] += 1
                                continue

                            try:
                                raw_departure_datetime = datetime.strptime(
                                    f"{dep_date} {dep_time}", "%Y-%m-%d %H:%M"
                                )
                            except:
                                try:
                                    raw_departure_datetime = datetime.strptime(
                                        f"{dep_date} {dep_time}", "%B %d, %Y %H:%M"
                                    )
                                except Exception as e:
                                    print(f"Error parsing departure time: {dep_date} {dep_time}")
                                    print(e)



                            if state is not None:
                                departure_time_zone = time_zone_dict[state]
                            else:
                                departure_time_zone = time_zone_dict[country]

                            tz_adjusted_departure_datetime = (
                                raw_departure_datetime.replace(
                                    tzinfo=ZoneInfo(departure_time_zone)
                                )
                            )

                            try:
                                raw_arrival_datetime = datetime.strptime(
                                    f"{arr_date} {arr_time}", "%Y-%m-%d %H:%M"
                                )
                            except:
                                try:
                                    raw_arrival_datetime = datetime.strptime(
                                        f"{arr_date} {arr_time}", "%B %d, %Y %H:%M"
                                    )
                                except Exception as e:
                                    print(f"Error parsing arrival time: {arr_date} {arr_time}")
                                    print(e)
                            tz_adjusted_arrival_datetime = (
                                raw_arrival_datetime.replace(
                                    tzinfo=ZoneInfo("America/New_York")
                                )
                            )

                            flight_time = (
                                tz_adjusted_arrival_datetime
                                - tz_adjusted_departure_datetime
                            )
                            raw_flight_time = (
                                raw_arrival_datetime - raw_departure_datetime
                            )
                            flight_hours = flight_time.total_seconds() / 3600
                            if flight_hours < 0:
                                flight_exclusions["Negative Flight Time"] += 1
                                continue
                            if flight_hours > 19:
                                flight_exclusions[
                                    "Flight Time longer than 19 hours"
                                ] += 1
                                continue
                            included_flights += 1
                            flight_hours = round(flight_hours)
                            flight_times[flight_hours] += 1
                            writer.writerow(
                                [
                                    origin,
                                    airport_code,
                                    arr_date,
                                    terminal,
                                    equipment,
                                    flight,
                                    airline,
                                    country,
                                    state,
                                    flight_time,
                                ]
                            )
    print(f"\nExcluded flights: Total {sum(flight_exclusions.values())}")
    for key, value in flight_exclusions.items():
        print(f"{key}: {value}")
    print(f"\nIncluded flights: {included_flights}")

    flight_times = dict(
        sorted(
            flight_times.items(),
            key=lambda item: item[0],
        )
    )
    print("\nFlight Time Distribution:")
    for hour, count in flight_times.items():
        print(f"{hour}: {count}")

    missing_airport_codes = dict(
        sorted(
            missing_airport_codes.items(),
            key=lambda item: item[1],
            reverse=True,
        )
    )
    print("\nMissing Airport Codes:")
    for key, value in missing_airport_codes.items():
        print(f"{key}: {value}")

def start():
    create_all_flights_tsv()


if __name__ == "__main__":
    start()
