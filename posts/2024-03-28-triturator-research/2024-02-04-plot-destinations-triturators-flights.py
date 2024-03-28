#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import ast
import seaborn as sns


from collections import defaultdict

def time_to_float(time):
    time = time.strip()
    try:
        time_object = dt.datetime.strptime(time, "%H:%M:%S")
        hours = (
            time_object.hour
            + time_object.minute / 60
            + time_object.second / 3600
        )
        return hours
    except Exception as e:
        print(f"conversion of time {time} failed: {e}")
        return 0

def return_airline_trit_dict():
    airline_trit_dict = defaultdict(str)
    with open("2024-03-27_airlines_triturators.tsv") as infile:
        lines = infile.readlines()
        lines = lines[1:]
        for line in lines:

            airline, triturator, ground_handling, source = line.split("\t")
            prettier_triturator = {
                "Swissport": "Swissport\nTriturator",
                "American Airlines": "American Airlines\nTriturator",
                "Unknown": "Unknown",
                "None": "No Triturator used"
            }
            airline = airline.strip()
            triturator = triturator.strip()
            triturator = prettier_triturator[triturator]
            airline_trit_dict[airline] = triturator
    return airline_trit_dict


def return_flights():
    airline_trit_dict = return_airline_trit_dict()
    unmatched_airlines = defaultdict(int)
    with open("all_flights.tsv") as infile:
        lines = infile.readlines()
        us_flights = []
        lines = lines[1:]
        for line in lines:
            (
                origin,
                origin_code,
                date,
                terminal,
                equipment,
                flight,
                airline,
                nation,
                state,
                flight_time,
            ) = line.split("\t")
            airlines = ast.literal_eval(airline)
            prime_airline = airlines[0]
            prime_airline = prime_airline.strip()
            if nation == "Canada":
                triturator_status = "Unknown"
            elif prime_airline in airline_trit_dict:
                triturator_status = airline_trit_dict[prime_airline]
            elif prime_airline == "JetBlue": # Sometimes it's JetBlue, sometimes it's JetBlue Airways
                triturator_status = "Swissport\nTriturator"

            elif nation != "United States":
                triturator_status = "Swissport\nTriturator"
                unmatched_airlines[prime_airline] += 1


            else:
                triturator_status = "Unknown"
                unmatched_airlines[prime_airline] += 1
            us_flights.append(
                [
                    origin,
                    origin_code,
                    date,
                    terminal,
                    equipment,
                    flight,
                    airline,
                    nation,
                    state,
                    flight_time,
                    triturator_status,
                ]
            )

        us_flights_df = pd.DataFrame(us_flights)

        us_flights_df.columns = [
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
            "Triturator Status",
        ]
        print("Unmatched airlines:")
        for unmatched_airline, flights in sorted(unmatched_airlines.items(), key=lambda item: item[1], reverse=True):
            print(unmatched_airline, flights)
        return us_flights_df


def return_plotting_df():
    df = return_flights()

    df["Flight Hours"] = df["Flight Time"].apply(time_to_float)

    df["Plotting Origin"] = np.where(
        df["Nation"] == "United States",
        df["State"],
        df["Nation"],
    )

    df = (
        df.groupby(["Plotting Origin", "Triturator Status"])
        .agg({"Flight Hours": "sum"})
        .reset_index()
    )

    df = df.pivot(
        index="Plotting Origin",
        columns="Triturator Status",
        values="Flight Hours",
    )

    df["Total"] = df.sum(axis=1)
    df = df.sort_values(by="Total", ascending=False)
    df = df.drop("Total", axis=1)

    return df


def origin_to_nation_dict():
    df = return_flights()
    df["Plotting Origin"] = np.where(
        df["Nation"] == "United States",
        df["State"],
        df["Nation"],
    )
    nation_dict = df.set_index("Plotting Origin")["Nation"].to_dict()
    return nation_dict


def return_destination_trit_plot():
    df = return_plotting_df()
    nation_dict = origin_to_nation_dict()

    df = df.head(50)
    fig, ax = plt.subplots(figsize=(8, 7))
    df.plot.barh(stacked=True, ax=ax, width=0.8, color = sns.color_palette("tab10"))

    for i, label in enumerate(ax.get_yticklabels()):
        origin = label.get_text()
        if nation_dict[origin] != "United States":
            #label.set_color("darkgrey")
            label.set_weight("bold")


    ax.get_legend().set_title("")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tick_params(axis="y", which="both", left=False, right=False)
    plt.ylabel("")
    plt.xlabel("Total Flight Hours")
    plt.title("Total Flight Hours per Country/State and Triturator (Top 50)")

    plt.tight_layout()
    plt.savefig("triturator_destination_flight_hours.png", dpi=600)


def start():
    return_destination_trit_plot()


if __name__ == "__main__":
    start()
