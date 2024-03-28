#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import ast

def return_flights():
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
            if nation == "Canada":
                triturator_status = "Unknown"
            elif prime_airline in ["United Airlines", "American Airlines"]:
                triturator_status = "American Airlines\nTriturator"

            elif prime_airline in [
                "JetBlue Airways",
                "Delta Air Lines",
                "Southwest Airlines",
            ]:
                triturator_status = "Swissport\nTriturator"

            elif prime_airline in [  # Swissport Ground Handling
                "Porter Airlines",
                "LATAM Airlines",
                "Hawaiian Airlines",
                "BermudAir",
                "Korean Air",
                "Iberia",
                "Fly Play",
                "SAS Scandinavian Airlines",
                "Qatar Airways",  # Closest match to 'Qatar'
                "Qatar Executive",  # Also a match for 'Qatar'
                "TAP Air Portugal",
                "Turkish Airlines",
                "Hainan Airlines",
                "El Al Israel Airlines",
                "Aer Lingus",
                "ITA Airways",
                "Condor",
            ]:
                triturator_status = "Swissport\nTriturator"

            elif nation != "United States":
                triturator_status = "Swissport\nTriturator"

            else:
                triturator_status = "Unknown"

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
                    prime_airline,
                    triturator_status,
                ]
            )

        us_flights_df = pd.DataFrame(us_flights)
        # set headers
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
            "Prime Airline",
            "Triturator Status",
        ]

        return us_flights_df


def return_plotting_df():
    df = return_flights()

    nicer_trit_labels = {
        "aa_triturator": "American Airlines \nTriturator",
        "pr_triturator": "Swissport\nTriturator",
        "unknown": "Unknown",
    }
    df = df.replace({"Triturator Status": nicer_trit_labels})

    # create a df with number of flights per triturator
    df = df.groupby("Triturator Status").size().reset_index(name="Flights")
    # sort the df by number of flights 
    df = df.sort_values(by="Flights", ascending=False)
    return df


def return_destination_trit_plot():
    df = return_plotting_df()

    fig, ax = plt.subplots(figsize=(8, 3))
    df.plot(
        kind="barh",
        x="Triturator Status",
        y="Flights",
        ax=ax,
        legend=False,
    )

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tick_params(axis="y", which="both", left=False, right=False)
    plt.ylabel("")
    plt.xlabel("Total Flights")
    plt.title("Total Flighs per Triturator")

    plt.tight_layout()
    print("hello")
    plt.savefig("triturator_flights.png", dpi=600)


def start():
    return_destination_trit_plot()


if __name__ == "__main__":
    start()
