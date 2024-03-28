#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import ast

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

def return_flights():
    with open("all_flights.tsv") as infile:
        lines = infile.readlines()        
        us_flights = []
        lines = lines[1:]
        for line in lines:
            origin, origin_code, date, terminal, equipment, flight, airline, nation, state, flight_time = line.split("\t")
            airlines = ast.literal_eval(airline) 
            prime_airline = airlines[0]
            if nation == "Canada":
                triturator_status = "Unknown"
            elif prime_airline in [
                "United Airlines",
                    "American Airlines"]:
                triturator_status = "American Airlines\nTriturator"

            elif prime_airline in [
                "JetBlue Airways",
                "Delta Air Lines",
                    "Southwest Airlines"]:
                triturator_status = "Swissport\nTriturator"

            elif nation != "United States":
                triturator_status = "Swissport\nTriturator"
        
            else:
                triturator_status = "Unknown"
 
            
            us_flights.append([origin, origin_code, date, terminal, equipment, flight, airline, nation, state, flight_time, prime_airline, triturator_status])

        us_flights_df = pd.DataFrame(us_flights)
        # set headers
        us_flights_df.columns = ["Origin", "Origin Code", "Date", "Terminal", "Equipment", "Flight", "Airline", "Nation", "State", "Flight Time", "Prime Airline", "Triturator Status"]

        return us_flights_df 


def return_plotting_df():
    df = return_flights()

    df["Flight Hours"] = df["Flight Time"].apply(time_to_float)


    df = (
        df.groupby(["Prime Airline", "Triturator Status"])
        .agg({"Flight Hours": "sum"})
        .reset_index()

    )

    df = df.pivot(
        index="Prime Airline", columns="Triturator Status", values="Flight Hours"
    )

    df["Total"] = df.sum(axis=1)
    df = df.sort_values(by="Total", ascending=False)
    df = df.drop("Total", axis=1)

    return df



def return_destination_trit_plot():
    df = return_plotting_df()

    df = df.head(50)
    fig, ax = plt.subplots(figsize=(8, 7))
    df.plot.barh(stacked=True, ax=ax, width=0.8)

    
    # drop title of legend box
    ax.legend().set_title("")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tick_params(axis="y", which="both", left=False, right=False)
    plt.ylabel("")
    plt.xlabel("Total Flight Hours")
    plt.title("Total Flight Hours per Airline and Triturator (Top 50)")

    plt.tight_layout()
    plt.savefig("triturator_airline_flight_hours.png", dpi=600)





def start():
    return_destination_trit_plot()




if __name__ == "__main__":
    start()
