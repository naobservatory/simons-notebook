import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import random
from typing import Dict
import locale
from scipy import stats
from dataclasses import dataclass
from decimal import Decimal
from typing import Optional, List
from dataclasses import field

EPSILON = 0.000001


@dataclass(kw_only=True, eq=True)
class PathogenProperties:
    name: str = "Sars_CoV-2"
    doubling_time: int = 3
    cv_doubling_time: float = 10.0
    genome_length: int = 30000


@dataclass(kw_only=True, eq=True)
class SamplingParameters:
    """Parameters for sampling in pathogen detection"""

    shedding_values: List[float] = field(
        default_factory=lambda: [
            5e-06,
            1e-04,
            2e-04,
            3e-04,
            7e-04,
            8e-04,
            1e-03,
            6e-03,
            8e-03,
            1e-02,
            2e-02,
            3e-02,
            7e-02,
            4e-01,
            8e-01,
        ]
    )
    sigma_shedding_values: float = 0.05
    shedding_duration: int = 7
    sigma_shedding_duration: float = 0.05
    sample_population: int = 100
    sample_cost: float = 200
    low_quality: bool = False
    direct_flaggable: bool = False


@dataclass(kw_only=True, eq=True)
class SequencingParameters:
    read_length: int = 10000
    sample_depth: int = int(8e5)
    run_cost: float = 450
    processing_delay: int = 4


@dataclass(kw_only=True, eq=True)
class SamplingSequencingSchedule:
    # Sampling and sequencing schedules (True for active days)
    sampling_m: bool = True
    sampling_t: bool = True
    sampling_w: bool = True
    sampling_r: bool = True
    sampling_f: bool = True
    sampling_s: bool = False
    sampling_u: bool = False
    sequencing_m: bool = True
    sequencing_t: bool = True
    sequencing_w: bool = True
    sequencing_r: bool = True
    sequencing_f: bool = True
    sequencing_s: bool = False
    sequencing_u: bool = False


@dataclass(kw_only=True, eq=True)
class GlobalSettings:
    min_observations: int = 2
    sites: int = 1
    population_size: int = 1e10
    overhead: float = 50.0


number_of_simulations = 1000


@dataclass(kw_only=True, eq=True)
class Inputs:
    pathogen_props: PathogenProperties
    sampling_params: SamplingParameters
    sequencing_params: SequencingParameters
    schedule: SamplingSequencingSchedule
    global_settings: GlobalSettings
    number_of_simulations: int = 1000

    @classmethod
    def create_default(cls):
        return cls(
            pathogen_props=PathogenProperties(),
            sampling_params=SamplingParameters(),
            sequencing_params=SequencingParameters(),
            schedule=SamplingSequencingSchedule(),
            global_settings=GlobalSettings(),
            number_of_simulations=number_of_simulations,
        )


def get_input_cv(input_value, input_cv):
    mean = float(input_value)
    cv = float(input_cv) / 100

    if cv < EPSILON:
        return mean

    stdev = cv * mean
    return random.normalvariate(mean, stdev)


def get_input_sigma(input_value, input_sigma):
    geom_mean = float(input_value)
    sigma = float(input_sigma)

    if sigma < EPSILON:
        return geom_mean

    return math.exp(random.normalvariate(math.log(geom_mean), sigma))


def get_inputs_biased(input_values, input_sigma):
    empirical_values = [float(x) for x in input_values]
    sigma = float(input_sigma)
    if sigma < EPSILON:
        return empirical_values

    bias = math.exp(random.normalvariate(0, sigma))

    return [empirical_value * bias for empirical_value in empirical_values]


def simulate_one(inputs: Inputs):
    pathogen_props = inputs.pathogen_props
    sampling_params = inputs.sampling_params
    sequencing_params = inputs.sequencing_params
    schedule = inputs.schedule
    global_settings = inputs.global_settings
    number_of_simulations = inputs.number_of_simulations

    day = 0
    population = global_settings.population_size
    r = math.log(2) / get_input_cv(
        pathogen_props.doubling_time, pathogen_props.cv_doubling_time
    )
    growth_factor = math.exp(r)
    cumulative_incidence = 1 / population

    detectable_days = get_input_sigma(
        sampling_params.shedding_duration, sampling_params.sigma_shedding_duration
    )
    ra_sicks = get_inputs_biased(
        sampling_params.shedding_values, sampling_params.sigma_shedding_values
    )

    n_min_observations = int(global_settings.min_observations)
    observations = 0

    n_sites = int(global_settings.sites)
    site_infos = [
        {
            "sample_sick": 0,
            "sample_total": 0,
            "day_offset": random.randint(0, 6),
        }
        for _ in range(n_sites)
    ]

    bp_genome_length = int(pathogen_props.genome_length)
    n_sample_population = int(sampling_params.sample_population)
    read_length_usable = min(int(sequencing_params.read_length), bp_genome_length)

    if sampling_params.low_quality:
        read_length_usable = min(read_length_usable, 120)

    fraction_useful_reads = read_length_usable / bp_genome_length

    if sampling_params.direct_flaggable:
        fraction_useful_reads = 1

    v_processing_delay_factor = growth_factor ** float(
        sequencing_params.processing_delay
    )
    n_reads = int(sequencing_params.sample_depth)

    should_sample = [getattr(schedule, f"sampling_{day}") for day in "mtwrfsu"]
    should_sequence = [getattr(schedule, f"sequencing_{day}") for day in "mtwrfsu"]

    while True:
        day += 1
        cumulative_incidence *= growth_factor

        for site in range(n_sites):
            day_of_week = (day + site_infos[site]["day_offset"]) % 7
            if should_sample[day_of_week]:
                daily_incidence = cumulative_incidence * r
                individual_probability_sick = sum(
                    daily_incidence / (growth_factor**i)
                    for i in range(int(detectable_days))
                )
                n_sick = stats.poisson.rvs(
                    n_sample_population * individual_probability_sick
                )
                site_infos[site]["sample_sick"] += n_sick
                site_infos[site]["sample_total"] += n_sample_population

            if should_sequence[day_of_week]:
                ra_sick = 0
                if site_infos[site]["sample_sick"] == 0:
                    ra_sick = 0
                elif len(ra_sicks) == 1:
                    ra_sick = ra_sicks[0]
                elif site_infos[site]["sample_sick"] > len(ra_sicks) * 3:
                    ra_sick = sum(ra_sicks) / len(ra_sicks)
                else:
                    ra_sick = sum(
                        random.choice(ra_sicks)
                        for _ in range(site_infos[site]["sample_sick"])
                    )
                    ra_sick /= site_infos[site]["sample_sick"]

                probability_read_is_useful = (
                    site_infos[site]["sample_sick"]
                    / site_infos[site]["sample_total"]
                    * ra_sick
                    * fraction_useful_reads
                )

                site_infos[site]["sample_sick"] = 0
                site_infos[site]["sample_total"] = 0

                if probability_read_is_useful > 0:
                    observations += stats.poisson.rvs(
                        n_reads * probability_read_is_useful
                    )
                    if observations >= n_min_observations:
                        # Convert cumulative incidence to percentage
                        cumulative_incidence_percentage = cumulative_incidence * 100
                        return cumulative_incidence * v_processing_delay_factor

        if cumulative_incidence > 1 or day > 365 * 10:
            return 1


def calculate_cost(inputs: Inputs):
    sampling_params = inputs.sampling_params
    sequencing_params = inputs.sequencing_params
    schedule = inputs.schedule
    global_settings = inputs.global_settings
    locale.setlocale(locale.LC_ALL, "en_US.UTF-8")
    n_samples_weekly = sum(getattr(schedule, f"sampling_{day}") for day in "mtwrfsu")
    n_sequences_weekly = sum(
        getattr(schedule, f"sequencing_{day}") for day in "mtwrfsu"
    )

    total_cost = (
        global_settings.sites
        * (1 + (global_settings.overhead / 100))
        * 52
        * (
            n_samples_weekly * sampling_params.sample_cost
            + n_sequences_weekly * sequencing_params.run_cost
        )
    )

    formatted_cost = locale.currency(total_cost, grouping=True)
    return formatted_cost


def simulate_many(inputs: Inputs, n_simulations: int):
    results = []
    for _ in range(n_simulations):
        results.append(simulate_one(inputs))
    return results


def plot_simulation(results: List[float], label: str):
    index = np.arange(len(results)) / 10
    cumulative_incidence = sorted(results)
    median = np.median(cumulative_incidence)
    plt.plot(index, cumulative_incidence, label=label)


def run_simulation(inputs: Inputs):
    costs = calculate_cost(inputs)
    n_simulations = inputs.number_of_simulations
    doubling_time = inputs.pathogen_props.doubling_time
    sample_size = inputs.sampling_params.sample_population
    sites = inputs.global_settings.sites
    label = f"Doubling time: {doubling_time}, Sample size: {sample_size}, Sites: {sites}, Annual Cost: {costs}"
    results = simulate_many(inputs, n_simulations)
    plot_simulation(results, label)
    median = str(round(float(np.median(results)) * 100, 1)) + "%"
    # print(f"Doubling time: {doubling_time}, Sample size: {sample_size}, Sites: {sites}")
    # print(f"Cost: {costs}, Median cumulative incidence: {median}")
    return median


inputs = Inputs.create_default()


swab_ras = pd.read_csv("data/adjusted_composite_ras.tsv", sep="\t")[
    "relative_abundance"
].tolist()
inputs.sampling_params.shedding_values = swab_ras


plt.figure(figsize=(10, 6))
plt.xlabel("%")
plt.ylabel("Cumulative Incidence")
plt.ylim(0, 0.08)
plt.yticks(
    np.arange(0, 0.09, 0.01), [f"{x*100:.0f}%" for x in np.arange(0, 0.09, 0.01)]
)
plt.title(f"Cumulative Incidence across 1000 simulations")
plt.grid(True, linestyle="--", alpha=0.7)
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)

medians = []
for _ in range(10):
    median = run_simulation(inputs)
    medians.append(median)

# plt.legend()
# plt.show()
print(medians)
 1.5% ,  1.4% ,  1.4% ,  1.6% ,  1.4% ,  1.4% ,  1.5% ,  1.5% ,  1.4% ,  1.4% 