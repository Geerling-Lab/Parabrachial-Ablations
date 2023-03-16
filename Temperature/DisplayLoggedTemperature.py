import datetime
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import glob
import os
import collections


def normal(x, o):
    return np.exp(-(x/o) ** 2 / 2)


def temperatures_from_file(file_path):
    temperatures = []
    with open(file_path) as f:
        for line in f:
            dt, temperature = line.strip().split(",")
            try:
                dt = datetime.datetime.strptime(dt, "%Y-%m-%d %H:%M:%S.%f")
            except ValueError:
                continue
            temperature = float(temperature)
            if temperature > 1e-6:
                temperatures.append((dt, temperature))
    return temperatures


def average_temperatures(temperatures):
    averaged_temperatures = []
    max_index = 0
    reached_end = False
    max_time = datetime.timedelta(seconds=30)
    sigma = 30
    normals = {0: 1}
    for i in range(60):
        normals[i] = normal(i, sigma)
    small_temps = collections.deque()  # Contains temperatures up to one minute on each side
    small_dt = collections.deque()
    for dt, temperature in temperatures:
        if not small_dt:
            small_dt.append(temperatures[max_index][0])
            small_temps.append(temperatures[max_index][1])
            max_index += 1
        while dt - small_dt[0] > max_time:
            small_dt.popleft()
            small_temps.popleft()
            if not small_dt:
                break
        if not reached_end:
            while temperatures[max_index][0] - dt < max_time:
                small_dt.append(temperatures[max_index][0])
                small_temps.append(temperatures[max_index][1])
                max_index += 1
                if max_index >= len(temperatures):
                    reached_end = True
                    break
        small_weights = []
        for sdt in small_dt:
            if sdt < dt:
                small_weights.append(normals[(dt - sdt).seconds])
            else:
                small_weights.append(normals[(sdt - dt).seconds])
        averaged_temperatures.append((dt, np.average(small_temps, weights=small_weights)))
    return averaged_temperatures


if __name__ == "__main__":
    file_directory = r"R:\Fillan\Parabrachial Ablations\Temperature\5864 5865\3 mg kg 6 27 22"
    temperatures = []
    for file in glob.glob(os.path.join(file_directory, "*.csv")):
        core_temperatures = np.array(average_temperatures(temperatures_from_file(file)))
        plt.plot(core_temperatures[:, 0], core_temperatures[:, 1], label=file)
    plt.plot([datetime.datetime(year=2022, month=6, day=27, hour=13, minute=52), datetime.datetime(year=2022, month=6, day=27, hour=13, minute=52)], [20, 40], 'k', label="CNO 3mg/kg")
    plt.legend()
    plt.show()