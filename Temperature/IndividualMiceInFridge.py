import numpy as np
import matplotlib.pyplot as plt
import datetime
import DisplayLoggedTemperature


def find_cold_start(movement_file):
    with open(movement_file) as f:
        breaks = []
        last_time = None
        for line in f:
            time, _ = line.strip().split(",")
            time = datetime.datetime.strptime(time, "%m/%d/%y %H:%M:%S")
            if last_time:
                time_delta = time - last_time
                if time_delta.total_seconds() > 20:
                    breaks.append(time)
            last_time = time
    if len(breaks) != 1:
        print("Multiple time breaks: %s" % str(breaks))
    return breaks[0]


if __name__ == "__main__":
    core_file = r"R:\Fillan\Parabrachial Ablations\Temperature in Fridge\6194 6195 Food\6194 Core.csv"
    flir_file = r"R:\Fillan\Parabrachial Ablations\Temperature in Fridge\6194 6195 Food\6194\Output.csv"
    movement_file = r"R:\Fillan\Parabrachial Ablations\Temperature in Fridge\6194 6195 Food\6194 Movement.txt"
    environmental_file = r"C:\Users\Fillan\Downloads\New Sensor 1.csv"
    cold_start = find_cold_start(movement_file)
    temps = DisplayLoggedTemperature.temperatures_from_file(core_file)
    core_temperatures = DisplayLoggedTemperature.average_temperatures(temps)

    first_time = core_temperatures[0][0]
    new_core_temperatures = []
    for time, temperature in core_temperatures:
        new_core_temperatures.append(((time - first_time).total_seconds() / 3600, temperature))
    core_temperatures = np.array(new_core_temperatures)
    fig, ax = plt.subplots(nrows=5, sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0)
    ax[0].plot(core_temperatures[:, 0], core_temperatures[:, 1], "k-")
    ax[0].set_ylabel("Core\nTemperature")

    bat_temperatures = []
    tail_temperatures = []
    with open(flir_file) as f:
        for line in f:
            words = line.strip().split(",")
            time = " ".join(words[0].split(" ")[1:])
            time = (datetime.datetime.strptime(time, "%m %d %Y %H %M %S") - first_time).total_seconds() / 3600
            bat_temperature = float(words[1])
            bat_temperatures.append((time, bat_temperature))
            if len(words) > 2:
                if len(words[2]) > 1:
                    tail_temperature = float(words[2])
                    tail_temperatures.append((time, tail_temperature))
    bat_temperatures = np.array(bat_temperatures)
    tail_temperatures = np.array(tail_temperatures)
    ax[3].plot(bat_temperatures[:, 0], bat_temperatures[:, 1], "k-")
    ax[3].set_ylabel("BAT\nTemperature")
    ax[2].plot(tail_temperatures[:, 0], tail_temperatures[:, 1], "k-")
    ax[2].set_ylabel("Tail\nTemperature")

    with open(movement_file) as f:
        movements = []
        for line in f:
            time, movement = line.strip().split(",")
            time = (datetime.datetime.strptime(time, "%m/%d/%y %H:%M:%S") - first_time).total_seconds() / 3600
            movement = float(movement)
            movements.append((time, movement))
    movements = np.array(movements)
    ax[1].stem(movements[:, 0], movements[:, 1], use_line_collection=True, markerfmt="None", linefmt="k", basefmt="none")
    ax[1].set_yticks([])
    ax[1].set_ylabel("Movement")

    with open(environmental_file) as f:
        temperatures = []
        reading_header = True
        for line in f:
            if reading_header:
                reading_header = False
                continue
            time, temperature, humidity = line.split(",")
            time = datetime.datetime.strptime(time[1:-1], "%Y-%m-%d %H:%M")
            if time < first_time:
                continue
            if time < cold_start or time > datetime.datetime(year=2022, month=10, day=5, hour=3, minute=23):
                temperatures.append(((time - first_time).total_seconds() / 3600, 22))
            else:
                temperatures.append(((time - first_time).total_seconds() / 3600, float(temperature[1:-1])))
    temperatures = np.array(temperatures)
    ax[4].plot(temperatures[:, 0], temperatures[:, 1], "k-")
    ax[4].set_ylabel("Environmental\nTemperature")

    rewarm = datetime.datetime(year=2022, month=10, day=5, hour=3, minute=23)
    refood = datetime.datetime(year=2022, month=10, day=5, hour=9, minute=29)
    for i in range(5):
        y_min, y_max = ax[i].get_ylim()
        x = (cold_start - first_time).total_seconds() / 3600
        ax[i].plot((x, x), (y_min, y_max), "b-")
        ax[i].set_ylim([y_min, y_max])


        y_min, y_max = ax[i].get_ylim()
        x = (rewarm - first_time).total_seconds() / 3600
        ax[i].plot((x, x), (y_min, y_max), "r-")
        ax[i].set_ylim([y_min, y_max])


        y_min, y_max = ax[i].get_ylim()
        x = (refood - first_time).total_seconds() / 3600
        ax[i].plot((x, x), (y_min, y_max), "g-")
        ax[i].set_ylim([y_min, y_max])
    x_min, x_max = ax[0].get_xlim()
    ax[0].set_xlim([core_temperatures[0, 0], core_temperatures[-1, 0]])
    plt.show()