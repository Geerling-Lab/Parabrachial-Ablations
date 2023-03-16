import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
import os
from CompareStressTests import read_keyfile, TrialType
import glob
import datetime
import DisplayLoggedTemperature
"""
This program attempts to correlate mouse temperature change with the absolute temperature, and create a scatter plot"""

if __name__ == "__main__":

    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    root_folder = r"R:\Fillan\Parabrachial Ablations\Temperature"
    key_file = os.path.join(root_folder, "Keyfile.txt")
    mouse_numbers = read_keyfile(key_file)
    epoch_length = datetime.timedelta(minutes=1)

    ablation_epochs = []
    dead_epochs = []
    fig, ax = plt.subplots()
    for core_file in glob.glob(os.path.join(root_folder, "*", "Cold *", "Core *.csv")) + \
                     glob.glob(os.path.join(root_folder, "*", "Cervical Dislocation *", "Core *.csv")):
        mouse_number = os.path.basename(core_file).split(" ")[1][:-4]
        mouse_data = mouse_numbers[mouse_number]
        if mouse_data.genotype == "Vglut2":
            if mouse_data.vector == "Cas3":
                if mouse_number in ["5739", "5437", "5543", "5677"]:
                    pass
                else:
                    continue
            elif mouse_data.vector == "mCh":
                if mouse_number in ["4978", "5544", "5435", "5436"]:
                    pass
                else:
                    continue
        elif "Dislocation" in core_file:
            pass
        else:
            continue
        start_time = mouse_numbers[mouse_number].start[TrialType.Cold].time()
        epoch_end = None
        start_temperature = 0
        last_temperature = 0

        times, temperatures = [], []
        for line in open(core_file, 'r'):
            dt, temperature = line.strip().split(",")
            try:
                dt = datetime.datetime.strptime(dt, "%Y-%m-%d %H:%M:%S.%f")
            except ValueError:
                continue
            temperature = float(temperature)
            if "Dislocation" in core_file:
                if dt.time() < start_time:
                    continue
            else:
                if dt.time() < start_time:
                    continue
            if epoch_end is None:
                epoch_end = dt + epoch_length
                start_temperature = temperature
                continue
            if dt > epoch_end:
                if temperature < 1e-6:
                    continue
                if temperature > 17:
                    if "Dislocation" in core_file:
                        dead_epochs.append((temperature, (temperature - start_temperature) / (epoch_length.total_seconds() / 60)))
                    else:
                        ablation_epochs.append((temperature, (temperature - start_temperature) / (epoch_length.total_seconds() / 60)))
                epoch_end = epoch_end + epoch_length
                start_temperature = temperature
                continue
            if start_temperature < 1e-6:
                start_temperature = temperature
        temps = DisplayLoggedTemperature.temperatures_from_file(core_file)
        core_temperatures = DisplayLoggedTemperature.average_temperatures(temps)
        temps = []
        for time, temperature in core_temperatures:
            relative_time = (time - datetime.datetime.combine(time.date(), start_time)).total_seconds() / 3600
            if 6 > relative_time > -0.5:
                temps.append((relative_time, temperature))
        temps = np.array(temps)
        if "Dislocation" in core_file:
            ax.plot(temps[:, 0], temps[:, 1], c='k')
        else:
            if mouse_data.vector == "Cas3":
                ax.plot(temps[:, 0], temps[:, 1], c='r')
            else:
                ax.plot(temps[:, 0], temps[:, 1], c='lightgrey')
    ax.plot([0, 0], [20, 20], c='r', label='Cas3')
    ax.plot([0, 0], [20, 20], c='lightgrey', label='mCh')
    ax.plot([0, 0], [20, 20], c='k', label='Cervical Dislocation')
    ax.set_ylabel("Core Temperature (°C)")
    ax.set_xlabel("Time (Hours)")
    ax.legend(prop={'size': 8})
    plt.show()
    ablation_epochs = np.array(ablation_epochs)
    dead_epochs = np.array(dead_epochs)
    plt.scatter(ablation_epochs[:, 0], ablation_epochs[:, 1], facecolor='r', edgecolors='r', s=3, label="Cas3")
    plt.scatter(dead_epochs[:, 0], dead_epochs[:, 1], facecolor='k', edgecolors='k', s=6, label="Cervical Dislocation")
    plt.xlabel("Temperature (C)")
    plt.legend(prop={'size': 8})
    plt.ylabel("Temperature Change (ΔC/min)")
    plt.ylim([-1.5, 1.0])
    plt.xlim([17, 38])
    plt.gca().invert_xaxis()
    plt.savefig("Temperature vs Temperature Change.eps", format="eps", bbox_inches="tight")
    plt.show()