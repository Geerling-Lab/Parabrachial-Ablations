import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd
import datetime
import os
import DisplayLoggedTemperature
from CompareStressTests import create_relative_temperatures
from CompareStressTests import read_keyfile

"""
This program analyzes simultaneous readings of the mouse's temperature during cold exposure
The mouse is in the beaker, while being monitored in the following ways:
Brown Adipose Tissue - Flir C2, images taken with my C# program, analyzed as maximum temperature
Tail temperature - Flir C2, images taken with my C# program, manually analyzed
Core temperature - E-mitter, data taken with my TemperatureLogger program
Shivering? - Olympus E-M5 III, manually scored?
"""
if __name__ == "__main__":
    font_path = r"C:\Users\Filla\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 2

    parent_directory = r"R:\Fillan\Parabrachial Ablations\Temperature\Peipheral Temperature Cold"
    root_folder = r"R:\Fillan\Parabrachial Ablations\Temperature"
    key_file = os.path.join(root_folder, "Keyfile.txt")
    mouse_numbers = read_keyfile(key_file)
    fig, ax = plt.subplots(nrows=3, sharex=True)

    plt.subplots_adjust(wspace=0, hspace=0)
    for mouse_number in os.listdir(parent_directory):
        mouse_directory = os.path.join(parent_directory, mouse_number)
        if not os.path.isdir(mouse_directory):
            continue
        flir_file = os.path.join(mouse_directory, "Flir.xlsx")
        core_file = os.path.join(mouse_directory, "Core.csv")
        if not os.path.exists(flir_file) or not os.path.exists(core_file):
            continue
        core_temperatures = DisplayLoggedTemperature.average_temperatures(
            DisplayLoggedTemperature.temperatures_from_file(core_file))
        """My core temperature logger program sometimes reads one false temperature on startup. So, the next line uses
        the second read of a core temperature as the start time of the test"""
        core_temperatures = core_temperatures[1:]
        start_time = core_temperatures[0][0]
        times, core_temperatures, _ = create_relative_temperatures(core_temperatures, start_time, relative_temp=False)
        i = 0
        for i, time in enumerate(times):
            if time > .33:
                break
        core_change = core_temperatures[i] - core_temperatures[0]
        ax[0].plot(times, core_temperatures)
        dfs = pd.read_excel(flir_file, sheet_name=None)
        df = dfs[list(dfs.keys())[0]]
        times, tail_temperatures, _ = create_relative_temperatures(
            [(datetime.datetime.combine(start_time.date(), a), b) for a, b in zip(df["Time"], df["Tail Temperature"])], start_time, relative_temp=False)
        tail_AUC = 0
        for i in range(1, len(times)):
            tail_AUC += (times[i] - times[i - 1]) * (tail_temperatures[i] + tail_temperatures[i - 1]) / 2
            if times[i] > .33:
                break
        ax[1].plot(times, tail_temperatures)
        times, bat_temperatures, _ = create_relative_temperatures(
            [(datetime.datetime.combine(start_time.date(), a), b) for a, b in zip(df["Time"], df["BAT Temperature"])], start_time, relative_temp=False)
        bat_AUC = 0
        for i in range(1, len(times)):
            bat_AUC += (times[i] - times[i - 1]) * (bat_temperatures[i] + bat_temperatures[i - 1]) / 2
            if times[i] > .33:
                break
        print(mouse_number, "%.2f" % tail_AUC, "%.2f" % bat_AUC, "%.2f" % core_change)
        ax[2].plot(times, bat_temperatures, label=mouse_number)
        #ax[0, 1].legend()
    ax[0].set_ylabel("Core (°C)")
    ax[1].set_ylabel("Tail (°C)")
    ax[2].set_ylabel("BAT (°C)")
    ax[2].set_xlim([0, 0.5])
    ax[2].set_xlabel("Hours")
    plt.show()