import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
import os
from CompareStressTests import read_keyfile, TrialType
import glob
import datetime
import scipy.stats
import random
"""
This program attempts to correlate mouse movement with the temperature change. If the mouse is moving, it's likely
to have a lot of transitions between valid and invalid temperature readings. We divide the record into epochs,
and create a scatter plot of number of transitions versus temperature change."""

if __name__ == "__main__":

    font_path = r"C:\Users\Filla\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    root_folder = r"R:\Fillan\Parabrachial Ablations\Temperature"
    key_file = os.path.join(root_folder, "Keyfile.txt")
    mouse_numbers = read_keyfile(key_file)
    epoch_length = datetime.timedelta(minutes=10)

    control_mice = ["4803", "5089", "5090", "5093", "4978", "5175", "5176"]
    """
    Plot all cold mice as blue
    """
    control_epochs = []
    ablation_epochs = []
    for core_file in glob.glob(os.path.join(root_folder, "*", "Cold *", "Core *.csv")):
        mouse_number = os.path.basename(core_file).split(" ")[1][:-4]
        print(mouse_number)
        start_time = mouse_numbers[mouse_number].start[TrialType.Cold].time()
        epoch_end = None
        transitions = 0
        start_temperature = 0
        last_temperature = 0

        for line in open(core_file, 'r'):
            dt, temperature = line.strip().split(",")
            try:
                dt = datetime.datetime.strptime(dt, "%Y-%m-%d %H:%M:%S.%f")
            except ValueError:
                continue
            temperature = float(temperature)
            if dt.time() < start_time:
                continue
            if epoch_end is None:
                epoch_end = dt + epoch_length
                transitions = 0
                start_temperature = temperature
                continue
            if dt > epoch_end:
                if temperature < 1e-6:
                    continue
                if mouse_number in control_mice:
                    control_epochs.append((transitions + random.random(), temperature))
                else:
                    ablation_epochs.append((transitions + random.random(), temperature))
                epoch_end = epoch_end + epoch_length
                transitions = 0
                start_temperature = temperature
                continue
            if (last_temperature < 1e-6) != (temperature < 1e-6):
                transitions += 1
            if start_temperature < 1e-6:
                start_temperature = temperature
            last_temperature = temperature
    control_epochs = np.array(control_epochs)
    ablation_epochs = np.array(ablation_epochs)
    plt.scatter(control_epochs[:, 0], control_epochs[:, 1], facecolor='none', edgecolors='#547df0', s=4, alpha=0.3)
    plt.scatter(ablation_epochs[:, 0], ablation_epochs[:, 1], facecolor='#1546cf', edgecolors='#1546cf', s=4, alpha=0.8)

    "Plot all control mice as green"
    control_epochs = []
    ablation_epochs = []
    for core_file in glob.glob(os.path.join(root_folder, "*", "Control *", "Core *.csv")):
        mouse_number = os.path.basename(core_file).split(" ")[1][:-4]
        print(mouse_number)
        start_time = mouse_numbers[mouse_number].start[TrialType.Cold].time()
        epoch_end = None
        transitions = 0
        start_temperature = 0
        last_temperature = 0

        for line in open(core_file, 'r'):
            dt, temperature = line.strip().split(",")
            try:
                dt = datetime.datetime.strptime(dt, "%Y-%m-%d %H:%M:%S.%f")
            except ValueError:
                continue
            temperature = float(temperature)
            if dt.time() < start_time:
                continue
            if epoch_end is None:
                epoch_end = dt + epoch_length
                transitions = 0
                start_temperature = temperature
                continue
            if dt > epoch_end:
                if temperature < 1e-6:
                    continue
                if mouse_number in control_mice:
                    control_epochs.append((transitions + random.random(), temperature))
                else:
                    ablation_epochs.append((transitions + random.random(), temperature))
                epoch_end = epoch_end + epoch_length
                transitions = 0
                start_temperature = temperature
                continue
            if (last_temperature < 1e-6) != (temperature < 1e-6):
                transitions += 1
            if start_temperature < 1e-6:
                start_temperature = temperature
            last_temperature = temperature
    control_epochs = np.array(control_epochs)
    ablation_epochs = np.array(ablation_epochs)
    plt.scatter(control_epochs[:, 0], control_epochs[:, 1], facecolor='none', edgecolors='#41f076', s=4, alpha=0.3)
    plt.scatter(ablation_epochs[:, 0], ablation_epochs[:, 1], facecolor='#04d142', edgecolors='#04d142', s=4, alpha=0.8)

    control_epochs = []
    ablation_epochs = []
    for core_file in glob.glob(os.path.join(root_folder, "*", "Hot *", "Core *.csv")):
        mouse_number = os.path.basename(core_file).split(" ")[1][:-4]
        print(mouse_number)
        start_time = mouse_numbers[mouse_number].start[TrialType.Cold].time()
        epoch_end = None
        transitions = 0
        start_temperature = 0
        last_temperature = 0

        for line in open(core_file, 'r'):
            dt, temperature = line.strip().split(",")
            try:
                dt = datetime.datetime.strptime(dt, "%Y-%m-%d %H:%M:%S.%f")
            except ValueError:
                continue
            temperature = float(temperature)
            if dt.time() < start_time:
                continue
            if epoch_end is None:
                epoch_end = dt + epoch_length
                transitions = 0
                start_temperature = temperature
                continue
            if dt > epoch_end:
                if temperature < 1e-6:
                    continue
                if mouse_number in control_mice:
                    control_epochs.append((transitions + random.random(), temperature))
                else:
                    ablation_epochs.append((transitions + random.random(), temperature))
                epoch_end = epoch_end + epoch_length
                transitions = 0
                start_temperature = temperature
                continue
            if (last_temperature < 1e-6) != (temperature < 1e-6):
                transitions += 1
            if start_temperature < 1e-6:
                start_temperature = temperature
            last_temperature = temperature
    control_epochs = np.array(control_epochs)
    ablation_epochs = np.array(ablation_epochs)
    plt.scatter(control_epochs[:, 0], control_epochs[:, 1], facecolor='none', edgecolors='#eb6a3f', s=4, alpha=0.3)
    plt.scatter(ablation_epochs[:, 0], ablation_epochs[:, 1], facecolor='#de3700', edgecolors='#de3700', s=4, alpha=0.8)

    plt.xlabel("Number of transitions")
    plt.ylabel("Core Temperature (°C)")
    plt.savefig("Movement vs Temperature Correlation.eps", format="eps", bbox_inches="tight")
    plt.show()
