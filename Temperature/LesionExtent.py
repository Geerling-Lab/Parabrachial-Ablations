import matplotlib.pyplot as plt
import os
import glob
import numpy as np
import scipy.stats
import CompareStressTests
import DisplayLoggedTemperature
import datetime
from CompareStressTests import read_keyfile, TrialType
from matplotlib import font_manager
import pandas as pd

class Casp3s:
    def __init__(self):
        self.stereologies = []
        self.cores = []

    def graph(self, ax):
        self.stereologies = np.array(self.stereologies)
        self.cores = np.array(self.cores)
        _, p = scipy.stats.pearsonr(self.stereologies, self.cores)
        ax.plot(np.unique(self.stereologies), np.poly1d(np.polyfit(self.stereologies, self.cores, 1))(np.unique(self.stereologies)), label="p=%.3f" % p, c="k", zorder=0)
        ax.legend()

if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 2

    root_folder = r"R:\Fillan\Parabrachial Ablations\Temperature"
    key_file = os.path.join(root_folder, "Keyfile.txt")
    mouse_data = read_keyfile(key_file)

    temperatures = {}
    f = open(r"R:\Fillan\Parabrachial Ablations\Temperature\Temperature at Times.csv", 'w')

    lmx1b_df = pd.read_excel(r"R:\Fillan\Parabrachial Ablations\Experiment 3 Lmx1b\Lmx1b mice.xlsx", sheet_name="Overview")
    foxp2_df = pd.read_excel(r"R:\Fillan\Parabrachial Ablations\Experiment 3 Foxp2\FoxP2 mice.xlsx", sheet_name="Overview")
    vglut2_df = pd.read_excel(r"R:\Fillan\Parabrachial Ablations\Experiment 3 Vglut2 without EEG\Vglut2 without EEG mice.xlsx", sheet_name="Overview")
    vglut2_with_eeg_df = pd.read_excel(r"R:\Fillan\Parabrachial Ablations\Experiment 3 Vglut2 with EEG\Vglut2 with EEG mice.xlsx", sheet_name="Overview")

    lmx1b_casp3s = Casp3s()
    foxp2_caps3s = Casp3s()
    vglut2_cap3s = Casp3s()
    fig, axs = plt.subplots(nrows=4, sharey=True)
    for core_file in glob.glob(os.path.join(root_folder, "*", "Cold *", "Core *.csv")):
        folder = os.path.dirname(core_file)
        mouse_number = os.path.basename(core_file).split(" ")[1][:-4]
        if "repeat" in core_file.lower():
            continue
        mouse_datum = mouse_data[mouse_number]

        if mouse_datum.genotype == "FoxP2" or mouse_datum.genotype == "Vglut2" or mouse_datum.genotype == "Lmx1b":
            temps = DisplayLoggedTemperature.temperatures_from_file(core_file)
            core_temperatures = DisplayLoggedTemperature.average_temperatures(temps)
            date = core_temperatures[0][0].date()
            start_time = datetime.datetime.combine(date, mouse_datum.start[CompareStressTests.TrialType.Cold].time())
            start_temperature = None
            delta_temperature = None
            for time, temperature in core_temperatures:
                if time > start_time and start_temperature is None:
                    start_temperature = temperature
                if start_temperature is not None:
                    delta_temperature = temperature - start_temperature
                if time > start_time + datetime.timedelta(hours=4):
                    break
            PBrel = None
            Stereology = None
            facecolor = None
            edgecolor = None
            marker = None
            if mouse_datum.genotype == "Vglut2":
                row = vglut2_df.loc[vglut2_df['Mouse'] == int(mouse_number)]
                if row.size == 0:
                    row = vglut2_with_eeg_df.loc[vglut2_with_eeg_df['Mouse'] == int(mouse_number)]
                Stereology = row["Stereology"].values[0]
                vector = row["Vector"].values[0]
                lPBrel = row["L PBreL"].values[0]
                rPBrel = row["R PBreL"].values[0]
                edgecolor = 'red'
                marker = "o"
                if vector == "Cas3":
                    facecolor = edgecolor
                    if not np.isnan(Stereology):
                        vglut2_cap3s.stereologies.append(Stereology)
                        vglut2_cap3s.cores.append(delta_temperature)
                else:
                    facecolor = 'w'
                axs[0].scatter([Stereology], [delta_temperature], facecolor=facecolor, edgecolor=edgecolor,
                               marker=marker, s=50)
            elif mouse_datum.genotype == "FoxP2":
                Stereology = foxp2_df.loc[foxp2_df['Mouse'] == int(mouse_number)]["Remaining Neurons"].values[0]
                lPBrel = foxp2_df.loc[foxp2_df['Mouse'] == int(mouse_number)]["L PBreL"].values[0]
                rPBrel = foxp2_df.loc[foxp2_df['Mouse'] == int(mouse_number)]["R PBreL"].values[0]
                vector = foxp2_df.loc[foxp2_df['Mouse'] == int(mouse_number)]["Vector"].values[0]
                edgecolor = 'orangered'
                if vector == "Cas3":
                    facecolor = edgecolor
                    if not np.isnan(Stereology):
                        foxp2_caps3s.stereologies.append(Stereology)
                        foxp2_caps3s.cores.append(delta_temperature)
                else:
                    facecolor = 'w'
                axs[1].scatter([Stereology], [delta_temperature], facecolor=facecolor, edgecolor=edgecolor,
                               marker=marker, s=50)
            else:
                Stereology = lmx1b_df.loc[lmx1b_df['Mouse'] == int(mouse_number)]["YOLO"].values[0]
                lPBrel = -1
                rPBrel = -1
                vector = lmx1b_df.loc[lmx1b_df['Mouse'] == int(mouse_number)]["Vector"].values[0]
                edgecolor = 'violet'
                if vector == "Cas3":
                    facecolor = edgecolor
                    if not np.isnan(Stereology):
                        lmx1b_casp3s.stereologies.append(Stereology)
                        lmx1b_casp3s.cores.append(delta_temperature)
                else:
                    facecolor = 'w'
                axs[2].scatter([Stereology], [delta_temperature], facecolor=facecolor, edgecolor=edgecolor,
                               marker=marker, s=50)
            if lPBrel >= 0:
                axs[3].scatter([lPBrel + rPBrel], [delta_temperature], facecolor=facecolor, edgecolor=edgecolor, marker=marker, s=50)

            print("%s:%.2f:%.2f" % (mouse_number, lPBrel+rPBrel, delta_temperature))
    vglut2_cap3s.graph(axs[0])
    foxp2_caps3s.graph(axs[1])
    lmx1b_casp3s.graph(axs[2])
    axs[0].set_xlabel("Stereological Count (Vglut2)")
    axs[0].set_ylabel("ΔTemperature over 4 hours (C)")
    axs[0].set_xlim([0, axs[0].get_xlim()[1]])
    axs[1].set_xlabel("YOLO Count (Foxp2)")
    axs[1].set_xlim([0, axs[1].get_xlim()[1]])
    axs[1].set_ylabel("ΔTemperature over 4 hours (C)")
    axs[2].set_xlabel("YOLO Count (Lmx1b)")
    axs[2].set_xlim([0, axs[2].get_xlim()[1]])
    axs[2].set_ylabel("ΔTemperature over 4 hours (C)")
    axs[3].set_xlabel("PBrel")
    axs[3].set_ylabel("ΔTemperature over 4 hours (C)")
    plt.show()
