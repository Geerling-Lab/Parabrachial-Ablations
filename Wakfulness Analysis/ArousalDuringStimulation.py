import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib.patches as patches
import datetime
import glob
import os
import argparse
import math
from Helper import simple_beeswarm

"""
This program generates averaged plots showing wakefulness before, during, and after optogenetic stimulation
It expects that the stimulations were recorded with the TTL signal in Sirenia, and then exported to an annotations file.
Annotations file should be stored in a folder together, and exported score files should also be stores together"""


class Stimulation:
    def __init__(self, start_time, end_time, frequency, mouse, annotation_file, score_file, spectra_file):
        self.start_time = start_time
        self.end_time = end_time
        if frequency < 3:
            self.frequency = 2
        elif 3 < frequency < 7:
            self.frequency = 5
        elif 7 < frequency < 15:
            self.frequency = 10
        else:
            self.frequency = 20
        self.mouse = mouse
        if self.mouse in ["JR7", "JR8", "JR9", "JR10", "JR11", "JR12", "05960", "06555"]:
            self.Chr2 = True
        else:
            self.Chr2 = False
        self.annotation_file = annotation_file
        self.score_file = score_file
        self.spectra_file = spectra_file
        self.scores = None
        self.spectra = None
        self.awake_at_stimulation_start = False
        self.get_scores()
        self.get_spectra()

    def get_scores(self):
        self.scores = np.zeros(shape=(70,))
        with open(self.score_file) as f:
            reading_header = True
            for line in f:
                if reading_header:
                    if line.split("\t")[0] == "Date":
                        reading_header = False
                    continue
                date, time, _, _, score = line.strip().split("\t")
                score = int(float(score))
                time = datetime.datetime.strptime("%s %s" % (date, time), "%m/%d/%Y %H:%M:%S")
                relative_start_time = int(round((time - self.start_time).total_seconds()))
                if relative_start_time < -35:
                    continue
                elif relative_start_time < -30:
                    self.scores[0: relative_start_time - 25] = score
                elif relative_start_time > 35:
                    self.scores[relative_start_time + 30: 70] = score
                    break
                else:
                    self.scores[relative_start_time + 30: relative_start_time + 35] = score
        if np.any(self.scores[20:31] != 2):
            print("Mouse %s was awake at stimulation (%i Hz) starting at %s" % (self.mouse, self.frequency, self.start_time))
            self.awake_at_stimulation_start = True
        if np.any(self.scores == 0):
            print("Unscored: Mouse %s Datetime %s Frequency %s" % (self.mouse, self.start_time, self.frequency))

    def get_spectra(self):
        self.spectra = np.zeros(shape=(70, 50))
        with open(self.spectra_file) as f:
            reading_header = True
            for line in f:
                if reading_header:
                    if line.split("\t")[0] == "Date":
                        reading_header = False
                    continue
                date, time, _, _, *frequencies = line.strip().split("\t")
                frequencies = np.array([float(x) for x in frequencies])
                time = datetime.datetime.strptime("%s %s" % (date, time), "%m/%d/%Y %H:%M:%S")
                relative_start_time = int(round((time - self.start_time).total_seconds()))
                if relative_start_time < -35:
                    continue
                elif relative_start_time < -30:
                    self.spectra[0: relative_start_time - 25, :] = frequencies
                elif relative_start_time > 35:
                    self.spectra[relative_start_time + 30: 70, :] = frequencies
                    break
                else:
                    self.spectra[relative_start_time + 30: relative_start_time + 35, :] = frequencies
        #self.spectra /= np.sum(self.spectra)


class CreateWakeGraph:
    def __init__(self, title):
        self.title = title
        self.Chr2_dict = {}
        self.mCh_dict = {}

    def add_stimulation(self, Chr2, frequency, mouse_number, wake):
        dict = self.Chr2_dict if Chr2 else self.mCh_dict
        if frequency not in dict:
            dict[frequency] = {}
        if mouse_number not in dict[frequency]:
            dict[frequency][mouse_number] = []
        dict[frequency][mouse_number].append(wake)

    def graph(self):
        frequencies = [2, 5, 10, 20]
        colors = ['r', 'grey']
        bar_colors = [(1, .7, .7), (.7, .7, .7)]
        fig, ax = plt.subplots()
        ax.set_ylim([0, 1])
        x_values, y_values = [], []
        for i, dict in enumerate([self.Chr2_dict, self.mCh_dict]):
            for j, frequency in enumerate(frequencies):
                for mouse_number, stimulations in dict[frequency].items():
                    if len(stimulations) >= 3:
                        dict[frequency][mouse_number] = sum(stimulations) / len(stimulations)
                        print("Mouse Number: %s\tVector: %s\tFrequency: %i Hz\tValue: %.2f\tn=%i" % (mouse_number, ["Chr2", "mCh"][i], frequency, sum(stimulations) / len(stimulations), len(stimulations)))
                    else:
                        dict[frequency][mouse_number] = np.nan
                values = list(dict[frequency].values())
                xs = simple_beeswarm(values, base_x=j+.3 * i, width=.08, nbins=2)
                ax.plot(xs, values, color=colors[i], marker="o", linestyle='none')
                ax.bar(j + .3 * i, np.nanmean(np.array(values)), color=bar_colors[i], width=.3)
        ax.set_xticks([x+.15 for x in range(len(frequencies))])
        ax.set_xticklabels(["%i Hz" % x for x in frequencies])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.suptitle(self.title)
        plt.show()


def parse_annotation_folder(annotation_folder, score_folder, spectra_folder):
    stimulations = []
    for file in glob.glob(os.path.join(annotation_folder, "*.tsv")):
        prefix = os.path.basename(file).split("Annotation")[0].strip()
        mouse_number = prefix.split(" ")[0]
        score_file = os.path.join(score_folder, prefix + " scores.tsv")
        spectra_file = os.path.join(spectra_folder, prefix + " spectra.tsv")
        if not os.path.exists(score_file):
            raise IOError("%s does not exist" % score_file)
        with open(file) as f:
            reading_header = True
            last_pulse_time = None
            pulse_train_start_time = None
            pulses_in_pulse_train = 0
            for line in f:
                if reading_header:
                    words = line.split("\t")
                    if words[0] == "0":
                        reading_header = False
                        trial_start_time = datetime.datetime.strptime(line.split("\t")[1], "%m/%d/%y %H:%M:%S.%f")
                    continue
                words = line.split("\t")
                if len(words) != 6:
                    continue
                number, pulse_time, _, _, channel, annotation = line.split("\t")
                if channel != "EEG2" and channel != "ALL":
                    continue
                pulse_time = datetime.datetime.strptime(pulse_time, "%m/%d/%y %H:%M:%S.%f")
                if not pulse_train_start_time:
                    pulse_train_start_time = pulse_time
                else:
                    pulse_train_duration = last_pulse_time - pulse_train_start_time
                    if (pulse_time - last_pulse_time).total_seconds() > 10:
                        if 25 > pulse_train_duration.total_seconds() > 1:
                            stimulations.append(Stimulation(
                                pulse_train_start_time, pulse_time,
                                pulses_in_pulse_train / pulse_train_duration.total_seconds(),
                                mouse_number, file, score_file, spectra_file))
                        pulse_train_start_time = pulse_time
                        pulses_in_pulse_train = 0
                    else:
                        pulses_in_pulse_train += 1
                last_pulse_time = pulse_time
    return stimulations


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--annotations", help="Path to annotation folder")
    ap.add_argument("-s", "--scores", help="Path to scores folder")
    ap.add_argument("-f", "--spectra", help="Path to spectra folder")
    args = vars(ap.parse_args())
    stimulations = parse_annotation_folder(args["annotations"], args["scores"], args["spectra"])

    wake_during_stimulation = CreateWakeGraph("Wake During Stimulation")
    wake_after_stimulation = CreateWakeGraph("Wake After Stimulation")
    for Chr2 in [True, False]:
        for frequency in [2, 5, 10, 20]:
            score_array = []
            spectra_array = {}
            for stimulation in stimulations:
                if stimulation.frequency == frequency and stimulation.Chr2 == Chr2 and not stimulation.awake_at_stimulation_start:
                    score_array.append(stimulation.scores)
                    wake_during_stimulation.add_stimulation(Chr2, frequency, stimulation.mouse, np.any(stimulation.scores[30:40] == 1))
                    wake_after_stimulation.add_stimulation(Chr2, frequency, stimulation.mouse, np.any(stimulation.scores[30:] == 1))
                    if stimulation.mouse not in spectra_array:
                        spectra_array[stimulation.mouse] = []
                    spectra_array[stimulation.mouse].append(stimulation.spectra)
            print("Frequency: %i Hz" % frequency)
            for mouse, spectra in spectra_array.items():
                print("Mouse: %s\t%i stimulations" % (mouse, len(spectra)))

            score_array = np.array(score_array)
            x = np.arange(-30, 40, 1)
            wake_probability = np.mean(score_array == 1, axis=0)
            nrem_probability = np.mean(score_array == 2, axis=0)
            rem_probability = np.mean(score_array == 3, axis=0)
            fig, ax = plt.subplots(nrows=3)
            ax[0].plot(x, wake_probability, "r-", label="Wake", lw=2)
            ax[0].plot(x, nrem_probability, "b-", label="NREM", lw=2)
            ax[0].plot(x, rem_probability, "g-", label="REM", lw=2)
            ax[1].set_xlabel("Seconds")
            ax[0].set_ylabel("Probability")
            rect = patches.Rectangle((0, -.05), 10, 1.1, facecolor=(.1, .83, .92, 0.2))
            ax[0].add_patch(rect)
            ax[0].set_ylim([0, 1])
            ax[0].set_xlim([-30, 40])
            ax[0].spines['right'].set_visible(False)
            ax[0].spines['top'].set_visible(False)
            ax[0].legend()

            averaged_spectra_array = {}
            for mouse, spectra_arrays in spectra_array.items():
                averaged_spectra_array[mouse] = np.mean(np.array(spectra_arrays), axis=0)
            spectra_array = np.mean(np.array(list(averaged_spectra_array.values())), axis=0)
            spectra_array = np.flipud(spectra_array.T)
            ax[1].imshow(spectra_array, extent=[-30, 40, 0, 50], aspect='auto', cmap="nipy_spectral")
            rect = patches.Rectangle((0, -1), 10, 52, facecolor=(0, 0, 0, 0), edgecolor=(1, 1, 1, 1))
            ax[1].add_patch(rect)
            ax[1].set_ylim([0, 25])

            pre_stim_spectras = []
            stim_spectras = []
            post_stim_spectras = []
            individual_mouse_kwargs = {"alpha": 0.2, "lw": 0.5, "zorder": 1}
            for mouse, spectra_array in averaged_spectra_array.items():
                pre_stim_spectra = np.mean(spectra_array[20:30, :], axis=0)
                pre_stim_spectras.append(pre_stim_spectra)
                #ax[2].plot(np.arange(pre_stim_spectra.size), pre_stim_spectra, c="violet", **individual_mouse_kwargs)
                stim_spectra = np.mean(spectra_array[30:40, :], axis=0)
                stim_spectras.append(stim_spectra)
                #ax[2].plot(np.arange(stim_spectra.size), stim_spectra, c="blue", **individual_mouse_kwargs)
                post_stim_spectra = np.mean(spectra_array[40:50, :], axis=0)
                post_stim_spectras.append(post_stim_spectra)
                #ax[2].plot(np.arange(post_stim_spectra.size), post_stim_spectra, c="navajowhite", **individual_mouse_kwargs)
            pre_stim_spectras = np.mean(np.array(pre_stim_spectras), axis=0)
            stim_spectras = np.mean(np.array(stim_spectras), axis=0)
            post_stim_spectras = np.mean(np.array(post_stim_spectras), axis=0)
            ax[2].plot(np.arange(pre_stim_spectras.size), pre_stim_spectras, c="violet", zorder=2, lw=2, label="Pre-stimulation")
            ax[2].plot(np.arange(stim_spectras.size), stim_spectras, c="blue", zorder=3, lw=2, label="Stimulation")
            ax[2].plot(np.arange(post_stim_spectras.size), post_stim_spectras, c="navajowhite", zorder=2, lw=2, label="Post-stimulation")
            ax[2].legend()
            ax[2].set_xlim([0, 25])
            ax[2].set_ylim([0, 500])
            ax[2].spines['right'].set_visible(False)
            ax[2].spines['top'].set_visible(False)

            if Chr2:
                fig.suptitle("Chr2 Frequency: %s Hz" % frequency)
            else:
                fig.suptitle("mCh Frequency: %s Hz" % frequency)
            print("Frequency: %s Hz: %s" % (frequency, list(averaged_spectra_array.keys())))
            ax[0].set_title("%i stimulations, %i mice" % (score_array.shape[0], len(averaged_spectra_array.keys())), fontsize=10)
            plt.show()

    wake_during_stimulation.graph()
    wake_after_stimulation.graph()

