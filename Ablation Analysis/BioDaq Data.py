"""
This program is designed to take in a series of "Experiment Tab" files from BioDaq
As well as a key file that describes when mice were part of the experiment, and when certain tests were being
performed
"""
import datetime
import matplotlib.pyplot as plt
import os
import numpy as np
import functools
import glob

import pandas
import pandas as pd
from matplotlib import font_manager, colors
from WobblePlot import graph_wobble_plot
import scipy.stats


class MouseColor:
    def __init__(self):
        self.colors = {}

    def get_color(self, mouse):
        if mouse in self.colors:
            return self.colors[mouse]
        else:
            start = None
            if "Vglut2+Cas3 (Baseline)" in mouse:
                start = ac.colors["Vglut2+Cas3 (Baseline)"]
            elif "Vglut2+Cas3" in mouse:
                start = ac.colors["Vglut2+Cas3"]
            elif "Vglut2+mCh (Baseline)" in mouse:
                start = ac.colors["Vglut2+mCh (Baseline)"]
            elif "Vglut2+mCh" in mouse:
                start = ac.colors["Vglut2+mCh"]
            elif "Lmx1b+Cas3 (Baseline)" in mouse:
                start = ac.colors["Lmx1b+Cas3 (Baseline)"]
            elif "Lmx1b+Cas3" in mouse:
                start = ac.colors["Lmx1b+Cas3"]
            elif "TH+Cas3" in mouse:
                start = ac.colors["TH+Cas3"]
            if start is None:
                c = list(np.random.rand(3,))
            else:
                c = []
                for i in range(3):
                    c.append(sorted([0, np.random.normal(start[i], 0, (1,))[0], 1])[1])
            self.colors[mouse] = c
            return c


class AverageColor:
    def __init__(self) -> None:
        self.colors = {"Vglut2+Cas3": 'r', "Vglut2+Cas3 (Baseline)": 'r', "Vglut2+Cas3 (Saline)": [1, 0.8, 0.2],
                       "Vglut2+mCh": 'grey', "Vglut2+mCh (Baseline)": 'grey', "Vglut2+mCh (Saline)": [0.2, 0.8, 1],
                       "Lmx1b+Cas3": 'violet', "Lmx1b+Cas3 (Baseline)": 'violet', "Lmx1b+Cas3 (Saline)": [1, 0.8, 0.2],
                       "Lmx1b+mCh": 'grey', "Lmx1b+mCh (Baseline)": 'grey', "Lmx1b+mCh (Saline)": [0.2, 0.8, 1],
                       "Foxp2+Cas3": 'orangered', "Foxp2+Cas3 (Baseline)": 'orangered', "Foxp2+Cas3 (Saline)": [1, 0.8, 0.2],
                       "Foxp2+mCh": 'grey', "Foxp2+mCh (Baseline)": 'grey', "Foxp2+mCh (Saline)": [0.2, 0.8, 1],
                       "TH+Cas3": [.5, .5, 0], "TH+Cas3 (Baseline)": [.8, .8, .5], "TH+Cas3 (Saline)": [.8, .5, .2]}

    def get_color(self, genotype):
        if genotype in self.colors:
            return list(colors.to_rgb(self.colors[genotype]))
        else:
            raise IOError("%s is not a valid genotype" % genotype)


class Trial:
    def __init__(self, start_time, time_delta, hopper, mouse="", genotype="", addendum=""):
        self.start_time = start_time
        self.time_delta = time_delta
        self.end_time = start_time + time_delta
        self.hopper = hopper
        self.values = [(0, 0)]
        self.total = 0
        self.mouse = mouse
        self.genotype = genotype + addendum
        self.pp_genotype = genotype.replace("+", os.linesep) + os.linesep + addendum.strip()
        self.label = "%s (%s)" % (mouse, genotype)

    def add_line(self, date_time, change, hopper, bout_length):
        if hopper == self.hopper and self.start_time < date_time < self.end_time:
            bout_start_from_trial_start = (date_time - self.start_time).total_seconds() / 60
            bout_end_from_trial_start = (date_time + bout_length - self.start_time).total_seconds() / 60
            self.values.append((bout_start_from_trial_start, self.total))
            self.total += change
            self.values.append((bout_end_from_trial_start, self.total))

    def add_last_time_point(self):
        self.values.append((self.time_delta.total_seconds() / 60, self.total))
        self.values.sort(key=lambda x: x[0])

    def graph(self, ax):
        x = [i[0] for i in self.values]
        y = [i[1] for i in self.values]
        l, = ax.plot(x, y, "-", c=ac.get_color(self.genotype))
        l.set_label(self.label)

    def value_at_time(self, time_delta):
        if time_delta >= self.values[-1][0]:
            return self.values[-1][1]
        i = 0
        while self.values[i][0] < time_delta:
            i += 1
        slope = (self.values[i][1] - self.values[i - 1][1]) / (self.values[i][0] - self.values[i - 1][0])
        value = self.values[i][1] + (time_delta - self.values[i][0]) * slope
        return value


def error_bar(ax, x, y, yerr, color, label=None):
    ax.plot(x, y, color=color, label=label)
    lower = [i - j for i, j in zip(y, yerr)]
    upper = [i + j for i, j in zip(y, yerr)]
    ax.fill_between(x, lower, upper, color=color, alpha=0.2, edgecolor='w')


class Trial_Type:
    def __init__(self, trial_name, trial_time, pp_dict=None, max=None, meal_max=None):
        self.trial_name = trial_name
        self.trial_time = trial_time
        self.trials = []
        self.max = max
        self.meal_max = meal_max
        if self.trial_time.total_seconds() <= 7200:
            self.increment_for_average = 10
        elif self.trial_time.total_seconds() < 86400:
            self.increment_for_average = 60  # In minutes
        else:
            self.increment_for_average = 1440
        default_pp_dict = {"Consumption": "Fluid", "Experiment": "Experiment", "Baseline": "Baseline", "Title": None}
        self.pp_dict = {}
        for key, value in default_pp_dict.items():
            if key in pp_dict:
                self.pp_dict[key] = pp_dict[key]
            else:
                self.pp_dict[key] = value

    def graph_individual_time_series(self):
        fig, ax = plt.subplots()
        for trial in self.trials:
            if "baseline" in trial.genotype.lower() or "saline" in trial.genotype.lower():
                continue
            trial.add_last_time_point()
            trial.graph(ax)
        ax.set_xlim([0, self.trial_time.total_seconds() / 60])
        plt.title(self.trial_name)
        plt.ylabel("Consumption (g)")
        plt.xlabel("Time (min)")
        plt.tight_layout()
        plt.show()

    def graph_average_time_series(self, writer=None):
        genotypes = {}
        for i in range(self.increment_for_average, int(self.trial_time.total_seconds() / 60) + 1, self.increment_for_average):
            for trial in self.trials:
                if trial.genotype not in genotypes:
                    genotypes[trial.genotype] = {0: [0]}
                if i not in genotypes[trial.genotype]:
                    genotypes[trial.genotype][i] = []
                genotypes[trial.genotype][i].append(trial.value_at_time(i))

        df = pd.DataFrame()
        df['Minutes'] = range(0, int(self.trial_time.total_seconds() / 60) + 1, self.increment_for_average)

        fig, axs = plt.subplots(nrows=2, sharey=True)
        for genotype, dict in genotypes.items():
            n = 0
            xs = []
            ys = []
            yerr = []
            for time, values in dict.items():
                n = max(n, len(values))
                xs.append(time)
                ys.append(np.array(values).mean())
                yerr.append(np.array(values).std())
            if "Sucrose" in self.trial_name or "Quinine" in self.trial_name:
                n = n / 2
            error_bar(axs[int("Baseline" not in genotype)], xs, ys, yerr, label="%s (n=%i)" % (genotype.split(" (")[0], n), color=ac.get_color(genotype))
            df[genotype + " Mean (n=%i)" % n] = list(map(lambda s: "%.2f" % s, ys))
            df[genotype + " Stdev"] = list(map(lambda s: "%.2f" % s, yerr))
        if writer:
            df.to_excel(writer, sheet_name="%s time series" % self.trial_name)
        plt.suptitle(self.trial_name)
        for ax in axs:
            trial_time_minutes = self.trial_time.total_seconds() / 60
            ax.set_xlim(0, trial_time_minutes)
            if trial_time_minutes > 500:
                ax.set_xticks(np.arange(0, trial_time_minutes + 1, 120))
            elif trial_time_minutes > 200:
                ax.set_xticks(np.arange(0, trial_time_minutes + 1, 30))
            ax.set_ylim(0, ax.get_ylim()[1])
        axs[0].set_ylabel("%s (g)" % self.pp_dict["Consumption"])
        axs[1].set_ylabel("%s (g)" % self.pp_dict["Consumption"])
        ax_0 = axs[0].twinx()
        ax_0.set_ylabel(self.pp_dict["Baseline"])
        ax_0.set_yticks([])
        ax_1 = axs[1].twinx()
        ax_1.set_yticks([])
        ax_1.set_ylabel(self.pp_dict["Experiment"])
        axs[0].legend(prop={'size': 8})
        axs[1].legend(prop={'size': 8})
        axs[0].set_xticklabels([])
        axs[1].set_xlabel("Time (min)")
        plt.tight_layout()
        plt.subplots_adjust(hspace=.1, wspace=0)
        if self.max:
            axs[0].set_ylim([0, self.max])
        plt.show()

    def graph_dot_plot(self):
        genotypes = {}
        for trial in self.trials:
            if trial.genotype not in genotypes:
                genotypes[trial.genotype] = {}
            if trial.mouse not in genotypes[trial.genotype]:
                genotypes[trial.genotype][trial.mouse] = []
            genotypes[trial.genotype][trial.mouse].append(max(0, trial.value_at_time(self.trial_time.total_seconds() / 60)))
        labels = []
        last_genotype = None
        fig, ax = plt.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 6]})
        fig.set_size_inches(4, 6)
        genotype_list = sorted(list(genotypes), reverse=True)
        num_statistically_significant = 0
        x_values = sorted(list(range(5, int(15*len(genotype_list)/2), 15)) + list(range(15, int(15*len(genotype_list)/2)+1, 15))) # x values for where ticks go. (5, 15, 20, 30, 35, 45 ...)
        for i, genotype in enumerate(genotype_list):
            print(self.trial_name, genotype)
            if last_genotype:
                if genotype not in last_genotype:
                    last_genotype = None
            pretty_genotype = genotype.replace("+", "\n").replace(" (", "\n").replace(")", "")
            if "Baseline" in pretty_genotype:
                pretty_genotype = pretty_genotype.replace("Baseline", self.pp_dict["Baseline"])
            else:
                pretty_genotype += "\n" + self.pp_dict["Experiment"]
            if "baseline" in genotype.lower():
                kwargs = {'facecolors': 'w'}
            else:
                kwargs = {}
            labels.append(pretty_genotype)
            total, count = 0, 0
            for mouse, values in genotypes[genotype].items():
                value = sum(values) / len(values)
                print(mouse, "%.2f" % value)
                if last_genotype is not None:
                    if mouse in genotypes[last_genotype]:
                        last_values = genotypes[last_genotype][mouse]
                        last_value = sum(last_values) / len(last_values)
                        ax[1].plot([x_values[i-1], x_values[i]], [last_value, value], color=ac.get_color(genotype), zorder=2)
                        ax[1].scatter(x=x_values[i], y=value, color=ac.get_color(genotype), s=50, zorder=3, **kwargs)
                        total += value
                        count += 1
                else:
                    ax[1].scatter(x=x_values[i], y=value, color=ac.get_color(genotype), s=50, zorder=3, **kwargs)
                    total += value
                    count += 1
            ax[1].hlines(total / count, x_values[i] - 2, x_values[i] + 2, zorder=1)
            last_genotype = genotype
            for j in range(i+1, len(genotype_list)):
                if ("Baseline" in genotype) == ("Baseline" in genotype_list[j]):
                    first_arr = [sum(values) / len(values) for values in genotypes[genotype].values()]
                    second_arr = [sum(values) / len(values) for values in list(genotypes[genotype_list[j]].values())]
                    _, p = scipy.stats.ttest_ind(first_arr, second_arr, equal_var=False)
                    if p < 1:
                        ax[0].plot([x_values[i], x_values[i]], [num_statistically_significant, num_statistically_significant+.5], "k-")
                        ax[0].plot([x_values[j], x_values[j]], [num_statistically_significant, num_statistically_significant+.5], "k-")
                        ax[0].plot([x_values[i], x_values[j]], [num_statistically_significant+.5, num_statistically_significant+.5], "k-")
                        ax[0].text(x=(x_values[i] + x_values[j]) / 2, y=num_statistically_significant+.5, ha="center", va="bottom", s="p=%.4f" % p, fontsize=12)
                        num_statistically_significant += 2
            ax[0].set_ylim([0, num_statistically_significant+1])

        ax[1].set_xticks(x_values)
        ax[1].set_xticklabels(labels)
        ax[0].set_title(self.trial_name)
        ymin, ymax = ax[1].get_ylim()
        if self.max:
            ymax = self.max
        ax[1].set_ylim([max(0, ymin), ymax])
        trial_time = self.trial_time.total_seconds() // 3600
        string = "%s consumed (g) over %i hours" % (self.pp_dict["Consumption"], trial_time)
        ax[1].set_ylabel(string)
        ax[1].spines['top'].set_visible(False)
        ax[0].set_xlim(ax[1].get_xlim())
        ax[0].spines['bottom'].set_visible(False)
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.show()

    def graph_bout_analysis(self):
        genotypes = {}
        lines = []
        fig, ax = plt.subplots()
        for trial in self.trials:
            if 'saline' in trial.genotype:
                continue
            if trial.genotype not in genotypes:
                genotypes[trial.genotype] = (ac.get_color(trial.genotype), [])
            for before_bout, after_bout in zip(trial.values[1::2], trial.values[2::2]):
                genotypes[trial.genotype][1].append(max(after_bout[1] - before_bout[1], 0))
        graph_wobble_plot(ax, genotypes, scale=0.1)
        ax.set_title("Trial")
        ax.set_xlabel("Bout Size (g)")
        ax.set_ylabel("Bout Size (g)")
        ax.set_xlabel("Genotype")
        ax.set_title(self.trial_name)
        plt.tight_layout()
        plt.show()

    def graph_mouse_bout_analysis(self, writer):
        genotypes = {}
        for trial in self.trials:
            if trial.genotype not in genotypes:
                genotypes[trial.genotype] = {}
            if trial.mouse not in genotypes[trial.genotype]:
                genotypes[trial.genotype][trial.mouse] = []
            for before_bout, after_bout in zip(trial.values[1::2], trial.values[2::2]):
                bout_size = round(after_bout[1] - before_bout[1], 2)
                genotypes[trial.genotype][trial.mouse].append(bout_size)
        values_for_dataframe = {}
        fig, ax = plt.subplots(nrows=2, ncols=len(genotypes) // 2, sharey='all')
        for subplot_index, genotype in enumerate(genotypes):
            values_for_dataframe[genotype] = []
            mouse_values = genotypes[genotype]
            ticks = []
            labels = []
            for row_index, mouse in enumerate(mouse_values):
                values = mouse_values[mouse]
                values_for_dataframe["%s %s" % (mouse, genotype)] = list(map(lambda s: "%.2f" % s, values))
                ax.ravel()[subplot_index].scatter(np.random.rand(len(values)) * 10 + row_index * 20 + 15, values, color=ac.get_color(genotype))
                ticks.append(row_index * 20 + 20)
                labels.append(mouse)
            ax.ravel()[subplot_index].set_xticks(ticks)
            ax.ravel()[subplot_index].set_xticklabels(labels, rotation=45)
            ax.ravel()[subplot_index].set_title(genotype)
        if writer:
            df = pd.DataFrame(dict([(k, pd.Series(v, dtype="str")) for k, v in values_for_dataframe.items()]))
            df.to_excel(writer, sheet_name="%s Bouts" % self.trial_name)
        fig.suptitle(self.trial_name)
        plt.tight_layout()
        plt.show()

    def graph_meal_analysis(self, writer):
        def custom_sort(item1, item2):
            if "Baseline" in item1 and "Baseline" not in item2:
                return 1
            if "Baseline" in item2 and "Baseline" not in item1:
                return -1
            if "mCh" in item1 and "mCh" not in item2:
                return 1
            if "mCh" in item2 and "mCh" not in item1:
                return -1
            return 0

        def fillan_bee_swarm(y, x_mean, x_var, x_marker_size=1.0, y_marker_size=None):
            def create_active_xs(active_markers):
                x = []
                if len(active_markers) % 2 == 0:
                    current_x = -x_marker_size / 2
                    for _ in active_markers:
                        x.append(current_x)
                        current_x = -current_x
                        if current_x < 0:
                            current_x -= x_marker_size
                else:
                    current_x = 0
                    for _ in active_markers:
                        x.append(current_x)
                        current_x = -current_x
                        if current_x <= 0:
                            current_x -= x_marker_size
                return x

            x = []
            if len(y) == 0:
                return []
            sorting_order = [b[0] for b in sorted(enumerate(y), key=lambda i: i[1])]
            y = sorted(y)
            if y_marker_size is None:
                range = y[-1] - y[0]
                y_marker_size = 4 * range / len(y)
            active_markers = [y[0]]
            for s in y[1:]:
                if s - active_markers[0] > y_marker_size or len(active_markers) * x_marker_size > x_var * 2:
                    x.extend(create_active_xs(active_markers))
                    active_markers = []
                active_markers.append(s)
            x.extend(create_active_xs(active_markers))

            unsorted_xs = [0] * len(y)
            for s, o in zip(sorting_order, x):
                unsorted_xs[s] = o + x_mean
            return unsorted_xs

        inter_meal_interval = 5  # In minutes. From Campos,...,Palmiter 2016
        minimum_meal_size = 0.06  # In grams
        fig, ax = plt.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 6]})
        fig.set_size_inches(6, 8)
        genotypes = {}
        for trial in self.trials:
            if trial.genotype not in genotypes:
                genotypes[trial.genotype] = []
            meal_begin_index = 0
            meal_end_index = 0
            while meal_end_index < len(trial.values) - 1:
                while True:
                    meal_end_index += 2
                    if meal_end_index == len(trial.values) - 1:
                        meal_size = trial.values[meal_end_index][1] - trial.values[meal_begin_index][1]
                        if meal_size > minimum_meal_size:
                            genotypes[trial.genotype].append(meal_size)
                        break
                    if trial.values[meal_end_index - 1][0] - trial.values[meal_end_index - 2][0] > inter_meal_interval:
                        meal_size = trial.values[meal_end_index - 2][1] - trial.values[meal_begin_index][1]
                        if meal_size > minimum_meal_size:
                            genotypes[trial.genotype].append(meal_size)
                        break
                meal_begin_index = meal_end_index - 1
        ticks = []
        labels = []
        genotype_list = sorted(list(genotypes), reverse=True, key=functools.cmp_to_key(custom_sort))
        num_statistically_significant = 0
        df = pd.DataFrame()
        for row_index, genotype in enumerate(genotype_list):
            values = genotypes[genotype]
            plt.scatter(fillan_bee_swarm(values, x_mean=row_index*10, x_var=4, x_marker_size=.5, y_marker_size=.05), values,
                                                  color=ac.get_color(genotype), s=30)
            ticks.append(row_index * 10)
            pretty_genotype = genotype.replace("+", "\n").replace(" (", "\n").replace(")", "")
            if "Baseline" in pretty_genotype:
                pretty_genotype = pretty_genotype.replace("Baseline", self.pp_dict["Baseline"])
            else:
                pretty_genotype += "\n" + self.pp_dict["Experiment"]
            labels.append(pretty_genotype)
            new_df = pandas.DataFrame({genotype: values})
            df = pandas.concat([df, new_df], axis=1)
            print("%s: %.3f ± %.3f (n=%i)" % (genotype, np.mean(values), np.std(values), len(values)))
            if len(values) > 0:
                plt.hlines(y=sum(values) / len(values), xmin=row_index * 10 - 2.5, xmax=row_index * 10 + 2.5, colors='k', zorder=2)
            for j in range(row_index + 1, len(genotype_list)):
                if ("Baseline" in genotype) == ("Baseline" in genotype_list[j]):
                    first_arr = genotypes[genotype]
                    second_arr = genotypes[genotype_list[j]]
                    _, p = scipy.stats.ttest_ind(first_arr, second_arr)
                    if p < 1:
                        ax[0].plot([row_index * 10, row_index * 10], [num_statistically_significant, num_statistically_significant + .5],
                                   "k-")
                        ax[0].plot([j * 10, j * 10], [num_statistically_significant, num_statistically_significant + .5],
                                   "k-")
                        ax[0].plot([row_index * 10, j * 10],
                                   [num_statistically_significant + .5, num_statistically_significant + .5], "k-")
                        ax[0].text(x=(row_index + j) * 5, y=num_statistically_significant + .5, ha="center", va="bottom",
                                   s="p=%.4f" % p, fontsize=12)
                        num_statistically_significant += 2
        if writer:
            df.to_excel(writer, sheet_name="%s meal sizes" % self.trial_name)
        ax[0].set_ylim([0, num_statistically_significant + 1])
        ax[1].set_xticks(ticks)
        ax[1].set_xticklabels(labels)
        ax[0].set_title(self.trial_name)
        ax[1].set_ylabel("%s Meal Size (g)" % self.pp_dict["Consumption"])
        ax[1].spines['top'].set_visible(False)
        ax[0].set_xlim(ax[1].get_xlim())
        ax[0].spines['bottom'].set_visible(False)
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        if self.meal_max:
            ax[1].set_ylim([0, self.meal_max])
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.show()

    def correlate_with_stereology(self, stereology_file_path):
        genotypes = {}
        df = pd.read_excel(stereology_file_path, sheet_name="Overview")
        for trial in self.trials:
            if trial.genotype not in genotypes:
                genotypes[trial.genotype] = {}
            if trial.mouse not in genotypes[trial.genotype]:
                genotypes[trial.genotype][trial.mouse] = []
            genotypes[trial.genotype][trial.mouse].append(max(0, trial.value_at_time(self.trial_time.total_seconds() / 60)))
        fig, axs = plt.subplots(ncols=2, sharey=True)
        genotype_list = sorted(list(genotypes), reverse=True)
        stereological_counts, consumptions, colors = [[], []], [[], []], [[], []]
        for genotype in genotypes:
            if "Baseline" in genotype:
                i = 0
            else:
                i = 1
            for mouse, values in genotypes[genotype].items():
                value = sum(values) / len(values)
                try:
                    stereological_count = df.loc[df['Mouse'] == int(mouse)]["Stereology"].iloc[0]
                    stereological_counts[i].append(stereological_count)
                    consumptions[i].append(value)
                    colors[i].append(ac.get_color(genotype))
                except IndexError as ie:
                    print("Mouse %s not found" % mouse)
        axs[0].scatter(stereological_counts[0], consumptions[0], c=colors[0])
        axs[1].scatter(stereological_counts[1], consumptions[1], c=colors[1])
        axs[0].set_xlabel("Stereological Count")
        axs[1].set_xlabel("Stereological Count")
        axs[0].set_ylabel("%s consumed (g) over %i hours" % (self.pp_dict["Consumption"], self.trial_time.total_seconds() // 3600))
        axs[0].set_title(self.pp_dict["Baseline"])
        axs[1].set_title(self.pp_dict["Experiment"])
        sc = np.array(stereological_counts[0])
        c = np.array(consumptions[0])
        nans = np.logical_or(np.isnan(sc), np.isnan(c))
        if np.sum(nans) > 0:
            _, _, _, p_0, _ = scipy.stats.linregress(sc[~nans], c[~nans])
            axs[0].annotate("p=%.4f" % p_0, xy=(.99, .99), xycoords='axes fraction', size=14, ha='right', va='top')
        sc = np.array(stereological_counts[1])
        c = np.array(consumptions[1])
        nans = np.logical_or(np.isnan(sc), np.isnan(c))
        if np.sum(nans) > 0:
            _, _, _, p_1, _ = scipy.stats.linregress(sc[~nans], c[~nans])
            axs[1].annotate("p=%.4f" % p_1, xy=(.99, .99), xycoords='axes fraction', size=14, ha='right', va='top')
        plt.tight_layout()
        plt.show()


def generate_trials_from_keyfile(key_file):
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    fast_refeed = Trial_Type(trial_name="Fast Refeed", trial_time=datetime.timedelta(hours=4), pp_dict={"Consumption": "Food", "Experiment": "Fasting"}, max=2.8, meal_max=2.3)
    dehydration = Trial_Type(trial_name="Dehydration", trial_time=datetime.timedelta(hours=4), pp_dict={"Consumption": "Water", "Experiment": "Dehydrated"}, max=4.5, meal_max=1)
    dehydration_saline = Trial_Type(trial_name="Dehydration Saline", trial_time=datetime.timedelta(hours=2), pp_dict={"Consumption": "Saline", "Experiment": "Dehydrated"}, max=4)
    furosemide_saline = Trial_Type(trial_name="Furosemide Saline", trial_time=datetime.timedelta(hours=6), pp_dict={"Consumption": "Saline", "Experiment": "Furosemide"})
    furosemide_water = Trial_Type(trial_name="Furosemide Water", trial_time=datetime.timedelta(hours=6), pp_dict={"Consumption": "Water", "Experiment": "Furosemide"})
    sucrose = Trial_Type(trial_name="Sucrose", trial_time=datetime.timedelta(hours=23), pp_dict={"Consumption": "Fluid", "Experiment": "Sucrose", "Baseline": "Water"}, max=36.5)
    quinine = Trial_Type(trial_name="Quinine", trial_time=datetime.timedelta(hours=23), pp_dict={"Consumption": "Fluid", "Experiment": "Quinine", "Baseline": "Water"}, max=7)
    hypertonic_water = Trial_Type(trial_name="Hypertonic Water", trial_time=datetime.timedelta(hours=1), pp_dict={"Consumption": "Water", "Experiment": "Hypertonic", "Baseline": "Isotonic"})
    hypertonic_saline = Trial_Type(trial_name="Hypertonic Saline", trial_time=datetime.timedelta(hours=1), pp_dict={"Consumption": "Saline", "Experiment": "Hypertonic", "Baseline": "Isotonic"})
    total_saline = Trial_Type(trial_name="Total Saline", trial_time=datetime.timedelta(days=7), pp_dict={"Consumption": "Saline"})

    trial_types = {"Fast Refeed": fast_refeed, "Dehydration": dehydration, "Quinine": quinine, "Sucrose": sucrose,
                   "Total Saline": total_saline,
                   "Furosemide Saline": furosemide_saline,
                   "Furosemide Water": furosemide_water,
                   "Hypertonic Water": hypertonic_water, "Hypertonic Saline": hypertonic_saline,
                   "Dehydration Saline": dehydration_saline}
    had_trial = False
    with open(key_file) as f:
        mice = []
        for line in f:
            line = line.strip()
            if len(line.split(",")) < 3: # If line doesn't have two commas in it
                continue
            if line[0] != ">":
                if had_trial:
                    mice = []
                had_trial = False
                mouse, start_hopper, genotype = line.split(",")
                mice.append((mouse, int(start_hopper), genotype))
            else:
                had_trial = True
                line = line[1:]
                try:
                    trial_name, date_and_time, hopper_indent = line.split(",")
                except ValueError as ve:
                    print(ve)
                    print(line)
                    raise ValueError
                date_time = datetime.datetime.strptime(date_and_time, "%m/%d/%Y %H:%M")
                hopper_indent = int(hopper_indent)
                if ":" in trial_name:
                    trial_name, genotype_addendum = trial_name.split(":")
                    genotype_addendum = " (%s)" % genotype_addendum
                else:
                    genotype_addendum = ""
                if trial_name not in trial_types:
                    raise IOError("%s has an unknown trial type" % line)
                for m in mice:
                    t = Trial(start_time=date_time, time_delta=trial_types[trial_name].trial_time,
                              hopper=m[1]+hopper_indent, mouse=m[0], genotype=m[2], addendum=genotype_addendum)
                    trial_types[trial_name].trials.append(t)
    for trial in trial_types["Furosemide Saline"].trials:
        if "baseline" in trial.genotype.lower():
            t = Trial(start_time=trial.start_time, time_delta=total_saline.trial_time, hopper=trial.hopper, mouse=trial.mouse, genotype=trial.genotype)
            trial_types["Total Saline"].trials.append(t)
    return trial_types.values()


if __name__ == "__main__":
    analysis_folder = r"R:\Fillan\Parabrachial Ablations\BioDaq"
    key_file = os.path.join(analysis_folder, "Keyfile Lmx1b.txt")
    output_file = os.path.join(analysis_folder, "BioDaq Summary Data Lmx1b.xlsx")
    mc = MouseColor()
    ac = AverageColor()
    amount_consumed = 0
    trial_types = generate_trials_from_keyfile(key_file)
    trials = {}
    """
    Trials is a list of all the trials, indexed by the hopper and the date of the trial
    When we're going through the experiment file, this allows us to only check if each bout needs to be added to the
    subset of trials with the correct hopper, speeding up the import
    """
    for trial_type in trial_types:
        for trial in trial_type.trials:
            if trial.hopper not in trials:
                trials[trial.hopper] = {}
            if trial.start_time.date() not in trials[trial.hopper]:
                trials[trial.hopper][trial.start_time.date()] = []
            trials[trial.hopper][trial.start_time.date()].append(trial)
            trial_end_date = (trial.start_time + trial.time_delta).date()
            if trial_end_date != trial.start_time.date():
                current_date = trial.start_time.date()
                while current_date < trial_end_date:
                    current_date = current_date + datetime.timedelta(days=1)
                    if current_date not in trials[trial.hopper]:
                        trials[trial.hopper][current_date] = []
                    trials[trial.hopper][current_date].append(trial)
    """
    This code imports each line of the experiment file, and puts it into the correct trial
    Trial_Types is the parent class, and has a type of "Furosemide" or "Fast Refeed"...
    Each Trial_Type has a bunch of Trials, for each of the individual mice
    """
    search_experiment_files = os.path.join(analysis_folder, "*_experiment.csv")
    for experiment_file in glob.glob(search_experiment_files):
        print(experiment_file)
        with open(experiment_file) as f:
            for line in f:
                line = line.strip()
                if len(line) < 8:
                    continue
                date, time, bout_length, start, change, hopper = line.split(",")
                date_and_time = date + " " + time
                try:
                    date_time = datetime.datetime.strptime(date_and_time, "%m/%d/%Y %H:%M:%S")
                except ValueError as ve:
                    print(ve)
                    print(line)
                    raise ValueError
                bout_length = datetime.timedelta(seconds=int(bout_length))
                start = float(start)
                change = float(change)
                hopper = int(hopper)
                if change < -0.05 or change > 9.0:  # Filtering as with LogMeIn
                    continue
                if hopper not in trials:
                    continue
                if date_time.date() not in trials[hopper]:
                    continue
                for trial in trials[hopper][date_time.date()]:
                    trial.add_line(date_time, change, hopper, bout_length)
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        for trial in trial_types:
            #trial.correlate_with_stereology(stereology_file_path=r"R:\Fillan\Parabrachial Ablations\Lesion Mice Overview.xlsx")
            trial.graph_meal_analysis(writer)
            #trial.graph_mouse_bout_analysis(writer)
            #trial.graph_bout_analysis()
            #trial.graph_individual_time_series()
            #trial.graph_average_time_series(writer)
            #trial.graph_dot_plot()
