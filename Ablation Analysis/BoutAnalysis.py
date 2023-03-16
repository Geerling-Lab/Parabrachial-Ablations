import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import font_manager


def plot_df(df, ax, baseline=True, minimum_bout_size=0.00):
    mCh_bouts = []
    Cas3_bouts = []
    for column in df:
        if "Cas3" in column:
            if baseline == ("Baseline" in column):
                Cas3_bouts.append(df[column])
        if "mCh" in column:
            if baseline == ("Baseline" in column):
                mCh_bouts.append(df[column])
    Cas3_n = len(Cas3_bouts) - 1
    mCh_n = len(mCh_bouts) - 1
    Cas3_bouts = np.array(Cas3_bouts).flatten()
    Cas3_bouts = Cas3_bouts[~np.isnan(Cas3_bouts)]
    mCh_bouts = np.array(mCh_bouts).flatten()
    mCh_bouts = mCh_bouts[~np.isnan(mCh_bouts)]
    print("%i Cas3 bouts" % Cas3_bouts.size)
    print("%i mCh bouts" % mCh_bouts.size)
    U, p = scipy.stats.ttest_ind(Cas3_bouts, mCh_bouts, equal_var=False)
    maximum_bout_size = max(np.nanmax(Cas3_bouts), np.nanmax(mCh_bouts))
    bins = np.arange(0, round(maximum_bout_size, 2), 0.01)
    cas3_bins, _ = np.histogram(Cas3_bouts[Cas3_bouts >= minimum_bout_size], bins=bins)
    mCh_bins, _ = np.histogram(mCh_bouts[mCh_bouts >= minimum_bout_size], bins=bins)
    cas3_bins = cas3_bins / Cas3_bouts.size
    mCh_bins = mCh_bins / mCh_bouts.size
    bins = bins[:-1]
    ax.bar(x=bins, height=cas3_bins, width=.01, color=[1, 0.3, 0.3], label="Cas3 (n=%i)" % Cas3_n, zorder=2)
    ax.bar(x=bins, height=-mCh_bins, width=.01, color=[0.5, 0.5, 0.5], label="mCh (n=%i)" % mCh_n, zorder=2)
    ax.plot([], [], c='w', label="p=%.2f" % p)
    ax.axhline(y=0, c='k')
    y_axis = max(np.max(cas3_bins), np.max(mCh_bins)) * 1.2
    ax.plot([np.mean(Cas3_bouts), np.mean(Cas3_bouts)], [0, y_axis], color=[0, 0, 1, 0.2], zorder=1)
    ax.plot([np.mean(mCh_bouts), np.mean(mCh_bouts)], [0, -y_axis], color=[0, 0, 1, 0.2], zorder=1)
    ax.set_ylim([-y_axis, y_axis])
    ax.set_xlim([0, maximum_bout_size])
    ax.set_yticks([])
    ax.legend(prop={'size': 8}, loc=1)


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 2

    genotypes = ["Vglut2", "Foxp2", "Lmx1b"]
    fig, axs = plt.subplots(9, 3, sharex=True)
    for i, genotype in enumerate(genotypes):
        x = i
        file_name = r"R:\Fillan\Parabrachial Ablations\BioDaq\BioDaq Summary Data %s.xlsx" % genotype
        axs[0, x].set_title(genotype)
        plt.subplots_adjust(wspace=0, hspace=0)
        df = pd.read_excel(file_name, sheet_name="Fast Refeed Bouts")
        plot_df(df, axs[0, x], baseline=True)
        plot_df(df, axs[1, x], baseline=False)
        df = pd.read_excel(file_name, sheet_name="Dehydration Bouts")
        plot_df(df, axs[2, x], baseline=True)
        plot_df(df, axs[3, x], baseline=False)
        df = pd.read_excel(file_name, sheet_name="Furosemide Water Bouts")
        plot_df(df, axs[4, x], baseline=False)
        df = pd.read_excel(file_name, sheet_name="Furosemide Saline Bouts")
        plot_df(df, axs[5, x], baseline=True)
        plot_df(df, axs[6, x], baseline=False)
        df = pd.read_excel(file_name, sheet_name="Dehydration Saline Bouts")
        plot_df(df, axs[7, x], baseline=False)
        df = pd.read_excel(file_name, sheet_name="Sucrose Bouts")
        plot_df(df, axs[8, x], baseline=False)

    axs[0, 0].set_ylabel("Food\nBaseline")
    axs[1, 0].set_ylabel("Food\nFast Refeed")
    axs[2, 0].set_ylabel("Water\nBaseline")
    axs[3, 0].set_ylabel("Water\nDehydration")
    axs[4, 0].set_ylabel("Water\nFurosemide")
    axs[5, 0].set_ylabel("Saline\nBaseline")
    axs[6, 0].set_ylabel("Saline\nFurosemide")
    axs[7, 0].set_ylabel("Saline\nDehydration")
    axs[8, 0].set_ylabel("Sucrose")
    plt.xlabel("Bout Size (g)")
    plt.xlim([0, .6])
    plt.show()