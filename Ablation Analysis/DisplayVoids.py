import matplotlib.pyplot as plt
import matplotlib.patches
import numpy as np
import glob
import os

if __name__ == "__main__":
    MAX_FRAMES = 3 * 60 * 60 * 2  # 2 hours at 3 frames/s
    directory = r"R:/Fillan/Parabrachial Arousal/MVT"
    os.chdir(directory)
    files = glob.glob("*Micturition.csv")
    num_files = len(files)
    height = (num_files + 1) * 100
    width = 120

    mCh_mice = ["02547", "02548", "02549", "02550", "02551", "02905", "02909", "02910", "02910", "03433", "03434", "04057", "04058", "04055", "04056", "04146", "04147", "04342", "04444"]
    # Genotype is False is mCh and True if Cas3
    files = [(x, False) if x.split(" ")[0] in mCh_mice else (x, True) for x in files]
    files.sort(key=lambda e: (e[1], e[0]), reverse=True)
    fig, ax1 = plt.subplots()
    ax1.set_xlim([0, width])
    ax1.set_ylim([0, height])
    ax1.xlabel = "Minute"
    y_labels = []
    y_sizes = []
    y = []
    for i, (file, genotype) in enumerate(files):
        line_y = i * 100 + 50
        name = file.split(" ")[0]
        if not genotype:
            c = [0, 0, 1]
        else:
            c = [1, 0, 0]
        if name in ["02547", "02548", "02549", "02550", "02551", "02883", "02885", "02886", "02905", "02910", "03433", "03434", "03558", "03559", "04055", "04056", "04146", "04147", "04149", "04148", "04339"]:
            name += " (F)"
        else:
            name += " (M)"
        y_labels.append(name)
        y.append(line_y)
        total_volume = 0
        with(open(os.path.join(directory, file))) as f:
            for line in f:
                words = line.split(",")
                if words[0].lower() == "frame":  # Remove header
                    continue
                else:
                    if len(words) > 2:
                        if len(words[2]) > 1:
                            size = float(words[2])
                        else:
                            size = 0.1
                    else:
                        size = 0.1
                    total_volume += size
                    line_height = size * 3
                    line_x = int(words[0]) / 180
                    line_height = 1
                    color = c + [size / (size + 3)]
                    ax1.plot([line_x, line_x], [line_y - line_height, line_y + line_height], color=color, linewidth=15, solid_capstyle='round')
        if i % 2 == 0:
            rect = matplotlib.patches.Rectangle((0, line_y - 50), 120, 100, facecolor=(0, 0, 0, 0.05))
            ax1.add_patch(rect)
        y_sizes.append("%.2f" % total_volume)
    ax1.set_xlabel("Minutes")
    ax1.set_ylabel("Mouse")
    ax1.set_yticks(y)
    ax1.set_yticklabels(y_labels)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Total Volume (cm^2)")
    ax2.set_ylim([0, height])
    ax2.set_yticks(y)
    ax2.set_yticklabels(y_sizes)
    plt.show()
