import matplotlib.pyplot as plt
import numpy as np
import math

def graph_wobble_plot(ax, input, scale=1):
    """
    :param input: dictionary with column titles and colors, lists; e.x.
    input['Cas3'] = ('r', [1, 2, 3, 4])
    input['mCh'] = ('g', [1, 1, 2, 2])
    :param scale: step size between dots on the same row
    :return:
    """
    column_spacing = 30
    dot_spacing = 0.1
    x = column_spacing
    columns_labels = []
    for column, items in input.items():
        color, values = items
        values = np.array(values)
        min = np.min(values) // scale * scale
        max = np.max(values) // scale * scale + scale + scale + scale
        bins = np.arange(min, max, step=scale)
        hist, bins = np.histogram(values, bins=bins)
        x_value_lookup = {}
        for i in range(hist.size):
            try:
                assert math.floor(bins[i] / scale + 0.5) not in x_value_lookup
            except:
                print(i, math.floor(bins[i] / scale + 0.5))
                raise AssertionError
            x_value_lookup[math.floor(bins[i] / scale + 0.5)] = x - dot_spacing * hist[i] / 2
        x_values = []
        for value in values:
            y_ = math.floor(value // scale)
            x_values.append(x_value_lookup[y_])
            x_value_lookup[y_] += dot_spacing
        ax.scatter(np.array(x_values), values, color=color)
        columns_labels.append(column.replace("+", "\n").replace(" (", "\n").replace(")", ""))
        x += column_spacing
    ax.set_xticks(np.arange(column_spacing, x, column_spacing))
    ax.set_xticklabels(columns_labels)


if __name__ == "__main__":
    input = {}
    input['Control mCh'] = ('b', np.array([-4, -3, -2, -1, -1, 0, 1, 2, 3, 4]))
    input['Experiment Cas3'] = ('r', list(np.random.normal(loc=5, scale=1, size=200)))
    input['Experiment mCh'] = ('c', list(np.random.normal(loc=1, scale=1, size=500)))
    input['Control Cas3'] = ('g', list(np.random.normal(loc=4, scale=1, size=50).astype(np.int8)))
    fig, ax = plt.subplots()
    graph_wobble_plot(ax, input, scale=0.1)
    plt.show()