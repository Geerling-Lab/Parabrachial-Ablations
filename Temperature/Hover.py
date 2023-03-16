import matplotlib.pyplot as plt
import numpy as np; np.random.seed(1)


class Axis:
    def __init__(self, ax, annot):
        self.ax = ax
        self.annot = annot
        self.lines = []

    def __eq__(self, other):
        return self.ax == other


def contains(line, event):
    print(line)
    print(event)
    return True


class Line_Event_Handler:
    def __init__(self, axs):
        self.lines = []
        self.names = []
        self.axs = []
        for ax in axs.ravel():
            annot = ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
            annot.set_visible(False)
            self.axs.append(Axis(ax, annot))

    def add_line(self, ax, line):
        #line.set_contains(contains)
        for a in self.axs:
            if a == ax:
                a.lines.append(line)

    def hover(self, event):
        for a in self.axs:
            vis = a.annot.get_visible()
            if event.inaxes == a.ax:
                for line in a.lines:
                    cont, _ = line.contains(event)
                    if cont:
                        a.annot.xy = event.xdata, event.ydata
                        a.annot.set_text(line.get_label())
                        a.annot.get_bbox_patch().set_alpha(0.4)
                        a.annot.set_visible(True)
                        fig.canvas.draw_idle()
                    else:
                        if vis:
                            #a.annot.set_visible(False)
                            fig.canvas.draw_idle()


fig, ax = plt.subplots(nrows=2, ncols=2)
line_event_handler = Line_Event_Handler(ax)
for i in range(ax.shape[0]):
    for j in range(ax.shape[1]):
        for k in range(3):
            x = np.sort(np.random.rand(15))
            y = np.sort(np.random.rand(15))
            name = str(k)

            line, = ax[i, j].plot(x, y, label=name)
            line_event_handler.add_line(ax[i, j], line)

fig.canvas.mpl_connect("motion_notify_event", line_event_handler.hover)

plt.show()