from PIL import Image, ImageFilter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.ndimage

im = Image.open(r"R:\Fillan\BoutonCount\Fig 1 Pipeline\Example Zoom\Zoom.png")
im = im.convert("L")
arr = np.array(im)
ar = arr.astype(np.float)
ar = scipy.ndimage.filters.gaussian_filter(ar, sigma=1)
X = np.arange(0, ar.shape[1], 1)
Y = np.arange(0, ar.shape[0], 1)
X, Y = np.meshgrid(X, Y)
fig = plt.figure()
ax = plt.gca(projection="3d")
ax.grid(b=None)
ax.set_zlim3d(50, 256)
ax.plot_surface(X, Y, ar, rstride=5, cstride=5, cmap="summer", antialiased=True, alpha=0.5, linewidth=0.5, edgecolors='k')
plt.axis("off")
plt.show()