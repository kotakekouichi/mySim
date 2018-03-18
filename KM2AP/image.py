import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-3, 3, 0.1);
y = np.sin(x);

plt.figure();
plt.plot(x, y);
plt.plot([1,2,3], [1,15., 4], 'r-')

x = np.random.randn(500);
y = np.random.randn(500);

heatmap, xedges, yedges = np.histogram2d(x, y, bins = 50);
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]];

print extent
plt.imshow(heatmap, extent = extent);
plt.show();
