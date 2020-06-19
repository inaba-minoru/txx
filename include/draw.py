from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure()
# ax = Axes3D(fig)
# X = np.arange(-4, 4, 0.25)
# Y = np.arange(-4, 4, 0.25)
# X, Y = np.meshgrid(X, Y)
# R = X ** 2 + Y ** 2
# Z = R

# fig = plt.figure()
# ax = Axes3D(fig)

# c = np.arange(0, 0.1, 0.001)
# eta = np.arange(99, 101, 0.05)

# c, eta = np.meshgrid(c, eta)

# g = np.sqrt(eta ** 2 - 1 + c ** 2)

# Z = (
#     (g - c) ** 2
#     * (1 + (c * (g + c) - 1) ** 2 / (c * (g - c) + 1) ** 2)
#     / (2 * (g + c) ** 2)
# )

# # # 具体函数方法可用 help(function) 查看，如：help(ax.plot_surface)
# ax.plot_surface(c, eta, Z, rstride=1, cstride=1, cmap="rainbow")

# plt.show()

# eta = 100
# c = 0.010002
# g = np.sqrt(eta ** 2 - 1 + c ** 2)

x = np.arange(0, 1, 0.01)
# y = np.arctan(0.01 * np.sqrt(x) / np.sqrt(1 - x))
y = np.arcsin(x)
plt.title("sine wave form")
# 使用 matplotlib 来绘制点
plt.plot(y, x)
plt.show()
