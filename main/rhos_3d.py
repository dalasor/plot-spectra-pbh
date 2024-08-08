import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
import matplotlib.ticker as ticker

# Определение точек для осей X и Y
y = np.log10([1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14])
x = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])  # [::-1]

# Создание сетки координат X и Y
X, Y = np.meshgrid(x, y)

# Расчет значений Z для каждой пары (X, Y)
# Здесь для демонстрации используем функцию Z = X^2 + Y^2, но вы можете использовать свои данные
z = np.array([[1.0e-20, 3.7e-14, 2.8e-9, 1.0e-7, 3.9e-6, 4.8e-4, 4.5e-3, 3.3e-3, 4.1e-4],
     [1.4e-7, 3.0e-8, 2.5e-7, 3.5e-6, 1.7e-4, 7.4e-3, 1.7e-2, 4.8e-3, 3.2e-4],
     [3.4e-8, 1.9e-7, 9.1e-7, 1.1e-5, 6.9e-4, 5.6e-3, 3.6e-3, 4.7e-4, 1.9e-5],
     [2.1e-9, 5.0e-9, 2.3e-8, 4.1e-7, 9.8e-6, 1.5e-5, 3.4e-6, 2.3e-7, 6.0e-9],
     [5.4e-7, 7.2e-7, 3.5e-7, 8.5e-6, 2.7e-5, 1.0e-5, 2.1e-6, 4.2e-8, 8.1e-10],
     [2.7e-7, 4.0e-6, 3.8e-5, 1.8e-4, 8.7e-5, 1.2e-5, 6.6e-7, 1.8e-8, 2.9e-10],
     [1.7e-9, 2.0e-8, 5.2e-8, 2.0e-8, 2.9e-9, 2.1e-10, 8.9e-12, 2.1e-13, 2.7e-15]])

# z = z.T

Z = np.log10(z)

# Преобразование сетки координат X и Y в одномерные массивы
X1 = X.ravel()
Y1 = Y.ravel()
Z1 = Z.ravel()
# Создание новой сетки для интерполяции
x_new = np.linspace(x.min(), x.max(), 500)
y_new = np.linspace(y.min(), y.max(), 500)
X_new, Y_new = np.meshgrid(x_new, y_new)
# Интерполяция данных
Zi = griddata((X1, Y1), Z1, (X_new, Y_new), method='cubic')


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# log_norm = LogNorm(vmin=Z.min(), vmax=Z.max())

# Построение графика
surf = ax.plot_surface(X_new, Y_new, Zi, cmap='BuPu')  #, norm=log_norm)

# Добавление цветовой шкалы
cbar = fig.colorbar(surf, ax=ax, shrink=0.9, aspect=10, pad=0.15)  #, norm=log_norm)

ax.set_ylabel(r'$\log_{10}(M_{\mathrm{PBH}}/\text{г})$')
ax.set_xlabel(r'$\sigma$')
ax.set_zlabel(r'$\log_{10}(\rho_{\nu}^{\mathrm{PBH}}/\rho_{\mathrm{C\nu B}})$')

# ax.set_xlim(ax.get_xlim()[::-1])
ax.set_ylim(ax.get_ylim()[::-1])
ax.set_zlim([-12.5, -1])
# ax.set_xscale('log')
# ax.set_zscale('log')

plt.show()
fig.clf()


# #### 2d plot #### #

f, ax = plt.subplots(1, 1)
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel(r'$\log_{10}(M\,/\,$г$)$')
ax.set_ylabel(r'$\sigma$')
ax.set_xlim([8.0, 14.0])
ax.set_ylim([0.0, 4.0])

vmin = -15.0
vmax = -1.8

# создаем маску для значений Zi > -2
mask = Zi > -2.0
# применяем маску к X_new и Y_new
x_values = Y_new[mask]
y_values = X_new[mask]
# создаем массив точек для ConvexHull
points = np.column_stack((x_values, y_values))

# Вычислите выпуклую оболочку этих точек
hull = ConvexHull(points)

# создаем полигон для отображения контура
hull_points = points[hull.vertices]
polygon = Polygon(hull_points, fill=True, facecolor='#520150', edgecolor='black', closed=True)

# создаем маску для значений Zi > -3
mask1 = Zi > -3.0
# применяем маску к X_new и Y_new
x_values1 = Y_new[mask1]
y_values1 = X_new[mask1]
# создаем массив точек для ConvexHull
points1 = np.column_stack((x_values1, y_values1))

# Вычислите выпуклую оболочку этих точек
hull1 = ConvexHull(points1)

# создаем полигон для отображения контура
hull_points1 = points1[hull1.vertices]
polygon1 = Polygon(hull_points1, fill=True, facecolor='#650762', edgecolor='black', closed=True)

color_mesh = ax.imshow(Zi.T, extent=(y.min(), y.max(), x.min(), x.max()),
                       origin='lower', cmap='BuPu', norm=Normalize(vmin=vmin, vmax=vmax))

f.colorbar(color_mesh, ax=ax, label=r'$\log_{10}(\rho_{\nu}^{\mathrm{PBH}}/\rho_{\mathrm{C\nu B}})$')

# plt.gca().add_patch(polygon1)
# plt.gca().add_patch(polygon)

# создаем контурный график
x_new = np.linspace(y.min(), y.max(), 500)
y_new = np.linspace(x.min(), x.max(), 500)
contour = plt.contour(x_new, y_new, Zi.T, levels=[-3.0, -2.0], colors='lightcyan', linewidths=2, linestyles='dashed')
plt.clabel(contour, inline=True, fontsize=10, colors='aqua')
ax.scatter(9.38, 2.78, color='lightcyan')
# plt.annotate(r'$-1.63$', (9.38, 2.78), fontsize=14, color='aqua')
ax.text(9.1, 2.9, r'$-1.63$', ha='center', fontsize=8, color='aqua')

# Находим максимум
print(np.max(Zi))
# Находим индекс максимального значения в Z
index = np.argmax(Zi)
# Преобразуем индекс обратно в координаты двумерной матрицы
i, j = np.unravel_index(index, Zi.shape)
print(i, j)
xm = X_new[i, j]
ym = Y_new[i, j]
print(xm)
print(ym)

plt.show()
f.clf()
