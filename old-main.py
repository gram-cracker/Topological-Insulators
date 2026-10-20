import numpy as np
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#import math as math

#Fourier–Bloch Spectral Method
#Variables Lattice

# pi = np.pi
# N=60

# t = np.linspace(0, 1, N)

# """
# l_1 = 1
# l_2 = 1
# """
# a =1
# e1 = np.array([1,0])
# e2 = np.array([0,1])

# # Reciprocal lattice vectors
# k1 = np.linalg.solve([e1, e2], [2*pi, 0])
# k2 = np.linalg.solve([e1, e2], [0, 2*pi])

# M_Point = 1/2 * k1 + 1/2 * k2



# # Symmetry matrices
# R = np.array([
#     [0, -1],
#     [1, 0]
# ])

# F = np.array([
#     [-1,0],
#     [0,1]
# ])

# R2 = np.array([
#     [0, 1],
#     [-1, 0]
# ])

# # Brillouin zone vertices
# vertices = np.array([
#     M_Point,
#     np.linalg.matmul(F, M_Point),
#     np.linalg.matmul(R, M_Point),
#     np.linalg.matmul(F, np.linalg.matmul(R, M_Point))
# ])



# # Small obstacle boundary 
# epsilon = 0.5
# # n_obstacle = 24
# obstacle = np.concatenate((
#     np.linspace([0, epsilon], [np.sqrt(1)/2 * epsilon, -epsilon/2], N // 2),
#     np.linspace([np.sqrt(1)/2 * epsilon, -epsilon/2], [-np.sqrt(1)/2 * epsilon, -epsilon/2], N // 2),
#     np.linspace([-np.sqrt(1)/2 * epsilon, -epsilon/2], [0, epsilon], N // 2)
# ))

# # Orthogonal basis
# num_phi = 10
# phi = []
# for k in range(-num_phi, num_phi+1):
#     phi.append(np.exp(t * complex(0,1) * k))
# phi = 1/np.sqrt(2*pi) * np.array(phi)




# def generate_lattice(b1, b2, l1_range, l2_range):
#     """
#     Generate lattice points from two basis vectors.

#     Parameters:
#         b1 (array-like): First basis vector (length 2).
#         b2 (array-like): Second basis vector (length 2).
#         m_range (tuple): (min_m, max_m) integer range for m.
#         n_range (tuple): (min_n, max_n) integer range for n.

#     Returns:
#         np.ndarray: Array of shape (num_points, 2) with lattice coordinates.
#     """
#     b1 = np.array(b1, dtype=float)
#     b2 = np.array(b2, dtype=float)

#     if b1.shape != (2,) or b2.shape != (2,):
#         raise ValueError("Basis vectors must be 2D.")

#     l1_vals = np.arange(l1_range[0], l1_range[1] + 1)
#     l2_vals = np.arange(l2_range[0], l2_range[1] + 1)

#     points = []
#     for m in l1_vals:
#         for n in l2_vals:
#             point = m * b1 + n * b2
#             points.append(point)

#     return np.array(points)



# lattice_points = generate_lattice(e1, e2, (-5, 5), (-5, 5))

# # Plotting
# plt.figure(figsize=(6, 6))
# plt.scatter(lattice_points[:, 0], lattice_points[:, 1], c='blue', s=20)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.title("2D Lattice from Basis Vectors")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.grid(True, linestyle='--', alpha=0.5)
# plt.show()




# #1D Case
# plot2 = plt.figure(1,clear=True).add_subplot(111)
# k = np.linspace(-pi/a, pi/a, N)
# for i in range(3):
#     G = 2 *i * pi/a
#     l1 = (k+G)**2
#     l2 = (k-G)**2
#     plot2.plot(k, l1, label=f'Band {i+1}', color='C'+str(i))
#     plot2.plot(k, l2, label=f'Band {i+1}', color='C'+str(i))
    
# plt.show()


#Use Finite Difference for 1D


#U2 =0 #second derivative
#u = x #something x

h = 1000
interval = (0,1)
v = np.linspace(-1,-1,h-2)
# v = np.zeros(h-2)
opp = np.zeros((h-2,h-2))
for i in range(h-2):
    for j in range(h-2):
        if i==j:
            opp[i,j] = -2/(h**2) + v[i]
        elif abs(i-j)==1:
            opp[i,j] = 1/(h**2)
        else:
            opp[i,j] = 0

print(opp)

vals, vecs = np.linalg.eig(opp)

for i in range(len(vecs)):
    vecs[i] = [0] + vecs[i] + [0]

print(vals)
print(vecs)

x = np.linspace(interval[0],interval[1],len(vecs[0]))
fig, ax = plt.subplots()
for i in range(min(len(vecs),5)):
    ax.scatter(x, vecs[i], label=f'Vector {i+1}')
plt.show()