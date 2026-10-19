import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from  Finite_Element import *
from Mesh import *

def geo_PEC_inclusions():
    pass





if __name__ == "__main__":
    plt.figure().clear()
    
    # --------------------
    # Geometry Description
    # --------------------
    
    # periodic cell
    
    lat_c = 1               # square lattice constant
    Lx = lat_c; Ly = lat_c 
    hx= 0.02;   hy= 0.02;   # "(smaller than the mesh size hmax)"
    Nx = Lx/hx; Ny = Ly/hy  # number of grid points in x and y directions
    xy_cell = np.array(
        np.zeros((1,Ny)) +
        np.array([i for i in range(hx,Lx,hx)]) +
        Lx*np.ones((1,Ny)) +
        np.array([i for i in range(Lx-hx,0,-hx)]) +
        np.array([i for i in range(Ly-hy,0,-hy)]) +
        np.zeros((1,Nx)) +
        np.array([i for i in range(hy,Ly,hy)]) +
        Ly*np.ones((1,Nx))
        ).T.conj()
    xy_box = np.array([0,Lx,0,Ly])          # computational domain
    xy_center = np.array([0.5*Lx,0.5*Ly])   # center of the inclusion
    
    # geomtery of the inclusions
    xy_inc = geo_PEC_inclusions()
    
    # -----------------
    # generate the mesh
    # -----------------
    
    p, t, t_num = mesh_generation_C4_PEC(xy_cell,xy_inc)
    
    # ------------------
    # build FEM matrices
    # ------------------
    
    M, Mx, My, Me, idx = FEM_matrices_C4_PEC_TE(p,t,xy_box,xy_inc)
    
    # -------------------------------------------------------------------------
    # 1. Solve the eigenvalues for the Bloch wave vectors along the M->X->T->M
    # -------------------------------------------------------------------------
    # Bloch wave vectors
    n_kappa = 100
    # kappa_vec1 = 1/lat_c * [ linspace(0,pi,n_kappa)',  zeros(n_kappa,1) ];          T -> X
    # kappa_vec2 = 1/lat_c * [ pi * ones(n_kappa,1),   linspace(0,pi,n_kappa)' ];     X -> M
    # kappa_vec3 = 1/lat_c * [ linspace(pi,0, n_kappa)', linspace(pi,0, n_kappa)' ];  M -> T
    
    # kappa_vec = [ kappa_vec1; kappa_vec2; kappa_vec3 ];
    # kappa_vec = kappa_vec3;

    kappa_vec1 = 1/lat_c * [ np.pi * np.ones((n_kappa,1)) + np.linspace(np.pi,0,n_kappa).T.conj() ]             # M -> X
    kappa_vec2 = 1/lat_c * [ np.linspace(np.pi,0,n_kappa).T.conj(),  np.zeros((n_kappa,1)) ]                    # X -> T
    kappa_vec3 = 1/lat_c * [ np.linspace(0,np.pi,n_kappa).T.conj(), np.linspace(0,np.pi,n_kappa).T.conj() ]     # T -> M

    kappa_vec = np.array([kappa_vec1, kappa_vec2, kappa_vec3])

    # kappa_vec1 = [ pi * ones(n_kappa,1),   linspace(pi,0,n_kappa)' ];   # X -> N
    # kappa_vec2 = [ linspace(pi,0,n_kappa)',  zeros(n_kappa,1) ];        # N -> T
    # kappa_vec3 = [ zeros(n_kappa,1),   linspace(0,pi,n_kappa)' ];       # T -> M
    # kappa_vec4 = [ linspace(0,pi,n_kappa)', pi * ones(n_kappa,1) ];     # M -> X
    
    # kappa_vec = 1/lat_c * [ kappa_vec1; kappa_vec2; kappa_vec3; kappa_vec4 ];

    # -------------------------------------------------------------------------
    # solve the first n_omega eigenvalues over each Bloch wave vector
    # -------------------------------------------------------------------------
    n_omega = 7
    omega = np.zeros((n_omega,np.size(kappa_vec,1)))
    for n in range(1, np.size(kappa_vec,1)):
        kappa = kappa_vec[n,:]
        A = M - 2j * ( kappa[1]*Mx + kappa[2]*My ) + (kappa*kappa.T.conj())*Me
        B = Me
        # FIX: lamb = np.linalg.eigvals(A,B,n_omega,'smallestabs')
        omega[:,n] = np.sqrt(lamb)*lat_c/(2*np.pi)

    # figure;
    # for n = 1:n_omega
    #     hold on; plot(real(omega(n,:)));
    # end