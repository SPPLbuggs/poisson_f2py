import main, sys
from setup_problem import setup_problem
import time
from mpi4py import MPI
import numpy as np

start = time.time()

comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

# grid spacing
nz,nr = 100,100

for i in range(1, len(sys.argv)/2+1):
        key = 2*i-1
        val = 2*i
        if sys.argv[key] == '-nz':
            nz = int(sys.argv[val])
        elif sys.argv[key] == '-nr':
            nr = int(sys.argv[val])

if rank == 0:
    print('Mesh size = {} by {}'.format(nz,nr))

if nz%size != 0:
    sys.exit('Error: Number of nodes in Z must be divisible by number of tasks')


# setup problem
neqn,z,type_z,r,type_r,phi,global_idx,local_idx,node_global = setup_problem(nz,nr,rank,size)

main.main(neqn,z,type_z,r,type_r,phi,local_idx,global_idx,node_global)

if rank == 0:
    print('Elapsed time = {:.2f} sec'.format(time.time()-start))


if rank == 0:
    Istart = 0
else:
    Istart = 1

if rank == size-1:
    Iend = len(z)
else:
    Iend = -1

#print phi[:,Istart:Iend].T
#print z[Istart:Iend]
phi_global = np.zeros([nz,nr],dtype='d')
z_global   = np.zeros(nz,dtype='d')
comm.Allgather([phi[:,Istart:Iend].T,MPI.DOUBLE], [phi_global, MPI.DOUBLE])
comm.Allgather([z[Istart:Iend].T,MPI.DOUBLE], [z_global, MPI.DOUBLE])

if rank == 0:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.contourf(r,z_global,phi_global)
    plt.colorbar()
    plt.xlabel('Radius (mm)')
    plt.ylabel('Z-axis (mm)')
    plt.title('Potential (V)')
    plt.show()
