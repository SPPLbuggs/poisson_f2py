import numpy as np
import sys

def setup_problem(nz,nr,rank,size):
    
    # Physical Params:
    length   = 3e-3
    radius   = 3e-3
    e_start  = 0.
    e_stop   = 5e-4
    
    # RF Params:
    v_l = -1000
    v_r = 0
    rf_fq   = 1e6
    final_time = 10e-6
    
    # Simulation Params:
    init = 1e14
    Te_init = 0.5
    len_0 = 1e-3
    phi_0 = 1e3
    tau_0 = 1e-6
    n_0 = 1e16
    
    if nz > 1 and nr > 1:
        stencil = 5
    else:
        stencil = 3
    
    # linear spacing
    z = np.linspace(0,length,nz)/len_0
    r = np.linspace(0,radius,nr)/len_0
    
    # hyperbolic tangent spacing
#    z = np.linspace(-1.5,1.5,nz)
#    z = np.tanh(z)*0.5+0.5
#    z = z-z[0]
#    z = z/z[-1]
#    z = z*length/len_0
    
#    r = np.linspace(-1.5,0,nz)
#    r = np.tanh(r)*0.5+0.5
#    r = r-r[0]
#    r = r/r[-1]
#    r = r*radius/len_0
    
    ne = np.ones([nr,nz])
    ni = ne
    nm = ne
    nt = ne*Te_init
    phi = np.zeros([nr,nz],dtype='d',order='F')

    type_z = np.zeros([nr,nz],'int')
    type_r = np.zeros([nr,nz],'int')
    
    if nz > 1:
        for i in range(nr):
            for j in range(nz):
                if i == nr - 1:
                    type_r[i,j] = 1
                elif i == 0:
                    type_r[i,j] = -1
            
                if j == nz-1:
                    if e_start/len_0 <= r[i] <= e_stop/len_0:
                        type_z[i,j] = 2
                        phi[i,j] = v_r/phi_0
                    else:
                        type_z[i,j] = 1
                elif j == 0:
                    if e_start/len_0 <= r[i] <= e_stop/len_0:
                        type_z[i,j] = -2
                        phi[i,j] = v_l/phi_0
                    else:
                        type_z[i,j] = -1
    
    global_idx  = np.zeros([nz*nr,2],'int',order='F')
    local_idx   = np.zeros([nz*nr,2],'int',order='F')
    node_global = np.zeros([nr,nz],'int',order='F')
    
    neqn = 0
    j_loc = 0
    temp = 0
    for j in range(nz):
        for i in range(nr):
            if abs(type_z[i,j]) != 2:
                global_idx[neqn,:] = [i+1,j+1]
                local_idx[neqn,:] = [i+1,j_loc+1]
                neqn += 1
                node_global[i,j] = neqn
        j_loc += 1
        if j_loc > nz/size-1+temp:
            temp = 1
            j_loc = 1
    
    z   = z[max(0,nz/size*rank-1):min(nz,nz/size*(rank+1))+1]
    phi = phi[:,max(0,nz/size*rank-1):min(nz,nz/size*(rank+1))+1]
    
    return neqn,z,type_z,r,type_r,phi,global_idx,local_idx,node_global

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit('Usage: python setup_problem.py -nx # -nr #')
    
    for i in range(1, len(sys.argv)/2+1):
        key = 2*i-1
        val = 2*i
        if sys.argv[key] == '-nx':
            nx = int(sys.argv[val])
        elif sys.argv[key] == '-nr':
            nr = int(sys.argv[val])
        else:
            sys.exit('Usage: python setup_problem.py -nx # -nr #')
    
    neqn,x,type_x,r,type_r,phi,node_idx,node_global = setup_problem(nx,nr)
    print 'x'
    print x
    print 'type_x'
    print type_x
    print 'r'
    print r
    print 'type_r'
    print type_r
    print 'node_idx'
    print node_idx
    print 'node_global'
    print node_global
