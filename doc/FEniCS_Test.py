from dolfin import *
from mpi4py import MPI

parameters["ghost_mode"] = "none"

mesh = UnitSquareMesh(20,20)
Q = FunctionSpace(mesh, "CG", 1)
F = Function(Q)


hdf = HDF5File(mesh.mpi_comm(), "a.h5", "w")

for i in range(10):
    hdf.write(F, "fun", i)

# Delete HDF5File object, closing file
del hdf

list_linear_solver_methods()
list_krylov_solver_preconditioners()
# Read functions back in

hdf = HDF5File(mesh.mpi_comm(), "a.h5", "r")
attr = hdf.attributes("fun")
nsteps = attr['count']
for i in range(nsteps):
    dataset = "fun/vector_%d"%i
    attr = hdf.attributes(dataset)
    print ('Retrieving time step:', attr['timestamp'])
    hdf.read(F, dataset)
