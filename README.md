# Power Flow Computation Using Newton-Raphson Method
## flow1
It is a unconstrained problem.
### source files
FLOW.H: project hearder
FLOW.C: five supportive functions in computation(PQ/PV node static)
main.cpp: main iteration algorithm and I/O
data.txt: 5 node power system network
dataFormat.txt: input data protocol
### output
exe cmd: iteration details & final results
iteration.txt: iteration time - error
## flow2
It is a constrained problem.
### source files
FLOW.H: project hearder
FLOW.C: five supportive functions in computation(based on PQ/PV node vector)
main.cpp: PQ/PV node vector, main iteration algorithm and I/O
data_2.txt: 8 node power system network with reactive power constraints
dataFormat.txt: input data protocol (same as "flow1/dataFormat.txt)
### output
same as "flow1"