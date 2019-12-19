# alpha-beta-core
Graph format:
a .meta file contains the number of edges and the number of nodes in each part. See data/example.meta
a .e file contains all the edges. See data/example.e

To build index with BasicDecom:
./abcore -BasicDecom path_to_graph (for example: ./abcore -BasicDecom ../data/example)

To build index with ComShrDecom:
./abcore -ComShrDecom path_to_graph (for example: ./abcore -ComShrDecom ../data/example)

To build index with ParallelDecom:
./abcore -ParallelDecom path_to_graph (for example: ./abcore -ParallelDecom ../data/example)

To query alpha beta core using BiCoreIndex:
./abcore -Query path_to_graph alpha beta (for example: ./abcore -Query ../data/example 2 3)

To insert edge with BiCore-Index-Ins:
./abcore -BiCore-Index-Ins path_to_graph vertex_1 vertex_2 (for example: ./abcore -BiCore-Index-Ins ../data/example 3 6)

To remove edge with BiCore-Index-Rem:
./abcore -BiCore-Index-Rem path_to_graph vertex_1 vertex_2 (for example: ./abcore -BiCore-Index-Rem ../data/example 3 7)

To insert edge with BiCore-Index-Ins*:
./abcore -BiCore-Index-Ins* path_to_graph vertex_1 vertex_2 (for example: ./abcore -BiCore-Index-Ins* ../data/example 3 6)

To remove edge with BiCore-Index-Rem*:
./abcore -BiCore-Index-Rem* path_to_graph vertex_1 vertex_2 (for example: ./abcore -BiCore-Index-Rem* ../data/example 3 7)
