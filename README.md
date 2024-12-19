# Friends-of-Friends

- Friends-of-Friends (FoF) is a percolation algorithm used to identify structures in the Universe based on the physical proximity of particles in the dark matter (detected only by their gravitational effects). It can be used, for example, to find the formation of new stars, galaxies, etc. 

- Astronomy accumulates unprecedented amount of data, which could benefit from High Performance Computing environments to be processed. I parallelized a version of the FoF clustering algorithm for CPU and NVIDIA GPUs using OpenACC
 
- [Poster](GHC-poster-ACMSRC2019.pdf) about this work won 1st place in the Student Research Competition (ACM SRC) at the Grace Hopper Celebration in 2019 in the undergraduate category. 
- Code and Jupyter notebook for the submission [here](./code/acmsrc/), note that it was using PGI compiler - it needs to be updated!
- [Presentation](GHC19-Solorzano-presentation.pptx) for GHC19

- The algorithm works by receiving as input a file with the particles n-body cosmological simulations, and a linking length parameter (R). It considers a sphere of radius R around each particle of the total input dataset, so all particles in this range belong to the same halo and are considered friends. The procedure continues by defining a sphere around each friend so "any friend of my friend is my friend", until no new friends are detected.


| Complexity | Serial | Parallel | Device | Libraries |
| ---- | ----- |----- |----- |----- |
| O (NlogN) | X | | CPU | | 
| O (N^2) | X | | CPU | | 
| O (N^2) | X | | CPU | | 
| O (N^2) | | X | CPU | OpenMP | 
| O (N^2) | | X | CPU + GPU | OpenACC | 

## Data

- Input files from Virgo Consortium:
http://www.mpa-garching.mpg.de/Virgo/data_download.html

- Controlled entries from "Comparativo de Resultados FoF N2 X FoF NlogN" file.
- It should expect: 4 groups found for file 1
		    1 group found for file 2

## Running

### NlogN

- Based on hirarchical oct-tree to represent the 3D space dimension, which facilitates relabeling particles belong to another group
Set the radius (R) in main.cpp.

```
g++ main.cpp grupo.cpp Corpo.cpp No.cpp Segmento.cpp Tupla.cpp -o teste
./teste file
```

### OpenMP

```
cd ./code/openmp/
gcc -fopenmp FoFOn2-openmp.cpp -o run
export OMP_NUM_THREADS=4
./run
```

### MPI

```
cd ./code/mpi/
mpicc FoFmpi.c -o run
mpiexec -np 2 ./run ./../../data/virgo/Virgo_1
# Enter the percolation radius (R)
```

### OpenACC

```
pgc++ -g FoF0n2-openacc.cpp -o teste -ta=nvidia:nvidia -Minfo
```

```
wget https://developer.download.nvidia.com/hpc-sdk/24.11/nvhpc_2024_2411_Linux_x86_64_cuda_12.6.tar.gz
```