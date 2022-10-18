# Setting docker up

There are Dockerfiles for images based on Ubuntu 22 ("impish") with OpenMPI (for non-Cray clusters) or MPICH (for Cray clusters like Daint), IntelMKL, VASP, wannier90 3.1, triqs 3.x.x, and Triqs MaxEnt included.

## Building the docker image
The Dockerfile is built with this command, where `<version name>` could be `3.0.0`:
```
docker build -t triqs_mpich:<version name> -f mpich_dockerfile ./
docker build -t triqs_openmpi:<version name> -f openmpi_dockerfile ./
```
Note that you need a vasp version as archive (`vasp.6.3.0.tgz`) in this directory to make the vasp CSC calculation work. Otherwise, you can just remove the lines in the docker file referring to the vasp installation.

## Running docker images with sarus in HPC clusters

To use your locally saved docker image, you first save it
```
docker save --output=docker-triqs.tar <image name>:<version name>
```
and then upload the tar to the HPC cluster and then load it into [sarus](https://user.cscs.ch/tools/containers/sarus/) via
```
sarus load docker-triqs.tar <image name>:<version name>
```
then you can run it within sarus.

# Running a docker container locally

You can start a docker container with either of these commands
```
docker run --rm -it -u $(id -u) -v ~$PWD:/work <image name>:<version name> bash
docker run --rm -it --shm-size=4g -e USER_ID=`id -u` -e GROUP_ID=`id -g` -p 8378:8378 -v $PWD:/work <image name>:<version name> bash
```
where the second command adds some important flags.
- The -e flags will translate your current user and group id into the container and make sure writing permissions are correct for the mounted volumes.
- The option --shm-size, which increases shared memory size.
This is hard coded in Docker to 64m and is often not sufficient and will produce SIBUS 7 errors when starting programs with mpirun! (see also https://github.com/moby/moby/issues/2606).
- The '-v' flags mounts a host directory as the docker directory given after the colon.
This way docker can permanently save data; otherwise, it will restart with clean directories each time.
Make sure you mount the directory where you save your data.
- All the flags are explained in 'docker run --help'.

Inside the docker, you can normally execute program. To run solid_dmft, for example, use
```
mpirun -n 4 solid_dmft
```
To start a jupyter-lab server from the current path, use
```
jupyter.sh
```
All these commands you can execute directly by just adding them to the `docker run ... bash` command with the `-c` flag, e.g.
```
docker run --rm -it --shm-size=4g -e USER_ID=`id -u` -e GROUP_ID=`id -g` -p 8378:8378 -v $PWD:/work <image name>:<version name> bash -c 'cd /work && mpirun -n 4 solid_dmft
```
