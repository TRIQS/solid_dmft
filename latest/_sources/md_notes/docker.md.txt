
# Docker

There are Dockerfiles for images based on Ubuntu 20 ("focal") with OpenMPI (for non-Cray clusters) or MPICH (for Cray clusters like Daint), IntelMKL, VASP, wannier90 2.1, triqs 3.x.x, and Triqs MaxEnt included.

## Building the docker image
The Dockerfile is built with this command, where `<version name>` could be `3.0.0`:
```
docker build -t triqs_mpich:<version name> -f mpich_dockerfile ./
docker build -t triqs_openmpi:<version name> -f openmpi_dockerfile ./
```
Note that you need a working, modified vasp version as archive (csc_vasp.tar.gz) in this directory to make the CSC calculation work.

## Pulling a docker image
Alternatively, you can pull an already-compiled image from the ETH gitlab container registry.
First [log in with a personal access token](https://gitlab.ethz.ch/help/user/packages/container_registry/index#authenticating-to-the-gitlab-container-registry).
This token you can save into a file and then log in into the registry with
```
cat <path to token> | docker login registry.ethz.ch -u <username gitlab> --password-stdin
```
and then run
```
docker pull registry.ethz.ch/d-matl-theory/uni-dmft/<image name>:<version name>
```
Just make sure that the version is the one that you want to have, it might not yet contain recent changes or bug fixes. Alternatively, there is the [official triqs docker image](https://hub.docker.com/r/flatironinstitute/triqs/), which however is not optimized for use on Daint.

## Getting docker images onto CSCS daint
First, you load the desired docker images with [sarus on daint](https://user.cscs.ch/tools/containers/sarus/).
Then there are two ways of getting the image on daint:

(1) For gitlab images (don't forget that you need the personal access token) or other, public image this can be done via:
```
sarus pull registry.ethz.ch/d-matl-theory/uni-dmft/<image name>:<version name>
sarus pull materialstheory/triqs
```
Pulling from the gitlab didn't work on daint when I tried, which leaves you with the second option.

(2) If you wish to use your locally saved docker image, you first have to save it
```
docker save --output=docker-triqs.tar <image name>:<version name>
```
and then upload the tar to daint and then load it via:
```
sarus load docker-triqs.tar <image name>:<version name>
```
then you can run it as shown in the example files.

### Running a docker container

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
Make sure you mount all the directories you need (where you save your data, where your uni-dmft directory is, ...)!
- All the flags are explained in 'docker run --help'.

Inside the docker, you can normally execute program. To run uni-dmft, for example, use
```
mpirun -n 4 python <path to uni-dmft>/run_dmft.py
```
To start a jupyter-lab server from the current path, use
```
jupyter.sh
```
All these commands you can execute directly by just adding them to the `docker run ... bash` command with the `-c` flag, e.g.
```
docker run --rm -it --shm-size=4g -e USER_ID=`id -u` -e GROUP_ID=`id -g` -p 8378:8378 -v $PWD:/work <image name>:<version name> bash -c 'cd /work && mpirun -n 4 python <path to uni-dmft>/run_dmft.py'
```
