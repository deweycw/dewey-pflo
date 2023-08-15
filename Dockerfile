FROM debian:latest

RUN apt update && apt upgrade -y
RUN apt install wget -y
RUN apt install python3 -y
RUN apt install gcc -y
RUN apt-get install gfortran -y
RUN apt-get install make -y
RUN apt-get install cmake -y
RUN apt-get install libhdf5-serial-dev -y 
RUN apt-get install metis -y
RUN apt-get install openmpi-bin openmpi-common libopenmpi-dev libgtk2.0-dev -y
RUN apt install git -y

RUN git clone https://gitlab.com/petsc/petsc.git petsc
WORKDIR /petsc
RUN git checkout v3.17.1
RUN ./configure --CFLAGS='-O3' --CXXFLAGS='-O3' --FFLAGS='-O3' --with-debugging=no --download-mpich=yes --download-hdf5=yes --download-hdf5-fortran-bindings=yes  --download-metis=yes --download-parmetis=yes --download-fblaslapack=1

RUN make PETSC_DIR=/petsc PETSC_ARCH=arch-linux-c-opt all

ENV PETSC_DIR='/petsc'
ENV PETSC_ARCH='arch-linux-c-opt'

WORKDIR /pflotran
RUN git clone https://github.com/deweycw/dewey-pflo.git

WORKDIR /pflotran/dewey-pflo/src/pflotran
RUN git checkout jin-bethke
RUN make pflotran
ENV PFLOTRAN_DIR='/pflotran/dewey-pflo/src/pflotran'

ARG USERNAME=pflotranUser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    #
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# ********************************************************
# * Anything else you want to do like clean up goes here *
# ********************************************************

# [Optional] Set the default user. Omit if you want to keep the default as root.
USER $USERNAME
WORKDIR /