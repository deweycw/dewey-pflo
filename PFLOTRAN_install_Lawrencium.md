# PFLOTRAN installation instructions on LBL Lawrencium cluster
### Christian Dewey, written and tested 2023-02-08


### **Step 1**: Load shared modules

        module load gcc/4.8.5
        module load lapack
        module load openmpi
        module load hdf5/1.8.18-gcc-p

Check that correct modules are loaded:

        module list

The gcc and hdf5 versions should appear as above. Lapack and openmpi should be versions 3.8.0-gcc and 2.0.2-gcc, respectively. 


### **Step 2**: Clone PETSc

        git clone https://gitlab.com/petsc/petsc.git petsc
        cd petsc
        git checkout v3.17.1


### **Step 3**: Configure and make PETSc

This line must be run from the petsc dir.

        ./configure --with-debugging=no --with-batch --with-shared-libraries=yes --known-mpi-shared-libraries=1 --with-mpi-dir=/global/software/sl-7.x86_64/modules/gcc/4.8.5/openmpi/2.0.2-gcc/ --known-64-bit-blas-indices --with-hdf5-dir=/global/software/sl-7.x86_64/modules/gcc/4.8.5/hdf5/1.8.18-gcc-p/ --download-hdf5-fortran-bindings=yes

After the configuration is complete, it will prompt you to make PETSc and provide a command for doing so. If PETSc compiles, it will prompt you test the compilation, and provide a command to do so. 

Check that PETSC_DIR and PETSC_ARCH are defined:

        echo $PETSC_DIR
        echo $PETSC_ARCH


### **Step 4**: Clone PFLOTRAN 

From $HOME, make a new directory for PFLOTRAN. Enter the directory and clone PFLOTRAN 

        git clone https://gitlab.com/pflotran/pflotran.git


### **Step 5**: Make PFLOTRAN 

        cd pflotran/src/pflotran
        git checkout maint/v4.0
        make pflotran 
        
This will take ~10 min.