# PFLOTRAN installation instructions on MacOS (arm64)
### Christian Dewey, written and tested 2023-02-08

### REQUIREMENTS: xcode command-line utils

### **Step 0**: If not already done, install Xcode command line tools and gcc, and make gcc aliases

At a new terminal window, type the following:

        xcode-select --install

This will prompt a window from MacOS asking if you want to install Xcode, and whether you agree to the terms. Click Yes/Agree, and proceed with installation. 

Once installation completes, restart computer. 

After restart, open a new terminal window and install gcc with homebrew:

        brew install gcc

Determine version of gcc installed:

        brew info gcc

This will display somehting like:
                                                                    
        ==> gcc: stable 12.2.0 (bottled), HEAD
        GNU compiler collection
        https://gcc.gnu.org/
        /opt/homebrew/Cellar/gcc/12.2.0 (1,470 files, 358.8MB) *
        Poured from bottle on 2023-02-07 at 17:42:04
        From: https://github.com/Homebrew/homebrew-core/blob/HEAD/Formula/gcc.rb
        License: GPL-3.0-or-later with GCC-exception-3.1
        ==> Dependencies
        Required: gmp ✔, isl ✔, libmpc ✔, mpfr ✔, zstd ✔
        ==> Options
        --HEAD
            Install HEAD version
        ==> Analytics
        install: 78,057 (30 days), 227,784 (90 days), 1,395,629 (365 days)
        install-on-request: 38,183 (30 days), 112,782 (90 days), 666,694 (365 days)
        build-error: 252 (30 days)

gcc version in shown in line 0; the location of the install is shown in line 3

Now create aliases for gcc. MacOS defaults to clang. The aliases force gcc to be used. 

Open ~/.zshrc:

        vi ~/.zshrc

Add the following lines to the file, where [VERSION] corresponds to the first number (e.g., '12') of the gcc version:

        alias gcc=gcc-[VERSION]
        alias g++=g++-[VERSION]

Write and save the changes in the ~/.zshrc, and then resource the file, or quit the terminal and reopen it. 


### **Step 1**: Install required packages with hombrew

        brew install open-mpi
        brew install hdf5-mpi
        brew install lapack
        brew install cmake

Determine installed versions and locations of install, as above with gcc:

        brew info open-mpi
        brew info hdf5-mpi
  

### **Step 2**: Clone PETSc

Move to the location where you wish to save PETSc. Then clone with:

        git clone https://gitlab.com/petsc/petsc.git petsc
        cd petsc
        git checkout v3.17.1


### **Step 3**: Configure and make PETSc

Enter the topmost PETSc dir and run configure:

        cd petsc

        ./configure --with-debugging=no --with-shared-libraries=0 --known-mpi-shared-libraries=1 --with-mpi-dir=/opt/homebrew/Cellar/open-mpi/4.1.4_2/ --known-64-bit-blas-indices --with-hdf5-dir=/opt/homebrew/Cellar/hdf5-mpi/1.12.2_1 / --download-hdf5-fortran-bindings=yes

Note that if the versions of hdf5-mpi and/or open-mpi differ from hdf5-mpi/1.12.2_1 and open-mpi/4.1.4_2, replace --with-mpi-dir= and --with-hdf5-dir= paths with correct values

After the configuration is complete, it will prompt you to make PETSc and provide a command for doing so. If PETSc compiles, it will prompt you test the compilation, and provide a command to do so. Note that the test may partially fail -- this is not necessarily a problem. 

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