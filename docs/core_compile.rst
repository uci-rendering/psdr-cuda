Compiling the system
=========================

Compiling PSDR-CUDA from scratch requires recent versions of CMake (at least 3.8.8), Python (at least 3.6), CUDA (at least 10.1), OptiX (at least 7.0), as well as a C++ compiler that supports the C++17 standard.


Linux
--------------------

The build process under Linux requires several external dependencies that are easily installed using the system-provided package manager (e.g., ``apt-get`` under Ubuntu).

To fetch all dependencies, enter the following commands on Ubuntu:

.. code-block:: bash

   # Install OpenEXR for image I/O
   sudo apt install libopenexr-dev

   # Install required Python packages
   sudo apt install python3-dev python3-distutils python3-setuptools

Additionally, the Enoki and pybind11 libraries are needed. To make pybind11 globally accessible, we recommend running ``sudo make install`` to copy all its header files to ``/usr/local/include/pybind11/``.

Lastly, a softlink to Nvidia OptiX needs to be created at ``ext/optix``. For example, if OptiX is installed under ``/home/user/NVIDIA-OptiX-SDK-7.1.0-linux64``, this can be done by running the following commands from inside PSDR-CUDA's root directory:

.. code-block:: bash

   mkdir -p ext
   ln -s /home/user/NVIDIA-OptiX-SDK-7.1.0-linux64 ext/optix

With all dependencies installed/compiled, building PSDR-CUDA requires defining a few Cmake variables including ``PYTHON_INCLUDE_PATH``, ``OpenEXR_ROOT``, and ``ENOKI_DIR``. Assuming the Enoki library is cloned to ``/home/user/enoki`` and built under ``/home/user/enoki/build``, this can be done by running the following commands from inside PSDR-CUDA's root directory:

.. code-block:: bash

   # Create a directory where build products are stored
   mkdir build
   cd build
   cmake -D CMAKE_C_COMPILER=gcc-9                       \
         -D CMAKE_CXX_COMPILER=g++-9                     \
         -D PYTHON_INCLUDE_PATH=/usr/include/python3.8/  \
         -D OpenEXR_ROOT=/usr/include/OpenEXR/           \
         -D ENOKI_DIR=/home/user/enoki/                  \
         ..
   make -j

Tested version
^^^^^^^^^^^^^^
* Ubuntu 20.04
* gcc 9.3.0


Windows
--------------------

To make the compilation under Windows more streamlined, we created a `package <https://www.dropbox.com/s/98rh15mswa5lxe5/ext_win64.7z?dl=0>`_ containing precompiled binaries of most of the dependencies. Utilizing this package requires having Visual Studio 2019, CUDA 11, and Python 3.8.

After decompressing this package to ``ext_win64`` under PSDR-CUDA's root directory, assuming Python to be installed under ``C:/Users/User/Anaconda3``, PSDR-CUDA can be built by running the following commands under the command prompt:

.. code-block:: batch

   rem Create a directory where build products are stored
   mkdir build
   cd build
   cmake -D PYTHON_ROOT="C:/Users/User/Anaconda3" -G "Visual Studio 16 2019" -A x64 ..
   cmake --build . --config Release
   copy /y lib\Release\*.pyd lib\

Tested version
^^^^^^^^^^^^^^
* Windows 10
* Visual Studio 2019 (Community Edition)
* Python 3.8.3 (Anaconda)
* CUDA 11
* OptiX 7.1


Mac OS
--------------------

Unfortunately, PSDR-CUDA does not work under the Mac OS due to the lack of CUDA support.
