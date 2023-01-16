@page build_system Build System
[TOC]
# Compile the code
For a more efficient compilation we are using ninja build for compilation.
Find detailed compilation instructions in the repository's README and user guide.

# Remcompile after significant changes
When recompiling with changed compilation flags, run

```bash
$ nj -t cleandead
$ nj -t clean
$ nj
```

in the target directory to clean the code of artefacts, that are only compiled with different compilation flags.

When recompiling after making adjustments to the file structure, such as renaming, adding or removing a file, remove the build folder in the root.

```bash
$ nj -t cleandead
$ nj -t clean
$ rm -r build
$ touch Make/MkFlags.ini
$ nj 
```

Changing MkFlags.ini is necessary, because this tells ninja to regenerate the build folder from the new directory structure.

# Cross-compilation setup
The CUDA library is a C++ library, while `HiRep` is written in C. This means that when compiling for GPU, we need to cross-compile. Make sure that your code is setup correctly, by following these instructions

* C source files always have the ending `.c`, C++ CUDA source files always have the `.cu`.
* C header files always have the ending `.h`, C++ CUDA header files always have the ending `.hpp`.
* A file needs to be a CUDA file only if it calls or contains a kernel. Pure references to CUDA objects should be done in pure C files. CUDA kernels can be declared either directly in the `.cu`-file or in a header file `.hpp` in the same directory that is included in only one corresponding `.cu` file with the line
```c
#include "./name_of_header.hpp"
```
There should never be any kernel header files in `Include`, because kernels should not be visible globally. What should be visible are higher level functions that do or do not call kernels internally and that in turn can be called equally from C or C++ code.
* Functions that are declared either in a C or C++ source file can be made visible by including it in a header in include. In order to make C++ functions visible to C code, they need to be declared with extern C linkage. However, declaring this linkage is only necessary if the header file is included in a C++ file. In order to make sure, that all functions are setup correctly, header files are generally wrapped in the following way
```c
#include "example1.h"
#include "example2.h"

#ifdef __cplusplus
    extern "C" {
#endif 

// Function declarations

#ifdef __cplusplus
    }
#endif
```


