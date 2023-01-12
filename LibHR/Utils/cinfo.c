char MACROS[] = "-DBC_T_PERIODIC -DBC_X_PERIODIC -DBC_Y_PERIODIC -DBC_Z_PERIODIC -DUPDATE_EO -DNDEBUG -DCHECK_SPINOR_MATCHING -DIO_FLUSH -DNG=2 -DGAUGE_SUN -DREPR_ADJOINT -DREPR_NAME=\"REPR_ADJOINT\"\n";
char CI_mkflags[] = "NG = 2\nREPR = REPR_ADJOINT\nGAUGE_GROUP = GAUGE_SUN\nMACRO += BC_T_PERIODIC\nMACRO += BC_X_PERIODIC\nMACRO += BC_Y_PERIODIC\nMACRO += BC_Z_PERIODIC\nMACRO += UPDATE_EO\nMACRO += NDEBUG\nMACRO += CHECK_SPINOR_MATCHING\nMACRO += IO_FLUSH\nCC = gcc\nMPICC = mpicc\nCFLAGS = -Wall -O1\nNVCC = nvcc\nGPUFLAGS = \nINCLUDE = \nLDFLAGS = \n\n";
char CI_cpuinfo[] = "No CPU info\n";
char CI_linux[] = "Darwin Claudios-MacBook-Pro.local 21.6.0 Darwin Kernel Version 21.6.0: Sat Jun 18 17:07:25 PDT 2022; root:xnu-8020.140.41~1/RELEASE_X86_64 x86_64\n";
char CI_gcc[] = "Apple clang version 14.0.0 (clang-1400.0.29.202)\nTarget: x86_64-apple-darwin21.6.0\nThread model: posix\nInstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin\n\n";
char CI_gitinfo[] = "";
char CI_gitrevision[] = "9d0038a8f3d3bcab57a424a20a233a5406858326";
