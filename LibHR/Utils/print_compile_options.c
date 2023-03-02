#include "Utils/print_compile_options.h"
#include "IO/logger.h"
#include <stdio.h>

void print_compiling_info() {
    printf("MACROS\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s", MACROS);
    printf("MkFlags\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s", CI_mkflags);
    printf("\n\n");
    printf("git branch\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s", CI_gitinfo);
    printf("\n\n");
    printf("git revision\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s", CI_gitrevision);
    printf("\n\n");
    printf("/proc/cpuinfo\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s", CI_cpuinfo);
    printf("\n\n");
    printf("/proc/version\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s", CI_linux);
    printf("\n\n");
    printf("gcc -v\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s", CI_gcc);
}

void print_compiling_info_short() {
    lprintf("SYSTEM", 0, "MACROS=%s", MACROS);
}
