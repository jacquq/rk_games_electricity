#ifndef UTILS_PATHS
#define UTILS_PATHS

/* 
========================================================================

Outils pour récupérer le path de l'exe, voir .hpp

========================================================================
*/

#if defined(_WIN32)
    #define NOMINMAX
    #include <windows.h>
    #include <Shlwapi.h>
    #include <io.h> 

    #define access _access_s
#endif

#ifdef __APPLE__
    #include <libgen.h>
    #include <limits.h>
    #include <mach-o/dyld.h>
    #include <unistd.h>
#endif

#ifdef __linux__
    #include <limits.h>
    #include <libgen.h>
    #include <unistd.h>

    #if defined(__sun)
        #define PROC_SELF_EXE "/proc/self/path/a.out"
    #else
        #define PROC_SELF_EXE "/proc/self/exe"
    #endif

#endif

#include <string>
#include <cstring>
#include <iostream>
using namespace std;

#if defined(_WIN32)
string getExecutablePath() {
   char rawPathName[MAX_PATH];
   GetModuleFileNameA(NULL, rawPathName, MAX_PATH);
   return string(rawPathName);
};
#endif

#ifdef __linux__
string getExecutablePath() {
   char rawPathName[PATH_MAX];
   realpath(PROC_SELF_EXE, rawPathName);
   return  string(rawPathName);
};
#endif

#ifdef __APPLE__
string getExecutablePath() {
    char rawPathName[PATH_MAX];
    char realPathName[PATH_MAX];
    uint32_t rawPathSize = (uint32_t)sizeof(rawPathName);

    if(!_NSGetExecutablePath(rawPathName, &rawPathSize)) {
        realpath(rawPathName, realPathName);
    }
    return  string(realPathName);
};
#endif

string getExecutableDir() {
    string path_exe = getExecutablePath();
    return path_exe.substr(0,path_exe.find_last_of("/\\"));
};

#endif