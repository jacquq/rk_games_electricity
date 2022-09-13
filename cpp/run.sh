#!/bin/sh

RED='\033[0;31m'
GREEN='\033[0;32m'
GRAY='\033[1;30m'
NC='\033[0m' # No Color

case "$(uname -s)" in

   Darwin|Linux)
    echo "${NC}[  0%]${GREEN} Construction 'build'${NC}"
    rm -rf "build"
    mkdir "build"
    cd "build"

    echo "${NC}[  0%]${GREEN} Compilation ...${GRAY}"
    export PATH=/opt/cmake/3.14.4/bin:$PATH
    cmake ..
    cmake --build . --config Release

    echo "${NC}[100%]${GREEN} Copy executable${NC}"
    cd ..
     cp "build/rk_games" "./rk_games";;
   CYGWIN*|MINGW32*|MSYS*|MINGW*)
    echo "[  0%] Construction 'build'"
    rm -rf "build"
    mkdir "build"
    cd "build"

    echo "[  0%] Compilation ..."
    cmake -DCMAKE_GENERATOR_PLATFORM=x64 .. -DBoost_NO_WARN_NEW_VERSIONS=1
    cmake --build . --config Release

    echo "[100%] Copy executable"
    cd ..
     cp "build/Release/rk_games.exe" "./rk_games.exe";;
   *)
     echo 'Other OS' 
     ;;
esac