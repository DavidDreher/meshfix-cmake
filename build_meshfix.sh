#!/bin/bash
set -euo pipefail

prev_path="$(pwd)";
cd "$(dirname "$(realpath "$0")")";
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
make -Cbuild
cd "$prev_path"
