# SurVClusterer

SurVClusterer reads SVs from multiple VCF files, each representing a different sample, and clusters SVs that are compatible, according to user-defined parameters.

Please note that it is still in active development.

## Installation

In order to compile the code, the following are required:
- A C and a C++ compiler are required. If the GCC suite is used, version 4.9.3 or above is required.
- A recent version of CMake (3.5 or above)

The following commands should be sufficient:

```
git clone https://github.com/Mesh89/SurVClusterer
cd SurVClusterer/
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```
