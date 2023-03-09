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

## Usage

SurVClusterer requires two files:

```
SurVClusterer filelist reference
```

where reference is the reference genome used to call the SVs, in FASTA format, and filelist is a list of the files to be clustered, in the format
```
SAMPLENAME1 /path/to/vcf/1
SAMPLENAME2 /path/to/vcf/2
```

For a full list of the parameters, run
```
SurVClusterer -h
```
