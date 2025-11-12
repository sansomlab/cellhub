# CellHub Reimplementation

## Prepare environment

### Development notes

```bash
$ module load Python/3.11 R
$ module list

Currently Loaded Modules:
  1) GCCcore/12.3.0                     33) LLVM/16.0.6-GCCcore-12.3.0            65) UCC/1.2.0-GCCcore-12.3.0             97) PROJ/9.2.0-GCCcore-12.3.0
  2) binutils/2.40-GCCcore-12.3.0       34) Mesa/23.1.4-GCCcore-12.3.0            66) OpenMPI/4.1.5-GCC-12.3.0             98) libgeotiff/1.7.1-GCCcore-12.3.0
  3) bzip2/1.0.8-GCCcore-12.3.0         35) libGLU/9.0.3-GCCcore-12.3.0           67) gompi/2023a                          99) cffi/1.15.1-GCCcore-12.3.0
  4) zlib/1.2.13-GCCcore-12.3.0         36) pixman/0.42.2-GCCcore-12.3.0          68) Szip/2.1.1-GCCcore-12.3.0           100) cryptography/41.0.1-GCCcore-12.3.0
  5) libreadline/8.2-GCCcore-12.3.0     37) gettext/0.21.1-GCCcore-12.3.0         69) HDF5/1.14.0-gompi-2023a             101) virtualenv/20.23.1-GCCcore-12.3.0
  6) ncurses/6.4-GCCcore-12.3.0         38) libxml2/2.11.4-GCCcore-12.3.0         70) UDUNITS/2.2.28-GCCcore-12.3.0       102) Python-bundle-PyPI/2023.06-GCCcore-12.3.0
  7) Tcl/8.6.13-GCCcore-12.3.0          39) PCRE2/10.42-GCCcore-12.3.0            71) GSL/2.7-GCC-12.3.0                  103) pybind11/2.11.1-GCCcore-12.3.0
  8) SQLite/3.42.0-GCCcore-12.3.0       40) GLib/2.77.1-GCCcore-12.3.0            72) ATK/2.38.0-GCCcore-12.3.0           104) SciPy-bundle/2023.07-gfbf-2023a
  9) XZ/5.4.2-GCCcore-12.3.0            41) cairo/1.17.8-GCCcore-12.3.0           73) DBus/1.15.4-GCCcore-12.3.0          105) libtirpc/1.3.3-GCCcore-12.3.0
 10) libffi/3.4.4-GCCcore-12.3.0        42) NASM/2.16.01-GCCcore-12.3.0           74) at-spi2-core/2.49.91-GCCcore-12.3.0 106) HDF/4.2.16-2-GCCcore-12.3.0
 11) OpenSSL/1.1                        43) libjpeg-turbo/2.1.5.1-GCCcore-12.3.0  75) at-spi2-atk/2.38.0-GCCcore-12.3.0   107) Boost/1.82.0-GCC-12.3.0
 12) Python/3.11.3-GCCcore-12.3.0       44) jbigkit/2.1-GCCcore-12.3.0            76) Gdk-Pixbuf/2.42.10-GCCcore-12.3.0   108) arpack-ng/3.9.0-foss-2023a
 13) GCC/12.3.0                         45) libdeflate/1.18-GCCcore-12.3.0        77) HarfBuzz/5.3.1-GCCcore-12.3.0       109) Armadillo/12.6.2-foss-2023a
 14) OpenBLAS/0.3.23-GCC-12.3.0         46) LibTIFF/4.5.0-GCCcore-12.3.0          78) libepoxy/1.5.10-GCCcore-12.3.0      110) CFITSIO/4.3.0-GCCcore-12.3.0
 15) FlexiBLAS/3.3.1-GCC-12.3.0         47) Java/11 -> Java/11.0.20               79) Wayland/1.22.0-GCCcore-12.3.0       111) giflib/5.2.1-GCCcore-12.3.0
 16) FFTW/3.3.10-GCC-12.3.0             48) Tk/8.6.13-GCCcore-12.3.0              80) GTK3/3.24.37-GCCcore-12.3.0         112) json-c/0.16-GCCcore-12.3.0
 17) gfbf/2023a                         49) cURL/8.0.1-GCCcore-12.3.0             81) Ghostscript/10.01.2-GCCcore-12.3.0  113) Xerces-C++/3.2.4-GCCcore-12.3.0
 18) expat/2.5.0-GCCcore-12.3.0         50) GMP/6.2.1-GCCcore-12.3.0              82) JasPer/4.0.0-GCCcore-12.3.0         114) Imath/3.1.7-GCCcore-12.3.0
 19) libpng/1.6.39-GCCcore-12.3.0       51) NLopt/2.7.1-GCCcore-12.3.0            83) LittleCMS/2.15-GCCcore-12.3.0       115) OpenEXR/3.1.7-GCCcore-12.3.0
 20) Brotli/1.0.9-GCCcore-12.3.0        52) libogg/1.3.5-GCCcore-12.3.0           84) Pango/1.50.14-GCCcore-12.3.0        116) Highway/1.0.4-GCCcore-12.3.0
 21) util-linux/2.39-GCCcore-12.3.0     53) FLAC/1.4.2-GCCcore-12.3.0             85) FriBidi/1.0.12-GCCcore-12.3.0       117) Brunsli/0.1-GCCcore-12.3.0
 22) fontconfig/2.14.2-GCCcore-12.3.0   54) libvorbis/1.3.7-GCCcore-12.3.0        86) ImageMagick/7.1.1-15-GCCcore-12.3.0 118) Qhull/2020.2-GCCcore-12.3.0
 23) freetype/2.13.0-GCCcore-12.3.0     55) libopus/1.4-GCCcore-12.3.0            87) GLPK/5.0-GCCcore-12.3.0             119) LERC/4.0.0-GCCcore-12.3.0
 24) xorg-macros/1.20.0-GCCcore-12.3.0  56) LAME/3.100-GCCcore-12.3.0             88) nodejs/18.17.1-GCCcore-12.3.0       120) OpenJPEG/2.5.0-GCCcore-12.3.0
 25) libpciaccess/0.17-GCCcore-12.3.0   57) libsndfile/1.2.2-GCCcore-12.3.0       89) FFTW.MPI/3.3.10-gompi-2023a         121) SWIG/4.1.1-GCCcore-12.3.0
 26) X11/20230603-GCCcore-12.3.0        58) ICU/73.2-GCCcore-12.3.0               90) ScaLAPACK/2.2.0-gompi-2023a-fb      122) GDAL/3.7.1-foss-2023a
 27) gzip/1.12-GCCcore-12.3.0           59) numactl/2.0.16-GCCcore-12.3.0         91) foss/2023a                          123) MPFR/4.2.0-GCCcore-12.3.0
 28) lz4/1.9.4-GCCcore-12.3.0           60) hwloc/2.9.1-GCCcore-12.3.0            92) netCDF/4.9.2-gompi-2023a            124) libgit2/1.7.1-GCCcore-12.3.0
 29) zstd/1.5.5-GCCcore-12.3.0          61) libevent/2.1.12-GCCcore-12.3.0        93) GEOS/3.12.0-GCC-12.3.0              125) R/4.5.1-gfbf-2023a-bare-noSciPy
 30) libdrm/2.4.115-GCCcore-12.3.0      62) UCX/1.14.1-GCCcore-12.3.0             94) libarchive/3.6.2-GCCcore-12.3.0
 31) libglvnd/1.6.0-GCCcore-12.3.0      63) libfabric/1.18.0-GCCcore-12.3.0       95) PCRE/8.45-GCCcore-12.3.0
 32) libunwind/1.6.2-GCCcore-12.3.0     64) PMIx/4.2.4-GCCcore-12.3.0             96) nlohmann_json/3.11.2-GCCcore-12.3.0
```

```bash
python -m venv cellhub-new-py311
source cellhub-new-py311/bin/activate
python -m pip install --upgrade pip
```

```bash
cd $CELLHUB_DIR
python setup.py develop
pip install snakemake snakemake-executor-plugin-slurm
pip install -r python/requirements.txt
Rscript R/install.packages.R # may need firstly install devtools, BiocManager, Cairo, stringi
R CMD INSTALL R/cellhub
```

```bash
$ snakemake --version
9.13.6
```

### Usage

To generate a `YAML` template:

```bash
snakemake -s ${PATH_TO_SNAKEFILE} --config target=${MODULE_TO_RUN} mode=genconfig --cores 1
```

To make real execution:

```bash
snakemake -s ${PATH_TO_SNAKEFILE} --config target=${MODULE_TO_RUN} --cores ${NCORES} -p --jobs ${NJOBS} \
    --executor ${EXECUTOR} --default-resources mem_mb=16000
    # ${EXECUTOR} can be `slurm`
```
