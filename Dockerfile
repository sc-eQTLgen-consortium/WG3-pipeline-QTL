
################## BASE IMAGE ######################

FROM ubuntu:22.04

################## METADATA ######################

LABEL base_image="ubuntu:22.04"
LABEL version="1.0.0"
LABEL software="WG3 Pipeline"
LABEL about.summary="WG3 sceQTLGen Consortium QTL Pipeline"
LABEL about.documentation="https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL.git"
LABEL about.tags="Genomics"

# Build syntax: docker build ./ -t wg3-pipeline-qtl:2024.02.05.0 --progress=plain > build.log 2>&1
# Total build takes ? minutes and has a size of 1.49 GB.
# Use dive wg3-pipeline-qtl:2024.02.05.0 to investigate memory usage.

################## MAINTAINER ######################

MAINTAINER Marc Jan Bonder <m.j.bonder@umcg.nl>, Martijn Vochteloo <m.vochteloo@umcg.nl>

################## INSTALLATION ######################

ADD . /tmp/repo
WORKDIR /tmp/repo

ENV PATH=/opt:/usr/games:/opt/conda/envs/py311/bin:/opt/conda/bin:/opt/GenotypeHarmonizer-1.4.27:/opt/ldstore_v2.0_x86_64:/opt/plink2:$PATH
ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN echo 'alias python=python3' >> ~/.bashrc

# Needed to prevent asking for geographic location when installing things.
RUN export TZ=Europe/Amsterdam \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# libssl3
RUN apt-get update -y \
    # libc-bin libc6 libsystemd0 libudev1
    && apt-get upgrade -y \
    # binutils binutils-common binutils-x86-64-linux-gnu build-essential bzip2 cpp
    # cpp-11 dpkg-dev g++ g++-11 gcc gcc-11 gcc-11-base libasan6 libatomic1
    # libbinutils libc-dev-bin libc6-dev libcc1-0 libcrypt-dev libctf-nobfd0
    # libctf0 libdpkg-perl libgcc-11-dev libgdbm-compat4 libgdbm6 libgomp1
    # libisl23 libitm1 liblsan0 libmpc3 libmpfr6 libnsl-dev libperl5.34
    # libquadmath0 libstdc++-11-dev libtirpc-dev libtsan0 libubsan1 linux-libc-dev
    # lto-disabled-list make patch perl perl-modules-5.34 rpcsvc-proto xz-utils
    && apt-get install -y --no-install-recommends build-essential \
    # ca-certificates openssl
    && apt-get install -y --no-install-recommends ca-certificates \
    # libpsl5 wget
    && apt-get install -y --no-install-recommends wget \
    # dirmngr gnupg gnupg-l10n gnupg-utils gnupg2 gpg gpg-agent gpg-wks-client
    # gpg-wks-server gpgconf gpgsm libassuan0 libksba8 libldap-2.5-0 libnpth0
    # libreadline8 libsasl2-2 libsasl2-modules-db libsqlite3-0 pinentry-curses
    # readline-common
    && apt-get install -y --no-install-recommends gnupg2 \
    # libdeflate0 libhtscodecs2 tabix
    && apt-get install -y --no-install-recommends tabix

##################################
############# PYTHON #############
##################################

# Reduce conda size by preventing Python from recreating a corresponding bytecode cache file (*.pyc) at runtime.
ENV PYTHONDONTWRITEBYTECODE=true

# Install Python.
# libexpat1 libmpdec3 libpython3-stdlib libpython3.10-minimal
# libpython3.10-stdlib media-types python3 python3-minimal python3.10
# python3.10-minimal
RUN apt-get install -y --no-install-recommends python3

# Install miniconda for the virtual environment.
# https://github.com/ContinuumIO/docker-images/blob/main/miniconda3/debian/Dockerfile
RUN cd /opt \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -O miniconda.sh -q \
    && mkdir -p /opt \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy

# Create and activate virtual environment
RUN eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" \
    && conda create -n py311 python=3.11.5 \
    # None
    && /opt/conda/envs/py311/bin/pip install numpy==1.26.0 \
    # python_dateutil-2.8.2 pytz-2023.3.post1 six-1.16.0 tzdata-2023.3
    && /opt/conda/envs/py311/bin/pip install pandas==2.1.1 \
    # h5py-3.10.0-cp311-cp311
    && /opt/conda/envs/py311/bin/pip install h5py==3.10.0 \
    # gdown-4.6.0 filelock-3.13.1 requests-2.31.0 tqdm-4.66.1 beautifulsoup4-4.12.3
    # soupsieve-2.5 charset_normalizer-3.3.2 idna-3.6 urllib3-2.2.0 certifi-2024.2.2
    # PySocks-1.7.1
    && /opt/conda/envs/py311/bin/pip install gdown==4.6.0 # v4.7.1 gives AttributeError: 'NoneType' object has no attribute 'groups'

# Creating a conda environment for limix and ld tools
RUN eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" \
    && conda create -n py38 python=3.8 \
    && conda install -n py38 blas=*=*mkl \
    && echo "blas=*=*mkl" > /opt/conda/envs/py38/conda-meta/pinned \
    # None
    && conda install -n py38 -c conda-forge bgen=4 \
    # None
    && /opt/conda/envs/py38/bin/pip install numpy==1.21 \
    # exceptiongroup-1.2.0 iniconfig-2.0.0 packaging-23.2 pluggy-1.4.0 tomli-2.0.1
    && /opt/conda/envs/py38/bin/pip install pytest==7.4.3 \
    # blosc2-2.0.0 cython-3.0.8 msgpack-1.0.7 numexpr-2.8.6 py-cpuinfo-9.0.0
    && /opt/conda/envs/py38/bin/pip install tables==3.8.0 \
    # joblib-1.3.2 scipy-1.10.1 threadpoolctl-3.2.0
    && /opt/conda/envs/py38/bin/pip install scikit-learn==1.3.2 \
    # contourpy-1.1.1 cycler-0.12.1 fonttools-4.48.1 importlib-resources-6.1.1 kiwisolver-1.4.5
    # matplotlib-3.7.4 matplotlib-venn-0.11.9 pillow-10.2.0 pyparsing-3.1.1 python-dateutil-2.8.2
    # six-1.16.0 zipp-3.17.0
    && /opt/conda/envs/py38/bin/pip install matplotlib-venn==0.11.9 \
    # Deprecated-1.2.14 cffi-1.16.0 click-8.1.7 cloudpickle-3.0.0 dask-2023.5.0 fsspec-2024.2.0
    # importlib-metadata-7.0.1 locket-1.0.0 pandas-2.0.3 partd-1.4.1 pycparser-2.21
    # pytz-2024.1 pyyaml-6.0.1 toolz-0.12.1 tqdm-4.66.1 tzdata-2023.4 wrapt-1.16.0 xarray-2023.1.0
    # zstandard-0.22.0
    && /opt/conda/envs/py38/bin/pip install pandas-plink==2.2.9 \
    # None
    && /opt/conda/envs/py38/bin/pip install h5py==3.10.0 \
    # appdirs-1.4.4 cachetools-5.3.2 cbgen-1.0.4 certifi-2024.2.2 charset-normalizer-3.3.2
    # idna-3.6 platformdirs-4.2.0 pooch-1.8.0 requests-2.31.0 texttable-1.7.0 urllib3-2.2.0
    && /opt/conda/envs/py38/bin/pip install bgen-reader==4.0.8 \
    # attrs-23.2.0 brent-search-2.0.1 liknorm-1.2.10 ndarray-listener-2.0.1 numpy-sugar-1.5.4 optimix-3.0.4
    # py-1.11.0 pytest-6.2.5 pytest-doctestplus-1.1.0 toml-0.10.2
    && /opt/conda/envs/py38/bin/pip install glimix-core==3.1.13 \
    # zstd-1.5.5.1
    && /opt/conda/envs/py38/bin/pip install https://files.pythonhosted.org/packages/a8/fd/f98ab7dea176f42cb61b80450b795ef19b329e8eb715b87b0d13c2a0854d/ldstore-0.1.9.tar.gz

RUN conda clean -y --all

##############################
############## R #############
##############################

# install the helper packages we need: cmake, dirmngr, and software-properties-common
# cmake cmake-data dh-elpa-helper emacsen-common libarchive13 libbrotli1
# libcurl4 libicu70 libjsoncpp25 libnghttp2-14 librhash0 librtmp1 libssh-4
# libuv1 libxml2
RUN apt-get install -y --no-install-recommends cmake \
    # dbus distro-info-data gir1.2-glib-2.0 gir1.2-packagekitglib-1.0 gpg gpgconf
    # iso-codes libapparmor1 libappstream4 libargon2-1 libassuan0 libcap2-bin
    # libcryptsetup12 libcurl3-gnutls libdbus-1-3 libdevmapper1.02.1 libdw1
    # libelf1 libgirepository-1.0-1 libglib2.0-0 libglib2.0-bin libglib2.0-data
    # libgstreamer1.0-0 libip4tc2 libjson-c5 libkmod2 libmpdec3
    # libpackagekit-glib2-18 libpam-systemd libpolkit-agent-1-0
    # libpolkit-gobject-1-0 libpython3-stdlib libpython3.10-minimal
    # libpython3.10-stdlib libreadline8 libsqlite3-0 libstemmer0d libunwind8
    # libxmlb2 libyaml-0-2 lsb-release media-types packagekit pkexec policykit-1
    # polkitd python-apt-common python3 python3-apt python3-blinker
    # python3-cffi-backend python3-cryptography python3-dbus python3-distro
    # python3-gi python3-httplib2 python3-importlib-metadata python3-jeepney
    # python3-jwt python3-keyring python3-launchpadlib python3-lazr.restfulclient
    # python3-lazr.uri python3-minimal python3-more-itertools python3-oauthlib
    # python3-pkg-resources python3-pyparsing python3-secretstorage python3-six
    # python3-software-properties python3-wadllib python3-zipp python3.10
    # python3.10-minimal readline-common software-properties-common systemd
    # systemd-sysv
    && apt-get install -y --no-install-recommends software-properties-common \
    # add the signing key (by Michael Rutter) for these repos
    # To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    # Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
    && wget -qO- "https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc" \
    && tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    # add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    # fontconfig fontconfig-config fonts-dejavu-core libblas3 libbsd0 libcairo2
    # libdatrie1 libfontconfig1 libfreetype6 libfribidi0 libgfortran5
    # libgraphite2-3 libharfbuzz0b libice6 libjbig0 libjpeg-turbo8 libjpeg8
    # liblapack3 libmd0 libpango-1.0-0 libpangocairo-1.0-0 libpangoft2-1.0-0
    # libpaper-utils libpaper1 libpixman-1-0 libpng16-16 libsm6 libtcl8.6
    # libthai-data libthai0 libtiff5 libtk8.6 libwebp7 libx11-6 libx11-data
    # libxau6 libxcb-render0 libxcb-shm0 libxcb1 libxdmcp6 libxext6 libxft2
    # libxrender1 libxss1 libxt6 r-base-core r-cran-boot r-cran-class
    # r-cran-cluster r-cran-codetools r-cran-foreign r-cran-kernsmooth
    # r-cran-lattice r-cran-mass r-cran-matrix r-cran-mgcv r-cran-nlme r-cran-nnet
    # r-cran-rpart r-cran-spatial r-cran-survival r-recommended tzdata ucf unzip
    # x11-common xdg-utils zip
    && apt-get install -y --no-install-recommends r-base \
    # gfortran gfortran-11 icu-devtools libblas-dev libbz2-dev libgfortran-11-dev
    # libicu-dev libjpeg-dev libjpeg-turbo8-dev libjpeg8-dev liblapack-dev
    # liblzma-dev libncurses-dev libncurses5-dev libpcre16-3 libpcre2-16-0
    # libpcre2-32-0 libpcre2-dev libpcre2-posix3 libpcre3-dev libpcre32-3
    # libpcrecpp0v5 libpng-dev libreadline-dev libxmuu1 pkg-config r-base-dev
    # xauth zlib1g-dev
    && apt-get install -y --no-install-recommends r-base-dev \
    # libbrotli-dev libexpat1-dev libfontconfig-dev libfontconfig1-dev
    # libfreetype-dev libfreetype6-dev uuid-dev
    && apt-get install -y --no-install-recommends libfontconfig1-dev \
    # libcurl4-openssl-dev
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev \
    # libxml2-dev
    && apt-get install -y --no-install-recommends libxml2-dev \
    # libssl-dev
    && apt-get install -y --no-install-recommends libssl-dev \
    # gir1.2-harfbuzz-0.0 libblkid-dev libffi-dev libglib2.0-dev
    # libglib2.0-dev-bin libgraphite2-dev libharfbuzz-dev libharfbuzz-gobject0
    && apt-get install -y --no-install-recommends libharfbuzz-dev \
    # libfribidi-dev
    && apt-get install -y --no-install-recommends libfribidi-dev \
    # libdeflate-dev libjbig-dev libtiff-dev libtiff5-dev libtiffxx5
    && apt-get install -y --no-install-recommends libtiff5-dev \
    # Required for hdf5r \
    # hdf5-helpers libaec-dev libaec0 libhdf5-103-1 libhdf5-cpp-103-1 libhdf5-dev
    # libhdf5-fortran-102 libhdf5-hl-100 libhdf5-hl-cpp-100 libhdf5-hl-fortran-100
    # libsz2
    && apt-get install -y --no-install-recommends libhdf5-dev \
    # libgsl-dev libgsl27 libgslcblas0
    && apt-get install -y --no-install-recommends libgsl-dev \
    # libcairo-gobject2 libcairo-script-interpreter2 libcairo2-dev libice-dev
    # liblzo2-2 libpixman-1-dev libpthread-stubs0-dev libsm-dev libx11-dev
    # libxau-dev libxcb-render0-dev libxcb-shm0-dev libxcb1-dev libxdmcp-dev
    # libxext-dev libxrender-dev x11proto-dev xorg-sgml-doctools xtrans-dev
    && apt-get install -y --no-install-recommends libcairo2-dev \
    # libxt-dev
    && apt-get install -y --no-install-recommends libxt-dev \
    # libgmp-dev libgmpxx4ldbl
    && apt-get install -y --no-install-recommends libgmp-dev \
    # libmpfr-dev
    && apt-get install -y --no-install-recommends libmpfr-dev \
    # libamd2 libbtf1 libcamd2 libccolamd2 libcholmod3 libcolamd2
    # libcxsparse3 libglpk-dev libglpk40 libklu1 libldl2 libltdl7
    # libmetis5 libmongoose2 librbio2 libsliplu1 libspqr2
    # libsuitesparse-dev libsuitesparseconfig5 libumfpack5
    && apt-get install -y --no-install-recommends libglpk-dev

# remotes_2.4.2.1
RUN R --slave -e 'install.packages("remotes")' \
    # getopt_1.20.4
    && R --slave -e 'remotes::install_version("optparse", version = "1.7.3", upgrade=FALSE)' \
    # R.methodsS3_1.8.2 R.oo_1.26.0
    && R --slave -e 'remotes::install_version("R.utils", version = "2.12.2", upgrade=FALSE)' \
    # None
    && R --slave -e 'remotes::install_version("data.table", version = "1.14.8", upgrade=FALSE)' \
    # Matrix_1.6-5 rematch2_2.1.2 diffobj_0.3.5 rprojroot_2.0.4
    # pkgbuild_1.4.3 fs_1.6.3 crayon_1.5.2 cpp11_0.4.7 pkgconfig_2.0.3
    # withr_3.0.0 waldo_0.5.2 ps_1.7.6 processx_3.8.3 praise_1.0.0
    # pkgload_1.3.4 jsonlite_1.8.8 evaluate_0.23 digest_0.6.34 desc_1.4.3
    # callr_3.7.3 brio_1.1.4 stringi_1.8.3 utf8_1.2.4 fansi_1.0.6
    # testthat_3.2.1 colorspace_2.1-0 vctrs_0.6.5 tidyselect_1.2.0
    # pillar_1.9.0 magrittr_2.0.3 tidyr_1.3.1 tibble_3.2.1 stringr_1.5.1
    # purrr_1.0.2 generics_0.1.3 ellipsis_0.3.2 backports_1.4.1
    # viridisLite_0.4.2 rlang_1.1.3 RColorBrewer_1.1-3 R6_2.5.1
    # munsell_0.5.0 lifecycle_1.0.4 labeling_0.4.3 glue_1.7.0
    # farver_2.1.1 cli_3.6.2 MatrixModels_0.5-3 SparseM_1.81 \
    #  numDeriv_2016.8-1.1 dplyr_1.1.4 broom_1.0.5 RcppEigen_0.3.3.9.4
    # Rcpp_1.0.12 nloptr_2.0.3 minqa_1.2.6 scales_1.3.0 lme4_1.1-35.1
    # quantreg_5.97 pbkrtest_0.5.2 abind_1.4-5 carData_3.0-5 isoband_0.2.7
    # gtable_0.3.4 car_3.1-2 corrplot_0.92 ggplot2_3.4.4 rstatix_0.7.2
    # polynom_1.4-1 gridExtra_2.3 ggsignif_0.6.4 cowplot_1.1.3 ggsci_3.0.0
    # ggrepel_0.9.5
    && R --slave -e 'remotes::install_version("matrixStats", version = "1.0.0", upgrade=FALSE)' \
    # rstatix -> 0.7.2, car -> 3.1-2, nloptr -> 2.0.3, lme4 -> 1.1-34, ggrepel -> 0.9.3, Matrix -> 1.6-4
    && R --slave -e 'remotes::install_version("ggpubr", version = "0.6.0", upgrade=FALSE)' \
    # prettyunits_1.2.0 hms_1.1.3 bit_4.0.5 progress_1.2.3 tzdb_0.4.0
    # bit64_4.0.5 vroom_1.6.5 clipr_0.8.0
    && R --slave -e 'remotes::install_version("readr", version = "2.1.4", upgrade=FALSE)' \
    # rappdirs_0.3.3 highr_0.10 fastmap_1.1.1 sass_0.4.8 mime_0.12
    # memoise_2.0.1 cachem_1.0.8 base64enc_0.1-3 yaml_2.3.8 xfun_0.41
    # tinytex_0.49 knitr_1.45 jquerylib_0.1.4 htmltools_0.5.7
    # fontawesome_0.5.2 bslib_0.6.1 systemfonts_1.0.5 sys_3.4.2
    # askpass_1.2.0 uuid_1.2-0 openssl_2.1.1 curl_5.2.0 httr_1.4.7
    # gargle_1.5.2 rematch_2.0.0 xml2_1.3.6 selectr_0.4-2 rstudioapi_0.15.0
    # rmarkdown_2.25 cellranger_1.1.0 textshaping_0.3.7 timechange_0.3.0
    # forcats_1.0.0 ids_1.0.1 googledrive_2.1.1 DBI_1.2.1 blob_1.2.4
    # rvest_1.0.3 reprex_2.1.0 readxl_1.4.3 ragg_1.2.7 modelr_0.1.11
    # lubridate_1.9.3 haven_2.5.4 googlesheets4_1.1.1 dtplyr_1.3.1
    # dbplyr_2.4.0 conflicted_1.2.0
    && R --slave -e 'remotes::install_version("tidyverse", version = "2.0.0", upgrade=FALSE)' \
    # None
    && R --slave -e 'remotes::install_version("pbapply", version = "1.7.2", upgrade=FALSE)' \
    #
    && R --slave -e 'remotes::install_version("hdf5r", version = "1.3.8", upgrade=FALSE)' \
    # parallelly_1.36.0 listenv_0.9.1 globals_0.16.2 future_1.33.1
    # progressr_0.14.0 future.apply_1.11.1 sp_2.1-3
    && R --slave -e 'remotes::install_version("SeuratObject", version = "4.1.4", upgrade=FALSE)' \
    # This needs to be in front of Seurat and needs to be version 1.5.1 to prevent 'ibxml/globals.h: No such file or directory' error with version 2.0.1.1
    && R --slave -e 'remotes::install_version("igraph", version="1.5.1", upgrade=FALSE)' \
    # sitmo_2.0.2. BH_1.84.0-0. spatstat.utils_3.0-4. tensor_1.5. polyclip_1.10-6. deldir_2.0-2. 
    # spatstat.geom_3.2-8. spatstat.data_3.0-4. promises_1.2.1. later_1.3.2. dotCall64_1.1-1. 
    # plyr_1.8.9. bitops_1.0-7. caTools_1.18.2. gtools_3.9.5. lazyeval_0.2.2. commonmark_1.9.1. 
    # sourcetools_0.1.7-1. xtable_1.8-4. httpuv_1.6.14. png_0.1-8. here_1.0.1. RcppTOML_0.2.2. 
    # dqrng_0.3.2. RcppProgress_0.4.2. irlba_2.3.5.1. RcppAnnoy_0.0.22. FNN_1.1.4. goftest_1.2-3. 
    # spatstat.sparse_3.0-3. spatstat.random_3.2-2. spam_2.10-0. RcppArmadillo_0.12.8.0.0. 
    # reshape2_1.4.4. gplots_3.1.3.1. crosstalk_1.2.1. htmlwidgets_1.6.4. shiny_1.8.0. zoo_1.8-12. 
    # reticulate_1.35.0. uwot_0.1.16. spatstat.explore_3.2-6. sctransform_0.4.1. scattermore_1.2. 
    # Rtsne_0.17. ROCR_1.0-11. RANN_2.6.1. plotly_4.10.4. patchwork_1.2.0. miniUI_0.1.1.1. 
    # lmtest_0.9-40. leiden_0.4.3.1. ica_1.0-3. ggridges_0.5.6. fitdistrplus_1.1-11.
    && R --slave -e 'remotes::install_version("Seurat", version = "4.4.0", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("scCustomize", version = "1.1.3", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("textTinyR", version = "1.1.7", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("locfit", version = "1.5-9.8", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("limma", version = "2.9.9", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("Hmisc", version = "5.1-1", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("pheatmap", version = "1.0.12", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("ggforce", version = "0.4.1", upgrade=FALSE)' \
    && R --slave -e 'remotes::install_version("ggnewscale", version = "0.4.9", upgrade=FALSE)'

RUN R --slave -e 'install.packages("devtools")' \
    # credentials_2.0.1 zip_2.3.1 gitcreds_0.1.2 httr2_1.0.0 ini_0.3.1 gert_2.0.1 gh_1.4.0
    # whisker_0.4.1 downlit_0.4.3 xopen_1.0.0 brew_1.0-10 usethis_2.2.2 pkgdown_2.0.7
    # profvis_0.3.8 rcmdcheck_1.4.0 roxygen2_7.3.1 rversions_2.1.2 sessioninfo_1.2.2
    # urlchecker_1.0.1 devtools_2.4.5
    && R --slave -e 'devtools::install_github("jrs95/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = FALSE, upgrade_dependencies=FALSE)'

# BiocManager 1.30.22
RUN R --slave -e 'install.packages("BiocManager")' \
    # 3.14 is for R version 4.1
    # BiocVersion_3.14.0
    && R --slave -e 'BiocManager::install(version = "3.14")' \
    && R --slave -e 'BiocManager::install("edgeR", version = "3.14")' \
    && R --slave -e 'BiocManager::install("qvalue", version = "3.14")' \
    && R --slave -e 'BiocManager::install("multtest", version = "3.14")' \
    && R --slave -e 'BiocManager::install("rhdf5", version = "3.14")' \
    && R --slave -e 'BiocManager::install("variancePartition", version = "3.14")' \
    && R --slave -e 'BiocManager::install("dreamlet", version = "3.14")'

#################################
############## JAVA #############
#################################

# Section build takes 19 seconds and has a size of 0.526 GB.

# adwaita-icon-theme ca-certificates-java gtk-update-icon-cache
# hicolor-icon-theme humanity-icon-theme java-common libasound2
# libasound2-data libatk1.0-0 libatk1.0-data libavahi-client3
# libavahi-common-data libavahi-common3 libcups2 libdrm-amdgpu1 libdrm-common
# libdrm-intel1 libdrm-nouveau2 libdrm-radeon1 libdrm2 libedit2
# libgdk-pixbuf-2.0-0 libgdk-pixbuf2.0-common libgif7 libgl1 libgl1-mesa-dri
# libglapi-mesa libglvnd0 libglx-mesa0 libglx0 libgtk2.0-0 libgtk2.0-common
# liblcms2-2 libllvm15 libnspr4 libnss3 libpciaccess0 libpcsclite1
# libsensors-config libsensors5 libx11-xcb1 libxcb-dri2-0 libxcb-dri3-0
# libxcb-glx0 libxcb-present0 libxcb-randr0 libxcb-sync1 libxcb-xfixes0
# libxcomposite1 libxcursor1 libxdamage1 libxfixes3 libxi6 libxinerama1
# libxrandr2 libxshmfence1 libxtst6 libxxf86vm1 openjdk-17-jdk
# openjdk-17-jdk-headless openjdk-17-jre openjdk-17-jre-headless
# shared-mime-info ubuntu-mono
RUN apt-get install -y --no-install-recommends openjdk-17-jdk

##################################
############## OTHER #############
##################################

# Using a fixed version of GenotypeHarmonizer-1.4.27 from MJB.
RUN cd /opt \
    && gdown 16bjxeP0cb4BCVBVxjxKxGu_8r8F-Ppd3 \
    && unzip GenotypeHarmonizer-1.4.27-SNAPSHOT.zip \
    && rm GenotypeHarmonizer-1.4.27-SNAPSHOT.zip \
    && chmod 777 GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar

RUN cd /opt \
    && gdown 1rPm-n8Zteq5v0t__OjZdi_YjdO5Mwl81 \
    && mkdir Genotype-IO-1.0.6-SNAPSHOT \
    && mv Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar Genotype-IO-1.0.6-SNAPSHOT/GenotypeIO.jar \
    && chmod 777 Genotype-IO-1.0.6-SNAPSHOT/GenotypeIO.jar

RUN cd /opt \
    && wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz \
    && tar -xzf ldstore_v2.0_x86_64.tgz \
    && rm ldstore_v2.0_x86_64.tgz \
    && chown -R root:root ldstore_v2.0_x86_64 \
    && chmod -R o+rx ldstore_v2.0_x86_64

RUN cd /opt \
    && mkdir plink2 \
    && cd plink2 \
      && wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240205.zip \
      && unzip -q plink2_linux_x86_64_20240205.zip \
      && rm plink2_linux_x86_64_20240205.zip

####################################
################ CLEAN #############
####################################

RUN apt-get clean \
    && apt-get autoremove -y
