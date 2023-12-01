# syntax=docker/dockerfile:1
########################################################################################################################
# SurVClusterer build stage
########################################################################################################################

FROM alpine:3.18.4 AS build

WORKDIR /usr/src

COPY . .

RUN apk update
RUN apk add make cmake build-base zlib-dev bzip2-dev xz-dev curl-dev autoconf

RUN chmod +x build_htslib.sh

RUN ./build_htslib.sh && \
  cmake -DCMAKE_BUILD_TYPE=Release . && make

########################################################################################################################
# sv_merging image
########################################################################################################################

FROM alpine:3.18.4

WORKDIR /usr/src

RUN apk update
RUN apk add libstdc++ libbz2 libcurl xz-dev

COPY --from=build \
    /usr/src/htslib-1.18 \
    /usr/src/htslib-1.18

COPY --from=build \
    /usr/src/clusterer \
    /usr/src/clusterer

COPY --from=build \
    /usr/src/compare-del \
    /usr/src/compare-del

COPY --from=build \
    /usr/src/compare-ins \
    /usr/src/compare-ins

ENV PATH=${PATH}:/usr/src/htslib-1.18