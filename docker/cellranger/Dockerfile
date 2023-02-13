# golob/cellranger:7.10
FROM --platform=linux/amd64 ubuntu:18.04
RUN export DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
curl \
pigz \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/local/

RUN curl -o cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1676340738&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NzYzNDA3Mzh9fX1dfQ__&Signature=B4uFQ3sUSrHi~c~F8mq30GxcaXvDw5X~b8dfwb2~efzo2Tj0H5oxccCgVN9xpvo~ib4WLwjXrpiYKvkVfNHEV-0Cr-GXIZNWhvjFxHqt19O5NHOjkWzVqk4bMA5WaeXpnA-DCiqfL-N5lAo8MyfPabP401xyMs61ZcB0flajUWut4170Fz4K6f~9bi5U~h0DVjxByIYjX8khx7ayVed~hnbAe~cRFi7ZV7b7zUKTWDFw~UQPP66hS1Wbduhh57KUNUlVNfNg5UowDtLj2PJMsrCX9gvO~TZDGGsajrMxt0VfnfH5NCyJ8eORy12c0At1H8dzLYgnoysHHJv134unug__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" \
&& tar xzvf cellranger-7.1.0.tar.gz \
&& rm cellranger-7.1.0.tar.gz

RUN echo PATH=/usr/local/cellranger-7.1.0/bin/:$PATH >> /root/.bashrc
RUN ln -s /usr/local/cellranger-7.1.0/bin/cellranger /usr/local/bin/cellranger

WORKDIR /root/



