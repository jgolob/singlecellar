# golob/cellranger:6.0.2
FROM ubuntu:18.04
RUN export DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
curl \
pigz \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/local/

RUN curl -o cellranger-6.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.2.tar.gz?Expires=1625024401&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MjUwMjQ0MDF9fX1dfQ__&Signature=MDIwvA55z2o2LAWDYQPmmHmFTgjxM1H4ROMaPLfECrJEksp4FRvVaKX9VJw5uY0~4jotWCNGA1qfG-HYvbS-dIQb1ELC2RVEjmASmJ2EdmYIOaGBy8hQaB~YLal5JI8Op89yxlB-RLIThnU0Y6JSaweT9i0Q5EdrJEzrZp8J20eZhwjtbaKV4OOAM6GRPGSgjVXjvsHKk4EXs2lNC-dlO~ZftmidE3O23O~O7ccaZyc5wPl6XhXl0~DF99bPIgI9KWF~oUPeAPdPAM3mFoDf0GkWsDeJTCQeJn1tCcjUVbWMKkbkFoXci4PQh8bQBwFZjU23YoxNQka1DAyDL4a~Rw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" \
&& tar xzvf cellranger-6.0.2.tar.gz \
&& rm cellranger-6.0.2.tar.gz

RUN echo PATH=/usr/local/cellranger-6.0.2/bin/:$PATH >> /root/.bashrc
RUN ln -s /usr/local/cellranger-6.0.2/bin/cellranger /usr/local/bin/cellranger

WORKDIR /root/



