FROM rocker/r-ubuntu:20.04 as dep

RUN apt update && apt -y install r-cran-tidyverse && rm -r /var/lib/apt/lists/
RUN apt update && apt -y install r-cran-plyr && rm -r /var/lib/apt/lists/
RUN apt update && apt -y install r-cran-shiny r-cran-shinyjs r-cran-shinyfiles && rm -r /var/lib/apt/lists/



FROM dep as install

COPY . /root/
WORKDIR /root
RUN R CMD INSTALL --no-multiarch --with-keep.source hgnc
RUN R CMD INSTALL --no-multiarch --with-keep.source rGEO.data
RUN R CMD INSTALL --no-multiarch --with-keep.source rGEO
RUN R CMD INSTALL --no-multiarch --with-keep.source qGSEA



FROM dep

LABEL maintainer="Zhuoer Dong <dongzhuoer@mail.nankai.edu.cn>"

COPY --from=install /usr/local/lib/R/site-library /usr/local/lib/R/site-library
WORKDIR /root
