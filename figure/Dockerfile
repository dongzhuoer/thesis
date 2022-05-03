FROM dongzhuoer/thesis 

LABEL maintainer="Zhuoer Dong <dongzhuoer@mail.nankai.edu.cn>"

RUN apt update && apt -y install r-cran-rmarkdown r-cran-prettydoc && rm -r /var/lib/apt/lists/
RUN apt update && apt -y install pandoc && rm -r /var/lib/apt/lists/

COPY *deb /root/
RUN dpkg -i /root/*deb
