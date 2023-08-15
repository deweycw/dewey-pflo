FROM deweycw/pflotran-base:latest

RUN apt update && apt upgrade -y

WORKDIR /pflotran
RUN git clone https://github.com/deweycw/dewey-pflo.git

WORKDIR /pflotran/dewey-pflo/src/pflotran
RUN git checkout docker
RUN make pflotran
ENV PFLOTRAN_DIR='/pflotran/dewey-pflo/'

WORKDIR /