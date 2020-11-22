CC 	   = g++
FFTW_LOC = ${HOME}/Install/fftw3
EIGEN_LOC = ${HOME}/Install/eigen
CFLAGS     = -O3 -Wno-unused-result -Wno-write-strings -std=c++11
LIBS      = -lm -O3
#CFLAGS     = -g -I${FFTW_LOC}/include -I${EIGEN_LOC} -Wno-unused-result -Wno-write-strings
#LIBS      = -g -lm -lfftw3_mpi -lfftw3 -L${FFTW_LOC}/lib


#############################################################################
# nothing should be changed below here

SRCS = main.cpp read_lammpstrj.cpp pbc_utils.cpp lc_order.cpp \
       
       
			 


OBJS = ${SRCS:.cpp=.o}

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c  $<

postproc-lammpstrj:  ${OBJS}
	$(CC) ${CFLAGS} ${DFLAGS} -o $@ ${OBJS} $(LIBS)

clean:
	rm -f *.o
	rm -f postproc-lammpstrj
	rm -f *~

