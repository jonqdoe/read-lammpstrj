CC 	   = mpic++ -std=c++11
FFTW_LOC = ${HOME}/Install/fftw3
#FFTW_LOC = /opt/seas/pkg/gcc/fftw3/mpi/double/3.3.7

EIGEN_LOC = ${HOME}/Install/eigen
CFLAGS     = -O3 -I${FFTW_LOC}/include -Wno-unused-result -Wno-write-strings -std=c++11 -I${FFTW_HOME}/include
LIBS      = -lm -O3 -lfftw3_mpi -lfftw3 -L${FFTW_LOC}/lib -L${FFTW_HOME}/lib
#CFLAGS     = -g -I${FFTW_LOC}/include -I${EIGEN_LOC} -Wno-unused-result -Wno-write-strings
#LIBS      = -g -lm -lfftw3_mpi -lfftw3 -L${FFTW_LOC}/lib


#############################################################################
# nothing should be changed below here

SRCS = main.cpp read_lammpstrj.cpp pbc_utils.cpp log_space.cpp \
       msd.cpp rdf.cpp lc_order.cpp nl-utils.cpp fftw_mpi_wrappers.cpp \
	   van-hove.cpp io_utils.cpp mesh_utils.cpp trim_lammpstrj.cpp cluster_analysis.cpp

TAGS := .tags
       
			 

OBJECTFOLDER := objects

OBJS = ${SRCS:.cpp=.o}

OBJS_MAKE = $(addprefix $(OBJECTFOLDER)/, ${OBJS})

$(OBJECTFOLDER)/%.o: %.cpp
	@mkdir -p $(dir $@)
	${CC} -g ${CFLAGS} ${DFLAGS} -c $< -o $@

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c  $<

postproc-lammpstrj:  ${OBJS_MAKE}
	$(CC) ${CFLAGS} ${DFLAGS} -o $@ ${OBJS_MAKE} $(LIBS)
	ctags -R -f $(TAGS) .

clean:
	rm -fr *.o $(OBJECTFOLDER)
	rm -f postproc-msd
	rm -f *~

