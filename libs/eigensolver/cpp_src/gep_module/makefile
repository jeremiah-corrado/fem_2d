main: solve_gep.o
	-${CXXLINKER} -g -o solve_gep solve_gep.o ${SLEPC_SYS_LIB} -lrt
	${RM} solve_gep.o

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common