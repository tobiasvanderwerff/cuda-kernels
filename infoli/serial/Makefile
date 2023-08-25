
# --- MACROS
# define program name
MAIN= infoli_simple

# define the C compiler to use
CC= nvcc

# define any compile-time flags
CFLAGS= -O3 -pg

# define any libraries to link into executable
LIBS= -lm

# define C source files
SRCS= infoli_simple.cu

# define C header files
HDRS= infoli_simple.h

# --- TARGETS
all: ${MAIN}

#Builds the program
${MAIN}: ${SRCS} ${HDRS}
	@echo #
	@echo "-- BUILDING PROGRAM --"
	${CC} ${SRCS} ${CFLAGS} ${LIBS} -o ${MAIN}

clean:
	@echo #
	@echo "-- CLEANING PROJECT FILES --"
	$(RM) *.o ${MAIN}

