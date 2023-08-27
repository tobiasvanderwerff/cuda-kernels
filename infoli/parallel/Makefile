
# --- MACROS
# define program name
MAIN= infoli_simple

# define the C compiler to use
CC= nvcc

# define any compile-time flags
CFLAGS= -O3 -pg --relocatable-device-code=true

# define any libraries to link into executable
LIBS= -lm

# define C source files
SRCS= src/*.cu src/**/*.cu

# define C header files
HDRS= src/*.h src/**/*.h

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
