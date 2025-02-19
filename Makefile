# directories
SDIR = src
ODIR = build
IDIR = include
#LDIR = /usr/lib/x86_64-linux-gnu

TARGET = scf.x

ARGS ?= geometry/h2o.xyz basis/6311.basis

DEFS ?= NDEBUG
DEFS_ = $(addprefix -D, $(DEFS))

# libraries
# LIBS_ := lapack lapacke boost_math_c99 blas
LIBS_ := lapack lapacke boost_math_c99 blas
LIBS = $(addprefix -l, $(LIBS_))
# sources & target objects
SRCS_ := $(shell ls $(SDIR))
OBJS = $(addprefix $(ODIR)/, $(SRCS_:%.cpp=%.o))

# default C++ flags
CXXFLAGS ?= -O3
override CXXFLAGS += -I$(IDIR)

# defult compiler
ifeq ($(origin CXX), default)
	CXX = g++
endif

.PHONY: all clean run

all: $(OBJS) $(TARGET)

$(ODIR)/%.o: $(SDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(DEFS_) $(CXXFLAGS) $(LIBS) -g -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) $^ $(LIBS) -o $(TARGET)

clean:
	@rm -rf $(ODIR) $(TARGET)
	@echo "Project cleaned!"

run: $(TARGET)
	./$(TARGET) $(ARGS)

clean-log:
	@rm -f logs/*
	@echo "Logs cleaned!"
