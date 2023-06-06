OPTIM = -O3 -fno-loop-optimize -fno-aggressive-loop-optimizations
CONF_DEBUG = -O0 -g -DDEBUG -DAS_N_U=20 -DAS_N_V=6 -DAS_SINGLE_FLOAT -DAS_AS_RECORD_COST -DAS_AS_RECORD_COST_N=15 -DAS_COST_TRUNCATE
CONF_SINGLE = $(OPTIM) -DAS_N_U=20 -DAS_N_V=6 -DAS_SINGLE_FLOAT -DAS_COST_TRUNCATE
CONF_DOUBLE = $(OPTIM) -DAS_N_U=20 -DAS_N_V=6 -DAS_COST_TRUNCATE

CONF?=DOUBLE
ifeq ($(CONF), DOUBLE)
CONF_OPTIONS = $(CONF_DOUBLE)
else ifeq ($(CONF), SINGLE)
CONF_OPTIONS = $(CONF_SINGLE)
else ifeq ($(CONF), DEBUG)
CONF_OPTIONS = $(CONF_DEBUG)
else
$(error only allow DOUBLE, SINGLE, DEBUG as CONF choices)
endif

CC = gcc
WARN_FLAGS = -Wall -W -Wwrite-strings -Winline -Wstrict-prototypes -Wnested-externs -Wpointer-arith -Wcast-align -Wcast-qual -Wshadow
CC_FLAGS = -fstack-usage -fwrapv -fPIC $(WARN_FLAGS) $(CONF_OPTIONS)
LINK_FLAGS = -lm
AR = ar
ARFLAGS = ruvs
PRODUCT = solveActiveSet.so
BIN_DIR = ./bin
INCLUDES = -Icommon/ -Ilib/
SOURCES = common/solveActiveSet.c solveActiveSet_cg.c solveActiveSet_chol.c solveActiveSet_qr.c solveActiveSet_qr_naive.c common/setupWLS.c lib/chol_math.c lib/qr_updates.c lib/qr_wrapper.c lib/qr_solve/qr_solve.c lib/qr_solve/r8lib_min.c lib/sparse_math.c

all: $(BIN_DIR)/$(PRODUCT)
	$(CC) $(CC_FLAGS) -o $(BIN_DIR)/main.o main.c $(BIN_DIR)/$(PRODUCT) $(INCLUDES) $(LINK_FLAGS)

$(BIN_DIR)/$(PRODUCT) : $(SOURCES)
	@mkdir -p $(dir $@)
	$(CC) -shared $(CC_FLAGS) -o $@ $(WARN_FLAGS) $(INCLUDES) $^

clean:
	@rm -rf *.o
	@rm -rf *.so
	@rm -rf *.su
	@rm -rf bin
