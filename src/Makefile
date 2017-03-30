all: release
EXEC = ancestry
EXEC_SEQ = ancestrySeq
UTIL_EXEC = 
DIR_EXEC = ../executable
DIR_EXEC_DBG = ../executable/dbg

$(DIR_EXEC):
	mkdir -p $@
$(DIR_EXEC_DBG):
	mkdir -p $@

##################################################
# Third-party libs.
TABIX_INC = third/tabix
TABIX_LIB = third/tabix/libtabix.a

GSL_INC = third/gsl/include
GSL_LIB = third/gsl/lib/libgsl.a

BZIP2_INC = third/bzip2
BZIP2_LIB = third/bzip2/libbz2.a

BLAS_LIB = third/openblas/lib/libopenblas.a

$(TABIX_INC) $(TABIX_LIB):
	(cd third; make tabix)
$(GSL_INC) $(GSL_LIB):
	(cd third; make gsl)
$(BZIP2_INC) $(BZIP2_LIB):
	(cd third; make bzip2)

THIRD_INC = $(TABIX_INC) $(GSL_INC) $(BZIP2_INC)
THIRD_LIB = $(TABIX_LIB) $(GSL_LIB) $(BZIP2_LIB) $(BLAS_LIB)
##################################################
# Our libs.
BASE_INC = ./base
BASE_LIB = ./base/lib-base.a
BASE_LIB_DBG = ./base/lib-dbg-base.a

VCF_INC = ./libVcf
VCF_LIB = ./libVcf/lib-vcf.a
VCF_LIB_DBG = ./libVcf/lib-dbg-vcf.a

$(BASE_LIB):
	(cd base; make)
$(BASE_LIB_DBG):
	(cd base; make debug)
$(VCF_LIB):
	(cd libVcf; make)
$(VCF_LIB_DBG):
	(cd libVcf; make debug)

##################################################

INCLUDE = $(VCF_INC) $(THIRD_INC) $(BASE_INC)
LIB = $(VCF_LIB) $(BASE_LIB) $(THIRD_LIB) 
LIB_DBG = $(VCF_LIB_DBG) $(BASE_LIB_DBG) $(THIRD_LIB)
CXX_INCLUDE = $(addprefix -I, $(INCLUDE)) -I.
CXX_LIB = $(LIB) -lz -lm -lpthread
CXX_LIB_DBG = $(LIB_DBG) -lz -lm -lpthread


DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS -std=c++0x #-Wall

.PHONY: release debug ## lib lib-dbg

# to build lib, we will use reverse order
# so third party lib is build first.
reverse = $(if $(1),$(call reverse,$(wordlist 2,$(words $(1)),$(1)))) $(firstword $(1))
lib = $(call reverse,$(LIB))
lib-dbg = $(call reverse,$(LIB_DBG))

release: CXX_FLAGS = -O2 -DNDEBUG $(DEFAULT_CXXFLAGS) -static
release: $(DIR_EXEC)/$(EXEC) $(DIR_EXEC)/$(EXEC_SEQ)
$(DIR_EXEC)/$(EXEC): Main.o \
                     $(lib) \
                     |$(DIR_EXEC)
	g++ -o $@ $< $(CXX_FLAGS) $(CXX_LIB)
$(DIR_EXEC)/$(EXEC_SEQ): Main.Seq.o \
                         $(lib) \
                         |$(DIR_EXEC)
	g++ -o $@ $< $(CXX_FLAGS) $(CXX_LIB)

debug: CXX_FLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS) 
debug: $(DIR_EXEC_DBG)/$(EXEC) $(DIR_EXEC_DBG)/$(EXEC_SEQ) 
$(DIR_EXEC_DBG)/$(EXEC): Main.o \
                         $(lib-dbg) \
                         | $(DIR_EXEC_DBG)
	g++ -o $@ $< $(CXX_FLAGS) $(CXX_LIB_DBG) 
$(DIR_EXEC_DBG)/$(EXEC_SEQ): Main.Seq.o \
                             $(lib-dbg) \
                             | $(DIR_EXEC_DBG)
	g++ -o $@ $< $(CXX_FLAGS) $(CXX_LIB_DBG) 


##################################################
GitVersion.h: ../.git/HEAD ../.git/index
	echo "const char *gitVersion = \"$(shell git rev-parse HEAD)\";" > $@

-include Main.d
Main.o: Main.cpp GitVersion.h
	g++ -MMD -c $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include Main.Seq.d
Main.Seq.o: Main.Seq.cpp GitVersion.h
	g++ -MMD -c $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__


##################################################
# build utils
util: $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC))
define BUILD_util
  TAR := $(DIR_EXEC)/$(notdir $(basename $(1)))
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O2 $(DEFAULT_CXXFLAGS) -static
  $$(TAR): $$(SRC) $(LIB) | $(DIR_EXEC)
	g++ -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util, $(s))))

util-dbg: $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))
define BUILD_util_dbg
  TAR := $(DIR_EXEC_DBG)/$(notdir $(basename $(1)))
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O0 -ggdb $(DEFAULT_CXXFLAGS)
  $$(TAR): $$(SRC) $(LIB_DBG) | $(DIR_EXEC_DBG)
	g++ -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB_DBG)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util_dbg, $(s))))


clean: 
	rm -rf *.o *.d \
        $(addprefix $(DIR_EXEC)/, $(EXEC)) \
        $(addprefix $(DIR_EXEC_DBG)/, $(EXEC)) \
        $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC)) \
        $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))

deepclean: clean
	rm -rf *~
	(cd base; make clean)
	(cd libVcf; make clean)
	(cd third; make deepclean)
	(cd libsrc; make clean)

# archive 
DATE=$(shell date '+%m%d')
tar:
	tar zvchf rvtest.$(DATE).tgz *.h *.cpp Makefile .git/HEAD .git/index third base libVcf libsrc

README.md:README.wiki
doc: README.md
	java -jar third/wiki2html.jar README.wiki > README.html 
	pandoc -f html -t markdown README.html > README.md 
