BUILDDIR = bin
OBJDIR = build

EXECNAME_1 = gfa1_parser
IDIR_1 = include/$(EXECNAME_1)
SRC_1 = src/$(EXECNAME_1)

EXECNAME_2 = unitig_distance
IDIR_2 = include/$(EXECNAME_2)
SRC_2 = src/$(EXECNAME_2)

# CXX = g++-13
CXX ?= g++

CXXFLAGS_1 = -std=c++11 -march=native -O2 -pedantic -Wall -I$(IDIR_1) -I$(IDIR_1)/containers
LDFLAGS_1 =

CXXFLAGS_2 = -std=c++11 -pthread -march=native -O2 -pedantic -Wall -I$(IDIR_2)
LDFLAGS_2 = -pthread

SOURCES_1 = $(shell find $(SRC_1) -type f -name *.cpp)
OBJECTS_1 = $(patsubst $(SRC_1)/%,$(OBJDIR)/%,$(SOURCES_1:.cpp=.o))
DEPENDS_1 = $(patsubst $(SRC_1)/%,$(OBJDIR)/%,$(SOURCES_1:.cpp=.d))

SOURCES_2 = $(shell find $(SRC_2) -type f -name *.cpp)
OBJECTS_2 = $(patsubst $(SRC_2)/%,$(OBJDIR)/%,$(SOURCES_2:.cpp=.o))
DEPENDS_2 = $(patsubst $(SRC_2)/%,$(OBJDIR)/%,$(SOURCES_2:.cpp=.d))

.PHONY: all clean

all: $(EXECNAME_1) $(EXECNAME_2)

clean:
	\rm $(OBJDIR)/*.o $(OBJDIR)/*.d $(BUILDDIR)/$(EXECNAME_1) $(BUILDDIR)/$(EXECNAME_2)

-include $(DEPENDS_1)

$(EXECNAME_1): $(OBJECTS_1)
	mkdir -p $(BUILDDIR); $(CXX) $(LDFLAGS_1) $(OBJECTS_1) -o $(BUILDDIR)/$(EXECNAME_1)

$(OBJDIR)/%.o: $(SRC_1)/%.cpp
	mkdir -p $(OBJDIR); $(CXX) $(CXXFLAGS_1) -MMD -MP -c $< -o $@

-include $(DEPENDS_2)

$(EXECNAME_2): $(OBJECTS_2)
	mkdir -p $(BUILDDIR); $(CXX) $(LDFLAGS_2) $(OBJECTS_2) -o $(BUILDDIR)/$(EXECNAME_2)

$(OBJDIR)/%.o: $(SRC_2)/%.cpp
	mkdir -p $(OBJDIR); $(CXX) $(CXXFLAGS_2) -MMD -MP -c $< -o $@
