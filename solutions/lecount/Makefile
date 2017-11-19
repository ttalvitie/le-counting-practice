CXX ?= g++
CFLAGS ?= -std=c++11 -Wall -O2 -Wno-sign-compare -march=native -g
LDFLAGS ?= -lgmp -lgmpxx
SRCS := $(shell find . -name '*.cpp')
BINS := lecount

.PHONY: clean all

all: $(BINS)

lecount: lecount.o arguments.o tools.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -f $(SRCS:%.cpp=%.o) $(SRCS:%.cpp=%.d) $(BINS)

-include $(SRCS:%.cpp=%.d)
