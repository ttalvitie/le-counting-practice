.PHONY: all clean

all:
	$(MAKE) -C instances
	$(MAKE) -C solutions

clean:
	$(MAKE) -C solutions clean
