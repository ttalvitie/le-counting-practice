.PHONY: solutions clean

solutions:
	for dir in `ls -d */` ; do $(MAKE) -C $$dir ; done

clean:
	for dir in `ls -d */` ; do $(MAKE) -C $$dir clean ; done
