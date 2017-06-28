DIRS = lighthouse

all:
	@for d in $(DIRS); do $(MAKE) -C examples/$$d || exit $$?; done

clean:
	@for d in $(DIRS); do $(MAKE) -C examples/$$d GNU=$(GNU) clean || exit $$?; done
    