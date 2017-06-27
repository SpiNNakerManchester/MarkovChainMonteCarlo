all: mcmc.aplx mcmc_coordinator.aplx

%.aplx: %.Makefile %.c
	"$(MAKE)" -f $<

clean:
	for d in mcmc mcmc_coordinator; do ("$(MAKE)" -f $$d.Makefile clean) || exit $$?; done
    