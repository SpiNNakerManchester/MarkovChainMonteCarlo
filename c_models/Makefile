DIRS = mcmc_models mcmc_coordinator mcmc_root_finder mcmc_cholesky

all: $(DIRS)
	@for d in $(DIRS); do $(MAKE) -C $$d || exit $$?; done

clean: $(DIRS)
	@for d in $(DIRS); do $(MAKE) -C $$d GNU=$(GNU) clean || exit $$?; done
    