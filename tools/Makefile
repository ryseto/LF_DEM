# load Makefile inputs
makeconfig = config/Makefile_config.mk
gotconfig = 0
ifneq ("$(wildcard $(makeconfig))","")
        include $(makeconfig)
        gotconfig = 1
endif

install:
	cd file_reading; python setup.py build_ext --inplace; cd ..
	cp generateYaplotFile.pl $(INSTALL_DIR)
	cp lfdem_utils.py $(ANACONDA_DIR)
	cp pyaplot.py $(ANACONDA_DIR)
	cp yapgen.py $(INSTALL_DIR)
	cp yapgen.py $(ANACONDA_DIR)
	cp binconf_info.py $(INSTALL_DIR)
	cp wall_conf.py $(INSTALL_DIR)
	cp file_reading/clfdem_file*so $(ANACONDA_DIR)
	cp file_reading/clfdem_file.pyx $(ANACONDA_DIR)
