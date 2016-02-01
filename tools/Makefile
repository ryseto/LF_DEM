# load Makefile inputs
makeconfig = config/Makefile_config.mk
gotconfig = 0
ifneq ("$(wildcard $(makeconfig))","")
        include $(makeconfig)
        gotconfig = 1
endif

install:
	cp generateYaplotFile.pl $(INSTALL_DIR)
	cp lfdem_file.py $(INSTALL_DIR)
	cp lfdem_file.py $(ANACONDA_DIR)
	cp pyaplot.py $(INSTALL_DIR)
	cp pyaplot.py $(ANACONDA_DIR)
	cp yapgen.py $(INSTALL_DIR)
	cp binconf_info.py $(INSTALL_DIR)
	cp wall_conf.py $(INSTALL_DIR)
