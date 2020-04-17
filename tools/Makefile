# load Makefile inputs

makeconfig = config/Makefile_config.mk
gotconfig = 0
ifneq ("$(wildcard $(makeconfig))","")
        include $(makeconfig)
        gotconfig = 1
endif

PYTHON_SITE_PACKAGES = $(shell env python -c "import site;print(site.getsitepackages()[0])")

all:
	cd file_reading; python3 setup.py build_ext --inplace; cd ..
	cp lfdem_utils.py $(PYTHON_SITE_PACKAGES)
	cp pyaplot.py $(PYTHON_SITE_PACKAGES)
	cp yapgen.py $(PYTHON_SITE_PACKAGES)
	cp file_reading/clfdem_file*so $(PYTHON_SITE_PACKAGES)
	cp file_reading/clfdem_file.pyx $(PYTHON_SITE_PACKAGES)

install:
	cp generateYaplotFile.pl $(INSTALL_DIR)
	cp yapgen.py $(INSTALL_DIR)
	cp binconf_info.py $(INSTALL_DIR)
	cp wall_conf.py $(INSTALL_DIR)
