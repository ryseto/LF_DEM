# load Makefile inputs
makeconfig = Makefile.in
ifneq ("$(wildcard $(makeconfig))","")
	include $(makeconfig)
endif

all:
	python setup.py build_ext --inplace

install:
	cp _pyLFDEM*.so $(PYTHON_SITE_PACKAGES)
	cp pyLFDEM.py $(PYTHON_SITE_PACKAGES)

swig:
	$(SWIG) -python -c++ pyLFDEM.i
