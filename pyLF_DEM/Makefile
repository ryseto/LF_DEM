# load Makefile inputs

SWIG = $(shell which swig)
PYTHON_SITE_PACKAGES = $(shell env python -c "import site;print(site.getsitepackages()[0])")

# you can override SWIG and PYTHON_SITE_PACKAGES in a file called Makefile.in
makeconfig = Makefile.in
ifneq ("$(wildcard $(makeconfig))","")
	include $(makeconfig)
endif

all:
	python setup.py build_ext --inplace

install:
	cp _pyLFDEM*.so $(PYTHON_SITE_PACKAGES)
	cp pyLFDEM.py $(PYTHON_SITE_PACKAGES)
	cp pyLFDEM_helper.py $(PYTHON_SITE_PACKAGES)

swig:
	$(SWIG) -python -c++ pyLFDEM.i
