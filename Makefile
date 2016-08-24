#!/usr/bin/make -f

PREFIX ?= /usr
LIBDIR ?= lib
LV2DIR ?= $(PREFIX)/$(LIBDIR)/lv2

#Basic Flags
OPTIMIZATIONS ?= -msse -msse2 -mfpmath=sse -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only

LDFLAGS ?= -Wl,--as-needed -shared -Wl,-Bstatic -Wl,-Bdynamic `pkg-config fftw3 --libs`
CFLAGS ?= $(OPTIMIZATIONS) -Wall -fPIC -DPIC -lm `pkg-config fftw3 --cflags --libs`

BUNDLE = nrepel.lv2
LIB_EXT=.so
###############################################################################

#library detection
ifeq ($(shell pkg-config --exists lv2 || echo no), no)
  $(error "LV2 SDK was not found")
endif
ifeq ($(shell pkg-config --exists fftw3 || echo no), no)
  $(error "FFTW was not found")
endif

#directory creation
$(BUNDLE): manifest.ttl nrepel.ttl nrepel$(LIB_EXT)
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp manifest.ttl nrepel.ttl $(BUNDLE)
	mv nrepel$(LIB_EXT) $(BUNDLE)

#file compiling
nrepel$(LIB_EXT): nrepel.c
	$(CXX) -o nrepel$(LIB_EXT) \
		$(CFLAGS) \
		nrepel.c \
		$(LV2FLAGS) $(LDFLAGS)


#targets
test: #ttl files control
	sord_validate_lv2 $(BUNDLE)

install: $(BUNDLE)
	install -d $(DESTDIR)$(LV2DIR)/$(BUNDLE)
	install -t $(DESTDIR)$(LV2DIR)/$(BUNDLE) $(BUNDLE)/*

uninstall:
	rm -rf $(DESTDIR)$(LV2DIR)/$(BUNDLE)

clean:
	rm -rf $(BUNDLE)

.PHONY: test clean install uninstall
