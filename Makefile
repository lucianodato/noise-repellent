#!/usr/bin/make -f
OPTIMIZATIONS ?= -msse -msse2 -mfpmath=sse -ffast-math -fomit-frame-pointer -fno-finite-math-only
PREFIX ?= /usr/local
CFLAGS ?= $(OPTIMIZATIONS) -Wall
LOADLIBES=-lm
STRIP?=strip
STRIPFLAGS?=-s
DEBUG?=0

#for debug building
ifeq ($(DEBUG), 1)
  override CFLAGS += -O0 -g3 -DDEBUG
else
  override CFLAGS += -O3 -DNDEBUG
endif

LV2DIR ?= $(PREFIX)/lib/lv2
LV2NAME=nrepel
BUNDLE=nrepel.lv2
BUILDDIR=build/
SRCDIR=src/
TTLDIR=lv2ttl/
DOCDIR=doc/

###############################################################################
#get os and configure compiling flags and library extension
LIB_EXT=.so

targets=
UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
  LV2LDFLAGS=-dynamiclib
  LIB_EXT=.dylib
  EXTENDED_RE=-E
  STRIPFLAGS=-u -r -arch all
else
  LV2LDFLAGS=-Wl,-Bstatic -Wl,-Bdynamic
  LIB_EXT=.so
  EXTENDED_RE=-r
endif

ifneq ($(XWIN),)
  CC=$(XWIN)-gcc
  STRIP=$(XWIN)-strip
  LV2LDFLAGS=-Wl,-Bstatic -Wl,-Bdynamic -Wl,--as-needed
  LIB_EXT=.dll
  override LDFLAGS += -static-libgcc -static-libstdc++
endif

targets+=$(BUILDDIR)$(LV2NAME)$(LIB_EXT)

###############################################################################

# check for build-dependencies
ifeq ($(shell pkg-config --exists lv2 || echo no), no)
  $(error "LV2 SDK was not found")
endif
ifeq ($(shell pkg-config --exists fftw3f || echo no), no)
  $(error "FFTW was not found")
endif

override CFLAGS += -fPIC -std=c99
override LOADLIBES += `pkg-config --libs fftw3f`

# build target definitions
default: all

all: $(BUILDDIR)manifest.ttl $(BUILDDIR)$(LV2NAME).ttl $(targets)

$(BUILDDIR)manifest.ttl: $(TTLDIR)manifest.ttl.in
	@mkdir -p $(BUILDDIR)
	sed "s/@LIB_EXT@/$(LIB_EXT)/" \
	  $(TTLDIR)manifest.ttl.in > $(BUILDDIR)manifest.ttl

$(BUILDDIR)$(LV2NAME).ttl: $(TTLDIR)$(LV2NAME).ttl
	@mkdir -p $(BUILDDIR)
	cp $(TTLDIR)$(LV2NAME).ttl $(BUILDDIR)$(LV2NAME).ttl

$(BUILDDIR)$(LV2NAME)$(LIB_EXT): $(SRCDIR)$(LV2NAME).c
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) \
		-o $(BUILDDIR)$(LV2NAME)$(LIB_EXT) $(SRCDIR)$(LV2NAME).c \
		-shared $(LV2LDFLAGS) $(LDFLAGS) $(LOADLIBES)

ifeq ($(DEBUG), 0)
	$(STRIP) $(STRIPFLAGS) $(BUILDDIR)$(LV2NAME)$(LIB_EXT)
endif

# install/uninstall/clean/doc target definitions
install: all
	install -d $(DESTDIR)$(LV2DIR)/$(BUNDLE)
	install -m644 $(BUILDDIR)$(LV2NAME)$(LIB_EXT) $(DESTDIR)$(LV2DIR)/$(BUNDLE)
	install -m644 $(BUILDDIR)manifest.ttl $(BUILDDIR)$(LV2NAME).ttl $(DESTDIR)$(LV2DIR)/$(BUNDLE)

uninstall:
	rm -f $(DESTDIR)$(LV2DIR)/$(BUNDLE)/manifest.ttl
	rm -f $(DESTDIR)$(LV2DIR)/$(BUNDLE)/$(LV2NAME).ttl
	rm -f $(DESTDIR)$(LV2DIR)/$(BUNDLE)/$(LV2NAME)$(LIB_EXT)
	-rmdir $(DESTDIR)$(LV2DIR)/$(BUNDLE)

clean:
	rm -f $(BUILDDIR)manifest.ttl $(BUILDDIR)$(LV2NAME).ttl $(BUILDDIR)$(LV2NAME)$(LIB_EXT) lv2syms
	-test -d $(BUILDDIR) && rmdir $(BUILDDIR) || true

doc:
	doxygen -s $(DOCDIR)doxygen.conf

.PHONY: doc clean all install uninstall
