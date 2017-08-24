#!/usr/bin/make -f
OPTIMIZATIONS ?= -msse -msse2 -mfpmath=sse -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only
PREFIX ?= /usr/local
CFLAGS ?= $(OPTIMIZATIONS) -Wall

STRIP?=strip
STRIPFLAGS?=-s
DEBUG?=0

nrepel_VERSION?=$(shell git describe --tags HEAD 2>/dev/null | sed 's/-g.*$$//;s/^v//' || echo "LV2")
###############################################################################
LIB_EXT=.so

LV2DIR ?= $(PREFIX)/lib/lv2
LOADLIBES=-lm
LV2NAME=nrepel
BUNDLE=nrepel.lv2
BUILDDIR=build/
SRCDIR=src/
TTLDIR=lv2ttl/
targets=

UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
  LV2LDFLAGS=-dynamiclib
  LIB_EXT=.dylib
  EXTENDED_RE=-E
  STRIPFLAGS=-u -r -arch all -s lv2syms
  targets+=lv2syms
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
# extract versions
LV2VERSION=$(nrepel_VERSION)
include git2lv2.mk

# check for build-dependencies
ifeq ($(shell pkg-config --exists lv2 || echo no), no)
  $(error "LV2 SDK was not found")
endif
ifeq ($(shell pkg-config --exists fftw3f || echo no), no)
  $(error "FFTW was not found")
endif

override CFLAGS += -fPIC -std=c99
override CFLAGS += `pkg-config --cflags lv2`
override LOADLIBES += `pkg-config --cflags --libs fftw3f`

#for debug building
ifeq ($(DEBUG), 1)
  override CFLAGS += -g3 -DDEBUG
else
  override CFLAGS += -DNDEBUG
endif

# build target definitions
default: all

all: $(BUILDDIR)manifest.ttl $(BUILDDIR)$(LV2NAME).ttl $(targets)

lv2syms:
	echo "_lv2_descriptor" > lv2syms

$(BUILDDIR)manifest.ttl: $(TTLDIR)manifest.ttl.in
	@mkdir -p $(BUILDDIR)
	sed "s/@LV2NAME@/$(LV2NAME)/;s/@LIB_EXT@/$(LIB_EXT)/" \
	  $(TTLDIR)manifest.ttl.in > $(BUILDDIR)manifest.ttl

$(BUILDDIR)$(LV2NAME).ttl: $(TTLDIR)$(LV2NAME).ttl.in
	@mkdir -p $(BUILDDIR)
	sed "s/@VERSION@/lv2:microVersion $(LV2MIC) ;lv2:minorVersion $(LV2MIN) ;/g" \
		$(TTLDIR)$(LV2NAME).ttl.in > $(BUILDDIR)$(LV2NAME).ttl

$(BUILDDIR)$(LV2NAME)$(LIB_EXT): $(SRCDIR)$(LV2NAME).c
	@mkdir -p $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) \
		-o $(BUILDDIR)$(LV2NAME)$(LIB_EXT) $(SRCDIR)$(LV2NAME).c \
		-shared $(LV2LDFLAGS) $(LDFLAGS) $(LOADLIBES)

ifeq ($(DEBUG), 0)
	$(STRIP) $(STRIPFLAGS) $(BUILDDIR)$(LV2NAME)$(LIB_EXT)
endif

# install/uninstall/clean target definitions
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

.PHONY: clean all install uninstall
