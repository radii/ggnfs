VERSION=0.77.1-050812
THISDIR=branch_0
HOME=.

TARGET=x86
#TARGET=ppc
#TARGET=x86_64

EXCLUDE_FILES=*.tar *.fb spairs.out* core msdebug.txt maplesyntax rels.bin* \
              cols* lpindex* *.pdf rels.bin.* lpindex* *.afb.0 relations.jpg
EXCLUDE=$(foreach opt,$(EXCLUDE_FILES),--exclude=$(opt))

all : src ;
	echo "#define GGNFS_VERSION \"$(VERSION)-$(TARGET)\"" > include/version.h
	$(MAKE) -C src -f Makefile.$(TARGET)

clean : ;
	$(MAKE) -C src -f Makefile.$(TARGET) clean
	rm -f core ggnfs-*.tar.gz
	cd ./tests; sh cleanup.sh

backup : clean ;
	echo "#define GGNFS_VERSION \"$(VERSION)\"" > include/version.h
	echo $(VERSION) > Version; date >> Version
	tar cvf ggnfs-$(VERSION).tar -C ../ $(EXCLUDE) ./$(THISDIR)/ 
	gzip ggnfs-$(VERSION).tar



test : ;
	echo $(EXCLUDE)

