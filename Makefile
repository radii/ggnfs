VERSION=0.77.1-050720
THISDIR=ggnfs
HOME=.

EXCLUDE_FILES=*.tar *.fb spairs.out* core msdebug.txt maplesyntax rels.bin* \
              cols* lpindex* *.pdf rels.bin.* lpindex* *.afb.0 relations.jpg
EXCLUDE=$(foreach opt,$(EXCLUDE_FILES),--exclude $(opt))

all : src ;
	echo "#define GGNFS_VERSION \"$(VERSION)\"" > include/version.h
	$(MAKE) -C src

clean : ;
	$(MAKE) -C src clean
	rm -f core ggnfs-*.tar.gz
	cd ./tests; ./cleanup.sh

backup : clean ;
	echo "#define GGNFS_VERSION \"$(VERSION)\"" > include/version.h
	echo $(VERSION) > Version; date >> Version
	tar cvf ggnfs-$(VERSION).tar -C ../ ./$(THISDIR)/ $(EXCLUDE)
	gzip ggnfs-$(VERSION).tar



test : ;
	echo $(EXCLUDE)

