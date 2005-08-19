VERSION=0.77.1-050819
THISDIR=branch_0
HOME=.

choosetarget :
	@echo "Possible targets are:"
	@echo "	pentium3              	 Intel Pentium 3"
	@echo "	pentium4              	 Intel Pentium 4"
	@echo "	athlon			 AMD Athlon (k7)"
	@echo "	x86_64                   AMD Opteron/Athlon64 (k8)"
	@echo "	ppc                	 PowerPC"
	@echo "	doc			 Documentation"
	@echo "	snapshot       		 Sources snapshot"
	@echo "	install			 Installation"
	@echo "	clean         		 Clean up"

pentium3 :
	@ARCH="pentium3" $(MAKE) x86common

pentium4 :
	@ARCH="pentium4" $(MAKE) x86common
	
athlon :
	@ARCH="athlon" $(MAKE) x86common

x86_64 :
	@TARGET="x86_64" ARCH="k8" $(MAKE) common
	
x86_32 :
	@TARGET="x86_64" ARCH="athlon" $(MAKE) common
ppc :
	@TARGET="ppc" ARCH="970" $(MAKE) common

doc :
	@cd doc/ggnfs-doc && $(MAKE)

EXCLUDE_FILES=*.tar *.fb spairs.out* core msdebug.txt maplesyntax rels.bin* \
              cols* lpindex* *.pdf rels.bin.* lpindex* *.afb.0 relations.jpg
EXCLUDE=$(foreach opt,$(EXCLUDE_FILES),--exclude=$(opt))

x86common : src ;
	echo "#define GGNFS_VERSION \"$(VERSION)-$(ARCH)\"" > include/version.h
	@cd src/lasieve4 && rm -f -r asm && ln -s piii asm
	@ARCH=$(ARCH) $(MAKE) -C src -f Makefile.x86

common : src ;
	echo "#define GGNFS_VERSION \"$(VERSION)-$(ARCH)\"" > include/version.h
	@ARCH=$(ARCH) $(MAKE) -C src -f Makefile.$(TARGET)

clean : ;
	$(MAKE) -C src -f Makefile.x86 clean
	$(MAKE) -C src -f Makefile.x86_64 clean
	$(MAKE) -C src -f Makefile.ppc clean
	rm -rf core ggnfs-*.tar.gz
	cd ./tests; sh cleanup.sh

snapshot : clean ;
	echo "#define GGNFS_VERSION \"$(VERSION)\"" > include/version.h
	echo $(VERSION) > Version; date -u >> Version
	tar cvf ggnfs-$(VERSION).tar -C ../ $(EXCLUDE) ./$(THISDIR)/ 
	gzip ggnfs-$(VERSION).tar

test : ;
	echo $(EXCLUDE)

