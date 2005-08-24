VERSION=0.77.1-20050824
THISDIR=branch_0
HOME=.

.PHONY: choosetarget pentium3 pentium4 athlon x86_64 x86_32 ppc doc \
        x86common common clean snapshot test

choosetarget :
	@echo "Possible targets are:"
	@echo "	pentium3              	 Intel Pentium 3"
	@echo "	pentium4              	 Intel Pentium 4"
	@echo "	pentium-m              	 Intel Pentium M"
	@echo "	athlon			 AMD Athlon (k7)"
	@echo "	x86_64                   AMD Opteron/Athlon64 (k8)"
	@echo "	ppc_970             	 PowerPC 970"
	@echo "	ppc_7450             	 PowerPC 7450"
	@echo "	doc			 Documentation"
	@echo "	snapshot       		 Sources snapshot"
	@echo "	install			 Installation"
	@echo "	clean         		 Clean up"

pentium3 :
	@ARCH="pentium3" $(MAKE) x86common

pentium4 :
	@ARCH="pentium4" $(MAKE) x86common

pentium-m :
	@ARCH="pentium-m" $(MAKE) x86common
	
athlon :
	@ARCH="athlon" $(MAKE) x86common

x86_64 :
	ARCH="k8" $(MAKE) common

x86_32 :
	ARCH="athlon" $(MAKE) common

ppc_970 :
	ARCH="970" $(MAKE) common
 
ppc_7450 :
	ARCH="7450" $(MAKE) common

doc :
	$(MAKE) -C doc/ggnfs-doc

EXCLUDE_FILES=*.tar *.fb spairs.out* core msdebug.txt maplesyntax rels.bin* \
              cols* lpindex* *.pdf rels.bin.* lpindex* *.afb.0 relations.jpg
EXCLUDE=$(foreach opt,$(EXCLUDE_FILES),--exclude=$(opt))

x86common :
	echo "#define GGNFS_VERSION \"$(VERSION)-$(ARCH)\"" > include/version.h
	@cd src/lasieve4 && rm -f -r asm && ln -s piii asm
	@ARCH=$(ARCH) $(MAKE) -C src -f Makefile.x86

common :
	echo "#define GGNFS_VERSION \"$(VERSION)-$(ARCH)\"" > include/version.h
	@cd src/lasieve4 && rm -f -r asm && ln -s ppc32 asm
	@ARCH=$(ARCH) $(MAKE) -C src

clean :
	-rm -f include/version.h
	$(MAKE) -C src -f Makefile.x86 clean
	$(MAKE) -C src clean
	-rm -f -r asm
	-rm -rf core ggnfs-*.tar.gz
	cd ./tests; sh cleanup.sh

snapshot : clean
	echo $(VERSION) > Version; date -u >> Version
	tar cvf ggnfs-$(VERSION).tar -C ../ $(EXCLUDE) ./$(THISDIR)/
	gzip ggnfs-$(VERSION).tar

test :
	echo $(EXCLUDE)

