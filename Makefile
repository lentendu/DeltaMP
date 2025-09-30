# object variables
version := 0.6
SHELL := /bin/bash
commit := $(shell if [ -d .git ] ; then git log -1 --pretty=format:%h ; else echo $(version) ; fi)
batch := $(shell grep BATCH_QUEUE config.txt | cut -f 2)
modulefiles := $(shell grep MODULEFILES_DIR config.txt | cut -f 2)
module := $(shell grep MODULE_NAME config.txt | cut -f 2)
module_load := $(shell grep MODULES_TO_BE_LOADED config.txt | cut -f 2)
deltamp := $(patsubst src/%,bin/%,$(patsubst %.main,%,$(wildcard $(addsuffix *.main,src/))))
batch_spec := $(addprefix bin/,$(notdir $(shell ls lib/$(batch)/deltamp.*)))
steps := $(patsubst src/%,bin/%,$(patsubst %.step,%.sh,$(wildcard $(addsuffix *.step,src/))))
standard := $(patsubst %.options,%-std.options,$(shell ls lib/$(batch)/*.options | grep -v "_\|-"))
stdval := $(shell grep QUEUE_PARTITION_STANDARD config.txt | cut -f 2)
highmems := $(patsubst %-std.options,%-highmem.options,$(standard))
hmval := $(shell grep QUEUE_PARTITION_HIGHMEM config.txt | cut -f 2)
arrays := $(patsubst %.options,%_array.options,$(standard))
large_arrays := $(patsubst %.options,%_array_large.options,lib/$(batch)/mp-std.options)
options := $(arrays) $(large_arrays) $(standard) $(highmems)
heads := $(patsubst %.options,%.head,$(options))
test_config := $(patsubst src/%,test/%,$(patsubst %.config,%.tsv,$(wildcard $(addsuffix *.config,src/))))
tempdir := $(shell grep TEMPD config.txt | cut -f 2)

# search paths
vpath %.main src
vpath %.step src
vpath %.head lib/$(batch)
vpath %.options lib/$(batch)
vpath %.config src

# main rule
.PHONY: all clean
all: $(deltamp) dep $(module) $(steps) $(batch_spec) $(test_config)

# rule to build deltamp, pipeline_master, restart_from_step, delete_subproject and check_previous_step
$(deltamp): bin/% : %.main | lib/$(batch)/option_variables
	SEDMAIN=$$(sed 's/^/s\//;s/\t/\//;s/$$/\/g/' $| |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SEDMAIN $< | sed 's#^\(VERSION\[DELTAMP\]=\)$$#\1$(version)#;s#^\(DELTAMP_BUILD=\)$$#\1$(CURDIR)#;s#^\(DELTAMP_COMMIT=\)$$#\1$(commit)#' > $@ && chmod +x $@

# rule to check dependencies once depending modules loaded
dep:
	@echo control dependencies
	(module load $$(echo "$(module_load)" | sed 's/;/ /g') ;\
	while read soft ver; do unset OUT; OUT=$$($$soft --version 2>&1 | sed -n "s/[Vv]ersion[:= ]*/v/;s/^$$soft \([0-9][0-9\.-]*\)/$$soft v\1/;/v[0-9\.-][0-9\.-]*/{s/.*\(v[0-9\.-][0-9\.-]*\).*/$$soft \1/p;q}"); if [ -z "$$OUT" ]; then OUT=$$($$soft -h 2>&1 | sed -n "s/[Vv]ersion[:= ]*/v/;s/^$$soft \([0-9][0-9\.-]*\)/$$soft v\1/;/v[0-9\.-][0-9\.-]*/{s/.*\(v[0-9\.-][0-9\.-]*\).*/$$soft \1/p;q}"); fi; if [ -z "$$OUT" ]; then echo $$soft not_found; else echo $$OUT | sed 's/ v/ /;s/\([0-9][0-9]*\.*[0-9]*\).*/\1/'; fi | awk -v v="$$ver" '{if($$2=="not_found"){print} else {if($$2<v){print $$0,"below requested "v} else {print $$0,"ok"}}}'; done < dependencies | column -t | awk '{print; if($$0~"not_found" || $$0~"below"){err=1}}END{if(err==1){print "\ndependency issue" ; exit err}}')
	
# rule to build the module file of the current version
$(module): src/deltamp.module | $(modulefiles)
	 awk -v M="$(module_load)" '{if($$1=="LOAD"){split(M,a,";");for(i in a){print "module\t\tload\t"a[i]}} else print $$0}' $< |\
	  sed 's#VERSION#$(version)#;s#PATH\tbin#PATH\t$(CURDIR)/bin#' > $|/$@

# rule to build the modulefiles directory
$(modulefiles):
	mkdir -p $@ $@/$(dir $(module))

# rule to build step scripts
$(steps): bin/%.sh : | %.step
	SEDSTEP=$$(sed 's/^/s@/;s/\t/@/;s/$$/@g/' lib/$(batch)/option_variables |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SEDSTEP $| | sed 's#DeltaMP/DELTAMP_VERSION#$(module)#;s#$$TEMPD#$(tempdir)#g' | cat $< - | sed 's/log\/NAME/log\/'$*'/' > $@

# header (type of job) specific dependencies
bin/init.sh bin/454_quality.sh bin/Illumina_quality.sh bin/Nanopore_quality.sh bin/Pacbio_quality.sh bin/doc.sh bin/archiver.sh : serial-std.head
bin/Illumina_fastq.sh bin/Illumina_pair_end.sh bin/Illumina_opt.sh bin/Nanopore_opt.sh bin/Pacbio_opt.sh : serial-std_array.head
bin/end.sh : serial-highmem.head
bin/get.sh bin/OTU.sh : mp-std.head
bin/Illumina_demulti.sh bin/454_demulti.sh bin/454_sff.sh bin/454_opt.sh bin/Nanopore_fastq.sh bin/Pacbio_fastq.sh bin/trim.sh : mp-std_array.head
bin/cut_db.sh bin/id.sh : mp-highmem.head
bin/asv.sh : mp-std_array_large.head

# rule to set parameters of job headers
$(heads): lib/$(batch)/%.head : %.options
	SEDHEAD=$$(cat config.txt <(grep "^QUEUE_P" lib/$(batch)/option_variables) | sed '2d;s/^/s@/;s/\t/@/;s/$$/@g/' |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SEDHEAD $< > $@

# rules to build array job headers
$(arrays): lib/$(batch)/%_array.options : %.options
	sed 's/MAX_MEMORY/MAX_ARRAY_MEMORY/;s/MAX_CPU/MAX_ARRAY_CPU/' $< > $@
ifeq ($(batch),GridEngine)
	sed 's/NAME/NAME.$$TASK_ID/' $@ > $@.temp && mv $@.temp $@
else ifeq ($(batch),Slurm)
	sed 's/NAME/NAME.%a/' $@ | cat - <(echo 'sleep $$(( $$SLURM_ARRAY_TASK_ID / 2 ))') > $@.temp && mv $@.temp $@
endif

# rules to build large array job headers
$(large_arrays): lib/$(batch)/%_array_large.options : %.options
ifeq ($(batch),GridEngine)
	sed 's/NAME/NAME.$$TASK_ID/' $< > $@
else ifeq ($(batch),Slurm)
	sed 's/NAME/NAME.%a/' $< > $@
endif

# rules to build high memory job headers
$(highmems): lib/$(batch)/%-highmem.options : %.options
	sed 's/MAX_MEMORY/MAX_HIGH_MEM_MEMORY/;s/MAX_CPU/MAX_HIGH_MEM_CPU/' $< > $@
ifeq ($(batch),GridEngine)
	sed '$$s/^$$/#$$ -l highmem\n/' $@ > $@.temp && mv $@.temp $@
endif
ifneq ($(strip $(hmval)),)
	if [[ "$(<F)" == *"mp"* ]] ; then \
		echo QUEUE_PREFIX QUEUE_PARTITION $(strip $(hmval)) >> $@ ; \
	elif [ "$(stdval)" != "" ] ; then \
		echo QUEUE_PREFIX QUEUE_PARTITION $(strip $(stdval)) | cat $< - > $@ ; \
	fi
endif

# rules to build standard job headers
$(standard): lib/$(batch)/%-std.options : %.options
ifneq ($(strip $(stdval)),)
	echo QUEUE_PREFIX QUEUE_PARTITION $(strip $(stdval)) | cat $< - > $@
else
	cp $< $@
endif

# Copy batch queuing system specific executables to bin
$(batch_spec): bin/deltamp.% : lib/$(batch)/deltamp.%
	cp $^ $@ && chmod +x $@

# rule to build test configuration file
$(test_config): test/%.tsv : %.config
	sed "s#USER#$$USER#;s#CURDIR#$(CURDIR)#" $< > $@

# clean rule
clean :
	rm -rf $(deltamp) $(steps) $(heads) $(arrays) $(large_arrays) $(standard) $(highmems) $(batch_spec) $(test_config)
