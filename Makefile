# object variables
version := 0.3
SHELL := /bin/bash
batch := $(shell grep BATCH_QUEUE config.txt | cut -f 2)
module := modulefiles/DeltaMP/$(version)
module_load := $(shell grep MODULES_TO_BE_LOADED config.txt | cut -f 2)
deltamp := $(patsubst src/%,bin/%,$(patsubst %.main,%,$(wildcard $(addsuffix *.main,src/))))
batch_spec := $(addprefix bin/,$(notdir $(shell ls lib/$(batch)/deltamp.*)))
steps := $(patsubst src/%,bin/%,$(patsubst %.step,%.sh,$(wildcard $(addsuffix *.step,src/))))
arrays := $(patsubst %.options,%_array.options,$(shell ls lib/$(batch)/*.options | grep -v "_"))
highmems := $(patsubst %.options,%_highmem.options,$(shell ls lib/$(batch)/*.options | grep -v "_"))
options := $(arrays) $(highmems) $(shell ls lib/$(batch)/*.options | grep -v "_")
heads := $(patsubst %.options,%.head,$(options))
test_config := $(patsubst src/%,test/%,$(patsubst %.config,%.tsv,$(wildcard $(addsuffix *.config,src/))))

# search paths
vpath %.main src
vpath %.step src
vpath %.head lib/$(batch)
vpath %.options lib/$(batch)
vpath %.config src

# main rule
.PHONY: all clean
all: $(deltamp) $(module) $(steps) $(batch_spec) $(test_config)

# rule to build deltamp, pipeline_master, restart_from_step, delete_subproject and check_previous_step
$(deltamp): bin/% : %.main | lib/$(batch)/option_variables
	SEDMAIN=$$(sed 's/^/s\//;s/\t/\//;s/$$/\/g/' $| |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SEDMAIN $< | sed 's#^\(VERSION\[DELTAMP\]=\)$$#\1$(version)#;s#^\(DELTAMP_BUILD=\)$$#\1$(CURDIR)#' > $@ && chmod +x $@

# rule to build the module file of the current version
$(module): src/deltamp.module | modulefiles
	 awk -v M="$(module_load)" '{if($$1=="LOAD"){split(M,a,";");for(i in a){print "module\t\tload\t"a[i]}} else print $$0}' $< |\
	  sed 's#VERSION#$(version)#;s#PATH\tbin#PATH\t$(CURDIR)/bin#' > $@

# rule to build the modulefiles directory
modulefiles:
	mkdir -p $@ $@/DeltaMP

# rule to build step scripts
$(steps): bin/%.sh : | %.step
	SEDSTEP=$$(sed 's/^/s@/;s/\t/@/;s/$$/@g/' lib/$(batch)/option_variables |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SEDSTEP $| | sed 's#/DELTAMP_VERSION#/$(version)#' | cat $< - | sed 's/log\/NAME/log\/'$*'/' > $@

# header (type of job) specific dependencies
bin/init.sh bin/454_quality.sh bin/Illumina_quality.sh bin/doc.sh bin/archiver.sh : serial.head
bin/Illumina_demulti.sh bin/Illumina_fastq.sh bin/Illumina_pair_end.sh bin/Illumina_opt.sh bin/454_raw_stat.sh bin/Illumina_raw_stat.sh : serial_array.head
bin/end.sh : serial_highmem.head
bin/get.sh bin/OTU.sh : mp.head
bin/454_demulti.sh bin/454_sff.sh bin/454_opt.sh bin/trim.sh : mp_array.head
bin/cut_db.sh bin/id.sh : mp_highmem.head

# rule to set parameters of job headers
$(heads): lib/$(batch)/%.head :  %.options
	SEDHEAD=$$(sed '2d;s/^/s@/;s/\t/@/;s/$$/@g/' config.txt |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SEDHEAD $< > $@

# rules to build array job headers
$(arrays): lib/$(batch)/%_array.options : %.options
	sed 's/MAX_MEMORY/MAX_ARRAY_MEMORY/;s/MAX_CPU/MAX_ARRAY_CPU/' $< > $@
ifeq ($(batch),GridEngine)
	sed 's/NAME/NAME.$$TASK_ID/' $@ > $@.temp && mv $@.temp $@
else ifeq ($(batch),Slurm)
	sed 's/NAME/NAME.%a/' $@ | cat - <(echo 'sleep $$SLURM_ARRAY_TASK_ID') > $@.temp && mv $@.temp $@
endif

# rules to build high memory job headers
$(highmems): lib/$(batch)/%_highmem.options : %.options
	sed 's/MAX_MEMORY/MAX_HIGH_MEM_MEMORY/;s/MAX_CPU/MAX_HIGH_MEM_CPU/' $< > $@
ifeq ($(batch),GridEngine)
	sed '$$s/^$$/#$$ -l highmem\n/' $@ > $@.temp && mv $@.temp $@
endif

# Copy batch queuing system specific executables to bin
$(batch_spec): bin/deltamp.% : lib/$(batch)/deltamp.%
	cp $^ $@ && chmod +x $@

# rule to build test configuration file
$(test_config): test/%.tsv : %.config
	sed "s#USER#$$USER#;s#CURDIR#$(CURDIR)#" $< > $@

# clean rule
clean :
	rm -r $(deltamp) modulefiles $(steps) $(heads) $(arrays) $(highmems) $(batch_spec) $(test_config)
