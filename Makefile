# object variables
version := 0.1
config := config.txt
batch := $(shell grep BATCH $(config) | cut -f 2)
module := modulefiles/DeltaMP/$(version)
module_load := $(shell grep MODULES_TO_BE_LOADED $(config) | cut -f 2)
max_cpus := $(shell grep MAX_CPU_PER_NODES $(config) | cut -f 2)
max_mem := $(shell grep MAX_MEMORY_FOR_HIGH_MEM $(config) | cut -f 2)
deltamp := $(patsubst src/%,bin/%,$(patsubst %.main,%,$(wildcard $(addsuffix *.main,src/))))
batch_spec := $(addprefix bin/,$(notdir $(shell ls lib/$(batch)/deltamp.*)))
steps := $(patsubst src/%,bin/%,$(patsubst %.step,%.sh,$(wildcard $(addsuffix *.step,src/))))
arrays := $(patsubst %.head,%_array.head,$(shell ls lib/$(batch)/*.head | grep -v "_"))
highmems := $(patsubst %.head,%_highmem.head,$(shell ls lib/$(batch)/*.head | grep -v "_"))
test_config := $(patsubst src/%,test/%,$(patsubst %.config,%.tsv,$(wildcard $(addsuffix *.config,src/))))

# search paths
vpath %.main src
vpath %.step src
vpath %.head lib/$(batch)
vpath %.config src

# main rule
.PHONY: all clean
all: $(deltamp) $(module) $(steps) $(batch_spec) $(test_config)

# rule to build deltamp, pipeline_master, restart_from_step, delete_subproject and check_previous_step
$(deltamp): bin/% : %.main | lib/$(batch)/option_variables
	SED=$$(sed 's/^/s\//;s/\t/\//;s/$$/\/g/' $| |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SED $< | sed 's#^\(VERSION\[DELTAMP\]=\)$$#\1$(version)#;s#^\(DELTAMP_BUILD=\)$$#\1$(CURDIR)#' > $@ && chmod +x $@

# rule to build the module file of the current version
$(module): src/deltamp.module | modulefiles
	 awk -v M="$(module_load)" '{if($$1=="LOAD"){split(M,a,";");for(i in a){print "module\t\tload\t"a[i]}} else print $$0}' $< |\
	  sed 's#VERSION#$(version)#;s#PATH\tbin#PATH\t$(CURDIR)/bin#' > $@

# rule to build the modulefiles directory
modulefiles:
	mkdir -p $@ $@/DeltaMP

# rule to build step scripts
$(steps): bin/%.sh : | %.step
ifeq ($(batch),GridEngine)
	cat $< $| | sed 's/log\/NAME/log\/'$*'/;s/MAX_CPUS$$/$(max_cpus)/;s/NCPUS/NSLOTS/g;s/ARRAY_TASK/SGE_TASK_ID/g;s/QUEUE_JOBNAME/-N/;s/JOBNAME/JOB_NAME/;s/QUEUE_HOLD/-hold_jid /;s/QUEUE_SEP/,/' > $@
else ifeq ($(batch),Slurm)
	cat $< $| | sed 's/log\/NAME/log\/'$*'/;s/MAX_CPUS$$/$(max_cpus)/;s/NCPUS/SLURM_CPUS_PER_TASK/g;s/ARRAY_TASK/SLURM_ARRAY_TASK_ID/g;s/QUEUE_JOBNAME/-J/;s/JOBNAME/SLURM_JOB_NAME/;s/QUEUE_HOLD/-d afterok:/;s/QUEUE_SEP/:/' > $@
endif

# header (type of job) specific dependencies
bin/init.sh bin/get.sh bin/454_quality.sh bin/Illumina_quality.sh bin/doc.sh bin/archiver.sh : serial.head
bin/Illumina_demulti.sh bin/Illumina_fastq.sh bin/Illumina_pair_end.sh bin/Illumina_opt.sh bin/454_raw_stat.sh bin/Illumina_raw_stat.sh bin/trim.sh : serial_array.head
bin/merge.sh bin/end.sh : serial_highmem.head
bin/OTU.sh : mp.head
bin/454_demulti.sh bin/454_sff.sh bin/454_opt.sh : mp_array.head
bin/cut_db.sh bin/id.sh : mp_highmem.head

# rules to build array job headers
$(arrays): lib/$(batch)/%_array.head : %.head
ifeq ($(batch),GridEngine)
	sed 's/NAME/NAME.$$TASK_ID/' $< > $@
else ifeq ($(batch),Slurm)
	sed 's/NAME/NAME.%a/' $< > $@
endif

# rules to build high memory job headers
$(highmems): lib/$(batch)/%_highmem.head : %.head
ifeq ($(batch),GridEngine)
	sed 's/mem=6G/mem=$(max_mem)G/;$$s/^$$/#$$ -l highmem\n/' $< > $@
else ifeq ($(batch),Slurm)
	sed 's/mem=64G/mem=$(max_mem)G/' $< > $@
endif

# Copy batch queuing system specific executables to bin
$(batch_spec): bin/deltamp.% : lib/$(batch)/deltamp.%
	cp $^ $@ && chmod +x $@

# rule to build test configuration file
$(test_config): test/%.tsv : %.config
	sed "s#USER#$$USER#;s#CURDIR#$(CURDIR)#" $< > $@

# clean rule
clean :
	rm -r $(deltamp) modulefiles $(steps) $(arrays) $(highmems) $(batch_spec) $(test_config)
