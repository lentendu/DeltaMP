# object variables
version := 0.1
config := config.txt
batch := $(shell grep BATCH $(config) | cut -f 2)
module := modulefiles/DeltaMP/$(version)
module_load := $(shell grep MODULES_TO_BE_LOADED $(config) | cut -f 2)
deltamp := $(patsubst src/%,bin/%,$(patsubst %.main,%,$(wildcard $(addsuffix *.main,src/))))
batch_spec := $(addprefix bin/,$(notdir $(shell ls lib/$(batch)/deltamp.*)))
steps := $(patsubst src/%,bin/%,$(patsubst %.step,%.sh,$(wildcard $(addsuffix *.step,src/))))
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

# rule to build deltamp and pipeline_master
$(deltamp): bin/% : %.main | lib/$(batch)/option_variables
	SED=$$(sed 's/^/s\//;s/\t/\//;s/$$/\//' $| |  tr "\n" ";" | sed 's/^/sed "/;s/;$$/"/'); \
	eval $$SED $< | sed 's#^\(VERSION\[DELTAMP\]=\)$$#\1$(version)#;s#^\(AMP_BUILD=\)$$#\1$(CURDIR)#' > $@ && chmod +x $@

# rule to build the module file of the current version
$(module): src/deltamp.module | modulefiles
	 awk -v M="$(module_load)" '{if($$1=="LOAD"){split(M,a,";");for(i in a){print "module\t\tload\t"a[i]}} else print $$0}' $< |\
	  sed 's#VERSION#$(version)#;s#PATH\tbin#PATH\t$(CURDIR)/bin#' > $@

# rule to build the modulefiles directory
modulefiles:
	mkdir -p $@ $@/DeltaMP

# rule to build step scripts
$(steps): bin/%.sh : 
ifeq ($(batch),GridEngine)
	cat $^ | sed 's/log\/NAME/log\/'$*'/;s/NCPUS/NSLOTS/g;s/ARRAY_TASK/SGE_TASK_ID/g;s/QUEUE_JOBNAME/-N/;s/JOBNAME/JOB_NAME/;s/QUEUE_HOLD/-hold_jid /;s/QUEUE_SEP/,/' > $@
else ifeq ($(batch),Slurm)
	cat $^ | sed 's/log\/NAME/log\/'$*'/;s/NCPUS/SLURM_CPUS_PER_TASK/g;s/ARRAY_TASK/SLURM_ARRAY_TASK_ID/g;s/QUEUE_JOBNAME/-J/;s/JOBNAME/SLURM_JOB_NAME/;s/QUEUE_HOLD/-d afterok:/;s/QUEUE_SEP/:/' > $@
endif

# header (type of job) specific dependencies
bin/init.sh bin/get.sh bin/454_quality.sh bin/Illumina_quality.sh bin/doc.sh bin/454_raw_stat.sh bin/Illumina_fastq.sh bin/Illumina_pair_end.sh bin/Illumina_raw_stat.sh bin/Illumina_opt.sh bin/trim.sh bin/archiver.sh : serial.head
bin/merge.sh bin/end.sh : serial_highmem.head
bin/OTU.sh bin/Illumina_demulti.sh bin/454_demulti.sh bin/454_sff.sh bin/454_opt.sh : mp.head
bin/cut_db.sh bin/id.sh : mp_highmem.head

# general dependencies
$(steps): bin/%.sh : %.step

# rules to build high memory job headers
$(highmems): lib/$(batch)/%_highmem.head : %.head
ifeq ($(batch),GridEngine)
	sed 's/mem=6G/mem=20G/;$$s/^$$/#$$ -l highmem\n/' $< > $@
else ifeq ($(batch),Slurm)
	sed 's/mem=6G/mem=20G/' $< > $@
endif

# Copy batch queueing system specific executables to bin
$(batch_spec): bin/deltamp.% : lib/$(batch)/deltamp.%
	cp $^ $@ && chmod +x $@

# rule to build test configuration file
$(test_config): test/%.tsv : %.config
	sed "s#USER#$$USER#;s#CURDIR#$(CURDIR)#" $< > $@

# clean rule
clean :
	rm -r $(deltamp) modulefiles $(steps) $(highmems) $(batch_spec) $(test_config)
