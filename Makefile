# object variables
config := config.txt
batch := $(shell cat ${config})
batch_spec := $(addprefix bin/,$(notdir $(shell ls lib/$(batch)/deltamp*)))
steps := $(patsubst src/%,bin/%,$(patsubst %.step,%.sh,$(wildcard $(addsuffix *.step,src/))))
highmems := $(patsubst %.head,%_highmem.head,$(shell ls lib/$(batch)/*.head | grep -v "_"))

# search paths
vpath %.step src
vpath %.head lib/$(batch)

# main rule
all: $(steps) $(batch_spec)

# rule to build step scripts
$(steps): bin/%.sh : 
ifeq ($(batch),GridEngine)
	cat $^ | sed 's/log\/NAME/log\/'$*'/;s/NCPUS/NSLOTS/g;s/ARRAY_TASK/SGE_TASK_ID/g' > $@
else ifeq ($(batch),Slurm)
	cat $^ | sed 's/log\/NAME/log\/'$*'/;s/NCPUS/SLURM_CPUS_PER_TASK/g;s/ARRAY_TASK/SLURM_ARRAY_TASK_ID/g' > $@
endif

# header (type of job) specific dependencies
bin/init.sh bin/get.sh bin/454_quality.sh bin/Illumina_quality.sh bin/doc.sh bin/454_raw_stat.sh bin/Illumina_fastq.sh bin/Illumina_pair_end.sh bin/Illumina_raw_stat.sh bin/Illumina_opt.sh bin/trim.sh : serial.head
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
$(batch_spec): bin/deltamp% : lib/$(batch)/deltamp%
	cp $^ $@ && chmod +x $@

# clean rule
clean :
	rm $(steps) $(highmems) $(batch_spec)
