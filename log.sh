
#!/bin/bash

#################################################
## Helen Kang
## 210913
## Script to process NextSeq runs of FT010_Fresh ATAC (adapted from 210611_FT007)


#################################################
## Directories
QSUB=$GROUP_HOME/bin/quick-sub
QSUBPAR=$GROUP_HOME/bin/quick-sub-long
QSUBPARBIGMEM=$GROUP_HOME/bin/quick-sub-bigmem-parallel
PROJECTNAME=tutorial
PROJECT=/oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/${PROJECTNAME}
DATADIR=/oak/stanford/groups/engreitz/Projects/SequencingRuns/
GEXDATADIR=${DATADIR}/210107_NB551514_0533_AH3MVFBGXH
ATACDATADIR=${DATADIR}/210108_NB551514_0534_AH3LGHBGXH
GEXFLOWCELL=$(grep "\<Flowcell\>" ${GEXDATADIR}/RunInfo.xml | sed -e "s/<Flowcell>//g" | sed -e "s/<\/Flowcell>//g") # doesn't fully work, need to remove the space in the front
ATACFLOWCELL=$(grep "\<Flowcell\>" ${ATACDATADIR}/RunInfo.xml | sed -e "s/<Flowcell>//g" | sed -e "s/<\/Flowcell>//g")
# FLOWCELL=HHJGLBGXJ
LOGS=${PROJECT}/logs/
mkdir -p ${LOGS}
SCRATCHDIR=$SCRATCH/${PROJECTNAME}
mkdir -p ${SCRATCHDIR}


#################################################
## make fastq files
# scRNA-seq
# index reference: https://support.10xgenomics.com/sample_index_sets/Dual_Index_Kit_TT_Set_A.csv
## need gex-samplesheet.csv
JOBNAME=gex_mkfastq # this run takes 10 minutes if there is no error
cd ${PROJECT}
$QSUB -j ${JOBNAME} -s ${LOGS}/${JOBNAME}.qsh -o ${LOGS}/${JOBNAME}.qout -t 12:00:00 -m 128G -p normal,owners,engreitz "cellranger-arc mkfastq --id=${PROJECTNAME}_gex \
    --run=${GEXDATADIR} \
    --csv=${PROJECT}/gex-samplesheet.csv"


# scATAC-seq
# index reference: https://support.10xgenomics.com/sample_index_sets/Single_Index_Kit_N_Set_A.csv
JOBNAME=atac_mkfastq
$QSUB -j ${JOBNAME} -s ${LOGS}/${JOBNAME}.qsh -o ${LOGS}/${JOBNAME}.qout -t 12:00:00 -m 64G  "cellranger-arc mkfastq --id=${PROJECTNAME}_atac \
    --run=${ATACDATADIR} \
    --csv=${PROJECT}/atac-samplesheet.csv"


#################################################
## get count matrix and BAM files
# gex only portion (when atac is not available)
cd ${PROJECT}
SAMPLE=(fresh frozen)
GEXFASTQDIR=$(echo "${PROJECT}/${PROJECTNAME}_gex/outs/fastq_path/${FLOWCELL}" | sed "s/\ //g" )
sampleLength=${#SAMPLE[@]}
# use cellranger 6.0.0
PATH=$(echo "$PATH" | sed -e 's/:\/home\/groups\/engreitz\/Software\/cellranger-4.0.0\/bin/g/') ## to remove cellranger 4.0.0

MEM_GB=196
local_mem=$((9 * ${MEM_GB} / 10))
NUM_THREADS=24

for (( i = 0; i < sampleLength; i++ ))
do
    id=${SAMPLE[$i]}
    $QSUBPAR -j ${id} -s ${LOGS}/${id}.qsh -o ${LOGS}/${id}.qout -m ${MEM_GB}GB -c ${NUM_THREADS} -t 24:00:00 -n 1 "source ~/.bashrc; cd ${PROJECT}; cellranger count \
	--fastqs=${GEXFASTQDIR} \
	--transcriptome=/oak/stanford/groups/engreitz/Users/kangh/2010_process_sequencing_data/CellRangerRef/refdata-gex-GRCh38-2020-A \
	--id=gex_${id} \
        --include-introns \
	--sample=${id} \
        --chemistry=ARC-v1 \
        --jobmode=local \
	--localcores=${NUM_THREADS} \
	--localmem=${local_mem} "
done


#################################################
## get count matrix and BAM files (ATAC only)
ATACFASTQDIR=/oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/210912_FT010_fresh_Telo_sortedEC/210912_FT010_fresh_Telo_sortedEC_atac/outs/fastq_path/HHK7TBGXJ
MEM_GB=192
local_mem=$((9 * ${MEM_GB} / 10))
NUM_THREADS=24

for sample in "${SAMPLE[@]}"
do
    JOBNAME=${sample}_count
    echo $JOBNAME
    $QSUBPAR -j ${JOBNAME} -s ${LOGS}/${JOBNAME}.qsh -o ${LOGS}/${JOBNAME}.qout -t 48:00:00 -m ${MEM_GB}GB -c ${NUM_THREADS} -n 1 "source ~/.bashrc; cd ${PROJECT};    cellranger-atac count \
		    --id=atac_${sample} \
		    --reference=$GROUP_HOME/Software/cellranger-arc-1.0.1/refdata-cellranger-arc-GRCh38-2020-A \
		    --fastqs=${ATACFASTQDIR}/${sample} \
		    --jobmode=local \
		    --localcores=${NUM_THREADS} \
		    --localmem=${local_mem}  "
done


#################################################
## get count matrix and BAM files (Multiome)
## need $sample_libraries.csv beforehand
MEM_GB=256
local_mem=$((9 * ${MEM_GB} / 10))
NUM_THREADS=30
PATH=$(echo $PATH | sed -e 's/:\/home\/groups\/engreitz\/Software\/cellranger-arc-1.0.1/g/')
for sample in "${SAMPLE[@]}"
do
    JOBNAME=${sample}_count
    $QSUBPARBIGMEM -j ${JOBNAME} -s ${LOGS}/${JOBNAME}.qsh -o ${LOGS}/${JOBNAME}.qout -t 24:00:00 -m ${MEM_GB}GB -c ${NUM_THREADS} -n 1 "source ~/.bashrc; cd ${PROJECT}; cellranger-arc --help; cellranger-arc count --id=multiome_${sample} \
    	--reference=$GROUP_HOME/Software/cellranger-arc-1.0.1/refdata-cellranger-arc-GRCh38-2020-A \
    	--libraries=${PROJECT}/${sample}_libraries.csv \
        --jobmode=local \
        --localcores=${NUM_THREADS} \
        --localmem=${local_mem} "
done


## if bigmem partition is too crowded, use engreitz or normal: (this is much slower. it usually takes ~8 hours to complete)
MEM_GB=196
local_mem=$((9 * ${MEM_GB} / 10))
NUM_THREADS=24
PATH=$(echo $PATH | sed -e 's/:\/home\/groups\/engreitz\/Software\/cellranger-arc-1.0.1/g/')
for sample in "${SAMPLE[@]}"
do
    JOBNAME=${sample}_count
    $QSUBPAR -j ${JOBNAME} -s ${LOGS}/${JOBNAME}.qsh -o ${LOGS}/${JOBNAME}.qout -t 48:00:00 -m ${MEM_GB}GB -c ${NUM_THREADS} -n 1 "source ~/.bashrc; cd ${PROJECT}; cellranger-arc --help; cellranger-arc count --id=multiome_${sample} \
    	--reference=$GROUP_HOME/Software/cellranger-arc-1.0.1/refdata-cellranger-arc-GRCh38-2020-A \
    	--libraries=${PROJECT}/${sample}_libraries.csv \
        --jobmode=local \
        --localcores=${NUM_THREADS} \
        --localmem=${local_mem} "
done



# #################################################
# ## CMO: get whitelisted barcodes
# CMOOUT=${PROJECT}/outputs/
# WHITELIST=${PROJECT}/outputs/CBC_whitelist.txt

# mkdir -p ${CMOOUT}
# Rscript /oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/210320_cell_hashing/get_CBC_whitelist.R --datadir ${PROJECT}/ --sampleName gex_CM_CMO --whitelist.out ${WHITELIST}

# # copying directories from combining FASTQ step
# # scRNA-seq
# GEXFASTQDIR=${PROJECT}/${PROJECTNAME}_gex/outs/fastq_path/$FLOWCELL
# # GEXFASTQOUT=${PROJECT}/gex_fastq_combined/
# # mkdir -p ${GEXFASTQOUT}

# ## extract CMO from R1 and R2
# CMOCITESEQOUT=${PROJECT}/CMO_CITE-seq-Count_outputs/
# SAMPLENAME=CM_CMO
# EXPECTED_CELLS=$(wc -l ${WHITELIST} | awk '{print $1}')
# mkdir -p ${CMOCITESEQOUT}

# ## 210520
# ## extract CMO from R1 and R2
# # CMO_SAMPLE=CM_CMO
# # GEXFASTQDIR=${PROJECT}/${PROJECTNAME}_gex/outs/fastq_path/$FLOWCELL
# # for (( i = 1; i < 5; i++ ))
# # do
# #     CMOCITESEQOUT=${PROJECT}/CMO_CITE-seq-Count_outputs/Lane_${i}/
# #     EXPECTED_CELLS=$(wc -l ${WHITELIST} | awk '{print $1}')
# #     mkdir -p ${CMOCITESEQOUT}
# #     JOBNAME=CMO_CITE-seq-Count_${i} ## TAKES 4+ HOURS!
# #     $QSUB -j ${JOBNAME} -s ${LOGS}/${JOBNAME}.qsh -o ${LOGS}/${JOBNAME}.qout -t 12:00:00 -m 96G  "source ~/.bashrc; conda activate cnmf_env; CITE-seq-Count -R1 ${GEXFASTQDIR}/${CMO_SAMPLE}_S3_L00${i}_R1_001.fastq.gz \
# #     -R2 ${GEXFASTQDIR}/${CMO_SAMPLE}_S3_L00${i}_R2_001.fastq.gz \
# #     -t ${CMOOUT}/TAG_LIST.csv \
# #     -cbf 1 -cbl 16 \
# #     -umif 17 -umil 28 \
# #     -cells ${EXPECTED_CELLS} \
# #     -wl ${WHITELIST} \
# #     -o ${CMOCITESEQOUT} " # debug mode will output >1G of text after 5 minutes of running
# # done





# # #################################################
# # ## make TDF files for IGV

# # # make tdf files from BAM
# # DATASET=(fresh frozen)
# # OME=(atac gex)
# # for sample in "${DATASET[@]}"
# # do
# #     for libtype in "${OME[@]}"
# #     do
# # 	BAMDIR=${PROJECT}/multiome_${sample}/outs/${libtype}_possorted_bam.bam
# # 	IGVDIR=${PROJECT}/IGV/
# # 	mkdir -p ${IGVDIR}
# # 	JOBNAME=toTDF.${sample}.${libtype}

# # 	$QSUB -j $JOBNAME -s ${LOG}/$JOBNAME.qsh -o ${LOG}/$JOBNAME.qout -t 4:00:00 -m 64G "source ~/.bashrc; conda activate EngreitzLab; \
# # bedtools genomecov -ibam ${BAMDIR} -bg > ${SCRATCHDIR}/${JOBNAME}.bedgraph; \
# # igvtools toTDF ${SCRATCHDIR}/${JOBNAME}.bedgraph ${IGVDIR}/${sample}_${libtype}_possorted_bam.bam.tdf /oak/stanford/projects/genomics-refs/refs/hg38/hg38.chrom.sizes"
# #     done
# # done
