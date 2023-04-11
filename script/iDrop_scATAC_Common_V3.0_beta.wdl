version 1.0
workflow iDrop_scATAC_common_2 {
    input{
        Array[Array[String]] Data
        String Outdir
        String SampleID
        String ProjectID
        String refCode
        String? readStructure
        Int ?ForceFrag
    }
    String outdir=Outdir
    String ID=SampleID
    String runID= ID
    String PID = "sentieon.q -P "+ProjectID+"_sentieon"

    String Rscript="/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/mawen/src/tools/Rscript"
    String root="/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/mawen/pipeline/01.scATAC/01.scATAC_v3.0"
    String? defaultConfig="/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/songyumo/pipeline/scATAC/script/C4scATAClib_seqT1_R1_70_R2_50.json"
    String config=select_first([readStructure,defaultConfig])
    String lib="/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/songyumo/pipeline/scATAC/script/lib.sh"
    String whitelist="/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/songyumo/test/4_scATAC_chromap/whitelist.txt"
    Map[String,Array[Int]] mapResource = {
      "0" : [4,20],
      "1" : [4,180],
      "2" : [20,430]
    }


    call makedir {
        input:
        Dir=outdir,
    }
    String refdir=makedir.refConfig[refCode]["bwaIndex"]
    String ref_index=makedir.refConfig[refCode]["chromapIndex"]
    String tss=makedir.refConfig[refCode]["TSS"]
    String promo=makedir.refConfig[refCode]["promoter"]
    String blacklis=makedir.refConfig[refCode]["blacklist"]
    String chromeSize=makedir.refConfig[refCode]["chromeSize"]
    String chrmt = makedir.refConfig[refCode]["chrM"]
    String species=makedir.refConfig[refCode]["refName"]
    String macs2sp=makedir.refConfig[refCode]["genomesize"]
    String model=makedir.refConfig[refCode]["model"]
    String bigG = makedir.refConfig[refCode]["bigGenome"]
    String d2cpara = if model == "1" then "-r "+species else if blacklis != "None" then "--bg "+chromeSize+" --ts "+tss+" --bl "+blacklis else "--bg "+chromeSize+" --ts "+tss

    scatter(fq in Data){
        String fq1=fq[0]
        String fq2=fq[1]
    }


    call mapping {
    input:
        fastq1=fq1,
        fastq2=fq2,
        outdir=makedir.Outdir,
        root=root,
        ref_index=ref_index,
        ref=refdir,
        whitelist=whitelist,
        cpp=mapResource[bigG][0],
        mem=mapResource[bigG][1]
  }

  call deconvolution { 
	input: 
        lib=lib,
        alnbed=mapping.alnbed,
        outdir=outdir,
        root=root,	
        ID=ID,
        ForceFrag=ForceFrag, 
        chrmt=chrmt,
        d2cpara=d2cpara
  }   
  call qc {
	input:
        lib=lib,
        root=root,
        tss=tss,
        FragmentFile=deconvolution.FragmentFile,
        outdir=outdir,
        ID=ID,
        Rscript=Rscript
  }

  call peakcount {
	input:
        lib=lib,
        root=root,
        listtxt=deconvolution.listtxt,
        macs2sp= macs2sp,
        promo= promo,
        outdir=outdir,
        FragmentFile=deconvolution.FragmentFile,
        ID=ID
  }

  call report { 
	input: 
        ID=ID,
        fastq1=fq1,
        fastq2=fq2, 
        root=root,
        refdir=refdir, 
        outdir=outdir, 
        Qcfie = qc.Qcfie,
        Qcfie2 = peakcount.Qcfie2 
        }
}


###### script  ######
 
task makedir {
    input{
        String Dir
    }
    command {
        mkdir -p ${Dir}
        mkdir -p ${Dir}/d2cfile
        mkdir -p ${Dir}/out
        mkdir -p ${Dir}/out/Peak
        mkdir -p ${Dir}/out/Promoter
        mkdir -p ${Dir}/temp
        mkdir -p ${Dir}/report/div
        mkdir -p ${Dir}/report/base64
        mkdir -p ${Dir}/report/table
    }
    output {
        String Outdir="${Dir}"
        Map[String,Map[String,String]] refConfig = read_json("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/mawen/pipeline/01.scATAC/01.scATAC_v3.0/script/idrop_scatac_ref.json")
    }
}



task  mapping {
  input {
    Array[String] fastq1
    Array[String] fastq2
    String outdir
    String ref_index
    String ref
    String root
    String ?lib
    String whitelist
    Int cpp
    Int mem
    }
  command {
    export PATH=/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/bin:$PATH
    export LD_LIBRARY_PATH="/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib:/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib64:$LD_LIBRARY_PATH"
    ${root}/bin/chromap --preset atac --bc-error-threshold 0 --trim-adapters -x ${ref_index} -r ${ref} -1 ${sep=',' fastq1} -2 ${sep=',' fastq2} -o ${outdir}/out/aln.bed --barcode ${sep=',' fastq1} --barcode-whitelist ${whitelist} --read-format bc:0:19,r1:20:-1 -t ${cpp} 2> ${outdir}/out/alignment_report.tsv
	
  }
  runtime{
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String alnbed="${outdir}/out/aln.bed"
  }
}


task deconvolution {
    input {
        String alnbed
        String outdir
        String root
        String ?lib
        String d2cpara
        String chrmt
        String ID
        Int ?ForceFrag
        Int ?mapq
        Int cpp=1
        Int mem=7
    }
  
  command <<<
	
    export LD_LIBRARY_PATH="/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/songyumo/pipeline/scATAC/software/v1.3.7/gcclib/lib:/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/songyumo/pipeline/scATAC/software/v1.3.7/gcclib/lib64:$LD_LIBRARY_PATH" && export PATH="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/Python-3/bin:$PATH"
    ${root}/bin/d2c_v1.4.1/bin/d2c merge -i ${alnbed} --mapq 30 --bf ${default=0 ForceFrag} -o ${outdir}/d2cfile -c 10 -n ${ID} --mc ${chrmt} ${d2cpara} --sat --bt1 CB
    cp ${outdir}/d2cfile/${ID}.CorrelationBarcodes.tsv.gz ${outdir}/report/plot_input2_Jaccard_Overlap_Knee.csv.gz
    cp ${outdir}/d2cfile/${ID}.barcodeCount.tsv ${outdir}/report/plot_input1_Bead_Barcode_Knee.csv
    tail -n +2 ${outdir}/d2cfile/${ID}.Metadata.tsv |awk '{print $1}' > ${outdir}/d2cfile/list.${ID}.txt
	
  >>>
  runtime {
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String FragmentFile="${outdir}/d2cfile/${ID}.fragments.tsv.gz"
    String listtxt= "${outdir}/d2cfile/list.${ID}.txt"
    File FragmentFile1="${outdir}/d2cfile/${ID}.fragments.tsv.gz"
    File listtxt1= "${outdir}/d2cfile/list.${ID}.txt"
  }
}

task qc {
    input {
        String outdir
        String Rscript
        String root
        String tss
        String FragmentFile
        String ID
        String ?lib
        Int cpp=1
        Int mem=1
    }

  command <<<
        if [ -f ${default=abjdbashj lib} ]; then
            source ${lib}
        fi
        ${Rscript} ${root}/script/plot_TSSEnrichment_FragSize.R -T ${tss} -F ${FragmentFile} -G ${ID} -O ${outdir}
        
        nf=$(gunzip -c ${outdir}/d2cfile/${ID}.fragments.tsv.gz | awk '{if($3-$2<147) print $0}' | wc -l)
        mn=$(gunzip -c ${outdir}/d2cfile/${ID}.fragments.tsv.gz | awk '{if($3-$2>147 && $3-$2<294) print $0}' | wc -l)
        total=$(gunzip -c ${outdir}/d2cfile/${ID}.fragments.tsv.gz | wc -l)
        echo "$nf $total" | awk '{print "Fraction of nucleosome-free regions:"$1/$2*100"%"}' >> ${outdir}/report/5_2.library.QC.csv
        echo "$mn $total" | awk '{print "Fraction of fragments mono-nucleosome regions:"$1/$2*100"%"}' >> ${outdir}/report/5_2.library.QC.csv

  >>>
  runtime {
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String Qcfie = "${outdir}/report/5_2.library.QC.csv"
    File Qcfie1 = "${outdir}/report/5_2.library.QC.csv"
    
  }
}


task peakcount {
    input{
        String root
        String outdir
        String ?lib
        String macs2sp
        String promo
        String ID
        String FragmentFile
        String listtxt
        Int cpp=1
        Int mem=14
    }

  command <<<
    if [ -f ${default=abjdbashj lib} ]; then
        source ${lib}
    fi
  
    ${root}/bin/macs2 callpeak -t ${FragmentFile} -f BED -g ${macs2sp} -n ${ID} -B -q 0.001 --nomodel --outdir ${outdir}/out

    ${root}/bin/Rscript ${root}/script/C4scATAC_Cluster_AnnotationV3.R -I ${outdir}/out/${ID}_peaks.narrowPeak -F ${FragmentFile} -C ${listtxt} -G ${promo} -Q ${outdir}/d2cfile/${ID}.Metadata.tsv -O ${outdir}
       

    
  >>>
  runtime {
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String Qcfie2 = "${outdir}/out/${ID}_peaks.narrowPeak"
    File Qcfie21 = "${outdir}/out/${ID}_peaks.narrowPeak"
  }
}



task report {
    input{
        String ID
        String ?lib
        Array[String] fastq1
        Array[String] fastq2
        String outdir
        String Qcfie
        String Qcfie2
        String refdir
        String root
        Int cpp=1
        Int mem=5
    }

  command <<<
        if [ -f ${default=abjdbashj lib} ]; then
            source ${lib}
        fi
        export LC_ALL=en_US.UTF-8
        echo "Sample ID:${ID}" > ${outdir}/report/2.sample.csv
        echo "FASTQ path:${sep=',' fastq1},${sep=',' fastq2}" >>  ${outdir}/report/2.sample.csv
        echo "Pipeline version:v0.1" >> ${outdir}/report/2.sample.csv
        echo "Reference path:${refdir}" >> ${outdir}/report/2.sample.csv
        
        grep "Number of reads:" ${outdir}/out/alignment_report.tsv| sed 's/\.//g' | awk -F ":" '{printf "Total number of reads :%'"'"'0.0f\n",$2}' > ${outdir}/report/3.sequencing.csv
        grep "Number of barcodes in whitelist\|Number of reads:" ${outdir}/out/alignment_report.tsv |sed 's/\.//g'|awk -F ':' '{print $2}'|awk 'NR==1{tmp=$1}NR>1{printf "Reads pairs with a valid barcode:%'"'"'0.2f%\n",$1*2/tmp*100}'  >> ${outdir}/report/3.sequencing.csv
        # grep "Q30 bases in Reads" ${outdir}/temp/sequencing_report.json | awk -F "," '{print "Q30 bases in Reads:"$2}' >> ${outdir}/report/3.sequencing.csv
        # # grep "Q30 bases in Cell Barcode" ${outdir}/temp/sequencing_report.json | awk -F "," '{print "Q30 bases in Barcode:"$2}' >> ${outdir}/report/3.sequencing.csv
        # grep "Fragments pass QC" ${outdir}/temp/sequencing_report.json | awk -F "," '{printf "Reads Pass QC:%'"'"'0.0f\n",$2}' >> ${outdir}/report/3.sequencing.csv
        
        grep "Number of mapped reads:" ${outdir}/out/alignment_report.tsv | sed 's/\.//g' |awk -F ":" '{print "Reads Mapped to Genome:"$2}' > ${outdir}/report/3.mapping.csv
        grep "Number of uniquely mapped reads:" ${outdir}/out/alignment_report.tsv | sed 's/\.//g'| awk -F ":" '{printf "Uniquely mapped reads:%'"'"'0.0f\n",$2}' >> ${outdir}/report/3.mapping.csv
        #grep "Mitochondria ratio," ${outdir}/temp/alignment_report.json | awk -F "," '{print "Mitochondria ratio:"$2}' >> ${outdir}/report/3.mapping.csv
        echo "Mitochondria ratio:0" >> ${outdir}/report/3.mapping.csv

        grep "bead_cutoff" ${outdir}/d2cfile/${ID}.d2cCutoff.tsv | awk '{printf "Bead threshold:%'"'"'0.0f\n",$2}' > ${outdir}/report/4.cells.csv
        wc -l ${outdir}/d2cfile/${ID}.barcodeMerge.tsv > ${outdir}/report/wc_barcode.txt
        awk -F "," '{printf "Bead number:%'"'"'0.0f\n",$1}' ${outdir}/report/wc_barcode.txt >> ${outdir}/report/4.cells.csv
        grep "cor_cutoff" ${outdir}/d2cfile/${ID}.d2cCutoff.tsv | awk '{printf "cor threshold:%0.5f\n",$2}' >> ${outdir}/report/4.cells.csv
        rm ${outdir}/report/wc_barcode.txt

        ${root}/bin/python ${root}/script/barcode.py --outPath ${outdir}   
        ${root}/bin/python ${root}/script/jaccard.py --outPath ${outdir}
        ${root}/bin/python ${root}/script/svg_to_base64_string.py --outPath ${outdir}
        ${root}/bin/python ${root}/script/data-table.py --outPath ${outdir}
        ${root}/bin/python ${root}/script/st.py --outPath ${outdir} --ID ${ID}
        ${root}/bin/python ${root}/script/generate.py --outPath ${outdir} --htmlTemplate ${root}/script/template-test-1-new-1.html --ID ${ID}


  >>>
  runtime {
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    File html = "${outdir}/report/${ID}_scATAC_analysis_report.html"
  }
}