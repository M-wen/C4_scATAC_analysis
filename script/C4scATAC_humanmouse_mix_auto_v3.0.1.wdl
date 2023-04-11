version 1.0
workflow iDrop_scATAC_mix_2 {
  input{
    Array[Array[String]] Data
    String Outdir
    String SampleID
    String ProjectID
    String reference
    String? readStructure
    Int ?ForceFrag
  }

  String outdir=Outdir
  String ID=SampleID
  String PID = "sentieon.q -P "+ProjectID+"_sentieon"

##########    genome infor   ############

  Map[String,Array[String]] refConfig = {
    "hg19-mm10" : ["hg19-mm10","2403862610","/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/src/C4scATAC_humanmouse_mix_v2.3/ref/fasta/genome.fa","/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/src/C4scATAC_humanmouse_mix_v2.3/ref/tss.bed","/hwfssz5/CNGB_WRITE/USER/wangxin/01.database/13.scATAC/02.malu_new/Promoter_4k_malu.bed","hg19_chrM,mm10_chrM","/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/pipeline/01.scATAC/01.scATAC_v3.0/database/chromap_index/HumanMouseMix_genome.index"]
  }
########################################

  String refdir=refConfig[reference][2]
  String root="/jdfssz1/ST_SUPERCELLS/PUB/snATAC/pipeline/02.scATAC_humanmouse_v3.0"
  String tss=refConfig[reference][3]
  String promo=refConfig[reference][4]
  String chrmt = refConfig[reference][5]
  String ref_index = refConfig[reference][6]
  String runID= ID
  String? defaultConfig="newT1"
  String config=select_first([readStructure,defaultConfig])
  String lib="/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/src/C4scATAC_humanmouse_mix_v2.3/lib.sh"
  String species=refConfig[reference][0]
  String macs2sp=refConfig[reference][1]
  String whitelist="/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/test/4_scATAC_chromap/whitelist.txt"
  String read_format = if config == "oldT1" then "bc:50:-1,r2:0:49" else if config == "500" then "bc:50:59,bc:66:-1,r2:0:49" else "bc:0:19,r1:20:-1"

 

########################################### 


  call makedir {
    input:
    Dir=outdir
  }
  scatter(fq in Data){
    String fq1=fq[0]
    String fq2=fq[1]
  }
  Array[String] barcode = if config == "oldT1" then fq2 else if config == "500" then fq2 else fq1

  call mapping {
    input:
      fastq1=fq1,
      fastq2=fq2,
      outdir=makedir.Outdir,
      root=root,
      ref_index=ref_index,
      ref=refdir,
      whitelist=whitelist,
      barcode=barcode,
      read_format=read_format
  }
   
  call deconvolution { 
    input: 
      lib=lib,
      alnbed=mapping.alnbed,
      outdir=outdir,
      root=root,	
      ID=ID,
      species = species,
      ForceFrag=ForceFrag, 
      chrmt=chrmt 
  }   
  
  call callpeak as hg19callpeak{
     input:
     root=root,
     outdir=outdir,
     lib=lib,
     ID=ID,
     fragment=deconvolution.human_fragment,
     spices="hg19",
     Barcodelist=deconvolution.human_barcode_list,
     gsize="hs",
     QC_stat=deconvolution.human_QC_stat
  }
  call callpeak as mm10callpeak{
   input:
     root=root,
     outdir=outdir,
     lib=lib,
     fragment=deconvolution.mouse_fragment,
     ID=ID,
     spices="mm10",
     Barcodelist=deconvolution.mouse_barcode_list,
     gsize="mm",
     QC_stat=deconvolution.mouse_QC_stat
  }
  
  call QC as humanQC{
    input:
       root=root,
       outdir=outdir,
       ID=ID,
       lib=lib,
       spices="hg19",
       sp="Human",
       QC_stat=deconvolution.human_QC_stat,
       tss_bed="/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/USER/liuyang9/pipe/C4scATAC/database/refdata-cellranger-atac-hg19-1.1.0/regions/tss.bed",
       fragments=deconvolution.human_fragment,
       Peak=hg19callpeak.Peak,
       DB_list=deconvolution.human_barcode_list,
       color="\"#E41A1C\""
  }
  call QC as mouseQC{
    input:
      root=root,
      outdir=outdir,
      ID=ID,
      lib=lib,
      spices="mm10",
      sp="Mouse",
      QC_stat=deconvolution.mouse_QC_stat,
      tss_bed="/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/USER/liuyang9/pipe/C4scATAC/database/refdata-cellranger-atac-mm10-1.1.0/regions/tss.bed",
      fragments=deconvolution.mouse_fragment,
      Peak=mm10callpeak.Peak,
      DB_list=deconvolution.mouse_barcode_list,
      color="\"#377EB8\""
  }

  call report {
    input:
      lib=lib,
      ID=ID,
      root=root,
      outdir=outdir,
      lib=lib,
      plot1_input=deconvolution.barcodeQuantSimple,
      plot2_input=deconvolution.implicatedBarcodes,
      plot3=deconvolution.plot3,
      plot4=deconvolution.plot4,
      plot5=deconvolution.plot5,
      plot6=deconvolution.plot6,
      plot7_human=humanQC.plot7,
      plot7_mouse=mouseQC.plot7,
      plot8_human=humanQC.plot8,
      plot8_mouse=mouseQC.plot8,
      library_QC_human=humanQC.library_QC,
      library_QC_mouse=mouseQC.library_QC,
      human_cell_report=humanQC.cell_report,
      mouse_cell_report=mouseQC.cell_report,
      alignment_report=mapping.aln_report,
      bapParams=deconvolution.bapParams,
      barcodeTranslate=deconvolution.barcodeTranslate,
      collision_rate=deconvolution.collision_rate
  }

}

###### script  ######
 
task makedir {
  input{
    String Dir
  }
  
  command <<<
    mkdir -p ~{Dir}
    mkdir -p ~{Dir}/00.log
    mkdir -p ~{Dir}/01.Data
    mkdir -p ~{Dir}/02.chromap
    mkdir -p ~{Dir}/03.d2cfile
    mkdir -p ~{Dir}/04.split_bam
    mkdir -p ~{Dir}/05.Peak
    mkdir -p ~{Dir}/05.Peak/Peak
    mkdir -p ~{Dir}/05.Peak/Peak/Mouse
    mkdir -p ~{Dir}/05.Peak/Peak/Human
    mkdir -p ~{Dir}/06.QC
    mkdir -p ~{Dir}/07.report/div
    mkdir -p ~{Dir}/07.report/base64
    mkdir -p ~{Dir}/07.report/table
  >>>
  output {
    String Outdir="${Dir}"
  }
}


task  mapping{
  input{
    Array[String] fastq1
    Array[String] fastq2
    String outdir
    String ref_index
    String ref
    String root
    String ?lib
    String whitelist
    String read_format
    Array[String] barcode
    Int cpp=4
    Int mem=35
  }

  command <<<
    export PATH=/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/bin:$PATH
    export LD_LIBRARY_PATH="/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib:/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib64:$LD_LIBRARY_PATH"
   ~{root}/bin/chromap --preset atac --bc-error-threshold 0 --trim-adapters -x ~{ref_index} -r ~{ref} -1 ~{sep=',' fastq1} -2 ~{sep=',' fastq2} -o ~{outdir}/02.chromap/aln.bed --barcode ~{sep=',' barcode} --barcode-whitelist ~{whitelist} --read-format ~{read_format} -t ~{cpp} 2> ~{outdir}/02.chromap/alignment_report.tsv
	
 >>> 
  runtime{
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String alnbed="${outdir}/02.chromap/aln.bed"
    String aln_report="${outdir}/02.chromap/alignment_report.tsv"
  }
}

task deconvolution {
  input{
    String alnbed
    String outdir
    String root
    String ?lib
    String species
    String chrmt
    String ID
    Int ?ForceFrag
    Int ?mapq
    Int cpp=1
    Int mem=7
  }

  command <<<
	
    export LD_LIBRARY_PATH="/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scATAC/software/v1.3.7/gcclib/lib:/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scATAC/software/v1.3.7/gcclib/lib64:$LD_LIBRARY_PATH" && export PATH="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/Python-3/bin:$PATH"
    ~{root}/bin/d2c_v1.4.3/bin/d2c merge -i ~{alnbed} --mapq 30 --bf ~{default=0 ForceFrag} -r ~{species} --mix-species -o ~{outdir}/03.d2cfile -c 10 -n ~{ID} --mc ~{chrmt}  --sat --bt1 CB --bg /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/ref/hg38_mm10_chrom_size.tsv --ts /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/src/C4scATAC_humanmouse_mix_v2.3/ref/tss.bed
    ~{root}/bin/Rscript ~{root}/script/sta_collision_rate.R -I ~{outdir}/03.d2cfile/~{ID}.Metadata.tsv -O ~{outdir}

    perl ~{root}/script/fishInWinter.fq.gz.pl  -bf table -ff table -bc 1 -fc 4 ~{outdir}/03.d2cfile/Mouse_barcode_list.txt ~{outdir}/03.d2cfile/~{ID}.fragments.tsv.gz > ~{outdir}/03.d2cfile/~{ID}.Mouse.fragments.tsv
    perl ~{root}/script/fishInWinter.fq.gz.pl  -bf table -ff table -bc 1 -fc 4 ~{outdir}/03.d2cfile/Human_barcode_list.txt ~{outdir}/03.d2cfile/~{ID}.fragments.tsv.gz > ~{outdir}/03.d2cfile/~{ID}.Human.fragments.tsv
    ~{root}/bin/bgzip ~{outdir}/03.d2cfile/~{ID}.Mouse.fragments.tsv
    ~{root}/bin/bgzip ~{outdir}/03.d2cfile/~{ID}.Human.fragments.tsv
    ~{root}/bin/tabix -p bed ~{outdir}/03.d2cfile/~{ID}.Mouse.fragments.tsv.gz 
    ~{root}/bin/tabix -p bed ~{outdir}/03.d2cfile/~{ID}.Human.fragments.tsv.gz
    tail -n +2 ~{outdir}/03.d2cfile/~{ID}.Metadata.tsv |awk '{print $1}' > ~{outdir}/03.d2cfile/list.~{ID}.txt
	
  >>>
  runtime{
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String QCstat="${outdir}/03.d2cfile/${ID}.Metadata.tsv"
    String total_fragments="${outdir}/03.d2cfile/${ID}.fragments.tsv.gz"
    String baplist="${outdir}/03.d2cfile/list.${ID}.txt"
    String human_barcode_list="${outdir}/03.d2cfile/Human_barcode_list.txt"
    String mouse_barcode_list="${outdir}/03.d2cfile/Mouse_barcode_list.txt"
    String human_QC_stat="${outdir}/03.d2cfile/Human_QC_stat.csv"
    String mouse_QC_stat="${outdir}/03.d2cfile/Mouse_QC_stat.csv"
    String human_fragment="${outdir}/03.d2cfile/${ID}.Human.fragments.tsv.gz"
    String mouse_fragment="${outdir}/03.d2cfile/${ID}.Mouse.fragments.tsv.gz"
    String implicatedBarcodes="${outdir}/03.d2cfile/${ID}.CorrelationBarcodes.tsv.gz"
    String barcodeQuantSimple="${outdir}/03.d2cfile/${ID}.barcodeCount.tsv"
    String bapParams="${outdir}/03.d2cfile/${ID}.d2cCutoff.tsv"
    String barcodeTranslate="${outdir}/03.d2cfile/${ID}.barcodeMerge.tsv"
    String plot3="${outdir}/03.d2cfile/Plot3_Merged_Species_barcode_pollution.svg"
    String plot4="${outdir}/03.d2cfile/Plot4_Number_of_beads_per_droplet.svg"
    String plot5="${outdir}/03.d2cfile/Plot5_human_TSS_QC.svg"
    String plot6="${outdir}/03.d2cfile/Plot6_mouse_TSS_QC.svg"
    String collision_rate="${outdir}/03.d2cfile/1.collision_rate.tsv"
  }
}





task callpeak{
  input{
    String root
    String outdir
    String ID
    String spices
    String fragment
    String Barcodelist
    String gsize
    String QC_stat
    String lib
    Int cpp=2
    Int mem=10
  }

    
    command <<<
      source ~{lib}
      ~{root}/bin/macs2 callpeak -t ~{fragment} -f BED -g ~{gsize} -n ~{spices} -B -q 0.01 --nomodel --outdir ~{outdir}/05.Peak


  >>>
  runtime{
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String Peak = "${outdir}/05.Peak/${spices}_peaks.narrowPeak"
    File Peak1 = "${outdir}/05.Peak/${spices}_peaks.narrowPeak"
  }
}


  
task QC{
  input{
    String root
    String outdir
    String ID
    String spices
    String tss_bed
    String sp
    String QC_stat
    String fragments
    String color
    String Peak
    String lib 
    String DB_list
    Int cpp=1
    Int mem=1
  }

    command <<<
        source ~{lib}
       # ~{Rscript} ~{root}/script/plot_TSSEnrichment_FragSize.R -T ~{tss_bed} -F ~{fragments} -G ~{spices} -O ~{outdir} -SP ~{sp} -C ~{color}
  	~{root}/bin/Rscript} ~{root}/script/plot_TSSEnrichment_FragSize.R -T ~{tss_bed} -F ~{fragments} -G ~{spices} -O ~{outdir} -SP ~{sp} -C ~{color}
        ~{root}/bin/Rscript ~{root}/script/scATAC_humanmouse_mix_QC_report.2.R -I ~{Peak} -F ~{fragments} -Q ~{QC_stat} -O ~{outdir} -SP ~{sp} -C ~{DB_list}

        nf=$(less  ~{fragments} | awk '{if($3-$2<147) print $0}' | wc -l)
        mn=$(less  ~{fragments} | awk '{if($3-$2>147 && $3-$2<294) print $0}' | wc -l)
        total=$(less ~{fragments} |wc -l )
        nft=$(echo "$nf $total" |awk '{printf("%.2f\n",$1/$2*100)}')
        mnt=$(echo "$mn $total" |awk '{printf("%.2f\n",$1/$2*100)}')
        echo $nft |awk '{print "Fraction of nucleosome-free regions:"$1"%"}' >> ~{outdir}/06.QC/5.~{sp}.library.QC.csv
        echo $mnt | awk '{print "Fraction of fragments mono-nucleosome regions:"$1"%"}' >> ~{outdir}/06.QC/5.~{sp}.library.QC.csv

    >>>
    runtime{
        backend:"SGE"
        cpu:cpp
        memory:"${mem} GB"
    }
    output{
       String plot7="${outdir}/06.QC/Plot7_TSS_${spices}.svg"
       String plot8="${outdir}/06.QC/plot8_InterSize_${spices}.svg"
       String library_QC="${outdir}/06.QC/5.${sp}.library.QC.csv"
       String cell_report="${outdir}/06.QC/1.${sp}.cell_report.csv"
    }

}



task report {
  input{
    String ID
    String ?lib
    String outdir
    String root
    String plot1_input
    String plot2_input
    String plot3
    String plot4
    String plot5
    String plot6
    String plot7_human
    String plot7_mouse
    String plot8_human
    String plot8_mouse
    String human_cell_report
    String mouse_cell_report
    String alignment_report
    String bapParams
    String barcodeTranslate
    String library_QC_human
    String library_QC_mouse
    String collision_rate
    Int cpp=1
    Int mem=5
  }

  command <<<

    source ~{lib}
    export LC_ALL=en_US.UTF-8
    
    cp ~{plot1_input} ~{outdir}/07.report/plot_input1_Bead_Barcode_Knee.csv && cp ~{plot2_input} ~{outdir}/07.report/plot_input2_Jaccard_Overlap_Knee.csv.gz
    
    cp ~{plot3} ~{outdir}/07.report/Plot3_Merged_Species_barcode_pollution.svg
    cp ~{plot4} ~{outdir}/07.report/Plot4_Number_of_beads_per_droplet.svg
    cp ~{plot5} ~{outdir}/07.report/Plot5_human_TSS_QC.svg
    cp ~{plot6} ~{outdir}/07.report/Plot6_mouse_TSS_QC.svg
    cp ~{plot7_human} ~{outdir}/07.report/Plot7_TSS_hg19.svg
    cp ~{plot7_mouse} ~{outdir}/07.report/Plot7_TSS_mm10.svg
    cp ~{plot8_human} ~{outdir}/07.report/plot8_InterSize_hg19.svg
    cp ~{plot8_mouse} ~{outdir}/07.report/plot8_InterSize_mm10.svg
    cp ~{human_cell_report} ~{outdir}/07.report/1.hg19.cell_report.csv
    cp ~{mouse_cell_report} ~{outdir}/07.report/1.mm10.cell_report.csv
    cp ~{collision_rate} ~{outdir}/07.report/1.collision_rate.tsv
    cp ~{library_QC_human} ~{outdir}/07.report/5.hg19.library.QC.csv
    cp ~{library_QC_mouse} ~{outdir}/07.report/5.mm10.library.QC.csv


    grep "Number of reads:" ~{alignment_report}| sed 's/\.//g' | awk -F ":" '{printf "Total number of reads :%'"'"'0.0f\n",$2/2}' > ~{outdir}/07.report/2.sequencing.csv
    grep "Number of barcodes in whitelist\|Number of reads:" ~{alignment_report} |sed 's/\.//g'|awk -F ':' '{print $2}'|awk 'NR==1{tmp=$1}NR>1{printf "Reads pairs with a valid barcode:%'"'"'0.2f%\n",$1*2/tmp*100}'  >> ~{outdir}/07.report/2.sequencing.csv

    # head -n5 ~{alignment_report} > ~{outdir}/07.report/3.alignment.csv 
    grep "Number of mapped reads:\|Number of barcodes in whitelist" ~{alignment_report} |sed 's/\.//g'|awk -F ':' '{print $2}' |awk 'NR==1{tmp=$1}NR>1{printf "Reads Mapped to Genome:%'"'"'0.2f%\n",tmp/($1*2)*100}' > ~{outdir}/07.report/3.mapping.csv
    grep "Number of uniquely mapped reads:" ~{alignment_report} | sed 's/\.//g'| awk -F ":" '{printf "Uniquely mapped reads:%'"'"'0.0f\n",$2}' >> ~{outdir}/07.report/3.mapping.csv
    #grep "Mitochondria ratio," ~{outdir}/temp/alignment_report.json | awk -F "," '{print "Mitochondria ratio:"$2}' >> ~{outdir}/07.report/3.mapping.csv
    echo "Mitochondria ratio:-" >> ~{outdir}/07.report/3.mapping.csv

    grep "bead_cutoff" ~{bapParams} | awk '{printf "Bead threshold:%'"'"'0.0f\n",$2}' > ~{outdir}/07.report/4.cells.csv
    wc -l ~{barcodeTranslate} > ~{outdir}/07.report/wc_barcode.txt
    awk '{printf "Bead number:%'"'"'0.0f\n",$1}' ~{outdir}/07.report/wc_barcode.txt >> ~{outdir}/07.report/4.cells.csv
    grep "cor_cutoff" ~{bapParams} | awk '{printf "Jaccard threshold:%0.5f\n",$2}' >> ~{outdir}/07.report/4.cells.csv
    #rm ~{outdir}/07.report/wc_barcode.txt
    
    ~{root}/bin/python ~{root}/script/barcode.py --outPath ~{outdir}
    ~{root}/bin/python ~{root}/script/jaccard.py --outPath ~{outdir}
    ~{root}/bin/python ~{root}/script/svg_to_base64_string.py --outPath ~{outdir}
    ~{root}/bin/python ~{root}/script/generate.py --outPath ~{outdir} --htmlTemplate ~{root}/script/scATAC-report-demo-template-test-1-new-1-sym-cellline.html --ID ~{ID}
    >>>
  runtime{
    backend:"SGE"
    cpu:cpp
    memory:"${mem} GB"
  }
  output{
    File html = "${outdir}/07.report/${ID}_C4-scATAC_HumanMouseMix_report.html"
  }
}
