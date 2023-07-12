import pandas as pd
import os
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
        '''input the outpath''',)
    parser.add_argument('--htmlTemplate', type=str, help=
        '''input the html template''',)
    parser.add_argument('--ID', type=str, help=
        '''input the ID''',)
    args = parser.parse_args()
    return [args.outPath, args.htmlTemplate, args.ID,]

def get_args_from_file():
    path=get_args()[0]
    csv = [path+'/report/1.cell_report.csv',\
    path+'/report/2.sample.csv',\
    path+'/report/3.sequencing.csv',\
    path+'/report/4.cells.csv',\
    path+'/report/5.library.QC.csv',path+'/report/5_2.library.QC.csv',path+'/report/3.mapping.csv']

    stat = dict()
    for i in range(len(csv)):
        if i==0:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['estimate_num_of_cell'] = df[1][0]
            stat['median_frag_per_cell'] = df[1][1]
            stat['median_frac_peaks'] = df[1][2]
            stat['median_frac_tss'] = df[1][3]
        if i==1:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['sample_id'] = df[1][0]
            stat['fastq_path'] = df[1][1]
            stat['pipeline_version'] = df[1][2]
            stat['ref_path'] = df[1][3]
        if i==2:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['read_pairs'] = df[1][0]
            stat['frac_valid_barcode'] = df[1][1]
        if i==3:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['bead_thres'] = df[1][0]
            stat['bead_number'] = df[1][1]
            stat['jaccard_thres'] = df[1][2]
        if i==4:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['frac_frag_overlap'] = df[1][0]
            stat['call_peak_number'] = df[1][1]
            stat['overlap_call_peak'] = df[1][2]
            stat['percent_dup'] = df[1][3]
        if i==5:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['nc_free_region'] = df[1][0]
            stat['mono_nc_region'] = df[1][1]
        if i==6:
            df = pd.read_csv(open(csv[i]),encoding="utf_8",dtype=str,header=None,sep=":")
            stat['map_rate'] = df[1][1]
            stat['properly_reads'] = df[1][2]
            stat['mit_rate'] = df[1][0]
    plot_file = [
    path+'/report/div/barcode_rank.div',\
    path+'/report/div/jaccard_rank.div',\

    path+'/report/base64/plot3_DropBeadsnum.base64',\
    path+'/report/base64/plot4_QC.base64',\
    path+'/report/base64/plot5_InterSize.base64',\
    path+'/report/base64/plot6_TSS.base64',\
    path+'/report/base64/plot7_Cluster_peak.base64',\
    path+'/report/base64/plot8_Cluster_depth.base64',\
    path+'/report/base64/plot9_Cluster_annotation.base64',\
    path+'/report/div/saturation.div',\
]
    plot_base64 = []
    for f in plot_file:
        if re.search('plot7_Cluster_peak.base64',f) or re.search('plot8_Cluster_depth.base64',f) or re.search('plot9_Cluster_annotation.base64',f):
            if os.path.exists(f):
                base64 = open(f).read()
                #img = ("<img src=%s height=500px width=100\%>" %base64)
                img = "<img src=\"data:image/svg+xml;base64,"+base64+"\">"
                plot_base64.append(img)
            else:
                plot_base64.append('''
                <p style="font_family=DIN Next LT Pro;font_size=18px;font_weight=400">
                The cluster plot has not been generated because the data quality is too low.
                <p>
                ''')
        else:
            plot_base64.append(open(f).read())
    plot_order = ['plot1','plot2','plot3','plot4','plot5','plot6','plot7','plot8','plot9','plot10']
    plot_dict = dict(zip(plot_order, plot_base64))

    data_tables_file = [path+'/report/table/cell-table.txt']
    data_tables = []
    data_tables_dict = {}
    for f_ in data_tables_file:
        if os.path.exists(f_):
            if re.search('peak-table.txt',f_):
                table='''<table id=\"table_id_example\" class=\"display\" style=\"font-family:DIN Next LT Pro;\">            <thead style=\"font-size:11px\"><tr>
                    <th>peak.ID</th>
                    <th>width</th>
                    <th>gene.symbol</th>
                    <th>diatance</th>
                    <th>p_val</th>
                    <th>avg_log2FC</th>
                    <th>pct.1</th>

                    <th>pct.2</th>
                    <th>p_val_adj</th>
                    <th>cluster</th>

                </tr>
            </thead>
                <tbody style=\"font-size:11px;\">
                '''+open(f_).read()+"</tbody></table>"
                data_tables.append(table)
            elif re.search('cell-table.txt',f_):
                table='''<table id=\"table_id_example1\" class=\"display\" style=\"height:auto;padding-bottom:0px;font-family:DIN Next LT Pro\">
                    <thead style=\"font-size:11px\">
                    <tr>

                    <th>predicated.cell.type</th>
                    <th>number</th>
                    <th>ratio</th>

                    </tr>
            </thead>
                <tbody style=\"font-size:11px\">'''+open(f_).read()+"</tbody></table>"
                cell_table = open(f_).read()
                if cell_table.count("<tr>") < 10:
                    stat['cell-table-height'] = 31 * int(cell_table.count("<tr>"))
                else:
                    stat['cell-table-height'] = 310
                data_tables.append(table)
                #data_tables.append(open(f_).read())
            #with open(file,'r') as fr:
        else:
            data_tables.append('''
            <p style="font_family=DIN Next LT Pro;font_size=18px;font_weight=400">
            The table has not been generated because the data quality is too low.
            <p>
            ''')
    #print(data_tables)
    table_order = ['table2']
    data_tables_dict = dict(zip(table_order, data_tables))

    return stat, plot_dict, data_tables_dict

def write_param_to_template():
    stat, plot_dict, data_tables_dict = get_args_from_file()
    template = open(get_args()[1]).read()
    ID = get_args()[2]
    path=get_args()[0]
    from string import Template

    html=Template(template)
    if os.path.exists(path+'/report/base64/plot9_Cluster_annotation.base64'):
        report=html.safe_substitute(sample_info=ID, estimate_num_of_cell=stat['estimate_num_of_cell'],\
                    median_frag_per_cell=stat['median_frag_per_cell'], median_frac_peaks=stat['median_frac_peaks'],\
                    median_frac_tss=stat['median_frac_tss'],sample_id=stat['sample_id'],\
                    fastq_path=stat['fastq_path'], pipeline_version=stat['pipeline_version'],\
                    ref_path=stat['ref_path'],\
                    read_pairs=stat['read_pairs'],\
                    frac_valid_barcode=stat['frac_valid_barcode'],\
                    bead_thres=stat['bead_thres'],\
                    bead_number=stat['bead_number'],\
                    jaccard_thres=stat['jaccard_thres'],\
                    frac_frag_overlap=stat['frac_frag_overlap'],\
                    nc_free_region=stat['nc_free_region'],mono_nc_region=stat['mono_nc_region'],\
                    call_peak_number=stat['call_peak_number'],overlap_call_peak=stat['overlap_call_peak'],\
                    cell_table_height = stat['cell-table-height'],\
                    map_rate = stat['map_rate'],properly_reads = stat['properly_reads'],mit_rate = stat['mit_rate'],\
                    percent_dup=stat['percent_dup'],plot1=plot_dict['plot1'],\
                    plot2=plot_dict['plot2'],plot3=plot_dict['plot3'],\
                    plot4=plot_dict['plot4'],plot5=plot_dict['plot5'],\
                    plot6=plot_dict['plot6'],plot7=plot_dict['plot7'],\
                    plot8=plot_dict['plot8'],plot9=plot_dict['plot9'],\

                    plot10=plot_dict['plot10'],
                    table2=data_tables_dict['table2'],
                    )
    else:
        report=html.safe_substitute(sample_info=ID, estimate_num_of_cell=stat['estimate_num_of_cell'],\
                    median_frag_per_cell=stat['median_frag_per_cell'], median_frac_peaks=stat['median_frac_peaks'],\
                    median_frac_tss=stat['median_frac_tss'],sample_id=stat['sample_id'],\
                    fastq_path=stat['fastq_path'], pipeline_version=stat['pipeline_version'],\
                    ref_path=stat['ref_path'],\
                    read_pairs=stat['read_pairs'],\
                    frac_valid_barcode=stat['frac_valid_barcode'],\
                    bead_thres=stat['bead_thres'],\
                    bead_number=stat['bead_number'],\
                    jaccard_thres=stat['jaccard_thres'],\
                    frac_frag_overlap=stat['frac_frag_overlap'],\
                    nc_free_region=stat['nc_free_region'],mono_nc_region=stat['mono_nc_region'],\
                    call_peak_number=stat['call_peak_number'],overlap_call_peak=stat['overlap_call_peak'],\
                    cell_table_height = stat['cell-table-height'],\
                    map_rate = stat['map_rate'],properly_reads = stat['properly_reads'],mit_rate = stat['mit_rate'],\
                    percent_dup=stat['percent_dup'],plot1=plot_dict['plot1'],\
                    plot2=plot_dict['plot2'],plot3=plot_dict['plot3'],\
                    plot4=plot_dict['plot4'],plot5=plot_dict['plot5'],\
                    plot6=plot_dict['plot6'],plot7=plot_dict['plot7'],\
                    plot8=plot_dict['plot8'],plot9="There is no such species reference for annnotation.",\
                    #plot9=plot_dict['plot9'],
                    plot10=plot_dict['plot10'],
                    table2=data_tables_dict['table2'],
                    )
    return report

if __name__ == '__main__':
    outpath=get_args()[0]
    ID = get_args()[2]
    get_args_from_file()
    report=write_param_to_template()
    fw = open(outpath+'/report/'+ID+'_scATAC_analysis_report.html','w')
    fw.write(report)

