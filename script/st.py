import plotly as py
from plotly.graph_objs import Scatter, Layout, Data, Scattergl
import pandas as pd
import plotly.graph_objs as go
from plotly.io import *
import os
import argparse
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
        '''input the outpath''',)
    parser.add_argument('--ID', type=str, help=
        '''input the outpath''',)
    args = parser.parse_args()
    return args.outPath,args.ID

if __name__ == '__main__':
    path , ID= get_args()    
    df=pd.read_csv(path+"/d2cfile/"+ID+".sequenceSaturation.tsv",sep="\t")
    data=Scattergl(x=df["mean_frags_per_cell"],y=df["saturation"], mode='lines')
    
    layout = Layout(xaxis=dict(
                        gridcolor="lightgrey",
                        title="mean_frags_per_cell",
                        color="black",
                        showline=True,
                        zeroline=True,
                        linewidth=1,fixedrange= True,
                        linecolor="black"
                        ),
                    yaxis = dict(
                        title="saturation",
                        gridcolor="lightgrey",
                        linewidth=1,fixedrange= True,
                        color="black",
                        linecolor="black"
                        ),
                        height=360,width=450,
                        plot_bgcolor='rgba(0,0,0,0)',
                        hovermode='closest',
                        paper_bgcolor='white',
                    title="saturation : "+str(df['saturation'].tolist()[-1])
    
                            
        )
                                
    #layout=
    fig=dict(data=data, layout=layout)
    #fig.show()
    config={'displayModeBar': False}
    py.offline.plot(fig, filename =path+"/report/div/saturation.html")
    fig2=py.offline.plot(fig, include_plotlyjs=False,show_link=False,output_type='div',config=config)
    fw = open(path+"/report/div/saturation.div",'w')
    fw.write(fig2)

