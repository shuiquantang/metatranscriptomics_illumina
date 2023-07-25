#from plots import plots
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import sys, getopt
from collections import OrderedDict

def get_df(txt_filename):
    names = ['Name', 'Score', 'Group_id', 'LDA(log10)', 'p-value']
    df = pd.read_csv(txt_filename, sep='\t', names=names)
    df = df.dropna()
    df = df.sort(['Group_id', 'LDA(log10)'], ascending=[True, True])
    df2 = df[['Name', 'Group_id', 'LDA(log10)']]
    df2 = df2.set_index(['Name'])
    return df2

def draw_lefse(df, pdf_filename, title_str):
    sns.set(style='white')
    fig = plt.figure()
    rows = len(np.unique(df['LDA(log10)']))
    categories = list(OrderedDict.fromkeys(df['Group_id']))
    print categories
    color_value = cm.rainbow(np.linspace(0, 1, len(categories)))
    colors = dict(zip(categories, color_value))
    ax = df['LDA(log10)'].plot(kind='barh', x='LDA(log10)', y='Name', color=[colors[i] for i in df['Group_id']], figsize = (6,rows*0.3))
    #fig.suptitle(title_str, fontsize=20)
    #plt.xlabel('LDA score(log10)', fontsize=18)
    ax.set_title(title_str,fontsize= 16)
    ax.set_xlabel("LDA score(log10)", fontsize=14)
    ax.set_ylabel("")
    legend = []
    for k in categories:
        patch = mpatches.Patch(color=colors[k], label=k)
        legend.insert(0,patch)
        
    plt.legend(handles=legend, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.savefig('test.png', dpi=300, bbox_inches='tight')
    plt.savefig(pdf_filename, bbox_inches='tight')


def main(args):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:",
        ["--abundance_table", "--output"]
        )
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        sys.exit(2)
    for o, a in opts:
        if o in ('-i', '--abundance_table'):# work directory has to be a full path
            txt_filename = a
	elif o in ('-o', '--output'):
            pdf_filename = a
        else:
            assert False, "unhandled option"
            
    title_str = 'Biomarkers Ordered by Effect Size (LDA Score)'
    lefse_df = get_df(txt_filename)
    draw_lefse(lefse_df, pdf_filename, title_str)
    
if __name__ == "__main__":
    main(sys.argv)





