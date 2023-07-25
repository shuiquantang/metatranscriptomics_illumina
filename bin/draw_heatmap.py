#from plots import plots
import pandas as pd
import numpy as np
import seaborn as sns
import os, argparse
import matplotlib.pyplot as plt
import sys, getopt

### get a dictionary of sample-to-category from mapping file
def get_cat(map_file, txt_filename):
    sample_to_cat = {}
    sample_to_color = {}
    categories = []
    with open(map_file, 'r') as fh1:
	for line in fh1:
	    line = line.rstrip()
	    if line:
		if line.startswith('#'):
		    pass
		else:
		    coln = line.split('\t')
		    try:
			sample_to_cat[coln[0]] = coln[3]
			if not coln[3] in categories:
			    categories.append(coln[3])
		    except IndexError:
			return [], {}, []

    #colormap = sns.diverging_palette(240, 10, n=len(categories))
    
    colormap = sns.color_palette("hls", len(categories))
    cat_to_color = dict(zip(categories, colormap))
    
    for i in sample_to_cat:
	cat = sample_to_cat[i]
	color_array = cat_to_color[cat]
	sample_to_color[i] = []
	sample_to_color[i] = color_array
    
    with open(txt_filename, 'r') as fh:
	header = fh.next().rstrip()
	samples = header.split('\t')
	
    samples = samples[1:]
    
    col_colors = []
    for i in samples:
	color = sample_to_color[i]
	col_colors.append(color)
	
    return col_colors, cat_to_color, categories

###Still to be modified to include additional filtering criteria 
def get_df(txt_filename, top_num):
    df = pd.read_csv(txt_filename, sep='\t')
    if top_num > 0:
    	df['abundance_sum'] = df.sum(axis=1)
    	df = df.sort(['abundance_sum'], ascending=False)[:top_num]
    	del df['abundance_sum']
    ## drop the columns with all zeros
    df = df.loc[:, (df != 0).any(axis=0)]
    df = df.set_index(['Taxon'])
    colns = len(df.columns)
    return df, colns
### row_cluster and col_cluster are turned on by default
### Modify to change the cluster pattern 
###plt.savefig can be modified to save as a png or pdf, png requires 300 to 600 dpi minimum
def draw_heatmap(df, out_filename, clustering, x_axis_font, y_axis_font, figure_width, colors, cat_to_colors, categories, row_clustering):
    #cg = sns.clustermap(df, method='average', metric='braycurtis', linewidths=0.1, row_cluster=True, col_cluster=clustering)
    sns.set(font_scale=1)
    if colors:
	cg = sns.clustermap(df, figsize=(figure_width, 10), method='average', metric='braycurtis', row_cluster=row_clustering, col_cluster=clustering, col_colors=colors, cmap='coolwarm')
	for cat in categories:
	    cg.ax_col_dendrogram.bar(0, 0, color=cat_to_colors[cat],label=cat, linewidth=0)
	cg.ax_col_dendrogram.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=4)
	
    else:
	cg = sns.clustermap(df, figsize=(figure_width, 10), method='average', metric='braycurtis', row_cluster=row_clustering, col_cluster=clustering, cmap='coolwarm')
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=y_axis_font)
    plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=x_axis_font)
    #plt.savefig(png_filename, dpi=600, bbox_inches='tight')
    plt.savefig(out_filename, bbox_inches='tight')


def main(args):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:m:t:a:b:r:",
        ["--abundance_table", "--mapping_file", "--taxatoshow", "--output_a", "--output_b", "--row_clustering"]
        )
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        sys.exit(2)
    for o, a in opts:
        if o in ('-i', '--abundance_table'):# work directory has to be a full path
            txt_filename = a
	elif o in ('-m', '--mapping_file'):
            map_file = a
        elif o in ('-t', '--taxatoshow'):
            top_num = int(a)
	elif o in ('-a', '--output_a'):
            out_filename1 = a
	elif o in ('-b', '--output_b'):
            out_filename2 = a
	elif o in ('-r', '--row_clustering'):
            row_clustering = a
        else:
            assert False, "unhandled option"
    
    (col_colors, cat_to_colors, categories) = get_cat(map_file, txt_filename)
    (df1, colns) = get_df(txt_filename, top_num)
    x_axis_font = int(20 - colns/12)
    figure_width = int(10 + (colns*2/15))
    y_axis_font = int(9 - (top_num-50)*2/50)
    if row_clustering == 'True':
	row_clustering_controller = True
    else:
	row_clustering_controller = False
    
    draw_heatmap(df1, out_filename1, True, x_axis_font, y_axis_font, figure_width, col_colors, cat_to_colors, categories, row_clustering_controller)
    draw_heatmap(df1, out_filename2, False, x_axis_font, y_axis_font, figure_width, col_colors, cat_to_colors, categories, row_clustering_controller)

if __name__ == "__main__":
    main(sys.argv)
