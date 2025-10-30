from collections import defaultdict
import polars as pl

def common_figure_style(fig, y_italic=False, publish=False, small=True):
	font_base = 2 if publish and small else 4
	font_fam = 'Helvetica' if publish else 'Avenir'
	font = font_fam + ('' if publish else ' Black')
	ytick_font = font
	if y_italic and not publish: ytick_font += ' Oblique'
	tick_font = dict(family=font, size=font_base*5)
	ytick = dict(family=ytick_font, size=font_base*5)
	if y_italic and publish: ytick['style'] = 'italic'
	legend_font = dict(family=font, size=font_base*4)

	axis_style = dict(
		title_text='',
		showgrid=False,
		showline=True,
		linewidth=1 if publish else 2,
		linecolor='black',
		mirror=True,
		ticks='outside',
		ticklen=6,
		tickwidth=1 if publish else 2,
		color='black'
	)

	# mirror=True mirrors the axis on both sides (box), versus 'ticks', or 'allticks' to also mirror tick marks
	fig.update_xaxes(tickfont=tick_font, **axis_style)
	fig.update_yaxes(tickfont=ytick, **axis_style)
	fig.update_layout(legend=dict(font=legend_font), yaxis_title_font=legend_font, xaxis_title_font=legend_font)


# standardize chosen colors for IS fams
IS_fam_color_map = {
	'ISL3': '#ff0000',
	'IS256': '#1f77b4',
	'IS6': '#2ca02c',
	'IS30': '#9467bd',
	'IS3': '#8c564b',
	'IS200/IS605': '#f0934e',
	'IS110': '#f1ce1e',
	'other': '#7f7f7f',
	# May want to add IS1182 and IS21 for Staph spp
	'IS1182': '#dd4477',
	'IS21': '#b82e2e'

}
IS_fam_color_map = defaultdict(lambda: '#7f7f7f', IS_fam_color_map)


# Mostly from MoMA Klein https://github.com/BlakeRMills/MoMAColors
MLST_color_map = {
	'Lactis': "#579ea4",
	'A2': "#995041", # this one doesn't have parity actually
	'117': "#df7713",
	'80': "#f9c000",
	'17': "#86ad34",
	'18': "#5d7298",
	'203': "#81b28d",
	'1421': "#7e1a2f",
	'16': "#bd777a",
	'other': "#858794",
	# 'other': "#515260",
	# '-': "#a5a6ae",
	'-': "#e2d8d6"
}

species_color_map = {
	'EF': '#e4c9e1',
	'FM': '#a3cf55'
}

# TODO: standardized color palette across genera & species
# This can be paired with eskape_order, which is why it's not necessary to define as a map.
taxon_colors = {
	'Enterobacter spp.': "#c8350d",
	'Staphylococcus': "#5d7298", ### Genus, matching S. aureus
	'S. aureus': "#5d7298", 
	'K. pneumoniae': "#df7713",
	'A. baumannii': "#f9c000",
	'P. aeruginosa': "#7e1a2f", 
	'E. coli': "#bd777a",
	'Streptococcus': "#86ad34", ### Genus, not in ESKAPE
	'Enterococcus': "#579ea4", ### Genus, matching E. faecium 
	'E. faecium': "#579ea4",
	'E. faecalis': "#81b28d" # related color to E. faecium but not quite. 
}


# "#53362e", "#744940", "#9f7064", "#c99582", "#e6bcac", "#e2d8d6", "#a5a6ae", "#666879", "#515260", "#3d3d47"


def add_family_name(fig, fam, color, y):
	fig.add_annotation(
		x=1.01, # Slightly to the right of the plot area
		y=y,
		xref='paper', # uses the paper coordinate system for the x-axis, rather than within the plot
		yref='y', # keep the y-axis reference in data coordinates
		text=fam,
		font=dict(family='Helvetica', size=20, color=color), # taken from display_clusters#common_figure_style -- should standardize
		showarrow=False,
		xanchor='left',
		yanchor='middle'
	)



# We frequently want to have box with jittered points, and to set the color
# of the box border to black without affecting the color of the inside of the boxes,
# which needs to be done in kind of a custom way.
def common_box_style(fig, optimize_smaller=False):
	scale = 2 if optimize_smaller else 1
	fig.update_traces(marker=dict(size=2*scale), jitter=0.3)
	for trace in fig.data:
		if trace.type == 'box':
			base_color = trace.marker.color
			trace.update(fillcolor=base_color, line=dict(color='black', width=1*scale))


# Function used several times to take a group_by and aggregate on mean, stdev, and disallow negative error-y
def mean_std_error(counts_df_grouped, count_column):
	return (
		counts_df_grouped
		.agg(
			mean=pl.col(count_column).mean(),
			std=pl.col(count_column).std()
		)
		# Prevent lower error bar from being negative in the plot
		.with_columns(
			error_y_minus = (
				pl.when(pl.col('std') > pl.col('mean')) # mean-std is negative
				.then('mean')
				.otherwise('std')
			)
		)
	)


# provide opacity between 0-1. If the hex already has opacity it is overridden.
# Caution: this helper function does not do any validation.
def hex2rgb(x, a=None):
	x = x.lstrip('#')
	r,g,b = tuple(int(x[i:i+2],16) for i in (0,2,4))
	if len(x) > 6: a = int(x[6:],16)/255 # opacity->fraction
	return f'rgb({r},{g},{b})' if a is None else f'rgba({r},{g},{b},{a:.2f})'