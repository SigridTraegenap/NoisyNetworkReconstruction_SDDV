import matplotlib as mpl
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Palatino']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
# Say, "the default sans-serif font is COMIC SANS"
#matplotlib.rcParams['font.sans-serif'] = "Comic Sans MS"
# Then, "ALWAYS use sans-serif fonts"
#matplotlib.rcParams['font.family'] = "sans-serif"
params = {
    'font.family': "sans-serif",
    'font.sans-serif': "CMU Sans Serif",
   'axes.labelsize': 18,
   'axes.spines.top'  : False ,
   'axes.spines.right'  : False ,
   'axes.linewidth' : 1.3,
   'font.size': 18,
   'legend.fontsize': 18,
   'xtick.labelsize': 16,
   'ytick.labelsize': 16,
   'xtick.major.top'      : False,
   'ytick.major.right'      : False,
   'figure.figsize': [4.5, 4.5],
   'lines.linewidth' : 2,
   #'errorbar.capsize' : 10,
'mathtext.fontset' : 'cm',
"figure.subplot.left"    : 0.2 , # the left side of the subplots of the figure
"figure.subplot.right"   : 0.9   , # the right side of the subplots of the figure
"figure.subplot.bottom"   : 0.17  ,  # the bottom of the subplots of the figure
"figure.subplot.top"     : 0.88 
}
mpl.rcParams.update(params)
