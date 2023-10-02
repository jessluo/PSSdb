import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from wcmatch.pathlib import Path # Handling of path object
import os
import yaml
import math as m
path_to_config=Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg= yaml.safe_load(config_file)


from plotnine import *
import seaborn as sns
import matplotlib.pyplot as plt

path_to_datafile=(Path(cfg['git_dir']).expanduser()/ cfg['dataset_subdir']) / 'NBSS_data' / 'Sensitivity_analysis'
path_files=list(path_to_datafile.rglob('*_1b_*.csv'))
path_files_1a=list(path_to_datafile.rglob('*_1a_*.csv'))
df=pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name.split('_')[0], Biovol_metric = path.name.split('_')[2]),path_files)).drop_duplicates().reset_index(drop=True)
colors = {'Biovolume-ellipsoid': 'teal', 'Biovolume-area': 'black', 'Biovolume-orig': 'slateblue'}
plot = (ggplot(data=df)+
        geom_violin(aes(x='Instrument', y='NBSS_slope_mean', color='Biovol_metric'),  position = position_dodge(width=1))+
        geom_point(aes(x='Instrument', y='NBSS_slope_mean', color='Biovol_metric'), position = position_dodge(width=1),size = 1, alpha=0.4, shape = 'o')+
        stat_summary(aes(x='Instrument', y='NBSS_slope_mean', color='Biovol_metric'),geom='point', fun_y=np.mean, shape=0, size = 10,  position = position_dodge(width=1))+
        labs(y=r'Mean slope ( dm$^{-3}$ $\mu$m$^{-3}$)', x='')+
        scale_color_manual(values = colors)+
        scale_fill_manual(values=['#00000000', '#00000000', '#00000000']) +
        theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              panel_border = element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              #panel_border=element_rect(color='#222222'),
              legend_title=element_text(family="serif", size=10),
              legend_position=[0.5, 0.95],
              legend_text=element_text(family="serif", size=10),
              axis_title=element_text(family="serif", size=15),
              axis_text_x=element_text(family="serif", size=12),
              axis_text_y=element_text(family="serif", size=12, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(5,5)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Supp_fig_3.1_Biovol_sensitivity.pdf'.format(str(Path.home())), dpi=600)

plot = (ggplot(data=df)+
        geom_violin(aes(x='Instrument', y='NBSS_intercept_mean', color='Biovol_metric'),  position = position_dodge(width=1))+
        geom_point(aes(x='Instrument', y='NBSS_intercept_mean', color='Biovol_metric'), position = position_dodge(width=1),size = 1, alpha=0.4, shape = 'o')+
        stat_summary(aes(x='Instrument', y='NBSS_intercept_mean', color='Biovol_metric'),geom='point', fun_y=np.mean, shape=0, size = 10,  position = position_dodge(width=1))+
        labs(y=r'Mean Intercept ( $\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)', x='')+
        scale_color_manual(values = colors)+
        scale_fill_manual(values=['#00000000', '#00000000', '#00000000']) +
        theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              panel_border = element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              #panel_border=element_rect(color='#222222'),
              legend_title=element_text(family="serif", size=10),
              legend_position=[0.5, 0.95],
              legend_text=element_text(family="serif", size=10),
              axis_title=element_text(family="serif", size=15),
              axis_text_x=element_text(family="serif", size=12),
              axis_text_y=element_text(family="serif", size=12, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(5,5)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Supp_fig_3.2_Biovol_sensitivity.pdf'.format(str(Path.home())), dpi=600)

plot = (ggplot(data=df)+
        geom_violin(aes(x='Instrument', y='NBSS_r2_mean', color='Biovol_metric'),  position = position_dodge(width=1))+
        geom_point(aes(x='Instrument', y='NBSS_r2_mean', color='Biovol_metric'), position = position_dodge(width=1),size = 1, alpha=0.4, shape = 'o')+
        stat_summary(aes(x='Instrument', y='NBSS_r2_mean', color='Biovol_metric'),geom='point', fun_y=np.mean, shape=0, size = 10,  position = position_dodge(width=1))+
        labs(y=r'Mean R$^{2}$', x='')+
        scale_color_manual(values = colors)+
        scale_fill_manual(values=['#00000000', '#00000000', '#00000000']) +
        theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              panel_border = element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              #panel_border=element_rect(color='#222222'),
              legend_title=element_text(family="serif", size=10),
              legend_position=[0.5, 0.95],
              legend_text=element_text(family="serif", size=10),
              axis_title=element_text(family="serif", size=15),
              axis_text_x=element_text(family="serif", size=12),
              axis_text_y=element_text(family="serif", size=12, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(5,5)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Supp_fig_3.3_Biovol_sensitivity.pdf'.format(str(Path.home())), dpi=600)

import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

df_UVP = df[df.Instrument == 'UVP'].reset_index(drop=True)
df_IFCB = df[df.Instrument == 'IFCB'].reset_index(drop=True)

stats.f_oneway(df_UVP['NBSS_slope_mean'][df_UVP['Biovol_metric'] == 'Biovolume-area'],
               df_UVP['NBSS_slope_mean'][df_UVP['Biovol_metric'] == 'Biovolume-ellipsoid'])

model_slope = ols('NBSS_slope_mean ~ C(Biovol_metric)', data=df_UVP).fit()
aov_table_slope = sm.stats.anova_lm(model_slope, typ=2)
aov_table_slope

model_intercept = ols('NBSS_slope_mean ~ C(Biovol_metric)', data=df_IFCB).fit()
aov_table_intercept  = sm.stats.anova_lm(model_intercept, typ=2)
aov_table_intercept

