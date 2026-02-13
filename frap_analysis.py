"""
Frap_anysis-Copy2_clean.py
Auto-generated from the cleaned notebook for publication.
Add usage notes and author information as needed.
"""

#!/usr/bin/env python
# coding: utf-8

# # Frap Analysis — cleaned for publication
# 
# This notebook contains the Frap analysis code. It has been cleaned for public release:
# 
# - Outputs and execution counts removed
# - Top-level description and usage added
# - requirements extracted (best-effort) in `requirements_Frap_anysis.txt`
# 
# If you publish on GitHub, include the notebook, the generated script, and the requirements file.
# 
# **Author:** please add your name
# 
# **License:** MIT (see LICENSE_MIT.txt)

# In[ ]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import scipy.optimize as opt 
import glob
from scipy.stats import sem
import matplotlib
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#2e4185ff","#E7872B"]) ##FF9B00","#2E4185"]
from sklearn.metrics import r2_score
from matplotlib.colors import LogNorm, Normalize


# In[ ]:


# Meta data conferims that it is 5sec between messuments


# In[ ]:


def csv2data(path,ele='-',cell=0):
    t = pd.read_csv(path)
    t = t.loc[0:70]

    pre_bleach_reff  = np.mean(t.loc[0:2]['total'] - t.loc[0:2]['background'])
    pre_bleach_ROI  = np.mean(t.loc[0:2]['bleach'] - t.loc[0:2]['background'])

    norm_reff   = pre_bleach_reff/(t['total'] - t['background'])
    norm_ROI = (t['bleach'] - t['background'])/pre_bleach_ROI
    t['norm_reff'] = norm_reff*norm_ROI#t['total']- t['background']
    t['norm_full'] = 100*(t['norm_reff']-t['norm_reff'][3])/(t['norm_reff'][0] -t['norm_reff'][3])
    t['mdSAM']  = ele
    t['cell']  = cell
    t['time_n'] = np.NaN
    c = 0
    for i in t.index:
        t.loc[i,'time_n'] = c 
        c = c+5
#    t['time_n'] = [i+5 for i in range(len(t))]
    return t.loc[3:,['time','time_n','norm_full','mdSAM','cell']]


# In[ ]:


def prisim(t,Y0,Plateau,K):
    '''fitting recovery 
     - I_zero: Bottom
     - a: 
     - b: beta '''
    return Y0 + (Plateau-Y0)*(1-np.exp(-K*t))
def tau_function(t,I_zero,a,b):
    '''fitting recovery 
     - I_zero: Bottom
     - a: 
     - b: beta '''
    return I_zero - a*(np.exp(-b*t))

def tau_function_v2(t,A,tau):#y(t)=A⋅(1−e−τ⋅t)
    '''fitting recovery 
     - I_zero: Bottom
     - a: 
     - b: beta '''
    return A*(1-np.exp(-tau*t))
#Y=Y0 + (Plateau-Y0)*(1-exp(-K*x))

def prisim_v2(t,Plateau,K):
    '''fitting recovery 
     - I_zero: Bottom
     - a: 
     - b: beta '''
    return  (Plateau)*(1-np.exp(-K*t))


# In[ ]:


df = pd.DataFrame()
for c in enumerate(glob.glob('FRAP/fiji_data_prossess/csvs/*.csv')):
#    print(c)
    df = pd.concat([df,csv2data(c[1],ele= '+',cell =c[0])],ignore_index=True)
for c in enumerate(glob.glob('FRAP/fiji_data_prossess/csvs/dAPP/*.csv')):
#    print(c)
    df = pd.concat([df,csv2data(c[1],ele= '-',cell =c[0])],ignore_index=True)
for c in enumerate(glob.glob('FRAP/fiji_data_prossess/csvs/APP2P/*.csv')):
#    print(c)
    df = pd.concat([df,csv2data(c[1],ele= '2P',cell =c[0])],ignore_index=True)
df['time_post_bleach'] =df['time_n']-15


# In[ ]:


def fit_prisim(mdSAM='-'):
    x_data = df[df.mdSAM==mdSAM].loc[:,['norm_full','time_post_bleach']].groupby(['time_post_bleach']).mean().index
    y_data = df[df.mdSAM==mdSAM].loc[:,['norm_full','time_post_bleach']].groupby(['time_post_bleach']).mean().norm_full
    y_std = df[df.mdSAM==mdSAM].loc[:,['norm_full','time_post_bleach']].groupby(['time_post_bleach']).std().norm_full
    fit_, covar_ = opt.curve_fit(prisim, x_data[1:], y_data[1:],
                               sigma = y_std[1:],
                               absolute_sigma = True,maxfev=100000)
#print(best_fit_ab)
    sigma_abs = np.sqrt(np.diagonal(covar_))
    y_pred = prisim(x_data, *fit_)
    return fit_,sigma_abs,r2_score(y_data, y_pred)



# In[ ]:


fit_APP,sigma_APP,APP_r2 = fit_prisim(mdSAM='+')
fit_dAPP,sigma_dAPP,dAPP_r2 = fit_prisim(mdSAM='-')
fit_APP2P,sigma_APP2P,APP2P_r2 = fit_prisim(mdSAM='2P')


# In[ ]:


APP2P_r2


# In[ ]:


print('dAPP',1/fit_dAPP[2],sigma_dAPP[2])
print('APP',1/fit_APP[2],sigma_APP[2])
print('APP2P',1/fit_APP2P[2],sigma_APP2P[2])


# In[ ]:


sigma_APP[1]


# In[ ]:


fig, ax = plt.subplots(figsize=(2,3)) 

ax.errorbar(x='dAPP', y=fit_dAPP[1], yerr=sigma_dAPP[1],
                fmt='none', capsize=5, capthick=0, color='#76787A',zorder=0,rasterized=True)
ax.bar(x='dAPP', height=fit_dAPP[1],color='#76787A',lw=2)
ax.errorbar(x='APP', y=fit_APP[1], yerr=sigma_APP[1],
                fmt='none', capsize=5, capthick=0, color='#C03830',zorder=0,rasterized=True)
ax.bar(x='APP', height=fit_APP[1],fill=True,color='#C03830',lw=2)
ax.errorbar(x='2P', y=fit_APP2P[1], yerr=sigma_APP2P[1],
                fmt='none', capsize=5, capthick=0, color='#4E4F51',zorder=0,rasterized=True)
ax.bar(x='2P', height=fit_APP2P[1],fill=True,color='#4E4F51',lw=2)
plt.locator_params(axis='y',nbins=4)
plt.savefig('FRAP/frap_v3_AMP.svg',bbox_inches='tight',transparent=True)


# In[ ]:


fig, ax = plt.subplots(figsize=(2,3)) 

ax.errorbar(x='dAPP', y=fit_dAPP[2], yerr=sigma_dAPP[2],
                fmt='none', capsize=5, capthick=0, color='#76787A',zorder=0,rasterized=True)
ax.bar(x='dAPP', height=fit_dAPP[2],color='#76787A',lw=2)
ax.errorbar(x='APP', y=fit_APP[2], yerr=sigma_APP[2],
                fmt='none', capsize=5, capthick=0, color='#C03830',zorder=0,rasterized=True)
ax.bar(x='APP', height=fit_APP[2],fill=True,color='#C03830',lw=2)
ax.errorbar(x='2P', y=fit_APP2P[2], yerr=sigma_APP2P[2],
                fmt='none', capsize=5, capthick=0, color='#4E4F51',zorder=0,rasterized=True)
ax.bar(x='2P', height=fit_APP2P[2],fill=True,color='#4E4F51',lw=2)
plt.ylim([0,0.035])
plt.locator_params(axis='y',nbins=4)
plt.ylabel('Rate constant')
plt.savefig('FRAP/frap_v3_K.svg',bbox_inches='tight',transparent=True)


# In[ ]:


fig, ax = plt.subplots(figsize=(3,3)) 
sns.lineplot(x='time_post_bleach',y='norm_full',data=df,errorbar='se',hue='mdSAM',hue_order=['+','-','2P']
             ,palette=['#C03830','#808080','#58595B'],ax=ax,
             err_style='band',estimator=np.mean,lw=0,marker="o", markersize=4, markeredgewidth=.0, markeredgecolor="black")

plt.ylabel('Releative recovery (%)')
#ax.grid(zorder=0)

ax.set_ylim([-5,70])
ax.set_xlim([-5,200])
ax.locator_params(axis='y', nbins=4)
ax.set_xticklabels(ax.get_xticklabels(),fontsize=16,rotation=0)
ax.set_yticklabels(ax.get_yticklabels(),fontsize=16,rotation=0)

refDose = np.linspace(0,200,100)
plt.plot([ i for i in refDose],[prisim(i,*fit_dAPP) for i in refDose],lw=3,c='#808080')
plt.plot([ i for i in refDose],[prisim(i,*fit_APP) for i in refDose],lw=3,c='#C03830',ls='-')
plt.plot([ i for i in refDose],[prisim(i,*fit_APP2P) for i in refDose],lw=3,c='#58595B')
# plt.plot([ i for i in refDose],[prisim_dAPP_kata(i) for i in refDose],lw=2,c='k')
# plt.plot([ i for i in refDose],[prisim_APP_kata(i) for i in refDose],lw=2,c='k')
ax.set_ylabel('')
ax.set_xlabel('')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.set_ylabel('Releative recovery (%)',fontsize=12)
ax.set_xlabel('Time (sec)',fontsize=12)
# ax.set_title('$y=A\\cdot (1-e^{-\\tau \\cdot t})$',fontsize=20)
# ax.text(62,75,"A =%.0f±%.0f\n$\\tau$ =%.3f±%.3f"%(A_dAPP,A_dAPP_error,tau_dAPP,tau_dAPP_error),fontsize=18)
# ax.text(62,-8,"A =%.0f±%.0f\n$\\tau$ =%.3f±%.3f"%(A_APP,A_APP_error,tau_APP,tau_APP_error),fontsize=18)


#ax.text(62,-8,"A = {} \n$\\tau$ = {}".format(APP_param_platu, APP_param_tau),fontsize=20)
ax.legend(title='mdSAM',title_fontsize=16,fontsize=12)
#ax.set_title('Structure',fontsize=16)
plt.legend([],frameon=False)
plt.savefig('FRAP/frap_v3.svg',bbox_inches='tight',transparent=True)


# In[ ]:


from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()
df.drop_duplicates


# In[ ]:





# In[ ]:


def flatten(xss):
    return [x for xs in xss for x in xs]
fig, ax = plt.subplots(figsize=(3,3)) 
sns.scatterplot(x=df[(df.mdSAM=='-')&(df.time_post_bleach.between(0,200))].time_post_bleach.drop_duplicates().to_list(),
                y=flatten(min_max_scaler.fit_transform(df[(df.mdSAM=='-')&(df.time_post_bleach.between(0,200))].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).mean())),color='#808080')
sns.scatterplot(x=df[(df.mdSAM=='+')&(df.time_post_bleach.between(0,200))].time_post_bleach.drop_duplicates().to_list(),
                y=flatten(min_max_scaler.fit_transform(df[(df.mdSAM=='+')&(df.time_post_bleach.between(0,200))].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).mean())),color='#C03830')
sns.scatterplot(x=df[(df.mdSAM=='2P')&(df.time_post_bleach.between(0,200))].time_post_bleach.drop_duplicates().to_list(),y=flatten(min_max_scaler.fit_transform(df[(df.mdSAM=='2P')&(df.time_post_bleach.between(0,200))].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).mean())),color='#58595B')
plt.xlim([0,200])
plt.xlabel('time (sec)')
plt.ylabel('normilized recovery')


# In[ ]:


df[df.mdSAM=='-'].pivot(columns='cell',index='time_post_bleach',values='norm_full').to_csv('Kata_fit_FRAP_dAPP.csv')
df[df.mdSAM=='+'].pivot(columns='cell',index='time_post_bleach',values='norm_full').to_csv('Kata_fit_FRAP_APP.csv')
df[df.mdSAM=='2P'].pivot(columns='cell',index='time_post_bleach',values='norm_full').to_csv('Kata_fit_FRAP_APP2P.csv')


# In[ ]:


sns.scatterplot(x='time_post_bleach',y='norm_full',data=df[df.mdSAM=='-'],hue='cell',palette='viridis')


# In[ ]:


import warnings # To suppress some warnings

# Suppress the specific FutureWarning
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")


# In[ ]:


l = []

y_std_app = df[df.mdSAM=='+'].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).std().norm_full/df[df.mdSAM=='+'].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).mean().norm_full
y_std_dapp = df[df.mdSAM=='-'].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).std().norm_full/df[df.mdSAM=='-'].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).mean().norm_full
y_std_app2p = df[df.mdSAM=='2P'].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).std().norm_full/df[df.mdSAM=='-'].loc[:,['time_post_bleach','norm_full']].groupby(['time_post_bleach']).mean().norm_full

for i in y_std_app.iloc[1:].index:
    l.append({'t':i,'val':100*y_std_app.loc[i],'mdSAM':'+'})
for i in y_std_dapp.iloc[1:].index:
    l.append({'t':i,'val':100*y_std_dapp.loc[i],'mdSAM':'-'})
for i in y_std_app2p.iloc[1:].index:
    l.append({'t':i,'val':100*y_std_app2p.loc[i],'mdSAM':'2P'})
sem_df = pd.DataFrame(l)


# In[ ]:


fig, ax = plt.subplots(figsize=(3,3)) 
sns.lineplot(x='t',y='val',data=sem_df,hue='mdSAM')
plt.yscale('log')
plt.xlim([0,200])
plt.savefig('frap_v3_CV.svg',bbox_inches='tight',transparent=True)


# In[ ]:


fig, ax = plt.subplots(figsize=(4,5)) 
g = sns.swarmplot(x='mdSAM',y='val',data=sem_df,hue='t',s=6,palette='binary_r',linewidth=.75,edgecolor='k',alpha=1)
plt.yscale('log')
plt.xlim([-1,2])
ax.grid()
plt.ylabel('Coefficient of Variation (%)',fontsize=16)
plt.xlabel('mdSAM',fontsize=16)
ax.set_xticklabels(ax.get_xticklabels(),fontsize=16,rotation=0)
ax.set_yticklabels(ax.get_yticklabels(),fontsize=16,rotation=0)
ax.set_title('FRAP CV',fontsize=16)
norm = Normalize(vmax=200,vmin=3)
ax.legend([],frameon=False)
sm = plt.cm.ScalarMappable(cmap='binary_r', norm=norm)
sm.set_array([])
cbar = ax.figure.colorbar(sm,shrink=.55)
cbar.ax.tick_params(labelsize=12)
cbar.set_label(label='Time (sec)', size='large')
plt.savefig('FRAP/frap_v2_CV.svg',bbox_inches='tight',transparent=True)


# In[ ]:


stats.ttest_ind(sem_df[sem_df.mdSAM=='+'].val,sem_df[sem_df.mdSAM=='-'].val)


# In[ ]:


sem_df[sem_df.mdSAM=='+'].val.describe()


# In[ ]:


1/(53.6133/235.2924)


# In[ ]:





# In[ ]:


fig, ax = plt.subplots(figsize=(3,4))
ytable  = ax.table(cellText=[['$66$','$43$'],['$56$','$70$'],['$39$','$48$'],['$.98$','$.90$']],
                      rowLabels=['$Plateau\,(\%)$','$Tau\,(\\frac{1}{sec})$','$Half‑time\,(s)$','$R^{2}$'],
                      colLabels=['$-$','$+$'],
                      loc='center',rowLoc='center',cellLoc='center',colColours=['gray','r'])
ax.axis('off')
plt.text(.33,.8,'$mdSAM$')
ytable.scale(.8, 2)
ytable.set_fontsize(14)
plt.savefig('FRAP/frap_table.svg',bbox_inches='tight',transparent=True)


# In[ ]:


std = []
for app in ['APP','dAPP']:
    for i in range(3,228+3,3):
        std.append({'time':i,'val':sem(df[(df.new_time==i)&(df.mdSAM==app)].Full_norm),'mdSAM':app})
std_df = pd.DataFrame(std)


# In[ ]:





# In[ ]:


m = []
for app in ['APP','dAPP']:
    for i in range(3,228+3,3):
        m.append({'time':i,'val':df[(df.new_time==i)&(df.mdSAM==app)].Full_norm.mean(),'mdSAM':app})
m_df = pd.DataFrame(m)


# In[ ]:


m_df[m_df.mdSAM=='dAPP'].val.to_clipboard(index=0)


# In[ ]:





# In[ ]:


fig, ax = plt.subplots(figsize=(3,3)) 
sns.stripplot(x='mdSAM',y='val',data=std_df[std_df.time.between(3,200)],hue='time',palette='viridis',s=5)#cividis_r
sns.boxplot(x='mdSAM',y='val',data=std_df[std_df.time.between(3,200)],
            boxprops=dict(alpha=1),palette=['w'],zorder=3,showfliers = False)#cividis_r

sm = plt.cm.ScalarMappable(cmap='viridis',norm=matplotlib.colors.Normalize(vmin=3, vmax=200))
sm.set_array([])
cbar = ax.figure.colorbar(sm,shrink=.5)
cbar.ax.tick_params(labelsize=10)
cbar.set_label(label='$time\,(Sec)$', size='small')
ax.grid()
plt.legend([],frameon=False)
plt.ylabel('$SEM$',fontsize=12)
ax.set_xticklabels(['$+$','$-$'],fontsize=12)
plt.xlabel('$mdSAM$',fontsize=12)
pval = stats.ttest_ind(std_df[(std_df.time.between(3,200))&(std_df.mdSAM=='APP')].val,
               std_df[(std_df.time.between(3,200))&(std_df.mdSAM=='dAPP')].val)[-1]
plt.text(0.4,.095,'$p=%.e$'%pval,fontsize=12)
plt.savefig('FRAP/sem.svg',bbox_inches='tight',transparent=True)


# In[ ]:


tau = []
for mdSAM in ['APP','dAPP']:
    for ex in set(df[(df.new_time>=0)&(df.mdSAM==mdSAM)].exp):
        fit_df = df[(df.new_time>=0)&(df.mdSAM==mdSAM)&(df.exp==ex)&(df.new_time<250)]
        t = df[(df.new_time>=0)&(df.mdSAM==mdSAM)&(df.exp==ex)].new_time
        au = df[(df.new_time>=0)&(df.mdSAM==mdSAM)&(df.exp==ex)].Full_norm
        Coefs_one, covMatrix = opt.curve_fit(tau_function,t,au,maxfev=5000)
        Coefs_two, covMatrix = opt.curve_fit(tau_function,t,au,maxfev=5000)
        # residual sum of squares
        ss_res_1 = np.sum((au -t.apply(lambda x: tau_function(x,*Coefs_one)))**2)

        # total sum of squares
        ss_tot_1 = np.sum((au - np.mean(au)) ** 2)

        # r-squared
        r2_1 = 1 - (ss_res_1 / ss_tot_1)
        tau.append({'mdSAM':mdSAM,'Tau':1/Coefs_one[-1],'I0':Coefs_one[0],'a':Coefs_one[1],'b':Coefs_one[-1],'r2':r2_1,'exp':ex})
#     for ex in set(df[(df.new_time>=0)&(df.mdSAM==mdSAM)].exp):
#         fit_df = df[(df.new_time>=0)&(df.mdSAM==mdSAM)&(df.exp==ex)&(df.new_time<250)]
#         t = df[(df.new_time>=0)&(df.mdSAM==mdSAM)&(df.exp==ex)].new_time
#         au = df[(df.new_time>=0)&(df.mdSAM==mdSAM)&(df.exp==ex)].Full_norm

#         Coefs_two, covMatrix = opt.curve_fit(tau_function_two_param,t,au,maxfev=5000)
#         # residual sum of squares
#         ss_res_2 = np.sum((au -t.apply(lambda x: tau_function_two_param(x,*Coefs_two)))**2)

#         # total sum of squares
#         ss_tot_2 = np.sum((au - np.mean(au)) ** 2)

#         # r-squared
#         r2_2 = 1 - (ss_res_2 / ss_tot_2)
#         tau.append({'mdSAM':mdSAM,'Tau':np.NaN,'I0':Coefs_two[0],'a':Coefs_two[1],'b':Coefs_two[-1],'g':Coefs_two[1],'d':Coefs_two[-1],'r2':r2_2,'exp':ex})


