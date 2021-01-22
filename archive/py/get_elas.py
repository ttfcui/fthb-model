
from model_iterate import lifecycle_iterate
from numpy import log
baseFile = 'lifecycle_TC03242016.f90'

base = lifecycle_iterate({'elasticity': 2}, baseFile)
#base.execSh(model='base')
#
#high_gamma = lifecycle_iterate({'elasticity': 2.02}, baseFile)
#high_gamma.execSh(model='high_gamma')
#
#low_gamma = lifecycle_iterate({'elasticity': 1.98}, baseFile)
#low_gamma.execSh(model='low_gamma')
#
#high_s = lifecycle_iterate({'sigma_z': 0.052621}, baseFile)
#high_s.execSh(model='high_s')
#
#low_s = lifecycle_iterate({'sigma_z': 0.051579}, baseFile)
#low_s.execSh(model='low_s')
#
#no_s = lifecycle_iterate({'sigma_z': 0}, baseFile)
#no_s.execSh(model='no_s')
#
#high_alpha = lifecycle_iterate({'elasticity2O': 0.8514}, baseFile)
#high_alpha.execSh(model='high_alpha')
#
low_alpha = lifecycle_iterate({'elasticity2O': 0.8686}, baseFile)
low_alpha.execSh(model='low_alpha')
#
#low_theta = lifecycle_iterate({'thetamatlab': 0.2475}, baseFile)
#low_theta.execSh(model='low_theta')
#
#high_theta = lifecycle_iterate({'thetamatlab': 0.2525}, baseFile)
#high_theta.execSh(model='high_theta')
#
#high_rrental = lifecycle_iterate({'r_rentalO': 0.076255}, baseFile)
#high_rrental.execSh(model='high_rrental')
#
#low_rrental = lifecycle_iterate({'r_rentalO': 0.074745}, baseFile)
#low_rrental.execSh(model='low_rrental')
#
#med_assets = lifecycle_iterate({'currenthouseholdstate(:,1)': 0.4}, baseFile)
#med_assets.execSh(model='med_assets')
#
#high_assets = lifecycle_iterate({'currenthouseholdstate(:,1)': 0.404}, baseFile)
#high_assets.execSh(model='high_assets')
#
#med_durables = lifecycle_iterate({'currenthouseholdstate(:,2)': 1}, baseFile)
#med_durables.execSh(model='med_durables')
#
#high_durables = lifecycle_iterate({'currenthouseholdstate(:,2)': 1.01}, baseFile)
#high_durables.execSh(model='high_durables')


allM = ['base', 'low_s', 'high_s', 'low_gamma', 'high_gamma', 'low_alpha', 
        'high_alpha', 'low_theta', 'high_theta', 'low_rrental', 'high_rrental',
        'med_assets', 'high_assets', 'med_durables', 'high_durables', 'no_s']
base.appendModels('fthb_dist', model=allM)
out = base.summaryTable('ageBought')

for name in allM[1:-5]:
    print 'Model: ' + name
    print 100*(log(out.loc['mean',name]) - log(out.loc['mean', 'base']))
    print 100*(log(out.loc['50%',name]) - log(out.loc['50%', 'base']))

print 'Model: assets'
print 100*(log(out.loc['mean','high_assets']) - log(out.loc['mean', 'med_assets']))
print 100*(log(out.loc['50%','high_assets']) - log(out.loc['50%', 'med_assets']))

print 'Model: durables'
print 100*(log(out.loc['mean','high_durables']) - log(out.loc['mean', 'med_durables']))
print 100*(log(out.loc['50%','high_durables']) - log(out.loc['50%', 'med_durables']))
