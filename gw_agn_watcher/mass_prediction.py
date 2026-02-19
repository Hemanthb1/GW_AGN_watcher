import pickle
import numpy as np
import pandas as pd
import healpy as hp
from glob import glob
from tqdm import tqdm
from gwpy.table import Table
from sklearn.ensemble import RandomForestClassifier

#### LOAD SIMULATED BAYESTAR EVENTS ####

rewrite = True
prefix = 'pop_o4a_hlv'

outdirs = []
for i in range(10):
    outdirs.append(f'data/{prefix}_{i}')
outdirs.append(f'data/{prefix}_bns')

if rewrite:
    for i in tqdm(range(len(outdirs))):
        if not i:
            sim = Table.read(outdirs[i]+'/coinc.xml', tablename='sim_inspiral:table').to_pandas()[['distance', 'mass1', 'mass2']]
            inj = Table.read(outdirs[i]+'/inj.xml', tablename='sim_inspiral:table').to_pandas()[['mass1', 'mass2', 'distance']]
            stats = pd.read_csv(outdirs[i]+'/bayestar.tsv', sep='\t', skiprows=1).sort_values('simulation_id')[['snr', 'distmean', 'diststd', 'area(50)', 'area(90)']].reset_index(drop=True)
        else:
            sim = pd.concat([sim, Table.read(outdirs[i]+'/coinc.xml', tablename='sim_inspiral:table').to_pandas()[['distance', 'mass1', 'mass2']]], ignore_index=True)
            inj = pd.concat([inj, Table.read(outdirs[i]+'/inj.xml', tablename='sim_inspiral:table').to_pandas()[['mass1', 'mass2', 'distance']]], ignore_index=True)
            stats = pd.concat([stats, pd.read_csv(outdirs[i]+'/bayestar.tsv', sep='\t', skiprows=1).sort_values('simulation_id')[['snr', 'distmean', 'diststd', 'area(50)', 'area(90)']].reset_index(drop=True)], ignore_index=True)

    det = pd.concat([sim[['mass1', 'mass2', 'distance']], stats[['distmean', 'diststd', 'area(50)', 'area(90)', 'snr']]], axis=1).rename(columns={'distance':'dist_true', 'area(50)':'area_50', 'area(90)':'area_90'})
    det.to_csv(f'{prefix}_det.csv', index=False)
    inj.to_csv(f'{prefix}_inj.csv', index=False)
else:
    det = pd.read_csv(f'{prefix}_det.csv')
    inj = pd.read_csv(f'{prefix}_inj.csv')

gw_output = ((det['mass1'].to_numpy()*det['mass2'].to_numpy())**(0.6))*(det['mass1'].to_numpy()+det['mass2'].to_numpy())**(-0.2)
chirp_inject = ((inj['mass1'].to_numpy()*inj['mass2'].to_numpy())**(0.6))*(inj['mass1'].to_numpy()+inj['mass2'].to_numpy())**(-0.2)
print('Number of detected events:', len(gw_output))

npix = 20
gw_input_adv = []
pass_mask = []
index = -1

for path in outdirs:
    print(f'Loading {path}...', flush=True)
    for i in tqdm(range(len(glob(path+'/downsamples/*_downsample.fits.gz')))):
        index += 1
        prob, distnorm, distmu, distsig = hp.read_map(f'{path}/downsamples/{i}_downsample.fits.gz', (0, 1, 2, 3))
        finite_mask = np.isfinite(distmu)
        if np.sum(finite_mask) < npix:
            pass_mask.append(False)
            continue
        pix_id_masked = np.flip(np.argsort(prob[finite_mask]))[:npix]
        pix_id = np.arange(len(prob))[finite_mask][pix_id_masked]
        pix_ra, pix_dec = hp.pix2ang(8, pix_id, lonlat=True)
        pix_prob = prob[pix_id]
        pix_dist = distmu[pix_id]
        pix_sig = distsig[pix_id]
        data_vector = np.concatenate((pix_ra, pix_dec, pix_prob, pix_dist, pix_sig))
        if len(data_vector)!=5*npix:
            print(np.sum(finite_mask))
        gw_input_adv.append(data_vector)
        pass_mask.append(True)
pass_mask = np.array(pass_mask)
gw_input_adv = np.array(gw_input_adv)
gw_output = gw_output[pass_mask]
print('Number events after skymap loading:', len(gw_output))

snr_mask = (det['snr']>=8)[pass_mask]
chirp_inject = ((inj['mass1'].to_numpy()*inj['mass2'].to_numpy())**(0.6))*(inj['mass1'].to_numpy()+inj['mass2'].to_numpy())**(-0.2)
gw_output = ((det['mass1'].to_numpy()*det['mass2'].to_numpy())**(0.6))*(det['mass1'].to_numpy()+det['mass2'].to_numpy())**(-0.2)
gw_output = gw_output[pass_mask][snr_mask]
gw_input_adv = gw_input_adv[snr_mask]
print('Number events after SNR cut:', len(gw_output))

#### TRAIN MODEL ####

nbins = 240
mchirp_max = 120
grid = np.linspace(0, mchirp_max, num=nbins, endpoint=True)
midpoints = grid[:-1] + offset

train_bins = np.digitize(gw_output, grid)
missing_no = np.setdiff1d(np.arange(1, nbins), train_bins)
clf = RandomForestClassifier().fit(gw_input_adv, train_bins)
clf.missing_bins = missing_no
clf.midpoints = midpoints

with open(f'mass_prediction_adv_{prefix.replace('pop_', '')}.dat', 'wb') as f:
    pickle.dump(clf, f)