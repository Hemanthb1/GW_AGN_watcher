import pickle
import os
import numpy as np
import pandas as pd
import healpy as hp
from glob import glob
from tqdm import tqdm
from gwpy.table import Table
from sklearn.ensemble import RandomForestClassifier

"""
Written by Isaac McMahon, 2026-02-22
Based on the Ordinal Classification photo-z method from https://doi.org/10.1093/mnras/stv1567
"""

def run_simulation(n_events, outdir='simulation_output', min_mass=1.2, max_mass=75.0, max_distance=6000.0, snr_threshold=8.0, downsample_nside=8):
    """
    Run a set of BAYESTAR simulations to use as training data. More parameters in the injection script can be changed,
    but are not implemented here. More information at https://lscsoft.docs.ligo.org/ligo.skymap/quickstart/bayestar-injections.html

    Parameters
    ----------
    n_events : int
        Number of CBCs to simulate
    outdir : str
        Name of the output directory to create
    min_mass : float, optional
        Minimum CBC component mass
    max_mass : float, optional
        Maximum CBC component mass
    max_distance : float, optional
        Maximum CBC luminosity distance in Mpc
    snr_threshold : float, optional
        Network SNR threshold for IFO detection
    downsample_nside : int, optional
        Skymap resolution for the mass estimation model
    """
    gps_start = 1000000000
    try:
        os.mkdir(outdir)
    except OSError as e:
        print(f'Output directory already exists', e)
        return
    os.system(f'''lalapps_inspinj --output {outdir}/inj.xml --f-lower 25 --waveform TaylorF2threePointFivePN \
        --t-distr uniform --time-step 1 --gps-start-time {gps_start} --gps-end-time {int(gps_start+n_events)} \
    	--m-distr log --min-mass1 {min_mass} --max-mass1 {max_mass} --min-mass2 {min_mass} --max-mass2 {max_mass} \
    	--d-distr volume --min-distance 1 --max-distance {int(max_distance)}e3 \
    	--l-distr random --i-distr uniform --enable-spin --min-spin1 0.0 --max-spin1 1.0 --min-spin2 0.0 --max-spin2 1.0 --aligned --disable-milkyway''')
    os.system(f'bayestar-sample-model-psd -o {outdir}/psd.xml --H1=aLIGO175MpcT1800545 --L1=aLIGO175MpcT1800545 --V1=aLIGOAdVO4T1800545')
    os.system(f'''bayestar-realize-coincs -o {outdir}/coinc.xml {outdir}/inj.xml --reference-psd {outdir}/psd.xml --detector H1 L1 V1 \
    	--measurement-error gaussian-noise --snr-threshold 2.0 --net-snr-threshold {snr_threshold} --min-triggers 2''')
    os.system(f'bayestar-localize-coincs -o {outdir}/ {outdir}/coinc.xml')
    os.system(f'ligolw_sqlite --preserve-ids --replace --database {outdir}/coinc.sqlite {outdir}/coinc.xml')
    os.system(f'ligo-skymap-stats -o {outdir}/bayestar.tsv --database {outdir}/coinc.sqlite {outdir}/*.fits --contour 50 90 --area 10 100')

    ds_dir = f'{outdir}/downsamples_nside{downsample_nside}'
    os.mkdir(ds_dir)
    for fitsfile in tqdm(glob(outdir+'/*.fits')):
        id = fitsfile.split('/')[-1].split('.')[0]
        os.system(f'ligo-skymap-flatten --nside {downsample_nside} {fitsfile} {ds_dir}/{id}_downsample_nside{downsample_nside}.fits.gz')
        os.remove(fitsfile)
    print('Simulation Complete')
    return

class MassEstimator():
    """
    Class to handle the mass classifier

    Parameters
    ----------
    data_npix : int, optional
        Number of pixels to use in the input vector of the classifier training
    downsample_nside : int, optional
        Skymap resolution for the mass estimation model
    """
    def __init__(self, data_npix=20, downsample_nside=8):
        self.data_npix = data_npix
        self.downsample_nside = downsample_nside

    def load_simulation_data(self, sim_dirs, det_data_file=None, inj_data_file=None, snr_threshold=8, overwrite=True):
        """
        Load simulation data from directories made by run_simulation as training data
    
        Parameters
        ----------
        sim_dirs : list
            List of paths to simulation output directories
        det_data_file : str, optional
            Path to existing detected simulation data, or name of file to write data to
        inj_data_file : str, optional
            Path to existing injected simulation data, or name of file to write data to
        snr_threshold : float, optional
            Network SNR threshold for IFO detection
        overwrite : bool, optional
            Flag to overwrite existing files
        """
        if det_data_file != None and inj_data_file != None:
            det = pd.read_csv(det_data_file)
            inj = pd.read_csv(inj_data_file)
        else:
            print('Both output files not provided, loading simulations')
            for i in tqdm(range(len(sim_dirs))):
                if not i:
                    sim = Table.read(sim_dirs[i]+'/coinc.xml', tablename='sim_inspiral:table').to_pandas()[['distance', 'mass1', 'mass2']]
                    inj = Table.read(sim_dirs[i]+'/inj.xml', tablename='sim_inspiral:table').to_pandas()[['mass1', 'mass2', 'distance']]
                    stats = pd.read_csv(sim_dirs[i]+'/bayestar.tsv', sep='\t', skiprows=1).sort_values('simulation_id')[['snr', 'distmean', 'diststd', 'area(50)', 'area(90)']].reset_index(drop=True)
                else:
                    sim = pd.concat([sim, Table.read(sim_dirs[i]+'/coinc.xml', tablename='sim_inspiral:table').to_pandas()[['distance', 'mass1', 'mass2']]], ignore_index=True)
                    inj = pd.concat([inj, Table.read(sim_dirs[i]+'/inj.xml', tablename='sim_inspiral:table').to_pandas()[['mass1', 'mass2', 'distance']]], ignore_index=True)
                    stats = pd.concat([stats, pd.read_csv(sim_dirs[i]+'/bayestar.tsv', sep='\t', skiprows=1).sort_values('simulation_id')[['snr', 'distmean', 'diststd', 'area(50)', 'area(90)']].reset_index(drop=True)], ignore_index=True)
        
            det = pd.concat([sim[['mass1', 'mass2', 'distance']], stats[['distmean', 'diststd', 'area(50)', 'area(90)', 'snr']]], axis=1).rename(columns={'distance':'dist_true', 'area(50)':'area_50', 'area(90)':'area_90'})
            if overwrite:
                det.to_csv(det_data_file, index=False)
                inj.to_csv(inj_data_file, index=False)

        gw_input = []
        pass_mask = []
        index = -1
        for path in sim_dirs:
            print(f'Loading {path}...', flush=True)
            for fitsfile in tqdm(glob(path+'/downsamples/*_downsample_nside*.fits.gz')):
                index += 1
                prob, distnorm, distmu, distsig = hp.read_map(fitsfile, (0, 1, 2, 3))
                finite_mask = np.isfinite(distmu)
                if np.sum(finite_mask) < self.data_npix:
                    pass_mask.append(False)
                    continue
                pix_id_masked = np.flip(np.argsort(prob[finite_mask]))[:self.data_npix]
                pix_id = np.arange(len(prob))[finite_mask][pix_id_masked]
                pix_ra, pix_dec = hp.pix2ang(self.downsample_nside, pix_id, lonlat=True)
                pix_prob = prob[pix_id]
                pix_dist = distmu[pix_id]
                pix_sig = distsig[pix_id]
                data_vector = np.concatenate((pix_ra, pix_dec, pix_prob, pix_dist, pix_sig))
                if len(data_vector)!=5*self.data_npix:
                    print(np.sum(finite_mask))
                gw_input.append(data_vector)
                pass_mask.append(True)
        
        self.chirp_inject = ((inj['mass1'].to_numpy()*inj['mass2'].to_numpy())**(0.6))*(inj['mass1'].to_numpy()+inj['mass2'].to_numpy())**(-0.2)
        gw_output = ((det['mass1'].to_numpy()*det['mass2'].to_numpy())**(0.6))*(det['mass1'].to_numpy()+det['mass2'].to_numpy())**(-0.2)
        print('Number events after skymap loading:', len(self.gw_output))
        
        self.pass_mask = np.array(pass_mask)
        snr_mask = (det['snr']>=snr_threshold)[self.pass_mask]
        self.gw_output = gw_output[self.pass_mask][snr_mask]
        self.gw_input = np.array(gw_input)[snr_mask]
        print('Number events after SNR cut:', len(self.gw_output))

    def train(self, model_path, nbins=240, mchirp_max=120., save=True):
        """
        Train Random Forest Classifier using simulated training data
    
        Parameters
        ----------
        model_path : str
            Path to existing pre-trained model data, or name of file to write model to
        nbins : int, optional
            Number of bins (classes) used in the classifier
        mchirp_max : float, optional
            Maximum source chirp mass fitted by the classifier
        save : bool, optional
            Flag to save trained model
        """
        model_path_file = f'{model_path}.dat'
        if os.path.exists(model_path_file):
            with open(model_path_file, 'rb') as f:
                self.clf = pickle.load(f)
            print('Model Loaded')
        else:
            grid = np.linspace(0, mchirp_max, num=nbins, endpoint=True)
            midpoints = grid[:-1] + mchirp_max/nbins
            
            train_bins = np.digitize(self.gw_output, grid)
            missing_no = np.setdiff1d(np.arange(1, nbins), train_bins)
            self.clf = RandomForestClassifier().fit(self.gw_input, train_bins)
            self.clf.missing_bins = missing_no
            self.clf.m_array = midpoints
            print('Training Complete')
            if save:
                with open(model_path_file, 'wb') as f:
                    pickle.dump(self.clf, f)

    def predict_mass(self, downsample_skymap):
        """
        Predict the chirp mass of a single event using a trained classifier
    
        Parameters
        ----------
        dowmsaple_skymap : str
            Path to skymap of event to be estimated.
            Skymap must be flat resolution and downsampled to the correct NSIDE according to the trained model.
            Use ligo-skymap-flatten to do this first.

        Returns
        ----------
        pdf : numpy.ndarray
            Discrete probability distribution for the estimated event
        midpoints : numpy.ndarray
            Midpoints of the chirp mass bins in Solar Masses
        """
        prob, distnorm, distmu, distsig = hp.read_map(downsample_skymap, (0, 1, 2, 3))
        finite_mask = np.isfinite(distmu)
        pix_id_masked = np.flip(np.argsort(prob[finite_mask]))[:self.data_npix]
        pix_id = np.arange(len(prob))[finite_mask][pix_id_masked]
        pix_ra, pix_dec = hp.pix2ang(self.downsample_nside, pix_id, lonlat=True)
        pix_prob = prob[pix_id]
        pix_dist = distmu[pix_id]
        pix_sig = distsig[pix_id]
        data_vector = np.concatenate((pix_id, pix_prob, pix_dist, pix_sig))
        test_event = [data_vector]
        
        pdf = self.clf.predict_proba(test_event)[0]
        if self.clf.missing_bins.shape[0]:
            for c in self.clf.missing_bins:
                pdf = np.insert(pdf, c-1, 0.)
        return pdf, self.clf.midpoints

