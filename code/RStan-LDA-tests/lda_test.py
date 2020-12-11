#!/usr/bin/python

# TODO For VI methods, calculate posterior means using the variational params
# TODO Figure out how to do nested parallel loops with pyMC3 NUTS sampling
# TODO improve saving of VI model & trace results

import numpy as np
import h5py
import os
import sys
import ntpath
import MMLDA as mods
import multiprocessing
import functools


def get_nuts_results(data_fn, model_dir, results_dir, draws, tune, nchain,
                     ncore, seed):
    data_type = os.path.splitext(ntpath.basename(data_fn))[0]
    model_fn = os.path.join(model_dir, 'lda_'+data_type+'_pyMC3_NUTS.h5')
    eval_fn = os.path.join(results_dir,
                           'lda_'+data_type+'_pyMC3_NUTS_convstats.csv')
    trace_fn = os.path.join(model_dir, 'lda_' + data_type + '_pyMC3_NUTS.p')
    if not os.path.isfile(model_fn):
        # Read data
        data = mods.data.read_test_data(data_fn)
        # Fit LDA model
        print("Fitting the model %s...\n" % model_fn)
        out = mods.lda.fit_NUTS(K=data['K'], X=data['R'], alp=data['alp'],
                                sig=data['sig'], eval_fn=eval_fn, draws=draws,
                                tune=tune, nchain=nchain, ncore=ncore,
                                seed=seed, trace_fn=trace_fn)
        # Save output
        hf = h5py.File(model_fn, 'w')
        hf.create_dataset('theta', data=out['theta'].T)
        hf.create_dataset('beta', data=out['beta'].T)
        hf.create_dataset('rt', data=out['rt'])
        hf.create_dataset('data_fn', data=data_fn)
        hf.create_dataset('trace_fn', data=trace_fn)
        hf.close()


def get_results(data_fn, model_dir, results_dir, model_type, samp_size, seed,
                iters):
    data_type = os.path.splitext(ntpath.basename(data_fn))[0]
    model_fn = os.path.join(model_dir,
                            'lda_'+data_type+'_pyMC3_'+model_type+'.h5')
    trace_fn = os.path.join(model_dir,
                            'lda_'+data_type+'_pyMC3_'+model_type+'.p')
    elbo_fn = os.path.join(results_dir,
                           'lda_'+data_type+'_pyMC3_'+model_type+'_elbo.pdf')
    if not os.path.isfile(model_fn):
        # Read data
        data = mods.data.read_test_data(data_fn)
        # Fit LDA model
        print("Fitting the model %s...\n" % model_fn)
        out = mods.lda.fit_vi(model_type=model_type, K=data['K'], X=data['R'],
                              alp=data['alp'], sig=data['sig'],
                              samp_size=samp_size, elbo_fn=elbo_fn,
                              trace_fn=trace_fn, iters=iters, seed=seed)
        # Save output
        hf = h5py.File(model_fn, 'w')
        hf.create_dataset('theta', data=out['theta'].T)
        hf.create_dataset('beta', data=out['beta'].T)
        hf.create_dataset('rt', data=out['rt'])
        hf.create_dataset('data_fn', data=data_fn)
        hf.create_dataset('trace_fn', data=trace_fn)
        hf.close()


def main(argv):
    # Set the seed for reporducible results
    seed = 123
    np.random.seed(seed)
    # Load parameters
    script_name = argv[0]
    data_dir = argv[1]
    model_dir = argv[2]
    results_dir = argv[3]
    nReps = int(argv[4])
    ncores = int(argv[5])
    draws = int(argv[6])
    tune = int(argv[7])
    vi_samp_size = int(argv[8])
    nchain = int(argv[9])
    vi_iters = int(argv[10])
    print("We are now running the script %s.\n" % script_name)
    print("The data directory is %s.\n" % data_dir)
    print("The model directory is %s.\n" % model_dir)
    print("The results directory is %s.\n" % results_dir)
    print("We will perform %d replications.\n" % nReps)
    print("We will use %d cores.\n" % ncores)
    print("We will have %d MCMC draws.\n" % draws)
    print("We will tune the MCMC chains with %d draws.\n" % tune)
    print("We will have %d iterations for VI methods.\n" % vi_iters)
    print("We will sample %d MCMC chains.\n" % nchain)
    # Check directory
    if not os.path.isdir(data_dir):
        raise Exception("The directory %s doesn't exist.\n" % data_dir)
    if not os.path.isdir(model_dir):
        print("Creating the directory %s.\n" % model_dir)
        os.mkdir(model_dir)
    if not os.path.isdir(results_dir):
        print("Creating the directory %s.\n" % results_dir)
        os.mkdir(results_dir)
    # List data files to process
    data_symm_fns = list(map(lambda x: os.path.join(
        data_dir, 'symmdata' + str(x) + '.h5'), range(1, nReps+1)))
    data_spar_fns = list(map(lambda x: os.path.join(
        data_dir, 'spardata' + str(x) + '.h5'), range(1, nReps+1)))
    data_fns = [*data_symm_fns, *data_spar_fns]
    # create multiprocessing pool
    pool = multiprocessing.Pool(processes=ncores)
    # Generate results for mean-field ADVI
    get_results_part = functools.partial(
        get_results, model_dir=model_dir, model_type='mfADVI', iters=vi_iters,
        samp_size=vi_samp_size, results_dir=results_dir, seed=seed)
    pool.map(get_results_part, data_fns)
    # Generate results for full-rank ADVI
    get_results_part = functools.partial(
        get_results, model_dir=model_dir, model_type='frADVI', iters=vi_iters,
        samp_size=vi_samp_size, results_dir=results_dir, seed=seed)
    pool.map(get_results_part, data_fns)
    # Generate results for full-rank ADVI
    get_results_part = functools.partial(
        get_results, model_dir=model_dir, model_type='SVGD', iters=vi_iters,
        samp_size=vi_samp_size, results_dir=results_dir, seed=seed)
    pool.map(get_results_part, data_fns)
    # Generate results for NUTS sampling
    get_nuts_results_part = functools.partial(
        get_nuts_results, model_dir=model_dir, results_dir=results_dir,
        draws=draws, tune=tune, nchain=nchain, ncore=1, seed=seed)
    pool.map(get_nuts_results_part, data_fns)

    # for fn in data_fns:
    #     get_nuts_results(data_fn=fn, model_dir=model_dir,
    #                      results_dir=results_dir, samp_size=nuts_samp_size,
    #                      nchain=nchain, ncore=1, seed=seed)


if __name__ == "__main__":
    print("\nRunning the Python script: %s.\n" % sys.argv[0])
    main(sys.argv)
    print("\nCOMPLETED LDA TESTING SUCCESSFULLY.\n")
