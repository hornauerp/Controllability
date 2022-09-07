import numpy as np
from elephant.statistics import instantaneous_rate
from elephant.kernels import GaussianKernel
from neo.core.spiketrainlist import SpikeTrainList
import quantities as pq
import os


def smooth_fr(input_path, output_path):
    times_path = os.path.join(input_path, 'spike_times.npy')
    ids_path = os.path.join(input_path, 'spike_templates.npy')
    times = np.load(times_path)
    times /= 20000
    ids = np.load(ids_path)
    # get neotrains
    neo_t: SpikeTrainList = SpikeTrainList.from_spike_time_array(times, ids, np.unique(ids), units='s',
                                                                 t_stop=np.max(times))

    # get smoothed data
    smoothed = instantaneous_rate(neo_t, sampling_period=1 * pq.ms, kernel=GaussianKernel(3 * pq.ms))
    np_array = np.array(smoothed).T
    np.save(output_path, np_array)

    return np_array
