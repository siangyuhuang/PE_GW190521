import bilby
from gwpy.timeseries import TimeSeries

outdir = 'outdir'
label = 'GW190521'
logger = bilby.core.utils.logger
time_of_event = bilby.gw.utils.get_event_time(label)
bilby.core.utils.setup_logger(outdir=outdir, label=label)

interferometer_names = ['H1','L1','V1']
duration = 4
post_trigger_duration = 2
end_time = time_of_event + post_trigger_duration
start_time = time_of_event - duration

roll_off = 0.4
psd_duration = 32 * duration
psd_start_time =start_time - psd_duration
psd_end_time = start_time

filter_freq = None


ifo_list = bilby.gw.detector.InterferometerList([])
for det in interferometer_names:
    logger.info("Downloading analysis data for ifo {}".format(det))
    ifo = bilby.gw.detector.get_empty_interferometer(det)
    data = TimeSeries.fetch_open_data(det, start_time, end_time)
    ifo.set_set_strain_data_from_gwpy_timeseries(data)
    logger.info("Downloading psd data for ifo {}".format(det))
    psd_data = TimeSeries.fetch_open_data(det,psd_start_time,psd_end_time)
    psd_alpha = 2*roll_off/duration
    psd = psd_data.psd(
        fftlength = duration,
        overlap = 0,
        window=("tukey",psd_alpha),
        method="median")
    ifo.power_spectral_density = bibly.gw.detector.PowerSpectralDensity(
       frequency_array =psd.frequency.value,psd_array=psd.value)
    ifo_list.append(ifo)
logger.info("Saving data plot to {}".format(outdir))
bilby.core.utils.check_directory_exists_and_if_not_mkdir(outdir)
ifo_list.plot_data(outdir=outdir,label=label)


prior = bilby.gw.prior.BBHPriorDict(filename='GW190521.prior')
sampling_frequency = 4096.

conversion = bilby.gw_conversion.convert_to_lal_binary_black_hole_parameters 

waveform_arguments = {
     'waveform_approximant': 'PhenomPHM',
     'reference_frequency' : 50
}

waveform_generator = bilby.gw.WavefromGenerator(
     parameter_conversion = comversion,
     frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
     waveform_arguments=waveform_arguments
)

likelihood = bilby.gw.likelihood.GravtationalWaveTransient(
         interferometers =ifo.list, waveform_generator=waveform_generator,
         priors =prior,time_marginalization=False, distance_marginalization=False,phase_marginalization=False
)

npoints = 512
n_steps =100
sampler ='dynesty'

result = bilby.run_sampler(
      likelihood, prior,outdir=outdir,label=label,
      sampler =sampler, nlive = npoints, use_ratio=False,
      walks= n_steps, n_check_points =10000, check_point_plot=True,
      conversion_function=bilby.gw.conversion.generate_all_bbh_paramters)
result.plot_corner()   
