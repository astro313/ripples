'''

Prepare files for RIPPLES

1) a binary file called u.bin containing the u (in meters)
2) a binary file called v.bin containing the v (in meters)
3) a binary file called vis_chan_0.bin containing the visibilities. For each visibility first the real part and then the imaginary part. So if you have 4 visibilities, for example, this would be length 8 vector of real_1, image_1, real_2, imag_2, real_3 ... etc.
4) a binary file called sigma_squared_inv.bin similar in format to vis_chan_0.bin containing the inverse of errors squared (we'll have to sit together for me to tell you how to compute the errors though).
5) a binary file called ant1.bin containing the names of 1st antennas (they should be indexed, like 0, 1, 2 , etc.) ant1 has the same length as u.bin indicating which antenna was used for that visibility.
6) another file called ant2.bin for the second antenna.
7) a file called frequencies.bin containing the frequencies of the visibilities in Hz


Need updates based on recent changes made to ripples_prep.py


'''

try:
    from casa import table as tb
except:
    from taskinit import tb
import numpy as np
import math

RAD2DEGREE = 180. / math.pi
DEGREE2RAD = math.pi / 180.
np.set_printoptions(threshold='nan')


def make_line_cube_nobinning(lineVism, lineimagename):
    """
    Make a line cube without binning to decide on channel ranges to split
    """

    mode = 'velocity'
    imagermode = 'csclean'
    cell = '0.0300arcsec'
    imsize = [800, 800]
    restfreq = '151.80526GHz'     # from proposal

    niter = 0
    threshold = '0.0mJy'
    interactive = False

    mode = 'channel'
    start = ''
    nchan = -1
    weighting = 'briggs'
    robust = 0.5

    for ext in ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual', '.flux.pbcoverage']:
        rmtables(lineimagename + ext)

    clean(vis=lineVis,
          imagename=lineimagename,
          spw='',
          mode=mode,
          start=start,
          nchan=nchan,
          restfreq=restfreq,
          imsize=imsize,
          cell=cell,
          niter=niter,
          threshold=threshold,
          interactive=interactive,
          imagermode=imagermode,
          stokes='I',
          weighting=weighting,
          robust=robust)

    exportfits(lineimagename+'.image', lineimagename+'.fits', overwrite=True)

    return lineimagename+'.fits'


def shift_phs_center(vis, pcd, field):
    """

    Parameters
    ----------
    vis: str
        MS file
    pcd: str
        in format like: '10h27m51.6s -43d54m18s'
    field: str

    Returns
    -------
    outvis: str
        filename of the output vis, with the shifted phase center

    Example
    -------
    fixvis(vis=vis,
           outpuvis='ngc3256-fixed.ms',
           field='NGC3256',
           phasecenter='J2000 10h27m51.6s -43d54m18s')
    """

    outvis = vis.replace('.ms', 'phsShift.ms')
    fixvis(vis=vis, outputvis=outvis, field=field,
           phasecenter='J2000 ' + pcd)

    return outvis


def ms2ripples(vis, savepath, timebinsec, Nsam=None, timebin=False, overwrite=True, verbose=True, debug=False):

    import astropy.constants as c

    if timebin:
        print("Performing time binning, by {:}s:".format(timebinsec))
        prefix = vis.replace('.ms', '_timebin' + '{:}s.ms'.format(timebinsec))
        prefix = os.path.basename(prefix)
        timebinVis = os.path.join(savepath, prefix)
        if overwrite:
            rmtables(timebinVis)
        default(split2)
        split2(vis=vis,
               timebin=str(timebinsec)+'s',
               outputvis=timebinVis,
               datacolumn='data',
               keepflags=False
            )
        if overwrite:
            plotms(vis=timebinVis, xaxis='uvdist', yaxis='amp', coloraxis='spw')
            if debug and not os.path.exists(timebinVis[:timebinVis.find('.ms')] + '.image'):
                clean(vis=timebinVis,
                  imagename=timebinVis[:timebinVis.find('.ms')],
                  spw='',
                  mode='mfs',
                  nchan=-1,
                  imsize=800,
                  cell='0.0300arcsec',
                  niter=0,
                  interactive=False,
                  stokes='I')
        oldvis = vis
        vis = timebinVis
    else:
        raw_input("No time binning...proceed with caution. Press Enter to continue.")

    # check that the MS has only 1 science target
    tb.open(timebinVis + '/FIELD')
    src = tb.getcell('NAME')
    print("Source in {}: {}\n").format(timebinVis, src)
    tb.close()

    ms.open(timebinVis)
    spwInfo = ms.getspectralwindowinfo()
    nchan = spwInfo["0"]["NumChan"]
    npol = spwInfo["0"]["NumCorr"]
    ms.close()

    tb.open(oldvis)
    _dataShape = tb.getcol('DATA').shape
    tb.close()

    tb.open(timebinVis)
    # tb.colnames
    data = tb.getcol('DATA')
    uvw = tb.getcol('UVW')
    uvwShape = uvw.shape
    nvis = len(uvw[1, :])       # before dropping flagged data
    flagRow = tb.getcol('FLAG_ROW')
    assert (flagRow == True).any() == False
    data_desc_id = tb.getcol("DATA_DESC_ID")
    sigma = tb.getcol('SIGMA')
    # weight = tb.getcol('WEIGHT')
    weight = 1./sigma**2
    # del sigma

    if Nsam is not None:
        from ripples_utils import pick_highSNR
        idx = pick_highSNR(weight, Nsam, plothist=verbose)
        uvw = uvw[:, idx]

    with open(os.path.join(savepath, 'u.bin'), 'wb') as f:
        f.write(uvw[0, :])
    with open(os.path.join(savepath, 'v.bin'), 'wb') as f:
        f.write(uvw[1, :])

    if debug:
        print uvw[0, :].max()
        print uvw[0, :].min()
        print uvw[1, :].max()
        print uvw[1, :].min()
        print uvw.dtype        # float64
    # del uvw

    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    assert len(ant1) == nvis
    assert len(ant2) == nvis
    tb.done()

    if Nsam is not None:
        ant1 = ant1[idx]
        ant2 = ant2[idx]
        assert len(ant1) == len(uvw[0, :])
    with open(os.path.join(savepath, 'ant1.bin'), 'wb') as f:
        ant1 = np.asarray(ant1, dtype=np.float_)
        # print ant1.dtype
        f.write(ant1)
    with open(os.path.join(savepath, 'ant2.bin'), 'wb') as f:
        ant2 = np.asarray(ant2, dtype=np.float_)
        f.write(ant2)

    # ant1.tofile(os.path.join(savepath, 'ant1.bin'))

    if debug:
        print ant1.max()
        print ant1.min()
        print ant2.max()
        print ant2.min()
        xx = np.fromfile(os.path.join(savepath, 'ant1.bin'))
        if Nsam is None:
            assert len(xx) == nvis
        else:
            assert len(xx) == len(idx)
        assert (xx == ant1).all()
        xx = np.fromfile(os.path.join(savepath, 'ant2.bin'))
        if Nsam is None:
            assert len(xx) == nvis
        else:
            assert len(xx) == len(idx)
        assert (xx == ant2).all()
    # del ant1, ant2

    if verbose:
        print("Number of coorelation: ", npol)
        print("data shape", data.shape)
        print("data shape before time binning", _dataShape)
        print("uvw shape", uvwShape)
        print("weight shpae", weight.shape)     # (npol, nvis)

    tb.open(vis + '/SPECTRAL_WINDOW')
    SPWFreqs = np.squeeze(tb.getcol("CHAN_FREQ"))
    tb.done()

    freq_per_vis = np.array([SPWFreqs[fff] for fff in data_desc_id])
    # freqs = np.mean(SPWFreqs)
    assert len(freq_per_vis) == nvis

    if Nsam is not None:
        freq_per_vis = freq_per_vis[idx]
    with open(os.path.join(savepath, 'frequencies.bin'), 'wb') as f:
        f.write(freq_per_vis)
    del data_desc_id, SPWFreqs

    if Nsam is not None:
        data = data[:, :, idx]
        weight = weight[:, idx]

    if debug:
        # randomly sample 1000 to image
        from ripples_utils import calc_img_from_vis
        _idx = np.random.choice(len(weight[0, :]), size=3000)
        if npol == 2:
            if data.shape[1] == 1:
                # expand the channel axis for weight
                __weight = weight[:, np.newaxis, :]
            _data = np.average(data, weights=__weight, axis=0)
        _weight = np.average(weight, axis=0)    # npol, nvis

        print _weight.shape
        print _data.shape
        _real = _data[0, _idx].real    # nchan, nvis
        _imag = _data[0, _idx].imag
        visOut = np.array(zip(_real, _imag)).flatten()
        _weight = _weight[_idx]
        __weight_factor = len(_idx)/np.sum(_weight * _real**2 + _weight * _imag**2)
        print __weight_factor
        _weight *= __weight_factor
        test_img = calc_img_from_vis(uvw[0, _idx], uvw[1, _idx], _weight, visOut, freq_per_vis[_idx], 800, pixsize=0.01)
        import pdb; pdb.set_trace()

    if npol == 1:
        real = data.real
        imag = data.imag
    elif npol == 2:
        print "Real and Imag shape before averaging two hands (npol, nchan, nvis): ", data.real.shape
        if data.shape[1] == 1:
            # expand the channel axis for weight
            weight = weight[:, np.newaxis, :]
        else:
            print weight.shape
            print("We shouldn't have to enter this condition.")
            import pdb; pdb.set_trace()

        # average the two hands
        data = np.average(data, weights=weight, axis=0)
        # print data.shape      (nchan, nvis)
        weight = np.average(weight, axis=0)
        print weight.shape  # should be (nchan, nvis)

        real = data.real
        imag = data.imag
        del data
        print "Shape after averaging two hands: ", real.shape
    elif npol > 2:
        raise NotImplementedError("more than 2 hands..")

    # rescale weight
    # uvmcmcfit way
    if Nsam is None:
        _factor = nvis/np.sum(weight * real**2 + weight * imag**2)
    else:
        _factor = len(idx)/np.sum(weight * real**2 + weight * imag**2)
    _sigmas = (weight**-0.5) * _factor
    if verbose:
        print "simple rescale, factor of: ", _factor
        print _sigmas.min(), _sigmas.max(), _sigmas.std()
        print "New sigma in [mjy/beam]", (_sigmas**-2).sum()**-0.5*1.e3
        del _sigmas

    # Yashar way
    # first grouping the visibilities into bins that probe the same signal
    # take differences btw visibilities that null the sky
    # then, we simply assume that the variance in the subtracted visibilities is equal to twice the noise variance
    plotms(timebinVis, xaxis='V', yaxis='U', coloraxis='baseline')

    scaling = []
    for a1 in np.unique(ant1):
          for a2 in np.unique(ant2):
                if a1 < a2:
                    baselineX = (ant1 == a1) & (ant2 == a2)

                    if debug:
                        print a1, a2
                        print ant1, ant2
                        print ""
                        print a1 in ant1
                        print a2 in ant2
                        print ""

                        print np.where(a1 == ant1)
                        print np.where(a2 == ant2)
                        print np.where((ant1 == a1) & (ant2 == a2) == True)

                    if baselineX.any() == True:    # important line! if we are picking a subset of points with nsam since we may miss some baselines.
                        if nchan == 1:
                            print real.shape
                            reals = real[0, baselineX]
                            imags = imag[0, baselineX]
                            sigs = weight[0, baselineX]**-0.5
                        else:
                            raise NotImplementedError("Not implemented for MS files with more than 1 channel per spw...")

                        # randomly split points into "two sets"
                        # subtract from "two set"
                        # subtract from neighboring
                        diffrs = reals - np.roll(reals, -1)
                        diffis = imags - np.roll(imags, -1)
                        std = np.mean([diffrs.std(), diffis.std()])

                        if debug:
                            print diffrs.min(), diffis.min()
                            print diffrs.max(), diffis.max()
                            print diffrs.std(), diffis.std()
                            print std / sigs.mean() / np.sqrt(2)
                        scaling.append(std / sigs.mean() / np.sqrt(2))
    del ant1, ant2
    sigma = weight**-0.5
    scaling = np.asarray(scaling).mean()
    sigma *= scaling
    print 'Scaling factor: ', scaling
    print 'Sigma after scaling [mJy/beam]: ', ((sigma**-2).sum())**-0.5 * 1E+3

    # real_1, imag_1, real_2, imag_2, etc
    visOut = np.array(zip(real, imag)).flatten()
    if Nsam is None:
        assert len(visOut) == int(nvis*2)
    weight = sigma**-2
    weight = np.array(zip(weight, weight)).flatten()    # TODO: want to make sure this works for nchan > 1
    assert len(weight) == len(visOut)

    if Nsam is not None:
        assert len(visOut) == Nsam * 2
    if Nsam is not None:
        assert len(weight) == Nsam * 2
    with open(os.path.join(savepath, 'vis_chan_0.bin'), 'wb') as f:
        f.write(visOut)
    with open(os.path.join(savepath, 'sigma_squared_inv.bin'), 'wb') as f:
        f.write(weight)

    if Nsam is None:
        blah = np.zeros((nvis))
    else:
        blah = np.zeros((len(idx)))
    with open(os.path.join(savepath, 'chan.bin'), 'wb') as f:
        f.write(blah)

    if debug:
        # image the SAVED visibilities
        uu = np.fromfile(os.path.join(savepath, 'u.bin'), 'wb')
        vv = np.fromfile(os.path.join(savepath, 'v.bin'), 'wb')
        weight = np.fromfile(os.path.join(savepath, 'sigma_squared_inv.bin'), 'wb')
        visout = np.fromfile(os.path.join(savepath, 'vis_chan_0.bin'), 'wb')
        calc_img_from_vis(uu, vv, weight, visOut, 1000, pixsize=0.05)

    return None


if __name__ == '__main__':
    import sys
    sys.path.append("./")

    from ripples_prep import get_source_name, get_phs_center, ms2ripples
    import config
    from ripples_utils import get_nchan, vel2chan
    # import subprocess

    # For line, do it for different channels, averaged, etc etc
    ccc = config.Configurable('ripples_prep.yaml')
    setupParam = ccc.config_dict
    savepath = setupParam['line']['savepath']
    overwrite = setupParam['line']['overwrite']
    lineMS = setupParam['line']['data']['lineMS']    # cube, after subtracting continuum
    unbinnedIm = setupParam['line']['data']['unbinnedIm']
    sliceChan = setupParam['line']['slice']
    shiftpcd = setupParam['shiftpcd']['shiftpcd']
    timebin = setupParam['line']['timebin']
    timebinsec = setupParam['line']['timebinsec']
    nsam = setupParam['cont']['nsam']


    if not os.path.exists(savepath):
        os.mkdir(savepath)
    _current = os.getcwd()

    if shiftpcd:
        field = setupParam['shiftpcd']['field']
        pcd = setupParam['shiftpcd']['pcd']
        outvis = shift_phs_center(lineMS, pcd, field)
        outpcd = get_phs_center(outvis)
        print "Phase Center to be shifted to", pcd
        print "Phase Center after fixvis: ", outpcd
        lineMS = outvis

    # if already shifted phase center, then use the MS file, regardless of whether we are setting shiftpcd in .yaml
    if not shiftpcd:
        outvis = lineMS.replace('.ms', 'phsShift.ms')
    if os.path.exists(outvis):
        lineMS = outvis

    # find ranges of channels to split in unbinned data cube
    if not os.path.exists(unbinnedIm + '.fits'):
        lineCubeFITS_unbinned = make_line_cube_nobinning(lineMS)
    else:
        lineCubeFITS_unbinned = unbinnedIm + '.fits'

    # split into different channel slices & time bin
    chanbin, _ = get_nchan(lineCubeFITS_unbinned)

    for ns, chan in sliceChan.iteritems():
        # change path
        _p = savepath + ns + '/'
        if not os.path.exists(_p):
            os.mkdir(_p)
        os.chdir(_p)
        outvis = 'lineSlice.ms'

        # split channel
        default(split2)
        rmtables(outvis)
        split2(vis=lineMS,      # absolute path
               spw=chan,
               outputvis=outvis,
               datacolumn='data',
               keepflags=False,
               width=chanbin
            )

        listobs(outvis)
        ms2ripples(contMS, savepath, timebinsec, nsam, timebin, \
               overwrite=overwrite, verbose=True, debug=False)        # if debug=True, will also image the saved visibilities using python.

        os.chdir(_current)

    objname = get_source_name(lineMS)
    tarballName = objname + '_line_timebin' + str(timebinsec) +'s.tgz'

    # change to ripples dir
    os.chdir(savepath)
    os.system('tar -zcvf ' +  tarballName + ' ' + '$(find -name "*.bin"')
    os.chdir(_current)
    print tarballName
    # On Voms: scp dleung@serenity.astro.corell.edu:/data/dleung/DATA/ALMA/QSO_c5/RXJ0911/calibrated/ripples/RXJ0911_line_timebin60s.tgz .
#