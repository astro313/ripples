'''

Prepare files for RIPPLES

1) a binary file called u.bin containing the u (in meters)
2) a binary file called v.bin containing the v (in meters)
3) a binary file called vis_chan_0.bin containing the visibilities. For each visibility first the real part and then the imaginary part. So if you have 4 visibilities, for example, this would be length 8 vector of real_1, image_1, real_2, imag_2, real_3 ... etc.
4) a binary file called sigma_squared_inv.bin similar in format to vis_chan_0.bin containing the inverse of errors squared (we'll have to sit together for me to tell you how to compute the errors though).
5) a binary file called ant1.bin containing the names of 1st antennas (they should be indexed, like 0, 1, 2 , etc.) ant1 has the same length as u.bin indicating which antenna was used for that visibility.
6) another file called ant2.bin for the second antenna.
7) a file called frequencies.bin containing the frequencies of the visibilities in Hz



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


def get_phs_center(msfile):
    """ get the phase center of the MS file using CASA taskinit.ms

    Parameters
    ----------
    msfile: str
        measurement set filename
    Returns
    -------
    pcd:
        phase center

    Note
    ----
    may need further modification if nField > 1 unless it is the same for all fields (need further investigation)
    """

    from taskinit import ms
    ms.open(msfile)
    pc = ms.getfielddirmeas()
    if not isinstance(pc, dict) is True:
        pc()
    epoch = pc['refer']
    pcd_dec = pc['m1']['value'] * 180 / np.pi
    pcd_ra = pc['m0']['value'] * 180 / np.pi
    if pcd_ra < 0:
        pcd_ra += 360
    ms.done()
    pcd = [pcd_ra, pcd_dec]
    return pcd


def ms2ripples(vis, savepath, timebinsec, Nsam=None, timebin=False, overwrite=True, verbose=True, debug=False):

    import astropy.constants as c
    from ripples_utils import calc_img_from_vis

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
    tb.open(vis + '/FIELD')
    src = tb.getcell('NAME')
    print("Source in {}: {}\n").format(vis, src)
    tb.close()

    ms.open(vis)
    spwInfo = ms.getspectralwindowinfo()
    nchan = spwInfo["0"]["NumChan"]
    npol = spwInfo["0"]["NumCorr"]
    ms.close()

    tb.open(vis)
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

    if debug:
        print uvw[0, :].max()
        print uvw[0, :].min()
        print uvw[1, :].max()
        print uvw[1, :].min()
        print uvw.dtype        # float64

    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    assert len(ant1) == nvis
    assert len(ant2) == nvis
    tb.done()

    if debug:
        print ant1.max()
        print ant1.min()
        print ant2.max()
        print ant2.min()


    if verbose:
        print("Number of coorelation: ", npol)
        print("data shape", data.shape)
        print("uvw shape", uvwShape)
        print("weight shpae", weight.shape)     # (npol, nvis)

    tb.open(vis + '/SPECTRAL_WINDOW')
    SPWFreqs = np.squeeze(tb.getcol("CHAN_FREQ"))
    tb.done()

    freq_per_vis = np.array([SPWFreqs[fff] for fff in data_desc_id])
    # freqs = np.mean(SPWFreqs)
    assert len(freq_per_vis) == nvis
    del data_desc_id, SPWFreqs

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
        weight = np.average(weight, axis=0)

        real = data.real
        imag = data.imag
        del data
        print "Shape after averaging two hands: ", real.shape

        # before rescaling weights:
        # plot using all the points
        # vis = np.array(zip(real[0, :], imag[0, :])).flatten()
        # yyy = calc_img_from_vis(uvw[0, :], uvw[1, :], weight[0, :], vis, freq_per_vis, 1200, 0.05)
        # import pdb; pdb.set_trace()

        # # only a subset of points
        # for _ in range(5):       # make 5 images, each with Npointss
        #     from ripples_utils import pick_highSNR
        #     print weight.shape
        #     _idx = pick_highSNR(weight, bestNsam=8000, plothist=False)
        #     vis = np.array(zip(real[0, _idx], imag[0, _idx])).flatten()
        #     calc_img_from_vis(uvw[0, _idx], uvw[1, _idx], weight[0, _idx], vis, freq_per_vis[_idx], 1200, 0.05)
        #     import pdb; pdb.set_trace()
    elif npol > 2:
        raise NotImplementedError("more than 2 hands..")


    print 1./np.sqrt(np.sum(weight))

    # Yashar way
    # first grouping the visibilities into bins that probe the same signal
    # take differences btw visibilities that null the sky
    # then, we simply assume that the variance in the subtracted visibilities is equal to twice the noise variance

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
    sigma = weight**-0.5
    scaling = np.asarray(scaling).mean()
    sigma *= scaling
    print 'Scaling factor: ', scaling
    print 'Sigma after scaling [mJy/beam]: ', ((sigma**-2).sum())**-0.5 * 1E+3

    # debug, total noise after scaling
    print 1./np.sqrt(np.sum(weight))

    # real_1, imag_1, real_2, imag_2, etc
    visOut = np.array(zip(real, imag)).flatten()
    assert len(visOut) == int(nvis*2)
    weight = sigma**-2
    weight = np.array(zip(weight, weight)).flatten()    # TODO: want to make sure this works for nchan > 1
    assert len(weight) == len(visOut)

    blah = np.zeros((nvis))

    if not Nsam:
        with open(os.path.join(savepath, 'u.bin'), 'wb') as f:
            f.write(uvw[0, :])
        with open(os.path.join(savepath, 'v.bin'), 'wb') as f:
            f.write(uvw[1, :])
        with open(os.path.join(savepath, 'ant1.bin'), 'wb') as f:
            ant1 = np.asarray(ant1, dtype=np.float_)
            # print ant1.dtype
            f.write(ant1)
        with open(os.path.join(savepath, 'ant2.bin'), 'wb') as f:
            ant2 = np.asarray(ant2, dtype=np.float_)
            f.write(ant2)
        with open(os.path.join(savepath, 'frequencies.bin'), 'wb') as f:
            f.write(freq_per_vis)
        with open(os.path.join(savepath, 'vis_chan_0.bin'), 'wb') as f:
            f.write(visOut)
        with open(os.path.join(savepath, 'sigma_squared_inv.bin'), 'wb') as f:
            f.write(weight)
        with open(os.path.join(savepath, 'chan.bin'), 'wb') as f:
            f.write(blah)
    else:
        from ripples_utils import pick_highSNR
        print weight.shape
        _idx = pick_highSNR(weight[::2], bestNsam=Nsam, plothist=False)    # already zipped the weight

        with open(os.path.join(savepath, 'u.bin'), 'wb') as f:
            f.write(uvw[0, _idx])
        with open(os.path.join(savepath, 'v.bin'), 'wb') as f:
            f.write(uvw[1, _idx])
        with open(os.path.join(savepath, 'ant1.bin'), 'wb') as f:
            ant1 = np.asarray(ant1[_idx], dtype=np.float_)
            # print ant1.dtype
            f.write(ant1)
        with open(os.path.join(savepath, 'ant2.bin'), 'wb') as f:
            ant2 = np.asarray(ant2[_idx], dtype=np.float_)
            f.write(ant2)
        with open(os.path.join(savepath, 'frequencies.bin'), 'wb') as f:
            f.write(freq_per_vis[_idx])
        visOut = np.array(zip(real[0, _idx], imag[0, _idx])).flatten()
        with open(os.path.join(savepath, 'vis_chan_0.bin'), 'wb') as f:
            f.write(visOut)
        with open(os.path.join(savepath, 'sigma_squared_inv.bin'), 'wb') as f:
            f.write(weight[_idx])
        with open(os.path.join(savepath, 'chan.bin'), 'wb') as f:
            f.write(blah[_idx])

    # After rescaling weights:
    # only a subset of points
    for _ in range(5):       # make 5 images, each with Npointss
        from ripples_utils import pick_highSNR
        print weight.shape
        _idx = pick_highSNR(weight[::2], bestNsam=4000, plothist=False)    # already zipped the weight
        vis = np.array(zip(real[0, _idx], imag[0, _idx])).flatten()
        calc_img_from_vis(uvw[0, _idx], uvw[1, _idx], weight[_idx], vis, freq_per_vis[_idx], 1200, 0.01)
        import pdb; pdb.set_trace()

    return None




def check_vis_fixvis(msFile, msFileFixvis):
    """
    compare visibilities from msFile before and after running fixvis()

    According to documentation, fixbis will recalculate the u,v,w coordinates relative to the new phase center.
    https://casa.nrao.edu/casadocs/casa-5.1.1/uv-manipulation/recalculation-of-uvw-values-fixvis

    """

    tb.open(msFile)
    data1 = tb.getcol('DATA')
    uvw1 = tb.getcol('UVW')
    tb.close()

    tb.open(msFileFixvis)
    data0 = tb.getcol('DATA')
    uvw0 = tb.getcol('UVW')
    tb.close()

    print data1[:100]
    print data0[:100]

    print uvw1[0, :100]
    print uvw0[0, :100]

    print uvw1[1, :100]
    print uvw0[1, :100]

    return None


def get_source_name(msFile):
    tb.open(msFile+'/SOURCE')
    objname = tb.getcol("NAME")
    objname = objname[0]
    tb.done()
    return objname


if __name__ == '__main__':
    import os
    import config

    ccc = config.Configurable('ripples_prep.yaml')
    setupParam = ccc.config_dict
    savepath = setupParam['cont']['savepath']
    overwrite = setupParam['cont']['overwrite']
    contMS = setupParam['cont']['data']['contMS']
    shiftpcd = setupParam['shiftpcd']['shiftpcd']
    timebin = setupParam['cont']['timebin']
    timebinsec = setupParam['cont']['timebinsec']
    nsam = setupParam['cont']['nsam']
    if nsam == 'None':
        nsam = None

    if not os.path.exists(savepath):
        os.mkdir(savepath)

    listobs(contMS)

    if shiftpcd:
        field = setupParam['shiftpcd']['field']
        pcd = setupParam['shiftpcd']['pcd']
        outvis = shift_phs_center(contMS, pcd, field)
        outpcd = get_phs_center(outvis)
        print "Phase Center to be shifted to", pcd
        print "Phase Center after fixvis: ", outpcd
        contMS = outvis

    # if already shifted phase center, then use the MS file, regardless of whether we are setting shiftpcd in .yaml
    if not shiftpcd:
        outvis = contMS.replace('.ms', 'phsShift.ms')

    # set contMS to the phsShift.ms
    if os.path.exists(outvis):
        contMS = outvis

        # debug, image contMS with CASA
        if not os.path.exists(contMS[:contMS.find('.ms')] + '.image'):
            clean(vis=contMS,
                  imagename=contMS[:contMS.find('.ms')],
                  spw='',
                  mode='mfs',
                  nchan=-1,
                  imsize=800,
                  cell='0.0300arcsec',
                  niter=0,
                  interactive=False,
                  stokes='I')

    ms2ripples(contMS, savepath, timebinsec, nsam, timebin, \
               overwrite=overwrite, verbose=False, debug=False)        # if debug=True, will also image the saved visibilities using python.

    # check fixvis()
    # ms1 = 'calibrated.LSRK_cont_timebin120s.ms'    # created manually for debugging
    # ms2 = 'ripples/calibrated.LSRK_contphsShift_timebin120s.ms'
    # check_vis_fixvis(ms1, ms2)

    objname = get_source_name(contMS)
    tarballName = objname + '_cont_timebin' + str(timebinsec) +'s.tgz'

    # change to ripples dir
    _current = os.getcwd()
    os.chdir(savepath)
    os.system('tar -zcvf '+ tarballName + ' *bin ')
    # p = subprocess.Popen(["scp -p 61022", tarballName, uploadDestination])
    # sts = p.wait()
    os.chdir(_current)
    print tarballName
    # On Voms: scp dleung@serenity.astro.cornell.edu:/data/dleung/DATA/ALMA/QSO_c5/RXJ0911/calibrated/ripples/RXJ0911_cont_timebin60s.tgz .

    from ripples_utils import test_image_vis
    img = test_image_vis(savepath, imsize=1200, pixsize_arcsec=0.01)
