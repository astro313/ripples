"""

Script to prepare VIS for Ripples


"""

from ripples_utils import calc_img_from_vis, check_binFiles
import os

def timebinning(vis, savepath, timebinsec, timebin=False, overwrite=True, debug=False):

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


def ms2ripples_yashar_getFromSPW(spwlist, visname='test/calibrated.LSRK_contphsShift_timebin660s.ms', outdir='test/', debug=False):
    """
    Get data from given spw, using ms tool

    NOTE
    ----
    Output checked with Yashar

    """

    from array import array
    import numpy as np
    import os

    ms.open(visname, nomodify=False)

    for i_spw in spwlist:
        ms.selectinit(datadescid=i_spw)
        recD = ms.getdata(["data"])
        aD = recD["data"]
        print aD.shape

        flah_d = ms.getdata(["flag"])
        flah_d = flah_d["flag"]

        # any flags?
        print (flah_d.flatten() == True).any()

        recS = ms.getdata(["sigma"])
        aS = recS["sigma"]     # print aS.shape
        weight = aS**-2

        # average two hands
        weight = np.average(weight, axis=0)
        print 1. / np.sqrt(np.sum(weight))

        anD1 = ms.getdata(["antenna1"])
        anD1 = anD1["antenna1"]

        anD2 = ms.getdata(["antenna2"])
        anD2 = anD2["antenna2"]

        UVW = ms.getdata(["UVW"])
        uvpoints = UVW["uvw"]
        u = uvpoints[0]
        v = uvpoints[1]
        nvis = len(v)
        print "nvis: {:}".format(nvis)
        print "nchan in {:}: {:}".format(i_spw, len(aD[0]))

        for j in range(0, len(aD[0])):

            arr = (aD[0][j] + aD[1][j]) / 2        # average the two hands
            float_array_real = array('d', arr.real)
            float_array_imag = array('d', arr.imag)
            vis = np.array(zip(float_array_real, float_array_imag)).flatten()

        # rescale weights
        # Yashar way
        # first grouping the visibilities into bins that probe the same signal
        # take differences btw visibilities that null the sky
        # then, we simply assume that the variance in the subtracted visibilities
        # is equal to twice the noise variance

        scaling = []
        for a1 in np.unique(anD1):
            for a2 in np.unique(anD2):
                if a1 < a2:
                    baselineX = (anD1 == a1) & (anD2 == a2)

                    if debug:
                        print a1, a2
                        print anD1, anD2
                        print ""
                        print a1 in anD1
                        print a2 in anD2
                        print ""

                        print np.where(a1 == anD1)
                        print np.where(a2 == anD2)
                        print np.where((anD1 == a1) & (anD2 == a2) == True)

                    # important line! if we are picking a subset of points with nsam
                    # since we may miss some baselines.
                    if baselineX.any() == True:
                        if len(aD[0]) == 1:
                            reals = arr.real[baselineX]
                            imags = arr.imag[baselineX]
                            sigs = weight[baselineX]**-0.5
                        else:
                            raise NotImplementedError(
                                "Not implemented for MS files with more than 1 channel per spw...")

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
        print 1. / np.sqrt(np.sum(weight))
        weight = sigma**-2
        weight = np.array(zip(weight, weight)).flatten()         # 2 * nvis elements

        if i_spw == spwlist[0]:
            w = weight
            print w.shape
            assert len(w) == 2 * nvis
            ant1 = anD1
            ant2 = anD2
            uu = u
            vv = v
            vvis = vis

        elif i_spw != spwlist[0]:
            uu = np.hstack((uu, u))
            vv = np.hstack((vv, v))
            ant1 = np.hstack((ant1, anD1))
            ant2 = np.hstack((ant2, anD2))
            w = np.hstack((w, weight))
            vvis = np.hstack((vvis, vis))
            print uu.shape
            print vv.shape
            print ant1.shape
            print ant2.shape
            print w.shape
            print vvis.shape
            assert len(vvis)/2 == len(uu)

    print("shapes before saving to bin files")
    print uu.shape
    print vv.shape
    print ant1.shape
    print ant2.shape
    print w.shape
    print vvis.shape

    assert len(w) == len(vvis)
    f = open(os.path.join(outdir, 'sigma_squared_inv.bin'), 'wb')
    for i in range(len(w)):
        w[i].tofile(f)
    f.close()

    output_file = open(os.path.join(outdir, 'vis_chan_0.bin'), 'wb')
    vvis.tofile(output_file)
    output_file.close()

    f = open(os.path.join(outdir, 'u.bin'), 'wb')
    for i in range(0, len(uu)):
        uu[i].tofile(f)
    f.close()

    f = open(os.path.join(outdir, 'v.bin'), 'wb')
    for i in range(0, len(vv)):
        vv[i].tofile(f)
    f.close()

    f = open(os.path.join(outdir, 'ant1.bin'), 'wb')
    ant1 = np.asarray(ant1, dtype=np.float_)
    for i in range(0, len(ant1)):
        ant1[i].tofile(f)
    f.close()

    f = open(os.path.join(outdir, 'ant2.bin'), 'wb')
    ant2 = np.asarray(ant2, dtype=np.float_)
    for i in range(0, len(ant2)):
        ant2[i].tofile(f)
    f.close()

    total_nvis = len(uu)
    blah = np.zeros((total_nvis))
    f = open(os.path.join(outdir, 'chan.bin'), 'wb')
    for i in range(len(blah)):
        blah[i].tofile(f)
    f.close()

    # added
    rec = ms.getdata(["axis_info"])
    chan_freqs = np.squeeze(rec["axis_info"]["freq_axis"]["chan_freq"])    # [chan]
    nchan = len(rec['axis_info']['freq_axis']['chan_freq'][:, 0])
    assert len(aD[0]) == nchan
    print chan_freqs

    freq_per_vis = np.array([chan_freqs for _ in range(total_nvis)])

    f = open(os.path.join(outdir, 'frequencies.bin'), 'wb')
    for i in range(len(freq_per_vis)):
        freq_per_vis[i].tofile(f)
    f.close()


if __name__ == '__main__':

    from ripples_utils import test_image_vis

    # using ms tool to extract data from all spw
    ms2ripples_yashar_getFromSPW(range(12), 'test/calibrated.LSRK_contphsShift_timebin660s.ms', 'test/', debug=False)
    check_binFiles('test/')
    test_image_vis('test', Nsam=4000)
    test_image_vis('test')  # all points

    # after running mstransform --> 1 SPWs
    ms2ripples_yashar_getFromSPW([0], 'test/calibrated.LSRK_contphsShift_timebin660s_1spw_1chan.ms', 'test_mstransform/')
