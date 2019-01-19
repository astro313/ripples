"""
Doesn't make any sense why the extracted data look like noise...

"""

from ripples_utils import calc_img_from_vis, check_binFiles
import os


def ms2ripples_debug(vis='calibrated.LSRK_contphsShift.ms', savepath='test', verbose=True, debug=False):

    # doesn't work with CASA v5.1, didn't test with casa 4.7 because numpy version issues
    # test under casa v4.5

    if debug:
        plotms(vis, xaxis='real', yaxis='imag', timebin=100)
        plotms(vis, xaxis='UVdist', yaxis='amp', timebin=100)

    import astropy.constants as c
    import numpy as np

    timebinsec = 1200      # actually same as if timebinsec = 660...
    print("Performing time binning, by {:}s:".format(timebinsec))
    prefix = vis.replace('.ms', '_timebin' + '{:}s.ms'.format(timebinsec))
    prefix = os.path.basename(prefix)
    timebinVis = os.path.join(savepath, prefix)
    # rmtables(timebinVis)
    default(split2)
    split2(vis=vis,
           timebin=str(timebinsec) + 's',
           outputvis=timebinVis,
           datacolumn='data',
           keepflags=False
           )

    vis = timebinVis

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
    weight = 1. / sigma**2
    # del sigma

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

    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    assert len(ant1) == nvis
    assert len(ant2) == nvis
    tb.done()

    from sys import getsizeof
    with open(os.path.join(savepath, 'ant1.bin'), 'wb') as f:
        print(len(ant1))
        print getsizeof(ant1)
        ant1 = np.asarray(ant1, dtype=np.float_)
        print getsizeof(ant1)
        # print ant1.dtype
        f.write(ant1)
    with open(os.path.join(savepath, 'ant2.bin'), 'wb') as f:
        print(len(ant2))
        print getsizeof(ant2)
        ant2 = np.asarray(ant2, dtype=np.float_)
        print getsizeof(ant2)
        f.write(ant2)

    # import pdb; pdb.set_trace()

    # ant1.tofile(os.path.join(savepath, 'ant1.bin'))

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

    with open(os.path.join(savepath, 'frequencies.bin'), 'wb') as f:
        f.write(freq_per_vis)
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
            import pdb
            pdb.set_trace()

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

    print 1. / np.sqrt(np.sum(weight))

    # Yashar way
    # first grouping the visibilities into bins that probe the same signal
    # take differences btw visibilities that null the sky
    # then, we simply assume that the variance in the subtracted visibilities
    # is equal to twice the noise variance

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

                # important line! if we are picking a subset of points with nsam
                # since we may miss some baselines.
                if baselineX.any() == True:
                    if nchan == 1:
                        reals = real[0, baselineX]
                        imags = imag[0, baselineX]
                        sigs = weight[0, baselineX]**-0.5
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
    del ant1, ant2
    sigma = weight**-0.5
    scaling = np.asarray(scaling).mean()
    sigma *= scaling
    print 'Scaling factor: ', scaling
    print 'Sigma after scaling [mJy/beam]: ', ((sigma**-2).sum())**-0.5 * 1E+3

    # debug, total noise after scaling
    print 1. / np.sqrt(np.sum(weight))

    # real_1, imag_1, real_2, imag_2, etc
    visOut = np.array(zip(real, imag)).flatten()
    assert len(visOut) == int(nvis * 2)
    weight = sigma**-2
    # TODO: want to make sure this works for nchan > 1
    weight = np.array(zip(weight, weight)).flatten()
    assert len(weight) == len(visOut)

    with open(os.path.join(savepath, 'vis_chan_0.bin'), 'wb') as f:
        f.write(visOut)
    with open(os.path.join(savepath, 'sigma_squared_inv.bin'), 'wb') as f:
        f.write(weight)

    blah = np.zeros((nvis))
    with open(os.path.join(savepath, 'chan.bin'), 'wb') as f:
        f.write(blah)

    # After rescaling weights:
    # only a subset of points
    for _ in range(5):       # make 5 images, each with Npointss
        from ripples_utils import pick_highSNR
        print weight.shape
        # already zipped the weight
        _idx = pick_highSNR(weight[::2], bestNsam=4000, plothist=False)
        vis = np.array(zip(real[0, _idx], imag[0, _idx])).flatten()
        calc_img_from_vis(uvw[0, _idx], uvw[1, _idx], weight[
                          _idx], vis, freq_per_vis[_idx], 1200, 0.01)
        import pdb
        pdb.set_trace()

    return None


def image_vis(binpath='../', imsize=800, pixsize_arcsec=0.01, Nsam=None):
    import numpy as np
    import os
    from ripples_utils import calc_img_from_vis

    # binary
    uu = np.fromfile(os.path.join(binpath, 'u.bin'))
    vv = np.fromfile(os.path.join(binpath, 'v.bin'))
    weight = np.fromfile(os.path.join(binpath, 'sigma_squared_inv.bin'))
    freqs = np.fromfile(os.path.join(binpath, 'frequencies.bin'))
    dataVis = np.fromfile(os.path.join(binpath, 'vis_chan_0.bin'))

    if Nsam:
        from ripples_utils import pick_highSNR
        for _ in range(5):       # make 5 images, each with Npointss
            # already zipped the weight
            _idx = pick_highSNR(weight[::2], bestNsam=Nsam, plothist=False)
            _vis_real = [dataVis[iii*2: 2*iii+1][0] for iii in _idx]
            _vis_imag = [dataVis[iii*2+1: 2*iii+2][0] for iii in _idx]
            _vis = np.array(zip(_vis_real, _vis_imag)).flatten()

            im_jybeam = calc_img_from_vis(uu[_idx], vv[_idx], weight[
                          _idx * 2], _vis, freqs[_idx], imsize, pixsize_arcsec)
            import pdb; pdb.set_trace()
    else:
        im_jybeam = calc_img_from_vis(
            uu, vv, weight, dataVis, freqs, imsize, pixsize_arcsec)
    return im_jybeam

# --------------------------------- ms tool -----------------------------------

def ms2ripples_yashar(visname='test/calibrated.LSRK_contphsShift_timebin660s.ms', i_spw=0):
    """
        NOTE
        ----
        Output checked with Yashar

    """
    from array import array
    import numpy as np

    ms.open(visname, nomodify=False)

    ms.selectinit(datadescid=i_spw)
    recD = ms.getdata(["data"])
    aD = recD["data"]

    time = ms.getdata(["time"])
    time = time["time"]

    flah_d = ms.getdata(["flag"])
    flah_d = flah_d["flag"]

    recS = ms.getdata(["sigma"])
    aS = recS["sigma"]     # print aS.shape

    anD1 = ms.getdata(["antenna1"])
    anD1 = anD1["antenna1"]

    anD2 = ms.getdata(["antenna2"])
    anD2 = anD2["antenna2"]

    UVW = ms.getdata(["UVW"])
    uvpoints = UVW["uvw"]
    u = uvpoints[0]
    v = uvpoints[1]
    nvis = len(v)
    print "nchan: ", len(aD[0])

    for j in range(0, len(aD[0])):
        output_file = open('test/vis_chan_0.bin', 'wb')

        arr = (aD[0][j] + aD[1][j]) / 2        # average the two hands
        float_array_real = array('d', arr.real)
        float_array_imag = array('d', arr.imag)
        print float_array_real[:5], float_array_imag[:5]
        vis = np.array(zip(float_array_real, float_array_imag)).flatten()
        print vis[0:10]

        vis.tofile(output_file)
        output_file.close()

        # read back in
        xx = np.fromfile('test/vis_chan_0.bin')
        assert len(xx) == nvis * 2
        assert (xx == vis).all()


    f = open('test/u.bin', 'wb')
    for i in range(0, len(u)):
        u[i].tofile(f)
    f.close()

    xx = np.fromfile('test/u.bin')
    assert len(xx) == nvis
    assert (xx == u).all()


    f = open('test/v.bin', 'wb')
    for i in range(0, len(v)):
        v[i].tofile(f)
    f.close()

    xx = np.fromfile('test/v.bin')
    assert len(xx) == nvis
    assert (xx == v).all()

    f = open('test/ant1.bin', 'wb')
    anD1 = np.asarray(anD1, dtype=np.float_)
    for i in range(0, len(anD1)):
        anD1[i].tofile(f)
    f.close()

    xx = np.fromfile('test/ant1.bin')
    assert len(xx) == nvis
    assert (xx == anD1).all()


    f = open('test/ant2.bin', 'wb')
    anD2 = np.asarray(anD2, dtype=np.float_)
    for i in range(0, len(anD2)):
        anD2[i].tofile(f)
    f.close()

    xx = np.fromfile('test/ant2.bin')
    assert len(xx) == nvis
    assert (xx == anD2).all()


    blah = np.zeros((nvis))
    f = open('test/chan.bin', 'wb')
    for i in range(len(blah)):
        blah[i].tofile(f)
    f.close()

    xx = np.fromfile('test/chan.bin')
    assert len(xx) == nvis
    assert (xx == blah).all()


    # added
    rec = ms.getdata(["axis_info"])
    chan_freqs = np.squeeze(rec["axis_info"]["freq_axis"]["chan_freq"])    # [chan]
    nchan = len(rec['axis_info']['freq_axis']['chan_freq'][:, 0])
    assert len(aD[0]) == nchan
    print chan_freqs

    freq_per_vis = np.array([chan_freqs for _ in range(nvis)])
    assert len(freq_per_vis) == nvis

    f = open('test/frequencies.bin', 'wb')
    for i in range(len(freq_per_vis)):
        freq_per_vis[i].tofile(f)
    f.close()

    # w/o rescaling for now
    weight = aS**-2
    print weight.shape
    weight = np.average(weight, axis=0)
    print weight[:10]
    weight = np.array(zip(weight, weight)).flatten()
    print weight[:20]

    f = open('test/sigma_squared_inv.bin', 'wb')
    for i in range(len(weight)):
        weight[i].tofile(f)
    f.close()



def ms2ripples_yashar_getallspw(spwlist, visname='test/calibrated.LSRK_contphsShift_timebin660s.ms'):
    """
    Get data from all the spw

    NOTE
    ----
    Output checked with Yashar

    """

    from array import array
    import numpy as np

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
        weight = np.average(weight, axis=0)
        weight = np.array(zip(weight, weight)).flatten()         # 2 * nvis elements

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

    # w/o rescaling for now
    f = open('test/sigma_squared_inv.bin', 'wb')
    for i in range(len(w)):
        w[i].tofile(f)
    f.close()

    output_file = open('test/vis_chan_0.bin', 'wb')
    vvis.tofile(output_file)
    output_file.close()

    f = open('test/u.bin', 'wb')
    for i in range(0, len(uu)):
        uu[i].tofile(f)
    f.close()

    f = open('test/v.bin', 'wb')
    for i in range(0, len(vv)):
        vv[i].tofile(f)
    f.close()

    f = open('test/ant1.bin', 'wb')
    ant1 = np.asarray(ant1, dtype=np.float_)
    for i in range(0, len(ant1)):
        ant1[i].tofile(f)
    f.close()


    f = open('test/ant2.bin', 'wb')
    ant2 = np.asarray(ant2, dtype=np.float_)
    for i in range(0, len(ant2)):
        ant2[i].tofile(f)
    f.close()

    total_nvis = len(uu)
    blah = np.zeros((total_nvis))
    f = open('test/chan.bin', 'wb')
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

    f = open('test/frequencies.bin', 'wb')
    for i in range(len(freq_per_vis)):
        freq_per_vis[i].tofile(f)
    f.close()


def ms2ripples_yashar_getallspw_mstransform(spwlist=[0], visname='test/calibrated.LSRK_contphsShift_timebin660s_1spw_1chan.ms'):
    """
    Get data from all the spw, after running mstranform to combine all SPWs

    """

    from array import array
    import numpy as np

    ms.open(visname, nomodify=False)

    for i_spw in spwlist:
        ms.selectinit(datadescid=i_spw)
        recD = ms.getdata(["data"])
        aD = recD["data"]
        print aD.shape

        flah_d = ms.getdata(["flag"])
        flah_d = flah_d["flag"]

        recS = ms.getdata(["sigma"])
        aS = recS["sigma"]     # print aS.shape
        weight = aS**-2
        weight = np.average(weight, axis=0)
        weight = np.array(zip(weight, weight)).flatten()         # 2 * nvis elements

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

        print u.shape
        print v.shape
        print anD1.shape
        print anD2.shape
        print weight.shape
        print vvis.shape
        assert len(vvis)/2 == len(u)

        # (26730,)
        # (26730,)
        # (26730,)
        # (26730,)
        # (53460,)
        # (53460,)

    # w/o rescaling for now
    f = open('test_mstranform/sigma_squared_inv.bin', 'wb')
    for i in range(len(w)):
        w[i].tofile(f)
    f.close()

    output_file = open('test_mstranform/vis_chan_0.bin', 'wb')
    vvis.tofile(output_file)
    output_file.close()

    f = open('test_mstranform/u.bin', 'wb')
    for i in range(0, len(u)):
        u[i].tofile(f)
    f.close()

    f = open('test_mstranform/v.bin', 'wb')
    for i in range(0, len(v)):
        v[i].tofile(f)
    f.close()

    f = open('test_mstranform/ant1.bin', 'wb')
    anD1 = np.asarray(anD1, dtype=np.float_)
    for i in range(0, len(anD1)):
        anD1[i].tofile(f)
    f.close()


    f = open('test_mstranform/ant2.bin', 'wb')
    anD2 = np.asarray(anD2, dtype=np.float_)
    for i in range(0, len(anD2)):
        anD2[i].tofile(f)
    f.close()

    blah = np.zeros((nvis))
    f = open('test_mstranform/chan.bin', 'wb')
    for i in range(len(blah)):
        blah[i].tofile(f)
    f.close()

    # added
    rec = ms.getdata(["axis_info"])
    chan_freqs = np.squeeze(rec["axis_info"]["freq_axis"]["chan_freq"])    # [chan]
    nchan = len(rec['axis_info']['freq_axis']['chan_freq'][:, 0])
    assert len(aD[0]) == nchan
    print chan_freqs
    import pdb; pdb.set_trace()

    freq_per_vis = np.array([chan_freqs for _ in range(nvis)])

    f = open('test_mstranform/frequencies.bin', 'wb')
    for i in range(len(freq_per_vis)):
        freq_per_vis[i].tofile(f)
    f.close()

    import pdb; pdb.set_trace()


def compare_tb_ms_uvw_data(vis):
    """
    compare the data points extracted using ms tool versus tb tool

    """

    import numpy as np
    tb.open(vis)
    # tb.colnames
    data = tb.getcol('DATA')
    uvw = tb.getcol('UVW')
    nvis = len(uvw[1, :])
    tb.close()
    print data.shape, data.min(), data.max()
    print uvw.shape
    print uvw[0, :].min(), uvw[0, :].max()
    print uvw[1, :].min(), uvw[1, :].max()
    print nvis

    # ms
    i_spw = 0
    ms.open(vis, nomodify=False)

    ms.selectinit(datadescid=i_spw)
    recD = ms.getdata(["data"])
    aD = recD["data"]

    UVW = ms.getdata(["UVW"])
    uvpoints = UVW["uvw"]
    u = uvpoints[0]
    v = uvpoints[1]
    nvis = len(v)

    print aD.shape, aD.min(), aD.max()
    print u.min(), u.max()
    print v.min(), v.max()
    print nvis

    # (2, 1, 103374) (-0.0166221950203-0.00523133110255j) (0.0200791154057-0.0129184629768j)
    # (3, 103374)
    # -1632.71137651 1877.10144335
    # -1802.52746357 2074.22241246
    # 103374

    # (2, 1, 8514) (-0.0134982932359-0.00121319107711j) (0.0133786378428+6.72544483677e-05j)
    # -1380.49549771 1854.19960841
    # -1802.52746357 2064.90236136
    # 8514

    # for nspw = 12, total nvis ~ 9000 * 12, consistent with tb tool

    import matplotlib.pyplot as plt
    plt.hist(np.average(data, axis=0).flatten(), bins=100)
    plt.title('from tbtool')
    plt.show(block=False)
    plt.savefig('test/data_tb.pdf')

    plt.figure()
    plt.hist(np.average(aD, axis=0).flatten(), bins=100)
    plt.title('from mstool')
    plt.show(block=False)
    plt.savefig('test/data_ms.pdf')

    plt.figure()
    plt.plot(u, v, 'ro', alpha=0.7)
    plt.title('from mstool')
    plt.show(block=False)
    plt.savefig('test/uv_ms.pdf')

    plt.figure()
    plt.plot(uvw[0, :], uvw[1, :], 'bo', alpha=0.3)
    plt.title('from tbtool')
    plt.show(block=False)
    plt.savefig('test/uv_tb.pdf')


    # compare with .uvfits
    exportuvfits(vis, vis.replace('.ms', '.uvfits'))
    import astropy.io.fits as fits
    xx = fits.open(vis.replace('.ms', '.uvfits'))
    dd = xx[0].data
    uu = dd['UU']
    vv = dd['VV']
    data = dd['DATA']
    nvis = data.shape[0]
    nspw = data.shape[3]
    nfreq = data.shape[4]
    npol = data.shape[5]

    print uu.shape, vv.shape
    print data.shape
    print nvis

    # (102512,) (102512,)
    # (102512, 1, 1, 12, 1, 2, 3)     # last col (real, imag, weight)
    # 102512

    data = np.average(data, axis=5)
    real = data[:, :, :, :, :, 0]
    imag = data[:, :, :, :, :, 1]
    data_complex = real + 1j * imag
    print data_complex.shape        # (102512, 1, 1, 12, 1)
    print data_complex[:, :, :, 0, :].min(), data_complex[:, :, :, 0, :].max()
    print data_complex[:, :, :, 1, :].min(), data_complex[:, :, :, 1, :].max()
    print data_complex[:, :, :, 2, :].min(), data_complex[:, :, :, 2, :].max()
    import pdb; pdb.set_trace()

    # very puzzling why there are nspw * nvis = 12 * 102512 points? Even more than what tb tool and ms tool return?

    # are most of the zeros?
    print np.squeeze((data_complex[1000:3000, :, :, 0, :].real > 0)).any()
    print np.squeeze((data_complex[17000:20000, :, :, 0, :].real > 0)).any()



def combinespw_tbtool_uvfits(vis):

    """

    Note
    ----
    mstransform combinespws=True without regridding, it only combines spws
    that are identical with each other. You have 3 sets of 4 spws (one for each observation ID),
    so it creates 4 channels in the single output spw. If you want to combine everything, you have
    to define more parameters, with regridms.
    """

    listobs(vis, listfile=vis.replace('.ms', '.log'), overwrite=True)
    vis1spw = vis[:vis.find('.ms')] + '_1spw.ms'
    mstranform(vis, vis1spw, datacolumn='data', combinewspw=True)   # num of unique spw will become the number of channels
    # also combine "all the channels"
    vis1spw1chan = vis[:vis.find('.ms')] + '_1spw_1chan.ms'
    mstransform(vis1spw, vis1spw1chan, datacolumn='data', combinespws=True, \
                regridms=True, nchan=1, width=4)
    listobs(vis1spw1chan, listfile=vis1spw1chan.replace('.ms', '.log'), \
            overwrite=True)

    import numpy as np
    tb.open(vis1spw1chan)
    # tb.colnames
    data = tb.getcol('DATA')
    uvw = tb.getcol('UVW')
    nvis = len(uvw[1, :])
    tb.close()
    print data.shape, data.min(), data.max()
    print uvw.shape
    print uvw[0, :].min(), uvw[0, :].max()
    print uvw[1, :].min(), uvw[1, :].max()
    print nvis
    # (2, 1, 26730) (-0.0106927361339+0.00347894337028j) (0.0101391784847+0.00074283964932j)
    # (3, 26730)
    # -1632.71137651 1877.10135185
    # -1802.52746357 2074.22084676
    # 26730

    # ms
    i_spw = 0
    ms.open(vis1spw1chan, nomodify=False)

    ms.selectinit(datadescid=i_spw)
    recD = ms.getdata(["data"])
    aD = recD["data"]

    UVW = ms.getdata(["UVW"])
    uvpoints = UVW["uvw"]
    u = uvpoints[0]
    v = uvpoints[1]
    nvis = len(v)

    print aD.shape, aD.min(), aD.max()
    print u.min(), u.max()
    print v.min(), v.max()
    print len(aD[0])
    print nvis

    # (2, 1, 26730) (-0.0106927361339+0.00347894337028j) (0.0101391784847+0.00074283964932j)
    # -1632.71137651 1877.10135185
    # -1802.52746357 2074.22084676
    # 1
    # 26730


    # compare with .uvfits
    exportuvfits(vis1spw1chan, vis1spw1chan.replace('.ms', '.uvfits'))
    import astropy.io.fits as fits
    xx = fits.open(vis1spw1chan.replace('.ms', '.uvfits'))
    dd = xx[0].data
    uu = dd['UU']
    vv = dd['VV']
    data = dd['DATA']
    nvis = data.shape[0]
    nspw = data.shape[3]
    nfreq = data.shape[4]
    npol = data.shape[5]

    print uu.shape, vv.shape
    print data.shape
    print nvis

    # (26730,) (26730,)
    # (26730, 1, 1, 1, 1, 2, 3)
    # 26730

    data = np.average(data, axis=5)
    real = data[:, :, :, :, :, 0]
    imag = data[:, :, :, :, :, 1]
    data_complex = real + 1j * imag
    print data_complex.shape        # (26730, 1, 1, 1, 1)
    nspw = data_complex.shape[3]
    nchan = data_complex.shape[4]

    import pdb; pdb.set_trace()


    return None



# my way: using tb tool
# ms2ripples_debug(savepath='test', verbose=False)
# image_vis(binpath='test')

# using ms tool
ms2ripples_yashar()
check_binFiles('test/')
image_vis('test')

# compare ms tool and tb tool and uvfits
compare_tb_ms_uvw_data(vis='test/calibrated.LSRK_contphsShift_timebin660s.ms')
combinespw_tbtool_uvfits('test/calibrated.LSRK_contphsShift_timebin660s.ms')


# using ms tool to extract data from all spw
ms2ripples_yashar_getallspw(range(12))
check_binFiles('test/')
image_vis('test', Nsam=4000)



# after running mstransform --> 1 SPWs
ms2ripples_yashar_getallspw_mstransform()

