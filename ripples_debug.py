"""
Doesn't make any sense why the extracted data from ripples_prep.py look like noise... from ripples_prep.py

It's solved using ms tool in this script, which will replace ripples_prep.py soon...



"""

from ripples_utils import calc_img_from_vis, check_binFiles
import os


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



if __name__ == '__main__':

    from ripples_utils import test_image_vis

    # using ms tool to extract data from all spw
    ms2ripples_yashar_getFromSPW(range(12), 'test/calibrated.LSRK_contphsShift_timebin660s.ms', 'test/', debug=False)
    check_binFiles('test/')
    test_image_vis('test', Nsam=4000)
    test_image_vis('test')  # all points

    # compare ms tool and tb tool and uvfits
    compare_tb_ms_uvw_data(vis='test/calibrated.LSRK_contphsShift_timebin660s.ms')
    combinespw_tbtool_uvfits('test/calibrated.LSRK_contphsShift_timebin660s.ms')

    # after running mstransform --> 1 SPWs
    ms2ripples_yashar_getFromSPW([0], 'test/calibrated.LSRK_contphsShift_timebin660s_1spw_1chan.ms', 'test_mstransform/')
